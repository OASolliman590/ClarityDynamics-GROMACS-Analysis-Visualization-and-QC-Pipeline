"""RMSF-guided residue alignment and PCA index preparation."""

from __future__ import annotations

import csv
import json
import statistics
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .config import PCATerminalDetectionConfig


PROTEIN_CA_ATOMS = {"CA"}
NON_PROTEIN_RESNAMES = {
    "TIP3", "SOL", "WAT", "HOH", "SOD", "NA", "CLA", "CL", "K", "POT",
    "UNK", "LIG",
}


@dataclass
class ResidueRecord:
    """Residue identity in one numbering system."""
    order: int
    resid: int
    resname: str
    chain: str = ""
    segment: str = ""


@dataclass
class ResidueMapEntry:
    """Alignment between RMSF row, GRO residue, and native residue labels."""
    rmsf_row: int
    rmsf_resid: int
    rmsf_value: float
    gro_resid: int
    native_resid: int
    resname: str
    chain: str = ""
    segment: str = ""


@dataclass
class TerminalDetection:
    """Detected terminal trimming result."""
    core_start_row: int
    core_end_row: int
    gro_start_resid: int
    gro_end_resid: int
    native_start_resid: int
    native_end_resid: int
    n_terminal_trimmed: int
    c_terminal_trimmed: int
    threshold: float
    median: float
    mad: float
    n_residues: int
    native_segments: List[Dict[str, object]]


@dataclass
class PCAIndexPlan:
    """make_ndx plan for the RMSF-guided PCA C-alpha group."""
    core_group_id: int
    calpha_group_id: int
    core_group_name: str
    calpha_group_name: str
    gro_start_resid: int
    gro_end_resid: int
    commands: List[str]

    @property
    def stdin(self) -> str:
        return "\n".join(self.commands) + "\n"


def parse_rmsf_file(path: Path) -> List[Tuple[int, float]]:
    """Parse a two-column RMSF file as (residue_or_row, value)."""
    values: List[Tuple[int, float]] = []
    with Path(path).open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith(("#", "@")):
                continue
            parts = stripped.split()
            if len(parts) < 2:
                continue
            try:
                values.append((int(float(parts[0])), float(parts[1])))
            except ValueError:
                continue
    if not values:
        raise ValueError(f"No RMSF values found in {path}")
    return values


def parse_gro_residues(path: Path) -> List[ResidueRecord]:
    """Parse protein residues from a GRO file using C-alpha atoms."""
    records: List[ResidueRecord] = []
    seen = set()
    lines = Path(path).read_text().splitlines()
    for line in lines[2:-1]:
        if len(line) < 20:
            continue
        try:
            resid = int(line[0:5])
        except ValueError:
            continue
        resname = line[5:10].strip()
        atom = line[10:15].strip()
        if atom not in PROTEIN_CA_ATOMS or resname in NON_PROTEIN_RESNAMES:
            continue
        key = (resid, resname)
        if key in seen:
            continue
        seen.add(key)
        records.append(ResidueRecord(order=len(records) + 1, resid=resid, resname=resname))
    if not records:
        raise ValueError(f"No protein C-alpha residues found in {path}")
    return records


def parse_pdb_residues(path: Path) -> List[ResidueRecord]:
    """Parse native protein residue labels from a PDB file."""
    records: List[ResidueRecord] = []
    seen = set()
    with Path(path).open() as handle:
        for line in handle:
            record = line[0:6].strip()
            if record not in {"ATOM", "HETATM"} or len(line) < 26:
                continue
            atom = line[12:16].strip()
            if atom not in PROTEIN_CA_ATOMS:
                continue
            resname = line[17:20].strip()
            if resname in NON_PROTEIN_RESNAMES:
                continue
            try:
                resid = int(line[22:26])
            except ValueError:
                continue
            chain = line[21].strip()
            segment = line[72:76].strip() if len(line) >= 76 else ""
            key = (chain, segment, resid, resname)
            if key in seen:
                continue
            seen.add(key)
            records.append(
                ResidueRecord(
                    order=len(records) + 1,
                    resid=resid,
                    resname=resname,
                    chain=chain,
                    segment=segment,
                )
            )
    if not records:
        raise ValueError(f"No protein C-alpha residues found in {path}")
    return records


def build_residue_map(
    rmsf_values: Sequence[Tuple[int, float]],
    gro_residues: Sequence[ResidueRecord],
    pdb_residues: Optional[Sequence[ResidueRecord]] = None,
) -> List[ResidueMapEntry]:
    """Align RMSF rows to GRO and native residue labels by residue order."""
    if len(rmsf_values) != len(gro_residues):
        raise ValueError(
            "RMSF residue count does not match GRO protein C-alpha count: "
            f"{len(rmsf_values)} vs {len(gro_residues)}"
        )
    native_records = list(pdb_residues or [])
    if len(native_records) != len(gro_residues):
        native_records = list(gro_residues)

    entries: List[ResidueMapEntry] = []
    for idx, ((rmsf_resid, rmsf_value), gro, native) in enumerate(
        zip(rmsf_values, gro_residues, native_records),
        start=1,
    ):
        entries.append(
            ResidueMapEntry(
                rmsf_row=idx,
                rmsf_resid=rmsf_resid,
                rmsf_value=rmsf_value,
                gro_resid=gro.resid,
                native_resid=native.resid,
                resname=native.resname or gro.resname,
                chain=native.chain,
                segment=native.segment,
            )
        )
    return entries


def smooth_values(values: Sequence[float], window: int) -> List[float]:
    """Return centered moving-average values."""
    if window <= 1 or len(values) <= 2:
        return [float(v) for v in values]
    half = window // 2
    smoothed: List[float] = []
    for idx in range(len(values)):
        start = max(0, idx - half)
        end = min(len(values), idx + half + 1)
        smoothed.append(sum(values[start:end]) / (end - start))
    return smoothed


def robust_threshold(values: Sequence[float], multiplier: float) -> Tuple[float, float, float]:
    """Return threshold, median, and MAD scaled to standard-deviation units."""
    median = statistics.median(values)
    deviations = [abs(value - median) for value in values]
    raw_mad = statistics.median(deviations)
    mad = raw_mad * 1.4826
    if mad == 0:
        # Flat core: the standard deviation would be inflated by the very terminal
        # outliers we want to detect, so fall back to robust scales instead.
        ordered = sorted(values)
        n = len(ordered)
        iqr = ordered[(3 * (n - 1)) // 4] - ordered[(n - 1) // 4] if n > 1 else 0.0
        mad = 0.7413 * iqr
        if mad == 0:
            mad = max(0.05 * abs(median), 1e-6)
    threshold = median + multiplier * mad
    return threshold, median, mad


def _terminal_trim(values: Sequence[float], threshold: float, stable_run: int, from_end: bool, limit: int) -> int:
    scan = list(reversed(values)) if from_end else list(values)
    max_scan = min(len(scan), limit)
    for idx in range(max_scan):
        stable = scan[idx:idx + stable_run]
        if len(stable) == stable_run and all(value <= threshold for value in stable):
            return idx
    return max_scan


def detect_terminals(
    residue_map: Sequence[ResidueMapEntry],
    config: Optional[PCATerminalDetectionConfig] = None,
) -> TerminalDetection:
    """Detect terminal overhangs from RMSF values and return a core range."""
    cfg = config or PCATerminalDetectionConfig()
    if not residue_map:
        raise ValueError("Cannot detect terminals from an empty residue map")

    if cfg.core_residue_range:
        start_row, end_row = cfg.core_residue_range
        if start_row < 1 or end_row > len(residue_map) or start_row > end_row:
            raise ValueError(f"Invalid PCA core_residue_range: {cfg.core_residue_range}")
        values = [entry.rmsf_value for entry in residue_map]
        threshold, median, mad = robust_threshold(values, cfg.mad_multiplier)
    elif not cfg.enabled:
        start_row = 1
        end_row = len(residue_map)
        values = [entry.rmsf_value for entry in residue_map]
        threshold, median, mad = robust_threshold(values, cfg.mad_multiplier)
    else:
        values = smooth_values([entry.rmsf_value for entry in residue_map], cfg.smoothing_window)
        threshold, median, mad = robust_threshold(values, cfg.mad_multiplier)
        limit = max(1, int(len(values) * cfg.max_terminal_fraction))
        n_trim = _terminal_trim(values, threshold, cfg.stable_run, from_end=False, limit=limit)
        c_trim = _terminal_trim(values, threshold, cfg.stable_run, from_end=True, limit=limit)
        start_row = n_trim + 1
        end_row = len(values) - c_trim
        if start_row > end_row:
            start_row = 1
            end_row = len(values)

    start_entry = residue_map[start_row - 1]
    end_entry = residue_map[end_row - 1]
    return TerminalDetection(
        core_start_row=start_row,
        core_end_row=end_row,
        gro_start_resid=start_entry.gro_resid,
        gro_end_resid=end_entry.gro_resid,
        native_start_resid=start_entry.native_resid,
        native_end_resid=end_entry.native_resid,
        n_terminal_trimmed=start_row - 1,
        c_terminal_trimmed=len(residue_map) - end_row,
        threshold=threshold,
        median=median,
        mad=mad,
        n_residues=len(residue_map),
        native_segments=compress_native_segments(residue_map[start_row - 1:end_row]),
    )


def compress_native_segments(entries: Sequence[ResidueMapEntry]) -> List[Dict[str, object]]:
    """Compress contiguous native residue labels into display ranges."""
    if not entries:
        return []
    segments: List[Dict[str, object]] = []
    start = entries[0]
    previous = entries[0]
    for entry in entries[1:]:
        same_label = entry.chain == previous.chain and entry.segment == previous.segment
        contiguous = entry.native_resid == previous.native_resid + 1
        if same_label and contiguous:
            previous = entry
            continue
        segments.append(_segment_dict(start, previous))
        start = entry
        previous = entry
    segments.append(_segment_dict(start, previous))
    return segments


def _segment_dict(start: ResidueMapEntry, end: ResidueMapEntry) -> Dict[str, object]:
    return {
        "chain": start.chain,
        "segment": start.segment,
        "start": start.native_resid,
        "end": end.native_resid,
    }


def build_index_plan(
    detection: TerminalDetection,
    existing_group_count: int,
    core_group_name: str = "core_no_terminals",
    calpha_group_name: str = "core_no_terminals_Calpha",
) -> PCAIndexPlan:
    """Build make_ndx commands for the terminal-excluded C-alpha group."""
    core_group_id = existing_group_count
    calpha_group_id = existing_group_count + 1
    commands = [
        f"ri {detection.gro_start_resid}-{detection.gro_end_resid}",
        f"name {core_group_id} {core_group_name}",
        f"{core_group_id} & a CA",
        f"name {calpha_group_id} {calpha_group_name}",
        "q",
    ]
    return PCAIndexPlan(
        core_group_id=core_group_id,
        calpha_group_id=calpha_group_id,
        core_group_name=core_group_name,
        calpha_group_name=calpha_group_name,
        gro_start_resid=detection.gro_start_resid,
        gro_end_resid=detection.gro_end_resid,
        commands=commands,
    )


def count_index_groups(index_path: Path) -> int:
    """Count groups in a GROMACS .ndx file."""
    count = 0
    with Path(index_path).open() as handle:
        for line in handle:
            if line.lstrip().startswith("["):
                count += 1
    return count


def write_residue_map(entries: Iterable[ResidueMapEntry], output_path: Path) -> None:
    """Write residue-map CSV."""
    rows = list(entries)
    if not rows:
        return
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def write_terminal_detection(detection: TerminalDetection, output_path: Path) -> None:
    """Write terminal-detection metadata JSON."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(asdict(detection), indent=2))


def write_index_plan(plan: PCAIndexPlan, output_path: Path) -> None:
    """Write the make_ndx plan for review/reuse."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Generated PCA index plan",
        f"# Core GRO residue range: {plan.gro_start_resid}-{plan.gro_end_resid}",
        f"# Core group: {plan.core_group_id} {plan.core_group_name}",
        f"# C-alpha group: {plan.calpha_group_id} {plan.calpha_group_name}",
        "",
        *plan.commands,
        "",
    ]
    output_path.write_text("\n".join(lines))
