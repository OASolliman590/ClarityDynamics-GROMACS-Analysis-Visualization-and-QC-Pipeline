from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from gromacs_analysis.mmpbsa.parser import MMPBSAParser
from gromacs_analysis.mmpbsa.processor import MMPBSAProcessor
from gromacs_analysis.mmpbsa.plotter import MMPBSAPlotter
from gromacs_analysis.mmpbsa.analyzer import MMPBSAAnalyzer
from gromacs_analysis.mmpbsa.config import MMPBSAConfig, MMPBSASystemConfig
from gromacs_analysis.utils.publication import PublicationConfig
from gromacs_analysis.utils.palette import categorical


def test_csv_results_normalize_component_names_and_preserve_saved_samples(tmp_path):
    source = tmp_path / "results.csv"
    pd.DataFrame({
        "Frame #": [8, 3, 5],
        "ΔVDWAALS": [-2.0, -4.0, -6.0],
        "ΔEEL": [1.0, 2.0, 3.0],
        "ΔEGB": [-1.0, -2.0, -3.0],
        "ΔESURF": [-0.5, -0.5, -0.5],
        "ΔGGAS": [-1.0, -2.0, -3.0],
        "ΔGSOLV": [-1.5, -2.5, -3.5],
        "ΔTOTAL": [-2.5, -4.5, -6.5],
        "TOTAL": [-2.5, -4.5, -6.5],
    }).to_csv(source, index=False)

    parsed = MMPBSAParser().parse_results(source)

    canonical = parsed["canonical_components"]["VDWAALS"]
    assert (canonical["mean"], canonical["sample_sd"], canonical["n"]) == (-4.0, 2.0, 3)
    assert list(parsed["canonical_components"]) == ["VDWAALS", "EEL", "EGB", "ESURF", "GGAS", "GSOLV", "TOTAL"]
    assert list(parsed["frame_data"]["Frame #"]) == [8, 3, 5]
    assert list(parsed["frame_data"]["TOTAL"]) == [-2.5, -4.5, -6.5]


def test_publication_component_payload_separates_raw_grouped_and_total(tmp_path):
    processor = MMPBSAProcessor.__new__(MMPBSAProcessor)
    processor.config = None
    results = {
        "delta_components": {
            "ΔVDWAALS": {"mean": -2.0, "std": 1.0, "uncertainty": "frame_sd"},
            "ΔEEL": {"mean": 3.0, "std": 2.0, "uncertainty": "frame_sd"},
            "ΔEGB": {"mean": -4.0, "std": 1.5, "uncertainty": "frame_sd"},
            "ΔESURF": {"mean": -1.0, "std": 0.5, "uncertainty": "frame_sd"},
            "ΔGGAS": {"mean": 1.0, "std": 1.0, "uncertainty": "frame_sd"},
            "ΔGSOLV": {"mean": -5.0, "std": 2.0, "uncertainty": "frame_sd"},
            "ΔTOTAL": {"mean": -6.0, "std": 3.0, "uncertainty": "frame_sd"},
        },
        "mode": "single",
    }

    payload = processor.publication_component_payload({"Ligand": {"results": results}})

    assert list(payload["component"]) == ["VDWAALS", "EEL", "EGB", "ESURF", "GGAS", "GSOLV", "TOTAL"]
    assert list(payload["group"]) == ["raw"] * 4 + ["grouped"] * 2 + ["total"]
    assert payload.loc[payload.component == "TOTAL", "uncertainty"].item() == "saved_sample_sd"


def test_saved_energy_figure_uses_index_and_centered_mean():
    plotter = MMPBSAPlotter.__new__(MMPBSAPlotter)
    plotter._slots = {"Ligand": 0}
    payload = pd.DataFrame({"sample_index": [1, 2, 3, 4, 5], "total": [-1., -3., -5., -7., -9.]})
    payload["moving_mean"] = [-2., -3., -5., -7., -8.]
    payload["overall_mean"] = -5.0

    figure = plotter.create_publication_energy_figure(payload, "Ligand", window=3, x_label="Saved sample index")

    assert len(figure.data) == 3
    assert list(figure.data[0].x) == [1, 2, 3, 4, 5]
    assert list(figure.data[1].y) == [-2., -3., -5., -7., -8.]
    assert "Saved sample index" in figure.layout.xaxis.title.text
    assert figure.data[0].marker.color == categorical(0)


def test_saved_energy_payload_uses_only_validated_matching_time_mapping():
    frames = pd.DataFrame({"Frame #": [10, 20, 30], "TOTAL": [-2.0, -4.0, -6.0]})
    results = {"frame_data": frames}
    processor = MMPBSAProcessor.__new__(MMPBSAProcessor)

    index_payload = processor.publication_energy_payload(results, window=3)
    timed_payload = processor.publication_energy_payload(results, window=3, sampling_metadata={
        "validated": True, "frames": [10, 20, 30], "time": [0.0, 1.0, 2.0], "unit": "ns", "source": "frames.csv",
    })
    rejected_payload = processor.publication_energy_payload(results, window=3, sampling_metadata={
        "validated": True, "frames": [1, 2, 3], "time": [0.0, 1.0, 2.0], "unit": "ns",
    })

    assert index_payload["x_label"] == "Saved sample index"
    assert timed_payload["x_label"] == "Time (ns)"
    assert timed_payload["samples"]["moving_mean"].tolist() == [-3.0, -4.0, -5.0]
    assert rejected_payload["x_label"] == "Saved sample index"


def test_saved_energy_time_mapping_requires_unique_frames_and_provenance():
    processor = MMPBSAProcessor.__new__(MMPBSAProcessor)
    duplicate_frames = {"frame_data": pd.DataFrame({"Frame #": [10, 10, 30], "TOTAL": [-1., -2., -3.]})}
    duplicate_mapping = processor.publication_energy_payload(duplicate_frames, sampling_metadata={
        "validated": True, "frames": [10, 10, 30], "time": [0., 1., 2.], "unit": "ns", "source": "manifest.csv",
    })
    missing_source = processor.publication_energy_payload({
        "frame_data": pd.DataFrame({"Frame #": [10, 20, 30], "TOTAL": [-1., -2., -3.]})
    }, sampling_metadata={
        "validated": True, "frames": [10, 20, 30], "time": [0., 1., 2.], "unit": "ns",
    })

    assert duplicate_mapping["x_label"] == "Saved sample index"
    assert duplicate_mapping["physical_time"]["status"] == "unavailable"
    assert missing_source["physical_time"]["status"] == "unavailable"
    assert "source" in missing_source["physical_time"]["reason"]


def test_dat_results_publish_canonical_components_and_reported_dispersion(tmp_path):
    source = tmp_path / "FINAL_RESULTS_MMPBSA.dat"
    source.write_text("""Delta (Complex - Receptor - Ligand):
Energy Component       Average     SD(Prop.)         SD   SEM(Prop.)        SEM
-------------------------------------------------------------------------------------
ΔVDWAALS                -2.00          0.30       0.20         0.03       0.02
ΔGGAS                   -1.00          0.40       0.25         0.04       0.03
ΔGSOLV                  -3.00          0.50       0.40         0.05       0.04
ΔTOTAL                  -4.00          0.60       0.50         0.06       0.05
-------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------
""")

    parsed = MMPBSAParser().parse_results(source)

    assert parsed["canonical_components"]["TOTAL"]["mean"] == -4.0
    assert parsed["canonical_components"]["TOTAL"]["sample_sd"] == 0.5
    assert parsed["canonical_components"]["TOTAL"]["uncertainty"] == "reported_frame_sd"


def test_dat_derived_grouped_mean_does_not_claim_unreported_dispersion(tmp_path):
    source = tmp_path / "FINAL_RESULTS_MMPBSA.dat"
    source.write_text("""Delta (Complex - Receptor - Ligand):
Energy Component       Average     SD(Prop.)         SD   SEM(Prop.)        SEM
-------------------------------------------------------------------------------------
ΔVDWAALS                -2.00          0.30       0.20         0.03       0.02
ΔEEL                     1.00          0.30       0.20         0.03       0.02
ΔEGB                    -3.00          0.50       0.40         0.05       0.04
ΔESURF                  -1.00          0.50       0.40         0.05       0.04
-------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------
""")

    parsed = MMPBSAParser().parse_results(source)

    assert parsed["canonical_components"]["GGAS"]["mean"] == -1.0
    assert parsed["canonical_components"]["GGAS"]["sample_sd"] is None
    assert parsed["canonical_components"]["GGAS"]["uncertainty"] == "unavailable"


@pytest.mark.parametrize("window", [0, -1, 2, 2.5, True])
def test_saved_energy_payload_rejects_nonpositive_or_even_window(window):
    with pytest.raises(ValueError, match="positive odd"):
        MMPBSAProcessor.publication_energy_payload({"frame_data": pd.DataFrame({"TOTAL": [-1.]})}, window=window)


def test_residue_selection_keeps_most_favorable_first_and_supports_all():
    processor = MMPBSAProcessor.__new__(MMPBSAProcessor)
    source = pd.DataFrame({"Residue": ["A:1", "A:2", "A:3", "A:4"], "TOTAL Avg.": [1.0, -3.0, -1.0, -5.0]})

    selected = processor.publication_residue_payload(source, top_n=2)
    all_residues = processor.publication_residue_payload(source, top_n="all")

    assert selected["Residue"].tolist() == ["A:4", "A:2"]
    assert all_residues["Residue"].tolist() == ["A:4", "A:2", "A:3", "A:1"]


def test_residue_selection_defaults_to_fifteen_and_breaks_ties_stably():
    processor = MMPBSAProcessor.__new__(MMPBSAProcessor)
    residues = [f"R{i:02d}" for i in range(17, 0, -1)]
    energies = {f"R{i:02d}": -float(i) for i in range(1, 15)}
    energies.update({"R14": -14.0, "R15": -14.0, "R16": -15.0, "R17": 2.0})
    source = pd.DataFrame({"Residue": residues, "TOTAL Avg.": [energies[r] for r in residues]})

    default = processor.publication_residue_payload(source)
    custom = processor.publication_residue_payload(source, top_n=4)
    all_residues = processor.publication_residue_payload(source, top_n="all")

    assert len(default) == 15
    assert default["Residue"].tolist() == [
        "R16", "R14", "R15", "R13", "R12", "R11", "R10", "R09", "R08", "R07",
        "R06", "R05", "R04", "R03", "R02",
    ]
    assert custom["Residue"].tolist() == ["R16", "R14", "R15", "R13"]
    assert len(all_residues) == 17


def test_replicate_component_statistics_use_only_parsed_replicates(tmp_path):
    first = tmp_path / "rep1.csv"
    second = tmp_path / "rep2.csv"
    pd.DataFrame({
        "Frame #": [1, 2], "VDWAALS": [-2.0, -4.0], "EEL": [1.0, 3.0],
        "TOTAL": [-1.0, -3.0],
    }).to_csv(first, index=False)
    pd.DataFrame({
        "Frame #": [1, 2], "VDWAALS": [-6.0, -8.0], "TOTAL": [-5.0, -7.0],
    }).to_csv(second, index=False)
    config = MMPBSAConfig(
        base_dir=tmp_path, output_dir=tmp_path / "out", protein_name="P", ligand_name="L",
        systems=[MMPBSASystemConfig(name="Sys", dir_pattern=".", replicates=3,
                                    results_file_pattern="rep*.csv")], file_format="csv",
    )
    processor = MMPBSAProcessor(config)
    aggregated = processor.load_all_data()["Sys"]["results"]
    payload = processor.publication_component_payload({"System": {"results": aggregated}})
    vdw = payload.loc[payload.component == "VDWAALS"].iloc[0]
    eel = payload.loc[payload.component == "EEL"].iloc[0]

    assert aggregated["n_replicates"] == 2  # two files found from a requested three
    assert (vdw["mean"], vdw["sd"], vdw["n"]) == (-5.0, 2.8284271247461903, 2)
    assert eel["n"] == 1
    assert pd.isna(eel["sd"])
    assert eel["uncertainty"] == "unavailable_single_replicate"
    assert aggregated["replicate_sources"] == {"1": str(first), "2": str(second)}


def test_replicate_saved_energy_means_do_not_cross_trace_boundaries():
    processor = MMPBSAProcessor.__new__(MMPBSAProcessor)
    frame_data = pd.DataFrame({
        "Replicate": [1, 1, 1, 2, 2, 2],
        "Frame #": [1, 2, 3, 1, 2, 3],
        "TOTAL": [-1., -3., -5., -10., -20., -30.],
    })

    payload = processor.publication_energy_payload({"frame_data": frame_data}, window=3)
    samples = payload["samples"]

    assert samples.groupby("replicate")["moving_mean"].apply(list).to_dict() == {
        1: [-2.0, -3.0, -4.0], 2: [-15.0, -20.0, -25.0],
    }
    assert samples.groupby("replicate")["overall_mean"].first().to_dict() == {1: -3.0, 2: -20.0}
    plotter = MMPBSAPlotter.__new__(MMPBSAPlotter)
    plotter._slots = {"Sys": 0}
    figure = plotter.create_publication_energy_figure(payload["samples"], "Sys", window=3)
    assert len(figure.data) == 6
    assert "window 3" in figure.data[1].name


def test_component_figure_explains_groups_uncertainty_and_missing_terms():
    payload = pd.DataFrame([
        {"system": "Sys", "component": "VDWAALS", "group": "raw", "mean": -2., "sd": .2,
         "n": 2, "uncertainty": "replicate_sd", "status": "available", "reason": None},
        {"system": "Sys", "component": "EEL", "group": "raw", "mean": None, "sd": None,
         "n": 0, "uncertainty": None, "status": "unavailable", "reason": "EEL missing"},
        {"system": "Sys", "component": "GGAS", "group": "grouped", "mean": -1., "sd": .1,
         "n": 2, "uncertainty": "replicate_sd", "status": "available", "reason": None},
        {"system": "Sys", "component": "TOTAL", "group": "total", "mean": -3., "sd": .3,
         "n": 2, "uncertainty": "replicate_sd", "status": "available", "reason": None},
    ])

    figure = MMPBSAPlotter.create_publication_components_figure(payload)

    assert {trace.name for trace in figure.data if trace.showlegend} >= {
        "Raw energetic terms", "Grouped sums", "TOTAL endpoint estimate",
    }
    annotation_text = " ".join(annotation.text for annotation in figure.layout.annotations)
    assert "replicate SD" in annotation_text
    assert "EEL" in annotation_text
    assert figure.layout.showlegend is True


def test_nonfinite_saved_total_is_unavailable_with_count():
    processor = MMPBSAProcessor.__new__(MMPBSAProcessor)

    payload = processor.publication_energy_payload({
        "frame_data": pd.DataFrame({"Frame #": [1, 2, 3], "TOTAL": [-1.0, np.inf, -3.0]})
    })

    assert payload["status"] == "unavailable"
    assert payload["invalid_sample_count"] == 1
    assert "finite" in payload["reason"]


def test_start_frame_recomputes_canonical_values_from_kept_saved_rows(tmp_path):
    input_dir = tmp_path / "rep"
    input_dir.mkdir()
    pd.DataFrame({"Frame #": [1, 2, 3], "VDWAALS": [100., 4., 6.], "TOTAL": [100., 1., 3.]}).to_csv(
        input_dir / "results.csv", index=False
    )
    config = MMPBSAConfig(
        base_dir=tmp_path, output_dir=tmp_path / "out", protein_name="P", ligand_name="L",
        systems=[MMPBSASystemConfig(name="Sys", dir_pattern="rep", results_file_pattern="results.csv")],
        file_format="csv", start_frame=2,
    )
    filtered = MMPBSAProcessor(config).load_all_data()["Sys"]["results"]

    assert filtered["canonical_components"]["VDWAALS"]["mean"] == 5.0
    assert filtered["canonical_components"]["TOTAL"]["mean"] == 2.0
    assert filtered["canonical_components"]["VDWAALS"]["n"] == 2


def test_nonfinite_component_summary_is_explicitly_unavailable():
    processor = MMPBSAProcessor.__new__(MMPBSAProcessor)
    payload = processor.publication_component_payload({"Sys": {"results": {
        "mode": "single", "canonical_components": {
            "VDWAALS": {"mean": -1.0, "sample_sd": 0.1, "n": 2},
            "EEL": {"mean": np.inf, "sample_sd": 0.2, "n": 2},
        },
    }}})
    eel = payload.loc[payload.component == "EEL"].iloc[0]

    assert eel["status"] == "unavailable"
    assert pd.isna(eel["mean"])
    assert pd.isna(eel["sd"])
    assert "finite" in eel["reason"]


def test_sanitized_system_labels_get_distinct_energy_artifact_names(tmp_path):
    names = ["A B", "A_B"]
    for index, name in enumerate(names, start=1):
        input_dir = tmp_path / f"rep{index}"
        input_dir.mkdir()
        pd.DataFrame({"Frame #": [1, 2], "TOTAL": [-1.0, -2.0]}).to_csv(
            input_dir / "results.csv", index=False
        )
    config = MMPBSAConfig(
        base_dir=tmp_path, output_dir=tmp_path / "out", protein_name="P", ligand_name="L",
        systems=[MMPBSASystemConfig(name=name, dir_pattern=f"rep{index}", results_file_pattern="results.csv")
                 for index, name in enumerate(names, start=1)],
        file_format="csv", publication=PublicationConfig(formats=["html"]),
    )
    config.plot_config.save_formats = ["html"]
    analyzer = MMPBSAAnalyzer(config)
    result = analyzer.run_analysis()

    pub_dir = tmp_path / "out" / "publication"
    assert (pub_dir / "saved_total_1_A_B.html").is_file()
    assert (pub_dir / "saved_total_2_A_B.html").is_file()
    assert (pub_dir / "saved_total_1_A_B_source.csv").is_file()
    assert (pub_dir / "saved_total_2_A_B_source.csv").is_file()
    assert set(key for key in result["publication_outputs"]["availability"] if key.startswith("saved_total:")) == {
        "saved_total:A B", "saved_total:A_B",
    }


def test_residue_figure_uses_horizontal_bars_with_adaptive_height():
    plotter = MMPBSAPlotter.__new__(MMPBSAPlotter)
    plotter._slots = {"Sys": 0}
    residues = pd.DataFrame({"Residue": [f"A:{i}" for i in range(1, 26)], "TOTAL Avg.": [-float(i) for i in range(1, 26)]})

    figure = plotter.create_publication_residue_figure({"Sys": residues})

    assert figure.data[0].orientation == "h"
    assert figure.layout.height >= 650
    assert "Residue" in figure.layout.yaxis.title.text


def test_analyzer_folder_run_preserves_legacy_plot_and_writes_publication_outputs(tmp_path):
    input_dir = tmp_path / "rep"
    input_dir.mkdir()
    pd.DataFrame({
        "Frame #": [8, 3, 5], "VDWAALS": [-2.0, -3.0, -4.0], "EEL": [1.0, 1.0, 2.0],
        "EGB": [-1.0, -2.0, -1.0], "ESURF": [-0.5, -0.5, -0.5],
        "GGAS": [-1.0, -2.0, -2.0], "GSOLV": [-1.5, -2.5, -1.5],
        "TOTAL": [-2.5, -4.5, -3.5],
    }).to_csv(input_dir / "results.csv", index=False)
    pd.DataFrame({
        "Frame #": [1, 2, 1, 2], "Residue": ["A:1", "A:1", "A:2", "A:2"],
        "TOTAL": [-2.0, -1.0, 0.5, 1.0],
    }).to_csv(input_dir / "decomp.csv", index=False)
    config = MMPBSAConfig(
        base_dir=tmp_path, output_dir=tmp_path / "out", protein_name="P", ligand_name="L",
        systems=[MMPBSASystemConfig(name="Ligand A", dir_pattern="rep",
                                    results_file_pattern="results.csv", decomp_file_pattern="decomp.csv")],
        file_format="csv",
        publication=PublicationConfig(formats=["html"]),
    )
    config.plot_config.save_formats = ["html"]
    analyzer = MMPBSAAnalyzer(config)
    result = analyzer.run_analysis()

    pub_dir = tmp_path / "out" / "publication"
    assert result["success"] is True
    assert (tmp_path / "out" / "plots" / "components_Ligand A.html").is_file()
    assert (pub_dir / "components.html").is_file()
    assert (pub_dir / "saved_total_1_Ligand_A.html").is_file()
    assert (pub_dir / "residues.html").is_file()
    assert pd.read_csv(pub_dir / "components_source.csv")["component"].tolist() == [
        "VDWAALS", "EEL", "EGB", "ESURF", "GGAS", "GSOLV", "TOTAL",
    ]
    assert pd.read_csv(pub_dir / "saved_total_1_Ligand_A_source.csv")["sample_index"].tolist() == [1, 2, 3]
    assert Path(result["publication_outputs"]["metadata"]).is_file()


def test_unusable_optional_decomposition_records_reason_and_keeps_other_figures(tmp_path):
    input_dir = tmp_path / "rep"
    input_dir.mkdir()
    pd.DataFrame({"Frame #": [1, 2, 3], "VDWAALS": [-1.0, -1.0, -1.0], "TOTAL": [-1.0, -2.0, -3.0]}).to_csv(
        input_dir / "results.csv", index=False
    )
    pd.DataFrame({"Residue": ["A:1"], "EEL Avg.": [0.5]}).to_csv(input_dir / "decomp.csv", index=False)
    config = MMPBSAConfig(
        base_dir=tmp_path, output_dir=tmp_path / "out", protein_name="P", ligand_name="L",
        systems=[MMPBSASystemConfig(name="Sys", dir_pattern="rep", results_file_pattern="results.csv",
                                    decomp_file_pattern="decomp.csv")], file_format="csv",
        publication=PublicationConfig(formats=["html"]),
    )
    analyzer = MMPBSAAnalyzer(config)
    config.plot_config.save_formats = ["html"]
    result = analyzer.run_analysis()

    pub_dir = tmp_path / "out" / "publication"
    assert (pub_dir / "components.html").is_file()
    assert (pub_dir / "saved_total_1_Sys.html").is_file()
    assert not (pub_dir / "residues.html").exists()
    assert result["publication_outputs"]["availability"]["residues:Sys"]["status"] == "unavailable"
    assert result["publication_outputs"]["availability"]["components"]["missing_components"]["Sys"]

