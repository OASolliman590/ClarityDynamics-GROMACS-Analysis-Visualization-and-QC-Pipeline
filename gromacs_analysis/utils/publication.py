"""Additive publication configuration, numerical sources, and export records."""
from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path
from typing import Union

from .figure_export import write_static


@dataclass
class PublicationConfig:
    """Shared settings for PCA and endpoint-energy publication reports."""
    enabled: bool = True
    subdirectory: str = 'publication'
    formats: list = field(default_factory=lambda: ['html', 'svg', 'png', 'pdf'])
    moving_window: int = 11
    top_n: Union[int, str] = 15
    scree_components: int = 10
    pca_manifest: str = 'publication_inputs.json'
    sampling_metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        for name in ('moving_window', 'scree_components'):
            value = getattr(self, name)
            if type(value) is not int or value < 1:
                raise ValueError(f'{name} must be a positive integer')
        if self.moving_window % 2 == 0:
            raise ValueError('moving_window must be a positive odd integer')
        if self.top_n != 'all' and (type(self.top_n) is not int or self.top_n < 1):
            raise ValueError("top_n must be a positive integer or 'all'")
        directory = Path(self.subdirectory)
        if directory.is_absolute() or '..' in directory.parts or directory == Path('.'):
            raise ValueError('subdirectory must be a separate relative publication directory')
        if not self.formats or any(f not in {'html', 'svg', 'png', 'pdf'} for f in self.formats):
            raise ValueError('publication formats must contain html, svg, png, or pdf')


def export_publication_figure(fig, output_dir: Path, stem: str, config: PublicationConfig) -> dict:
    """Attempt every requested format and expose each actual export outcome."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    records = {}
    for fmt in config.formats:
        path = output_dir / f'{stem}.{fmt}'
        try:
            if fmt == 'html':
                fig.write_html(str(path), include_plotlyjs=True)
            else:
                write_static(fig, path, fmt=fmt, width=1800, height=650, scale=2)
            records[fmt] = {'status': 'available', 'path': str(path)}
        except Exception as exc:
            records[fmt] = {'status': 'unavailable', 'reason': str(exc)}
    return records


def write_source_table(df, output_dir: Path, stem: str) -> str:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f'{stem}.csv'
    df.to_csv(path, index=False)
    return str(path)


def write_metadata(payload: dict, output_dir: Path, stem: str) -> str:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f'{stem}.json'
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')
    return str(path)
