"""Exercise format reporting at the exporter boundary and a real static image."""
from pathlib import Path
import plotly.graph_objects as go
import pytest
from gromacs_analysis.utils.figure_export import reset_static_export_state
from gromacs_analysis.utils.publication import PublicationConfig, export_publication_figure


def test_export_records_html_even_when_static_backend_is_unavailable(tmp_path, monkeypatch):
    from gromacs_analysis.utils.figure_export import StaticExportUnavailable
    def unavailable(*args, **kwargs):
        raise StaticExportUnavailable('Test backend missing')
    monkeypatch.setattr('gromacs_analysis.utils.publication.write_static', unavailable)
    result = export_publication_figure(go.Figure(go.Bar(x=['a'], y=[2])), tmp_path, 'example', PublicationConfig())
    assert set(result) == {'html', 'svg', 'png', 'pdf'}
    assert Path(result['html']['path']).is_file()
    for fmt in ('svg', 'png', 'pdf'):
        assert result[fmt] == {'status': 'unavailable', 'reason': 'Test backend missing'}


def test_real_static_vector_and_raster_exports(tmp_path):
    reset_static_export_state()
    result = export_publication_figure(go.Figure(go.Bar(x=['a'], y=[2])), tmp_path, 'example', PublicationConfig())
    assert result['html']['status'] == 'available'
    if result['svg']['status'] == 'unavailable':
        pytest.skip(result['svg']['reason'])
    assert '<svg' in Path(result['svg']['path']).read_text()
    assert Path(result['png']['path']).read_bytes().startswith(b'\x89PNG')
    assert Path(result['pdf']['path']).read_bytes().startswith(b'%PDF')
