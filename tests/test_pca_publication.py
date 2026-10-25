"""Publication behavior through synthetic folders and returned Plotly figures."""
import json
import numpy as np
import pandas as pd
import pytest
from gromacs_analysis.pca.config import PCAConfig
from gromacs_analysis.pca.analyzer import PCAAnalyzer
from gromacs_analysis.pca.processor import PCAProcessor
from gromacs_analysis.pca.plotter import PCAPlotter
from gromacs_analysis.md.config import PlotConfig


def folder(tmp_path, bases=('shared', 'shared'), diagnostics=False):
    (tmp_path / 'eigenvals.xvg').write_text('1 5\n2 3\n3 1\n4 1\n')
    states = []
    for name, basis in zip(('Apo', 'Bound'), bases):
        pd.DataFrame({'PC1': [1, 2, 3], 'PC2': [4, 5, 6], 'PC3': [7, 8, 9]}).to_csv(tmp_path / f'{name}.csv', index=False)
        states.append({'name': name, 'projections': f'{name}.csv', 'eigenvalues': 'eigenvals.xvg',
                       'basis_id': basis, 'basis_verified': True, 'basis_source': 'joint fit log', 'unit': 'nm'})
    manifest = {'states': states}
    if diagnostics:
        pd.DataFrame([
            {'metric': 'rmsip', 'value': 0.8, 'state_pair': 'Apo|Bound', 'compatible': True, 'source': 'overlap log'},
            {'metric': 'window_sensitivity', 'value': 0.7, 'state_pair': 'Apo|Bound', 'window': '0-10 ns', 'compatible': True, 'source': 'window log'},
            {'metric': 'split_half_rmsip', 'value': 0.9, 'state_pair': 'Apo', 'compatible': True, 'source': 'split log'}
        ]).to_csv(tmp_path / 'diagnostics.csv', index=False)
        manifest['diagnostics'] = 'diagnostics.csv'
        (tmp_path / 'cc_all.xvg').write_text('1 0.2\n2 0.1\n3 0.05\n')
    (tmp_path / 'publication_inputs.json').write_text(json.dumps(manifest))
    return PCAConfig(base_dir=tmp_path, output_dir=tmp_path / 'out', ligand_name='Bound',
                     plot_config=PlotConfig(save_formats=['html']))


def payload(config):
    processor = PCAProcessor(config)
    detected = processor.auto_detect_files()
    return processor.build_publication_payload(processor.load_all_data(detected), detected)


def test_verified_shared_basis_composite_uses_full_variance(tmp_path):
    config = folder(tmp_path)
    result = payload(config)
    fig = PCAPlotter(config).plot_publication_composite(result)
    projections = [t for t in fig.data if t.type in ('scatter', 'scatter3d') and t.name in ('Apo', 'Bound')]
    assert [(t.type, t.name) for t in projections] == [('scatter', 'Apo'), ('scatter3d', 'Apo'), ('scatter', 'Bound'), ('scatter3d', 'Bound')]
    assert projections[0].marker.color == projections[1].marker.color
    assert '50.0%' in fig.layout.xaxis.title.text
    assert '10.0%' in fig.layout.scene.zaxis.title.text
    assert list(next(t for t in fig.data if t.name == 'Cumulative variance').y) == [50, 80, 90, 100]
    assert any('PC1–PC3: 90.0%' in a.text for a in fig.layout.annotations)
    assert any('sampled coordinate covariance' in a.text for a in fig.layout.annotations)


def test_mismatched_basis_blocks_shared_overlay_but_keeps_summary(tmp_path):
    config = folder(tmp_path, bases=('first', 'other'))
    config.publication.formats = ['html']
    result = PCAAnalyzer(config).run_analysis()
    pub = result['publication']
    assert pub['availability']['composite']['status'] == 'unavailable'
    assert 'basis' in pub['availability']['composite']['reason']
    assert not (config.output_dir / 'publication' / 'pca_composite.html').exists()
    summary = pd.read_csv(pub['diagnostics_csv'])
    assert 'pc1_fraction' in summary.metric.values
    assert 'unavailable' in summary[summary.metric == 'shared_overlay'].status.values
    assert (config.output_dir / 'PCA_analysis_summary.md').exists()


def test_analyzer_exports_auditable_sources_and_diagnostic_availability(tmp_path):
    config = folder(tmp_path, diagnostics=True)
    config.publication.formats = ['html']
    result = PCAAnalyzer(config).run_analysis()['publication']
    summary = pd.read_csv(result['diagnostics_csv'])
    for name in ('rmsip', 'window_sensitivity', 'split_half_rmsip', 'cosine_content'):
        assert set(summary[summary.metric == name].status) == {'available'}
    metadata = json.loads(open(result['metadata_json']).read())
    assert metadata['basis_id'] == 'shared'
    assert metadata['state_labels'] == ['Apo', 'Bound']
    assert metadata['exports']['html']['status'] == 'available'
    sources = pd.read_csv(result['projections_csv'])
    assert list(sources.state.unique()) == ['Apo', 'Bound']
    assert list(sources.PC3) == [7, 8, 9, 7, 8, 9]


def test_partial_diagnostics_never_imply_zero_and_single_state_is_valid(tmp_path):
    config = folder(tmp_path, bases=('unverified',))
    path = tmp_path / 'publication_inputs.json'
    manifest = json.loads(path.read_text())
    manifest['states'][0]['basis_verified'] = False
    path.write_text(json.dumps(manifest))
    result = payload(config)
    assert result['availability']['composite']['status'] == 'available'
    records = result['diagnostics']
    for metric in ('rmsip', 'window_sensitivity', 'split_half_rmsip', 'cosine_content'):
        row = next(r for r in records if r['metric'] == metric)
        assert row['status'] == 'unavailable' and row['value'] is None and row['reason']


@pytest.mark.parametrize('n_components', [1, 3])
def test_native_joint_result_retains_shared_basis_and_full_denominator(tmp_path, n_components):
    from gromacs_analysis.pca.config import PCAJointConfig, PCAJointSystem
    from gromacs_analysis.pca.joint_analyzer import PCAJointAnalyzer
    from gromacs_analysis.utils.publication import PublicationConfig
    rng = np.random.default_rng(18)
    coords = {name: rng.normal(size=(20, 5, 3)) for name in ('Apo', 'Bound')}
    config = PCAJointConfig(output_dir=tmp_path, protein_name='Protein', n_components=n_components,
        systems=[PCAJointSystem(name=n, topology='unused.pdb', trajectories=['unused.xtc']) for n in coords],
        plot_config=PlotConfig(save_formats=['html']), publication=PublicationConfig(formats=['html']))
    result = PCAJointAnalyzer(config).run_analysis(coords_by_system=coords)
    pub = result['publication']
    metadata = json.loads(open(pub['metadata_json']).read())
    assert metadata['basis_id']
    if n_components == 1:
        assert metadata['availability']['composite']['status'] == 'unavailable'
        return
    assert metadata['availability']['composite']['status'] == 'available'
    fractions = pd.read_csv(pub['variance_csv'])
    # The three displayed eigenvalues do not exhaust the 3N-dimensional covariance.
    assert fractions.groupby('state').fraction.sum().max() < 1.0


def test_missing_optional_diagnostic_file_preserves_projection_composite(tmp_path):
    config = folder(tmp_path)
    path = tmp_path / 'publication_inputs.json'
    manifest = json.loads(path.read_text())
    manifest['diagnostics'] = 'absent.csv'
    path.write_text(json.dumps(manifest))
    config.publication.formats = ['html']
    result = PCAAnalyzer(config).run_analysis()['publication']
    assert result['availability']['composite']['status'] == 'available'
    summary = pd.read_csv(result['diagnostics_csv'])
    assert (summary.loc[summary.metric == 'rmsip', 'status'] == 'unavailable').all()


def test_unverified_matching_basis_identifiers_cannot_authorize_overlay(tmp_path):
    config = folder(tmp_path)
    path = tmp_path / 'publication_inputs.json'
    manifest = json.loads(path.read_text())
    manifest['states'][1]['basis_verified'] = False
    path.write_text(json.dumps(manifest))
    result = payload(config)
    assert result['availability']['composite']['status'] == 'unavailable'
    with pytest.raises(ValueError, match='verified basis'):
        PCAPlotter(config).plot_publication_composite(result)


def test_disabled_publication_preserves_legacy_analyzer_outputs(tmp_path):
    config = folder(tmp_path)
    config.publication.enabled = False
    result = PCAAnalyzer(config).run_analysis()
    assert result['success'] and result['publication'] is None
    assert (config.output_dir / 'PCA_analysis_summary.md').exists()
    assert not (config.output_dir / 'publication').exists()
