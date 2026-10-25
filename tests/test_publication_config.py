"""Publication configuration is shared by direct and unified stage entry points."""
from dataclasses import asdict
import pytest
import yaml
from gromacs_analysis.pca.yaml_loader import load_yaml_config as load_pca
from gromacs_analysis.mmpbsa.yaml_loader import load_yaml_config as load_mmpbsa
from gromacs_analysis.pipeline.yaml_loader import load_pipeline_config, save_pipeline_config


@pytest.mark.parametrize('stage,loader', [('pca', load_pca), ('mmpbsa', load_mmpbsa)])
def test_publication_yaml_parity_and_roundtrip(tmp_path, stage, loader):
    data = {'base_dir': str(tmp_path), 'output_dir': str(tmp_path / 'out'),
            'protein_name': 'Protein', 'ligand_name': 'Bound',
            'publication': {'top_n': 'all', 'moving_window': 5, 'formats': ['html'],
                            'scree_components': 4, 'subdirectory': 'publication/custom',
                            'sampling_metadata': {'Bound': {'validated': True, 'frames': [2, 4],
                                                          'time': [0.2, 0.4], 'unit': 'ns', 'source': 'sampling log'}}}}
    if stage == 'mmpbsa':
        data['systems'] = [{'name': 'Bound', 'dir_pattern': 'bound'}]
    direct = tmp_path / 'direct.yaml'
    direct.write_text(yaml.safe_dump(data))
    pipeline = tmp_path / 'pipeline.yaml'
    pipeline.write_text(yaml.safe_dump({stage: data}))
    expected = asdict(loader(direct).publication)
    config = load_pipeline_config(pipeline)
    assert asdict(config.stages[stage].config.publication) == expected
    saved = tmp_path / 'saved.yaml'
    save_pipeline_config(config, saved)
    assert asdict(load_pipeline_config(saved).stages[stage].config.publication) == expected


@pytest.mark.parametrize('field,value', [('moving_window', 0), ('moving_window', 2),
    ('moving_window', True), ('top_n', 0), ('top_n', 1.2), ('scree_components', 0),
    ('subdirectory', '../escape'), ('subdirectory', '.'), ('formats', ['exe'])])
def test_invalid_publication_settings_are_rejected(tmp_path, field, value):
    from gromacs_analysis.utils.publication import PublicationConfig
    with pytest.raises(ValueError):
        PublicationConfig(**{field: value})


def test_batch_pca_inherits_publication_and_allows_system_override(tmp_path):
    pipeline = tmp_path / 'batch.yaml'
    pipeline.write_text(yaml.safe_dump({'pca_batch': {
        'base_dir': str(tmp_path), 'publication': {'formats': ['html'], 'scree_components': 5},
        'systems': [{'name': 'First', 'dir_name': 'first'},
                    {'name': 'Second', 'dir_name': 'second', 'publication': {'scree_components': 2}}]}}))
    configs = load_pipeline_config(pipeline).stages['pca_batch'].config
    assert [c.publication.scree_components for c in configs] == [5, 2]
    assert [c.publication.formats for c in configs] == [['html'], ['html']]


def test_joint_publication_settings_survive_pipeline_save(tmp_path):
    pipeline = tmp_path / 'joint.yaml'
    pipeline.write_text(yaml.safe_dump({'pca_joint': {
        'base_dir': str(tmp_path), 'output_dir': str(tmp_path / 'out'),
        'publication': {'formats': ['html'], 'scree_components': 5},
        'systems': [{'name': 'Bound', 'topology': 'top.pdb', 'trajectories': ['traj.xtc']}]}}))
    original = load_pipeline_config(pipeline)
    saved = tmp_path / 'saved_joint.yaml'
    save_pipeline_config(original, saved)
    restored = load_pipeline_config(saved)
    assert restored.stages['pca_joint'].config is not None
    assert asdict(restored.stages['pca_joint'].config.publication) == asdict(original.stages['pca_joint'].config.publication)
