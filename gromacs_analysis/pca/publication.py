"""Numerical payloads and rendering for additive PCA publication figures.

Imported shared bases are attested in publication_inputs.json. Legacy file names
never establish basis equivalence. The native joint analyzer supplies its own
verified result, including the full covariance denominator.
"""
from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from ..utils.palette import APO_HEX, categorical
from ..utils.publication import export_publication_figure, write_metadata, write_source_table

DIAGNOSTICS = ('rmsip', 'window_sensitivity', 'split_half_rmsip', 'cosine_content')
CAPTION = ('Projections describe sampled coordinate covariance; they do not establish '
           'equilibrium populations or a free-energy landscape.')


def diagnostic(metric, value=None, unit='fraction', reason='', state_pair='', window='', source=''):
    return dict(metric=metric, value=value, unit=unit,
                status='available' if value is not None else 'unavailable', reason=reason,
                state_pair=state_pair, window=window, source=source)


def prepare_payload(states, diagnostics=(), scree_components=10):
    """Validate scientific availability and account for the full eigenvalue sum."""
    records = list(diagnostics)
    prepared = []
    color_slot = 0
    for index, state in enumerate(states):
        state = dict(state)
        eigenvalues = np.asarray(state.get('eigenvalues', []), dtype=float)
        total = state.get('total_variance', float(eigenvalues.sum()))
        valid = (len(eigenvalues) > 0 and np.isfinite(eigenvalues).all()
                 and (eigenvalues >= 0).all() and np.isfinite(total) and total > 0
                 and total + 1e-10 >= eigenvalues.sum())
        fractions = eigenvalues / total if valid else np.array([])
        state['eigenvalues'] = [float(v) if np.isfinite(v) else None for v in eigenvalues]
        state['fractions'] = fractions.tolist()
        state['cumulative'] = np.cumsum(fractions).tolist()
        state['scree_components'] = scree_components
        is_apo = state['name'].lower() == 'apo'
        state['color'] = state.get('color') or (APO_HEX if is_apo else categorical(color_slot))
        state['symbol'] = ['circle', 'diamond', 'square', 'x'][(color_slot // 8) % 4]
        if not is_apo:
            color_slot += 1
        state.setdefault('unit', 'source coordinate units')
        state.setdefault('basis_id', None)
        state.setdefault('basis_verified', False)
        for pc in range(3):
            records.append(diagnostic(f'pc{pc + 1}_fraction',
                float(fractions[pc]) if pc < len(fractions) else None,
                reason='' if pc < len(fractions) else 'Valid eigenvalue spectrum unavailable',
                state_pair=state['name'], source=state.get('eigenvalue_source', '')))
        records.append(diagnostic('pc1_pc3_cumulative',
            float(sum(fractions[:3])) if len(fractions) >= 3 else None,
            reason='' if len(fractions) >= 3 else 'Three component fractions unavailable',
            state_pair=state['name'], source=state.get('eigenvalue_source', '')))
        prepared.append(state)
    reason = ''
    if not prepared:
        reason = 'No PCA projection inputs available'
    elif len(prepared) > 1:
        basis_ids = {s['basis_id'] for s in prepared}
        if (len(basis_ids) != 1 or not next(iter(basis_ids), None)
                or not all(s['basis_verified'] is True and s.get('basis_source') for s in prepared)):
            reason = 'Shared eigenbasis unavailable: all states require the same verified basis identifier and provenance'
        elif any(s['fractions'] != prepared[0]['fractions'] or s['unit'] != prepared[0]['unit'] for s in prepared[1:]):
            reason = 'Shared basis has inconsistent eigenvalue fractions or coordinate units'
    if not reason and any(not s['projections'].get('PC1') or not s['projections'].get('PC2') for s in prepared):
        reason = 'PC1 and PC2 projection samples are required for the composite'
    available = not reason
    records.append(diagnostic('shared_overlay', 1.0 if available else None, unit='boolean', reason=reason,
                              state_pair='|'.join(s['name'] for s in prepared)))
    for metric in DIAGNOSTICS:
        if not any(r['metric'] == metric for r in records):
            records.append(diagnostic(metric, reason='Dedicated compatible input not supplied'))
    return {'states': prepared, 'diagnostics': records,
            'availability': {'composite': {'status': 'available' if available else 'unavailable', 'reason': reason}}}


def imported_payload(processor, all_data, detected):
    """Use an explicit manifest when present, otherwise retain single-state inputs."""
    config = processor.config
    manifest_path = Path(config.base_dir) / config.publication.pca_manifest
    if manifest_path.exists():
        parsed = processor.parser.parse_publication_manifest(manifest_path)
        states, records = parsed['states'], parsed['diagnostics']
    else:
        states, records, seen = [], [], set()
        for role, name in [('holo', detected.get('ligand_name') or 'Bound'), ('apo', 'Apo')]:
            source = detected.get(f'proj_3d_{role}') or detected.get(f'proj_2d_{role}', {}).get('12')
            if source is None or str(source) in seen:
                continue
            seen.add(str(source))
            points = all_data.get(f'proj_3d_{role}')
            if points:
                projection = {f'PC{i+1}': list(col) for i, col in enumerate(points)}
            else:
                xy = all_data.get(f'proj_2d_{role}', {}).get('12', {})
                projection = {'PC1': xy.get('x', []), 'PC2': xy.get('y', [])}
            # The legacy spectrum in base_dir belongs to that directory's bound state.
            spectrum = all_data.get('eigenvals', {}).get('y', []) if role == 'holo' else []
            states.append(dict(name=name, projections=projection, eigenvalues=spectrum, color=config.color if role == 'holo' else APO_HEX,
                               source=str(source), eigenvalue_source=str(detected.get('eigenvals') or '')))
    for role, name in [('holo', detected.get('ligand_name') or 'Bound'), ('apo', 'Apo')]:
        cosine = all_data.get(f'cc_{role}', {})
        for pc, value in zip(cosine.get('x', []), cosine.get('y', [])):
            valid = np.isfinite(value) and 0 <= value <= 1
            records.append(diagnostic('cosine_content', float(value) if valid else None,
                reason='' if valid else 'Cosine content must be finite and between zero and one',
                state_pair=name, window=f'PC{int(pc)}', source=str(detected.get(f'cc_{role}') or '')))
    return prepare_payload(states, records, config.publication.scree_components)


def composite_figure(payload):
    """Render one verified basis in 2D, 3D, and a scree/cumulative panel."""
    if payload['availability']['composite']['status'] != 'available':
        raise ValueError(payload['availability']['composite']['reason'])
    states = payload['states']
    fig = make_subplots(rows=1, cols=3, specs=[[{}, {'type': 'scene'}, {}]],
                        subplot_titles=('PC1–PC2', 'PC1–PC2–PC3', 'Explained variance'),
                        horizontal_spacing=0.08, column_widths=[0.32, 0.36, 0.32])
    for state in states:
        projection = state['projections']
        color = state['color']
        if projection.get('PC1') and projection.get('PC2'):
            fig.add_trace(go.Scatter(x=projection['PC1'], y=projection['PC2'], mode='markers',
                marker=dict(color=color, symbol=state['symbol'], size=5, opacity=0.55), name=state['name'], legendgroup=state['name']), row=1, col=1)
        if projection.get('PC3'):
            fig.add_trace(go.Scatter3d(x=projection['PC1'], y=projection['PC2'], z=projection['PC3'], mode='markers',
                marker=dict(color=color, symbol=state['symbol'], size=3, opacity=0.55), name=state['name'], legendgroup=state['name'], showlegend=False), row=1, col=2)
    first = states[0]
    fractions = first['fractions']
    def label(pc):
        variance = f'{100 * fractions[pc - 1]:.1f}%' if pc <= len(fractions) else 'variance unavailable'
        return f'PC{pc} ({variance}; {first["unit"]})'
    fig.update_xaxes(title_text=label(1), row=1, col=1)
    fig.update_yaxes(title_text=label(2), row=1, col=1)
    fig.update_scenes(xaxis=dict(title=dict(text=label(1).split(';')[0] + ')', font=dict(size=12)), tickfont=dict(size=11)),
                      yaxis=dict(title=dict(text=label(2).split(';')[0] + ')', font=dict(size=12)), tickfont=dict(size=11)),
                      zaxis=dict(title=dict(text=label(3).split(';')[0] + ')', font=dict(size=12)), tickfont=dict(size=11)))
    fig.update_scenes(camera=dict(eye=dict(x=1.7, y=1.7, z=1.3)))
    k = min(len(fractions), first['scree_components'])
    x = list(range(1, k + 1))
    fig.add_trace(go.Bar(x=x, y=[100*f for f in fractions[:k]], name='Individual variance', marker_color=categorical(0)), row=1, col=3)
    fig.add_trace(go.Scatter(x=x, y=[100*f for f in first['cumulative'][:k]], name='Cumulative variance', mode='lines+markers',
                            line=dict(color=categorical(1))), row=1, col=3)
    fig.update_xaxes(title_text='Principal component', dtick=1, row=1, col=3)
    fig.update_yaxes(title_text='Variance (%)', range=[0, 100], row=1, col=3)
    cumulative = f'PC1–PC3: {100 * sum(fractions[:3]):.1f}%' if len(fractions) >= 3 else 'PC1–PC3 fraction unavailable'
    fig.add_annotation(text=cumulative, x=0.5, y=1.15, xref='paper', yref='paper', showarrow=False)
    if not fractions:
        fig.add_annotation(text='Scree unavailable: valid full spectrum required', x=0.9, y=0.5, xref='paper', yref='paper', showarrow=False)
    if not all(s['projections'].get('PC3') for s in states):
        fig.add_annotation(text='3D unavailable for states without PC3', x=0.5, y=0.5, xref='paper', yref='paper', showarrow=False)
    fig.add_annotation(text=CAPTION, x=0.5, y=-0.25, xref='paper', yref='paper', showarrow=False, font=dict(size=13))
    fig.update_layout(template='plotly_white', width=1800, height=650, font=dict(family='Arial', size=16),
                      margin=dict(t=110, b=160, l=70, r=70), legend=dict(orientation='h', y=-0.1))
    return fig


def save_publication(payload, output_dir, config, plotter=None):
    """Write exact plotted samples, variance fractions, diagnostics and export outcomes."""
    directory = Path(output_dir) / config.subdirectory
    points, spectrum = [], []
    for state in payload['states']:
        projection = state['projections']
        for i, pc1 in enumerate(projection.get('PC1', [])):
            points.append(dict(state=state['name'], sample_index=i+1, PC1=pc1, PC2=projection['PC2'][i] if projection.get('PC2') else None,
                PC3=projection['PC3'][i] if projection.get('PC3') else None, basis_id=state['basis_id'],
                basis_verified=state['basis_verified'], unit=state['unit'], source=state.get('source', '')))
        for pc, (fraction, cumulative) in enumerate(zip(state['fractions'], state['cumulative']), start=1):
            spectrum.append(dict(state=state['name'], component=pc, fraction=fraction, cumulative=cumulative,
                                 displayed=pc <= config.scree_components, source=state.get('eigenvalue_source', '')))
    diagnostics_csv = write_source_table(pd.DataFrame(payload['diagnostics']), directory, 'pca_diagnostics')
    projections_csv = write_source_table(pd.DataFrame(points, columns=['state', 'sample_index', 'PC1', 'PC2', 'PC3', 'basis_id', 'basis_verified', 'unit', 'source']), directory, 'pca_projections')
    variance_csv = write_source_table(pd.DataFrame(spectrum, columns=['state', 'component', 'fraction', 'cumulative', 'displayed', 'source']), directory, 'pca_variance')
    exports = {}
    if payload['availability']['composite']['status'] == 'available':
        fig = plotter.plot_publication_composite(payload) if plotter else composite_figure(payload)
        exports = export_publication_figure(fig, directory, 'pca_composite', config)
    else:
        # Independent state figures remain useful even when a shared overlay is unsafe.
        for index, state in enumerate(payload['states'], start=1):
            single = prepare_payload([state], scree_components=config.scree_components)
            if single['availability']['composite']['status'] != 'available':
                continue
            exports[f'state_{index}'] = export_publication_figure(composite_figure(single), directory, f'pca_state_{index}', config)
    metadata = dict(state_labels=[s['name'] for s in payload['states']],
        basis_id=payload['states'][0]['basis_id'] if payload['states'] and len({s['basis_id'] for s in payload['states']}) == 1 else None,
        states=[{k: v for k, v in s.items() if k != 'projections'} for s in payload['states']],
        interpretation=CAPTION, scree_components=config.scree_components, window=None, smoothing_width=None,
        availability=payload['availability'], diagnostics=payload['diagnostics'], exports=exports,
        source_files=[projections_csv, variance_csv, diagnostics_csv])
    metadata_json = write_metadata(metadata, directory, 'pca_metadata')
    return dict(output_dir=str(directory), availability=payload['availability'], exports=exports,
                diagnostics_csv=diagnostics_csv, projections_csv=projections_csv, variance_csv=variance_csv, metadata_json=metadata_json)


def joint_payload(result, config):
    """Adapt a computed joint result, whose shared projection basis is known."""
    import hashlib
    fingerprint = hashlib.sha256(result.eigenvectors.tobytes() + result.mean_structure.tobytes()).hexdigest()
    states = []
    system_configs = {system.name: system for system in config.systems}
    color_slot = 0
    for name in result.system_names:
        system = system_configs[name]
        color = system.color or (APO_HEX if system.is_apo else categorical(color_slot))
        if not system.is_apo:
            color_slot += 1
        projection = result.projections[name]
        states.append(dict(name=name, projections={f'PC{i+1}': projection[:, i].tolist() for i in range(min(3, projection.shape[1]))},
            eigenvalues=result.eigenvalues.tolist(), total_variance=result.total_variance,
            basis_id=fingerprint, basis_verified=True, basis_source='joint_pca: common covariance eigenvectors and mean',
            unit='Å', color=color, source=str(config.output_dir / 'data' / f'projections_{name}.csv'),
            eigenvalue_source=str(config.output_dir / 'data' / 'eigenvalues.csv'),
            input_provenance=dict(topology=str(system.topology), trajectories=[str(t) for t in system.trajectories],
                                  selection=config.selection, stride=config.stride, fit_iterations=config.fit_iterations)))
    records = []
    for i, j in itertools.combinations(range(len(result.system_names)), 2):
        records.append(diagnostic('rmsip', float(result.rmsip_pairwise[i, j]),
            state_pair=f'{result.system_names[i]}|{result.system_names[j]}',
            window=f'first {result.parameters["n_overlap"]} components', source='data/rmsip_pairwise.csv'))
    return prepare_payload(states, records, config.publication.scree_components)
