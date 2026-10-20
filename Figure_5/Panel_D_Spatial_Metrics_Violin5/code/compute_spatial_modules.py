import argparse
import re
import sys
from pathlib import Path
import pandas as pd
CODE = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(CODE / 'common'))
import panel_io as io

def read_sample(path, tag):
    if not path.is_file():
        raise FileNotFoundError('A required spatial input is missing')
    allowed = {'sample_id', 'cell_id', 'group', 'x', 'y', 'sRNHI', 'state'}
    d = pd.read_csv(path, usecols=lambda name: name in allowed, dtype={'sample_id': 'string', 'cell_id': 'string'})
    required = {'cell_id', 'x', 'y', 'sRNHI', 'state'}
    if not required.issubset(d.columns):
        raise ValueError('A spatial input is missing required analysis fields')
    d['cell_id'] = io.anonymous_ids(d['cell_id'], 'Cell')
    if d['cell_id'].duplicated().any():
        raise ValueError('Cell IDs must be unique within a sample')
    if 'sample_id' in d.columns:
        sample_ids = io.anonymous_ids(d['sample_id'], 'Sample')
        if any((value != tag for value in sample_ids)):
            raise ValueError('Sample IDs do not match the anonymous input label')
    groups = {str(i): 'AML' if i <= 4 else 'Healthy' for i in range(1, 8)}
    if tag not in groups:
        raise ValueError('Spatial sample IDs must be integers from 1 through 7')
    group = groups[tag]
    if 'group' in d.columns and (not d['group'].eq(group).all()):
        raise ValueError('The analysis group does not match the anonymous sample configuration')
    if not d['state'].isin(['preserved', 'impaired']).all():
        raise ValueError('Spatial states must be preserved or impaired')
    return (d[['cell_id', 'x', 'y', 'sRNHI', 'state']].dropna(subset=['x', 'y', 'sRNHI']).reset_index(drop=True), group)

def main():
    parser = argparse.ArgumentParser(prog='fragmentation_official_modules')
    parser.add_argument('input_dir', nargs='?', type=Path, default=CODE / 'data' / 'figure_5' / 'spatial')
    parser.add_argument('output_dir', nargs='?', type=Path, default=CODE / 'output' / 'figure_5' / 'fragmentation')
    args = parser.parse_args()
    import anndata as ad
    import squidpy as sq
    import cellcharter as cc
    args.output_dir.mkdir(parents=True, exist_ok=True)
    files = []
    for path in sorted(args.input_dir.glob('*.csv')):
        match = re.fullmatch('(?:moran_)?([1-9][0-9]*)\\.csv', path.name)
        if match:
            files.append((path, match.group(1)))
    if not files:
        raise FileNotFoundError('No anonymous per-sample spatial inputs are available')
    if len({tag for _, tag in files}) != len(files):
        raise ValueError('Each anonymous sample must have exactly one spatial input')
    moran_rows, comp_rows = ([], [])
    for path, tag in files:
        d, group = read_sample(path, tag)
        if len(d) < 20:
            continue
        a = ad.AnnData(X=d[['sRNHI']].to_numpy(float), obs=d[['cell_id']].copy(), var=pd.DataFrame(index=['sRNHI']))
        a.obsm['spatial'] = d[['x', 'y']].to_numpy(float)
        sq.gr.spatial_neighbors(a, coord_type='generic', delaunay=True)
        r = sq.gr.spatial_autocorr(a, mode='moran', genes=['sRNHI'], n_perms=999, seed=20260909, n_jobs=8, show_progress_bar=False, copy=True)
        if r is None:
            r = a.uns.get('moranI')
        if r is None:
            r = a.uns.get('gearyC')
        if r is None or r.empty:
            raise RuntimeError('The spatial autocorrelation module returned no result')
        rr = r.iloc[0]
        moran_rows.append({'sample_id': tag, 'group': group, 'n': len(d), 'moran_I_official': float(rr['I']), 'p_value_sim': float(rr['pval_sim']), 'p_value_norm': float(rr['pval_norm'])})
        a2 = ad.AnnData(X=d[['sRNHI']].to_numpy(float), obs=d[['cell_id', 'state']].copy())
        a2.obs['state'] = pd.Categorical(a2.obs['state'].astype(str))
        a2.obsm['spatial'] = d[['x', 'y']].to_numpy(float)
        sq.gr.spatial_neighbors(a2, coord_type='generic', delaunay=True)
        try:
            comp = cc.gr.connected_components(a2, cluster_key='state', min_cells=1, copy=True)
        except Exception:
            comp = cc.gr.connected_components(a2, cluster_key='state', copy=True)
        a2.obs['component'] = comp
        pres = a2.obs[a2.obs['state'] == 'preserved']
        if len(pres) == 0:
            comp_rows.append({'sample_id': tag, 'group': group, 'n_preserved': 0, 'n_components_official': 0, 'largest_component_official': 0, 'largest_cc_fraction_official': 0.0, 'mean_component_size_official': 0.0})
            continue
        vc = pres['component'].value_counts()
        comp_rows.append({'sample_id': tag, 'group': group, 'n_preserved': int(len(pres)), 'n_components_official': int(vc.size), 'largest_component_official': int(vc.max()), 'largest_cc_fraction_official': round(float(vc.max()) / len(pres), 4), 'mean_component_size_official': round(float(vc.mean()), 2)})
        print(f"Sample {tag}: Moran I={moran_rows[-1]['moran_I_official']:.4f}, components={vc.size}, largest={vc.max()}/{len(pres)}")
    moran_columns = ['sample_id', 'group', 'n', 'moran_I_official', 'p_value_sim', 'p_value_norm']
    comp_columns = ['sample_id', 'group', 'n_preserved', 'n_components_official', 'largest_component_official', 'largest_cc_fraction_official', 'mean_component_size_official']
    pd.DataFrame(moran_rows, columns=moran_columns).to_csv(args.output_dir / 'per_sample_moran_official.csv', index=False)
    pd.DataFrame(comp_rows, columns=comp_columns).to_csv(args.output_dir / 'per_sample_components_official.csv', index=False)
    print('Completed spatial module analysis')
if __name__ == '__main__':
    main()
