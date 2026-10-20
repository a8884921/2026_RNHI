import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
CODE = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(CODE / 'common'))
import panel_io as io
THRESH = -1.256

def delaunay_edges(xy):
    tri = Delaunay(xy)
    edges = set()
    for simp in tri.simplices:
        for i in range(3):
            a, b = (simp[i], simp[(i + 1) % 3])
            edges.add(tuple(sorted((a, b))))
    return np.array(list(edges))

def connected_components(n, edges):
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = (find(a), find(b))
        if ra != rb:
            parent[rb] = ra
    for a, b in edges:
        union(a, b)
    comp = {}
    for i in range(n):
        comp.setdefault(find(i), []).append(i)
    return list(comp.values())

def morans_i(values, edges):
    vals = np.asarray(values, float)
    n = len(vals)
    if n < 3:
        return np.nan
    W = np.zeros((n, n))
    for a, b in edges:
        W[a, b] = 1
        W[b, a] = 1
    rs = W.sum(1)
    rs[rs == 0] = 1
    W = W / rs[:, None]
    z = vals - vals.mean()
    if z.std() == 0:
        return np.nan
    num = n * np.sum(W * z[:, None] * z[None, :])
    den = (z ** 2).sum() * W.sum()
    return float(num / den) if den != 0 else np.nan

def run_sample(d, tag, group):
    d = d.dropna(subset=['x', 'y']).copy()
    d = d[d['sRNHI_z'].notna()]
    if len(d) < 20:
        return None
    xy = d[['x', 'y']].to_numpy(float)
    n = len(d)
    edges = delaunay_edges(xy)
    preserved = (d['sRNHI_z'].values >= THRESH).astype(int)
    comps = connected_components(n, [(a, b) for a, b in edges if preserved[a] and preserved[b]])
    sizes = sorted([len(c) for c in comps], reverse=True)
    n_pres = int(preserved.sum())
    frac = n_pres / n
    largest = sizes[0] if sizes else 0
    n_comp = len(sizes)
    mean_cc = np.mean(sizes) if sizes else 0
    boundary = sum((1 for a, b in edges if preserved[a] != preserved[b]))
    bd = boundary / max(1, n_pres)
    mi = morans_i(d['sRNHI_z'].values, edges)
    return {'sample_id': tag, 'group': group, 'n_spots': n, 'n_preserved': n_pres, 'preserved_area_fraction': round(frac, 4), 'n_connected_components': n_comp, 'largest_connected_component': int(largest), 'largest_cc_fraction': round(largest / max(1, n_pres), 4), 'mean_component_size': round(mean_cc, 2), 'boundary_density': round(bd, 4), 'morans_I_sRNHI': round(mi, 4)}

def read_spatial_input(path, columns, allowed_samples):
    if not path.is_file():
        raise FileNotFoundError('A required spatial input is missing')
    frame = pd.read_csv(path, usecols=columns, dtype={'sample_id': 'string', 'cell_id': 'string'})
    frame['sample_id'] = io.anonymous_ids(frame['sample_id'], 'Sample')
    frame['cell_id'] = io.anonymous_ids(frame['cell_id'], 'Cell')
    if not set(frame['sample_id']).issubset(allowed_samples):
        raise ValueError('Spatial sample IDs are outside the configured anonymous group')
    if frame.duplicated(['sample_id', 'cell_id']).any():
        raise ValueError('Spatial sample and cell ID pairs must be unique')
    return frame

def main():
    parser = argparse.ArgumentParser(prog='fragmentation_step1')
    parser.add_argument('output_dir', nargs='?', type=Path, default=CODE / 'output' / 'figure_5' / 'fragmentation')
    parser.add_argument('frozen_dir', nargs='?', type=Path, default=CODE / 'data' / 'figure_5' / 'spatial')
    parser.add_argument('niche_dir', nargs='?', type=Path, default=CODE / 'data' / 'figure_5' / 'spatial')
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    aml_ids = {'1', '2', '3', '4'}
    healthy_ids = {'5', '6', '7'}
    aml = read_spatial_input(args.frozen_dir / 'aml_niche_mapping_srnhi.csv', ['sample_id', 'cell_id', 'x', 'y', 'sRNHI_z'], aml_ids)
    hv = read_spatial_input(args.frozen_dir / 'healthy_heldout_validation.csv', ['sample_id', 'cell_id', 'sRNHI_z'], healthy_ids)
    hl = read_spatial_input(args.niche_dir / 'niche_labels.csv', ['sample_id', 'cell_id', 'x', 'y'], healthy_ids)
    hmerge = hv.merge(hl, on=['sample_id', 'cell_id'], how='left', validate='one_to_one')
    rows = []
    for tag, group in aml.groupby('sample_id'):
        result = run_sample(group, tag, 'AML')
        if result:
            rows.append(result)
    for tag, group in hmerge.groupby('sample_id'):
        result = run_sample(group, tag, 'Healthy')
        if result:
            rows.append(result)
    columns = ['sample_id', 'group', 'n_spots', 'n_preserved', 'preserved_area_fraction', 'n_connected_components', 'largest_connected_component', 'largest_cc_fraction', 'mean_component_size', 'boundary_density', 'morans_I_sRNHI']
    res = pd.DataFrame(rows, columns=columns)
    res.to_csv(args.output_dir / 'per_sample_fragmentation.csv', index=False)
    print(res.to_string(index=False))
    print('Completed spatial fragmentation analysis')
if __name__ == '__main__':
    main()
