from pathlib import Path as _ProjectPath
import sys as _project_sys
_PROJECT_ROOT = next((p for p in _ProjectPath(__file__).resolve().parents if (p / 'common' / 'panel_io.py').is_file()), None)
if _PROJECT_ROOT is None:
    raise RuntimeError('The common analysis module is required.')
_project_sys.path.insert(0, str(_PROJECT_ROOT / 'common'))
import panel_io as _project_io
import csv, os, re, sys, json, math, collections
HERE = os.path.dirname(os.path.abspath(__file__))
A4 = os.path.dirname(HERE)
BASIS = str(_PROJECT_ROOT / 'data' / 'figure_1' / 'marker_heatmap' / 'music_basis_long.csv')
OUT = str(_PROJECT_ROOT / 'data' / 'figure_1' / 'marker_heatmap')
os.makedirs(OUT, exist_ok=True)
ORDER = ['Normal_HSC_MPP', 'Normal_LMPP', 'Normal_MyeloidProgenitor', 'Normal_Monocyte', 'Normal_Granulocytic', 'Normal_Megakaryocytic', 'Normal_EarlyErythroid', 'Normal_LateErythroid', 'Normal_B', 'Normal_T', 'Normal_NK', 'Normal_Dendritic', 'AML_LSPC_Quiescent', 'AML_LSPC_Primed', 'AML_LSPC_Cycling', 'AML_GMP_like', 'AML_ProMono_like', 'AML_Mono_like', 'AML_cDC_like']
MARKERS = {'Normal_HSC_MPP': ['CD34', 'AVP', 'CRHBP', 'HLF', 'PROM1', 'MEIS1', 'MLLT3', 'HOPX', 'SPINK2', 'NPR3', 'MYCT1', 'BEX1', 'BEX2', 'GATA2'], 'Normal_LMPP': ['CD34', 'FLT3', 'SPINK2', 'DNTT', 'ADA', 'SPRY1', 'UHRF1', 'LAT2', 'HOXA7', 'IL7R', 'SATB1', 'PROM1', 'MEIS1'], 'Normal_MyeloidProgenitor': ['MPO', 'ELANE', 'PRTN3', 'AZU1', 'CTSG', 'MS4A3', 'PRSS57', 'CEBPA', 'SPI1', 'CD34'], 'Normal_Monocyte': ['LYZ', 'CD14', 'FCN1', 'VCAN', 'MNDA', 'S100A8', 'S100A9', 'S100A12', 'CSF1R', 'FCGR3A', 'LST1', 'IFITM3', 'CD68', 'CTSS', 'CTSB', 'TYROBP', 'MS4A7', 'S100A4', 'S100A6'], 'Normal_Granulocytic': ['FCGR3B', 'CSF3R', 'CAMP', 'LCN2', 'DEFA3', 'DEFA4', 'CEACAM8', 'MMP8', 'MMP9', 'S100A8', 'S100A9', 'MS4A3', 'CLC', 'MS4A2', 'HPGDS', 'CSF2RB', 'ELANE', 'MPO', 'PRTN3', 'AZU1', 'CTSG'], 'Normal_Megakaryocytic': ['PF4', 'PPBP', 'ITGA2B', 'GP9', 'GP1BB', 'VWF', 'TREML1', 'SELP', 'NFE2', 'THBS1', 'CLU', 'GNG11', 'NRGN', 'CAVIN2', 'MMD', 'LTBP1', 'SPARC', 'TUBB1', 'GP6'], 'Normal_EarlyErythroid': ['GATA1', 'KLF1', 'CA1', 'TFRC', 'KIT', 'GYPA', 'CD36', 'EPOR', 'ANK1', 'ALAS2', 'AHSP', 'HES6', 'ADD2'], 'Normal_LateErythroid': ['HBB', 'HBA1', 'HBA2', 'HBD', 'AHSP', 'ALAS2', 'SLC4A1', 'GYPA', 'BPGM', 'ANK1', 'EPB42', 'BNIP3L', 'SLC25A37', 'SLC2A1', 'TMCC2'], 'Normal_B': ['CD79A', 'CD79B', 'MS4A1', 'CD19', 'CD37', 'CD74', 'HLA-DRA', 'IGKC', 'VPREB3', 'TCL1A', 'BANK1', 'CD22', 'CD24', 'CD83'], 'Normal_T': ['CD3D', 'CD3E', 'CD3G', 'CD2', 'IL7R', 'LCK', 'CD8A', 'CD4', 'CCR7', 'SELL', 'TRBC2', 'CD247', 'CD27', 'SPOCK2'], 'Normal_NK': ['NKG7', 'GNLY', 'KLRD1', 'KLRF1', 'PRF1', 'GZMB', 'GZMA', 'GZMH', 'NCAM1', 'FCGR3A', 'TYROBP', 'CTSW', 'CCL5', 'CCL4', 'SPON2', 'FGFBP2', 'CLIC3', 'KLRB1', 'NCR1', 'XCL1'], 'Normal_Dendritic': ['FCER1A', 'CLEC10A', 'CD1C', 'CST3', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DRA', 'HLA-DQA1', 'IRF8', 'CD74', 'LILRA4', 'IRF7', 'CLEC4C', 'SCT', 'SERPINF1', 'PLD4', 'TCF4'], 'AML_LSPC_Quiescent': ['CD34', 'PROM1', 'HLF', 'AVP', 'HOPX', 'MEIS1', 'MECOM', 'SPINK2', 'NPR3', 'GATA2', 'MYCT1'], 'AML_LSPC_Primed': ['CD34', 'PROM1', 'SPINK2', 'GATA2', 'MEIS1', 'FLT3', 'SATB1', 'MPO', 'CEBPA', 'SPI1', 'MS4A3', 'ANGPT1', 'PRSS57'], 'AML_LSPC_Cycling': ['MKI67', 'TOP2A', 'PCNA', 'TYMS', 'CCNB1', 'CDK1', 'BIRC5', 'UBE2C', 'CENPF', 'NUSAP1', 'ASPM', 'CHEK1', 'CDC7', 'DHFR', 'RRM2'], 'AML_GMP_like': ['MPO', 'ELANE', 'PRTN3', 'AZU1', 'CTSG', 'CSF3R', 'SPI1', 'CEBPA', 'CEACAM8', 'MS4A3', 'KIT', 'CD34', 'PRSS57'], 'AML_ProMono_like': ['LYZ', 'FCN1', 'S100A8', 'S100A9', 'MNDA', 'VCAN', 'CTSG', 'CEBPE', 'SPI1', 'CFD', 'MS4A3', 'CEBPA', 'PRTN3'], 'AML_Mono_like': ['LYZ', 'CD14', 'FCN1', 'VCAN', 'MNDA', 'CSF1R', 'S100A12', 'FCGR3A', 'LST1', 'IFITM3', 'CD68', 'CTSS', 'CTSB', 'TYROBP', 'MS4A7', 'S100A4'], 'AML_cDC_like': ['FCER1A', 'CLEC10A', 'CD1C', 'CST3', 'HLA-DPA1', 'IRF8', 'CD74', 'LILRA4', 'CLEC4C', 'IRF7', 'PLD4', 'SCT', 'HLA-DRA']}
BAD = re.compile('^(MT-|RPL|RPS|MRPL|MRPS|EIF|EEF|FAU$|TMSB|B2M$|ACTB$|ACTG1$|GAPDH$|TUBB$|TUBA|HSP|DNAJ|SNR|SRSF|HNRNP|NPM1$|-AS|MALAT1$|NEAT1$|XIST$|LINC|MIR|SNORD|SNORA|RPPH|POLR|PSM|ATP5|COX[0-9]|NDUF|UQCR|SDH|CYCS$|ENO1$|PKM$|LDHA$|PGK1$|VIM$|EEF1)')
MIN_LABEL_PER_BLOCK = 3
MAX_PER_BLOCK = 12
TOP_ADDON_PER_STATE = 5

def num(x):
    try:
        return float(x)
    except:
        return 0.0
mt = collections.defaultdict(dict)
with open(BASIS, encoding='utf-8-sig') as fh:
    for r in csv.DictReader(fh):
        mt[r['Gene']][r['CellType']] = num(r['MTheta'])
genes_all = sorted(mt)
Z = {}
for g in genes_all:
    row = [mt[g].get(s, 0.0) for s in ORDER]
    m = sum(row) / len(row)
    sd = (sum(((v - m) ** 2 for v in row)) / (len(row) - 1)) ** 0.5 or 1.0
    Z[g] = [(v - m) / sd for v in row]
print('basis genes:', len(genes_all), 'states:', len(ORDER))
cand_states = collections.defaultdict(list)
for s in ORDER:
    for k, g in enumerate(MARKERS[s]):
        if g in Z:
            cand_states[g].append((ORDER.index(s), k + 1))
panel = sorted(cand_states)
print('panel markers present in basis:', len(panel))

def prio_in(g, j):
    ps = [p for a, p in cand_states[g] if a == j]
    return min(ps) if ps else min((p for _, p in cand_states[g]))
Z_GATE = float(os.environ.get('Z_GATE', '1.0'))
MIN_ROWS_PER_BLOCK = int(os.environ.get('MIN_ROWS_PER_BLOCK', '3'))
assign = {}
for g in panel:
    ok = [(j, Z[g][j]) for j, pr in cand_states[g] if Z[g][j] >= Z_GATE]
    if not ok:
        continue
    assign[g] = max(ok, key=lambda t: (t[1], -min((p for a, p in cand_states[g] if a == t[0]))))[0]
dropped = [g for g in panel if g not in assign]
counts = collections.Counter((ORDER[j] for j in assign.values()))
forced = set()
for s_ in ORDER:
    j = ORDER.index(s_)
    guard = 0
    while counts[s_] < MIN_ROWS_PER_BLOCK and guard < 30:
        guard += 1
        donors = [g for g in panel if g in assign and assign[g] != j and any((a == j for a, _ in cand_states[g])) and (counts[ORDER[assign[g]]] > MIN_ROWS_PER_BLOCK)]
        if not donors:
            donors = [g for g in panel if g in assign and assign[g] != j and any((a == j for a, _ in cand_states[g])) and (counts[ORDER[assign[g]]] > 1)]
        if not donors:
            break
        g = min(donors, key=lambda x: (prio_in(x, j), -Z[x][j]))
        counts[ORDER[assign[g]]] -= 1
        assign[g] = j
        counts[s_] += 1
        forced.add(g)
rows = []
for g in sorted(assign, key=lambda g: (ORDER.index(ORDER[assign[g]]),)):
    j = assign[g]
    pr = min((p for a, p in cand_states[g] if a == j)) if any((a == j for a, _ in cand_states[g])) else min((p for _, p in cand_states[g]))
    rows.append(dict(feature=g, block=ORDER[j], is_marker=1, priority=pr, forced=1 if g in forced else 0, z_in_block=round(Z[g][j], 4)))
lab = {}
for s in ORDER:
    ms = [r for r in rows if r['block'] == s and r['is_marker']]
    ms.sort(key=lambda r: (r['priority'], -r['z_in_block']))
    for r in ms[:MIN_LABEL_PER_BLOCK]:
        lab[r['feature']] = 1
for r in rows:
    r['label'] = lab.get(r['feature'], 0)
with open(os.path.join(OUT, 'heatmap_matrix_v4.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['feature'] + ORDER)
    for r in sorted(rows, key=lambda r: (ORDER.index(r['block']), -r['z_in_block'])):
        g = r['feature']
        w.writerow([g] + [mt[g].get(s, 0.0) for s in ORDER])
with open(os.path.join(OUT, 'gene_meta_v4.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['feature', 'block', 'is_marker', 'priority', 'forced', 'z_in_block', 'label'])
    w.writeheader()
    for r in sorted(rows, key=lambda r: (ORDER.index(r['block']), -r['z_in_block'])):
        w.writerow(r)
with open(os.path.join(OUT, 'col_annotation.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['sample', 'group'])
    for s in ORDER:
        w.writerow([s, 'Leukemia' if s.startswith('AML_') else 'Normal'])
print('\nZ_GATE=%.2f | minimum rows per block=%d | panel=%d | included markers=%d | excluded markers=%d | labels=%d | reassigned markers=%d' % (Z_GATE, MIN_ROWS_PER_BLOCK, len(panel), len(assign), len(dropped), sum((r['label'] for r in rows)), len(forced)))
print('\nper-block rows / labels')
for s in ORDER:
    ms = [r for r in rows if r['block'] == s]
    mk = [r for r in ms if r['is_marker']]
    lb = [r['feature'] for r in mk if r['label']]
    print('  %-26s markers=%2d addon=%2d' % (s, len(mk), len(ms) - len(mk)))
    print('      labels : %s' % (', '.join(lb) if lb else 'None'))
    print('      markers: %s' % ', '.join(sorted((r['feature'] for r in mk))))
print('\nExcluded canonical markers with within-state z < %.1f:' % Z_GATE)
for g in sorted(dropped):
    j = max(range(len(ORDER)), key=lambda k: Z[g][k])
    own = [ORDER[a] for a, _ in cand_states[g]]
    print('   %-10s own=%-34s z_own_max=%.2f  global_peak=%s' % (g, ','.join(own), max((Z[g][ORDER.index(x)] for x in own)), ORDER[j]))
print('Marker heatmap inputs written.')
