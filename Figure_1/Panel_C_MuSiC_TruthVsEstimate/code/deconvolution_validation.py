from pathlib import Path as _ProjectPath
import sys as _project_sys
_PROJECT_ROOT = next((p for p in _ProjectPath(__file__).resolve().parents if (p / 'common' / 'panel_io.py').is_file()), None)
if _PROJECT_ROOT is None:
    raise RuntimeError('The common analysis module is required.')
_project_sys.path.insert(0, str(_PROJECT_ROOT / 'common'))
import panel_io as _project_io
import csv, math, collections
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
OUT = str(_PROJECT_ROOT / 'output' / 'figure_1' / 'deconvolution')
_ProjectPath(OUT).mkdir(parents=True, exist_ok=True)

def rd(p):
    with open(p, encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))

def pearson(x, y):
    n = len(x)
    mx = sum(x) / n
    my = sum(y) / n
    d1 = math.sqrt(sum(((a - mx) ** 2 for a in x)))
    d2 = math.sqrt(sum(((b - my) ** 2 for b in y)))
    return sum(((a - mx) * (b - my) for a, b in zip(x, y))) / (d1 * d2) if d1 * d2 else float('nan')

def metrics(rows):
    T = [float(r['Truth']) for r in rows]
    E = [float(r['Estimate']) for r in rows]
    return (pearson(T, E), math.sqrt(sum(((a - b) ** 2 for a, b in zip(T, E))) / len(T)), sum((abs(a - b) for a, b in zip(T, E))) / len(T))
PAL = {'Healthy pure': '#0072B2', 'AML pure': '#D55E00', 'Mixed': '#009E73'}
GROUPS = ('Healthy pure', 'AML pure', 'Mixed')
R = str(_PROJECT_ROOT / 'data' / 'figure_1' / 'deconvolution')
pairs = rd(str(_PROJECT_ROOT / 'data' / 'figure_1' / 'deconvolution' / 'music_accuracy_pairs.csv'))
design = {r['pseudo_sample']: r for r in rd(str(_PROJECT_ROOT / 'data' / 'figure_1' / 'deconvolution' / 'sametech_design.csv'))}
qc = {r['Sample']: float(r['R2']) for r in rd(str(_PROJECT_ROOT / 'data' / 'figure_1' / 'deconvolution' / 'music_standard_qc.csv'))}
_project_io.anonymous_ids([r['Sample'] for r in pairs], 'Sample')
_project_io.anonymous_ids(design, 'Sample')
_project_io.anonymous_ids(qc, 'Sample')
if any((r.get('group') not in GROUPS for r in design.values())):
    raise ValueError('The sample design requires explicit analysis groups.')
samples_all = sorted({r['Sample'] for r in pairs})
keep = [s for s in samples_all if int(design[s]['n_cells']) >= 200]

def grp(s):
    return design[s]['group']
pairs_ok = [dict(r) for r in pairs if r['Sample'] in keep]
cnt = collections.Counter()
dcell = collections.Counter()
for s in keep:
    g = grp(s)
    dcell[g] += int(design[s]['n_cells'])
    cnt[g] += 1
Rok = metrics(pairs_ok)
FIG_W, FIG_H = (13.5, 8.6)
fig = plt.figure(figsize=(FIG_W, FIG_H))
fig.suptitle('Figure 1E-1  Same-platform hold-out pseudo-bulk benchmark  -  MuSiC, locked 19-state reference', fontsize=13, fontweight='bold', x=0.5, y=0.995, va='top', ha='center')
fig.text(0.5, 0.96, 'True vs estimated cell-type proportions  (QC: pseudo-bulk depth >= 200 cells)', fontsize=10.5, ha='center', va='top', color='0.2')
import numpy as np
from scipy import stats as sstats
G_AML = set()
for r in pairs_ok:
    G_AML.add(r['CellType'])

def _per_sample(sub):
    T = np.array([float(r['Truth']) for r in sub])
    E = np.array([float(r['Estimate']) for r in sub])
    return (float(np.corrcoef(T, E)[0, 1]) if T.std() > 0 and E.std() > 0 else float('nan'), float(np.sqrt(np.mean((T - E) ** 2))), float(np.mean(np.abs(T - E))))
by_sample = collections.defaultdict(list)
for r in pairs_ok:
    by_sample[r['Sample']].append(r)
samp = sorted(by_sample)
per = {s: _per_sample(by_sample[s]) for s in samp}
gof = {s: grp(s) for s in samp}
_T = np.array([float(r['Truth']) for r in pairs_ok])
_E = np.array([float(r['Estimate']) for r in pairs_ok])
obsR = float(np.corrcoef(_T, _E)[0, 1])
obsRMSE = float(np.sqrt(np.mean((_T - _E) ** 2)))
obsmAD = float(np.mean(np.abs(_T - _E)))
rng = np.random.default_rng(20260910)
boot = []
for _ in range(5000):
    pick = rng.choice(samp, len(samp), replace=True)
    tt = np.concatenate([[float(r['Truth']) for r in by_sample[s]] for s in pick])
    ee = np.concatenate([[float(r['Estimate']) for r in by_sample[s]] for s in pick])
    boot.append((float(np.corrcoef(tt, ee)[0, 1]), float(np.sqrt(np.mean((tt - ee) ** 2))), float(np.mean(np.abs(tt - ee)))))
boot = np.array(boot)
ciR, ciRMSE, cimAD = (np.percentile(boot[:, i], [2.5, 97.5]) for i in range(3))
G = {g: np.array([per[s][0] for s in samp if gof[s] == g]) for g in GROUPS}
kw = sstats.kruskal(*[G[g] for g in GROUPS])
pw = {}
for a, b in (('Healthy pure', 'AML pure'), ('Healthy pure', 'Mixed'), ('AML pure', 'Mixed')):
    pw[f'{a} vs {b}'] = float(sstats.mannwhitneyu(G[a], G[b], alternative='two-sided').pvalue)
_ks = sorted(pw, key=lambda k: pw[k])
_m = len(_ks)
pwa = {}
_prev = 1.0
for _i, _k in enumerate(reversed(_ks)):
    _v = min(_prev, pw[_k] * _m / (_m - _i))
    pwa[_k] = _v
    _prev = _v
carry = np.array([sum((float(r['Estimate']) for r in by_sample[s] if r['CellType'] in G_AML and r['CellType'].startswith('AML_'))) for s in samp if gof[s] == 'Healthy pure'])
wct = sstats.wilcoxon(carry - 0.01, alternative='greater')
cboot = np.array([np.median(rng.choice(carry, len(carry), True)) for _ in range(20000)])
ct_ci = (float(np.percentile(cboot, 2.5)), float(np.percentile(cboot, 97.5)))
with open(OUT + '/Figure1E-1_stats_tests.csv', 'w', newline='', encoding='utf-8') as _f:
    import csv as _csv
    _w = _csv.writer(_f)
    _w.writerow(['quantity', 'estimate', 'ci_low', 'ci_high', 'p_value', 'test'])
    _w.writerow(['Pearson R (all QC-passed pairs)', f'{obsR:.4f}', f'{ciR[0]:.4f}', f'{ciR[1]:.4f}', '', 'bootstrap over 59 pseudo-bulks (5000 resamples)'])
    _w.writerow(['RMSE', f'{obsRMSE:.4f}', f'{ciRMSE[0]:.4f}', f'{ciRMSE[1]:.4f}', '', 'bootstrap over pseudo-bulks'])
    _w.writerow(['mAD', f'{obsmAD:.4f}', f'{cimAD[0]:.4f}', f'{cimAD[1]:.4f}', '', 'bootstrap over pseudo-bulks'])
    _w.writerow(['per-sample R, 3 groups', 'medians ' + ' / '.join((f'{g}={np.median(G[g]):.3f}' for g in GROUPS)), '', '', f'{kw.pvalue:.4f}', 'Kruskal-Wallis'])
    for k, v in pwa.items():
        _w.writerow(['per-sample R pairwise', k, '', '', f'{v:.4f}', 'Mann-Whitney U (BH-adjusted)'])
    _w.writerow(['AML-state estimate inside pure-healthy bulks', f'median={np.median(carry):.4f}; mean={carry.mean():.4f}; max={carry.max():.4f}', f'{ct_ci[0]:.4f}', f'{ct_ci[1]:.4f}', f'{wct.pvalue:.4f}', 'one-sided Wilcoxon (H1: > 1%)'])
DATA = [0.05, 0.17, 0.4, 0.63]
ax = fig.add_axes(DATA)
for g in GROUPS:
    sub = [r for r in pairs_ok if grp(r['Sample']) == g]
    ax.scatter([float(r['Truth']) for r in sub], [float(r['Estimate']) for r in sub], s=12, alpha=0.6, color=PAL[g])
ax.plot([0, 1], [0, 1], ls='--', c='black', lw=1)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xlabel('True proportion (known from single-cell composition)', fontsize=10, labelpad=9)
ax.set_ylabel('MuSiC-estimated proportion', fontsize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(alpha=0.2)
leg_handles = [Line2D([], [], marker='o', linestyle='none', markersize=7, markerfacecolor=PAL[g], markeredgecolor='white', markeredgewidth=0.6, alpha=0.85, label=g) for g in GROUPS]
leg = ax.legend(handles=leg_handles, loc='upper left', fontsize=9.5, title='Sample group', title_fontsize=9.5, frameon=True, framealpha=0.92, edgecolor='0.75', borderpad=0.55, handletextpad=0.5, labelspacing=0.45)
leg.set_zorder(6)
PX = 0.505
PITCH_K = 1.6
y_frac = 0.925
artists = []

def emit(text, fs=9.5, color='0.12', weight='normal', indent=False, extra_gap_pt=0):
    global y_frac
    if text:
        t = fig.text(PX + (0.018 if indent else 0.0), y_frac, text, ha='left', va='top', fontsize=fs, color=color, fontweight=weight)
        artists.append(t)
    y_frac -= (fs * PITCH_K + extra_gap_pt) / (72.0 * FIG_H)
emit('Accuracy (QC-passed)', 11, '0.1', 'bold', extra_gap_pt=8)
emit('   Pearson R   =  %.3f' % Rok[0], 10)
emit('   RMSE        =  %.4f' % Rok[1], 10)
emit('   mAD         =  %.4f' % Rok[2], 10)
emit('   matched points = %d   ( %d pseudo-bulks )' % (len(pairs_ok), len(keep)), 9.2)
y_frac -= 12 / (72.0 * FIG_H)
emit('Per-group samples / cells', 11, '0.1', 'bold', extra_gap_pt=8)
group_txt = {}
for g in GROUPS:
    emit(g, 10, PAL[g], 'bold')
    group_txt[g] = artists[-1]
    emit('    n = %d pseudo-bulks   |   %d cells' % (cnt[g], dcell[g]), 9.2, '0.2')
emit('    total cells = %d' % sum(dcell.values()), 9.6, '0.1', 'bold')
y_frac -= 12 / (72.0 * FIG_H)
y_frac -= 12 / (72.0 * FIG_H)
fig.canvas.draw()
ren = fig.canvas.get_renderer()
for g, t in group_txt.items():
    bb = t.get_window_extent(ren)
    xc, yc = fig.transFigure.inverted().transform((bb.x0 - 9.0, (bb.y0 + bb.y1) / 2.0))
    fig.add_artist(Line2D([xc], [yc], marker='o', markersize=6.2, linestyle='none', markerfacecolor=PAL[g], markeredgecolor='white', markeredgewidth=0.5, alpha=0.85, transform=fig.transFigure, clip_on=False, zorder=5))
tl_txt = f'R = {obsR:.3f} (95% CI {ciR[0]:.3f}-{ciR[1]:.3f})\nRMSE = {obsRMSE:.4f} ({ciRMSE[0]:.4f}-{ciRMSE[1]:.4f}); mAD = {obsmAD:.4f} ({cimAD[0]:.4f}-{cimAD[1]:.4f})\nKruskal-Wallis P = {kw.pvalue:.3f}\nCross-talk: n = {len(carry)}, median = {float(np.median(carry)):.3f} (95% CI {ct_ci[0]:.3f}-{ct_ci[1]:.3f}), P = {wct.pvalue:.3f}'
tl = fig.text(0.02, 0.935, tl_txt, ha='left', va='top', fontsize=8.1, color='0.1', linespacing=1.5, bbox=dict(boxstyle='round,pad=0.45', facecolor='white', edgecolor='0.7', alpha=0.95))
fig.savefig(OUT + '/Figure1E-1_True_vs_Estimated_scatter.png', dpi=300)
fig.savefig(OUT + '/Figure1E-1_True_vs_Estimated_scatter.pdf')
plt.close(fig)
