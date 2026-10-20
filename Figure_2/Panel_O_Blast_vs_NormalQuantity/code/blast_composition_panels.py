from pathlib import Path as _ProjectPath
import sys as _project_sys
_PROJECT_ROOT = next((p for p in _ProjectPath(__file__).resolve().parents if (p / 'common' / 'panel_io.py').is_file()), None)
if _PROJECT_ROOT is None:
    raise RuntimeError('The common analysis module is required.')
_project_sys.path.insert(0, str(_PROJECT_ROOT / 'common'))
import panel_io as _project_io
import importlib.util as _ilu
import sys as _sys
from pathlib import Path as _Path
HERE = _Path(__file__).resolve().parent
ROOT = _PROJECT_ROOT
STATS = _project_io.STATS
FIGS = _PROJECT_ROOT / 'output' / 'figure_2'
FIGS.mkdir(parents=True, exist_ok=True)
_spec = _ilu.spec_from_file_location('figstyle', _PROJECT_ROOT / 'common' / 'figstyle.py')
_fs = _ilu.module_from_spec(_spec)
_sys.modules['figstyle'] = _fs
_spec.loader.exec_module(_fs)
from figstyle import *
import pandas as pd
from scipy.stats import spearmanr
m = pd.read_csv(STATS / 'cohort_1_clinical_merged.csv', dtype={'sample_id': str}, usecols=['sample_id', 'blast_pct', 'nt', 'RNHI', 'c_HSC_MPP', 'c_LMPP', 'c_MyeProg', 'c_Mono', 'c_Gran', 'c_Mega', 'c_EarlyEry', 'c_LateEry', 'c_B', 'c_T', 'c_NK', 'c_DC'])
m['sample_id'] = _project_io.anonymous_ids(m['sample_id'], 'Sample')
b = m.dropna(subset=['blast_pct']).copy()
L12 = ['HSC_MPP', 'LMPP', 'MyeProg', 'Mono', 'Gran', 'Mega', 'EarlyEry', 'LateEry', 'B', 'T', 'NK', 'DC']
GROUPS = [('HSPC', ['HSC_MPP', 'LMPP', 'MyeProg'], '#4477AA'), ('Myeloid', ['Mono', 'Gran', 'DC'], '#228833'), ('Ery/Mk', ['EarlyEry', 'LateEry', 'Mega'], '#EE6677'), ('Lymphoid', ['B', 'T', 'NK'], '#AA3377')]
example_rows = pd.read_csv(STATS / 'figure2_examples.csv', usecols=['sample_id', 'group'], dtype={'sample_id': str})
example_rows['sample_id'] = _project_io.anonymous_ids(example_rows['sample_id'], 'Example sample')
if len(example_rows) != 2 or set(example_rows['group']) != {'RNHI-low', 'RNHI-high'} or (not example_rows['sample_id'].is_unique):
    raise ValueError('Exactly one example from each RNHI group is required.')
if not set(example_rows['sample_id']).issubset(set(m['sample_id'])):
    raise ValueError('Example samples must match the clinical analysis table.')
EXAMPLES = list(example_rows[['sample_id', 'group']].itertuples(index=False, name=None))

def lowess(x, y, frac=0.75, grid=None):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    grid = np.linspace(x.min(), x.max(), 80) if grid is None else grid
    k = max(3, int(np.ceil(frac * len(x))))
    out = np.empty(len(grid))
    for i, g in enumerate(grid):
        dist = np.abs(x - g)
        h = np.sort(dist)[k - 1] + 1e-12
        w = np.clip(1 - (dist / h) ** 3, 0, None) ** 3
        W = np.sum(w)
        xm = np.sum(w * x) / W
        ym = np.sum(w * y) / W
        sxx = np.sum(w * (x - xm) ** 2)
        beta = np.sum(w * (x - xm) * (y - ym)) / sxx if sxx > 0 else 0.0
        out[i] = ym + beta * (g - xm)
    return (grid, out)

def smooth_band(ax, x, y, seed=11, B=200):
    grid, fit = lowess(x, y)
    rng = np.random.default_rng(seed)
    boot = np.array([lowess(x[idx], y[idx], grid=grid)[1] for idx in (rng.integers(0, len(x), len(x)) for _ in range(B))])
    lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
    ax.fill_between(grid, lo, hi, color='#999999', alpha=0.35, lw=0, zorder=2)
    ax.plot(grid, fit, color='#404040', lw=0.8, zorder=3)

def sci_p(p):
    if p >= 0.001:
        return f'P = {p:.3f}' if p < 0.1 else f'P = {p:.2f}'
    e = int(np.floor(np.log10(p)))
    mnt = round(p / 10 ** e, 1)
    return f'P = {mnt:.1f} $\\times$ 10$^{{{e}}}$'.replace('-', '−')
W_PT, H_PT = (595.276 - 150.0, 841.89 - 668.0)
fig = plt.figure(figsize=(W_PT / 72, H_PT / 72))
fx = lambda x_pt: x_pt / W_PT
fy = lambda y_pt: 1 - y_pt / H_PT
rows = []
for k, (col, ylab, xl) in enumerate([('nt', 'NormalQuantity\n(sum of normal state fractions)', 42), ('RNHI', 'RNHI', 170)]):
    letter = 'OP'[k]
    fig.text(fx(xl - 40), fy(4), letter, fontsize=18, color='#231815', va='top', ha='left')
    ax = fig.add_axes([fx(xl), fy(136), fx(88), (136 - 42) / H_PT])
    x = b.blast_pct.to_numpy(float)
    y = b[col].to_numpy(float)
    rho = spearmanr(x, y)
    rows.append(dict(panel=letter, x='blast_pct', y=col, n=len(b), rho=rho.correlation, P=rho.pvalue))
    ax.scatter(x, y, s=5.5, color='#0072B2', alpha=0.85, lw=0, zorder=4, rasterized=False)
    smooth_band(ax, x, y)
    ax.set_xlim(-3, 103)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.grid(color='#EBEBEB', lw=0.5)
    ax.set_axisbelow(True)
    ax.set_xlabel('Measured blast (%)', fontsize=6.1, labelpad=1.5)
    ax.set_ylabel(ylab, fontsize=6.1, labelpad=2, linespacing=1.15)
    ax.tick_params(labelsize=6.0, length=2)
    if col == 'RNHI':
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:g}'.replace('-', '−')))
    ax.set_title(f'Cohort 1 (n = {len(b)})\nρ = {rho.correlation:.2f}, '.replace('-', '−') + sci_p(rho.pvalue), fontsize=7.0, color='#333333', fontweight='normal', pad=3, linespacing=1.3)
pd.DataFrame(rows).to_csv(STATS / 'figure2_O_P_statistics.csv', index=False)
fig.text(fx(266), fy(4), 'Q', fontsize=18, color='#231815', va='top', ha='left')
hd = pd.read_csv(_PROJECT_ROOT / 'data' / 'reference' / 'healthy_donor_counts.csv', usecols=['sample_id', 'HSC_MPP', 'LMPP', 'MyeProg', 'Mono', 'Gran', 'Mega', 'EarlyEry', 'LateEry', 'B', 'T', 'NK', 'DC'], dtype={'sample_id': str}).set_index('sample_id')
_project_io.anonymous_ids(hd.index, 'Reference sample')
H = hd[L12].astype(float)
Hs = H.div(H.sum(axis=1), axis=0)
bars, labels = ([], [])
for sid, lab in EXAMPLES:
    r = m[m.sample_id == sid].iloc[0]
    v = r[['c_' + l for l in L12]].astype(float).to_numpy()
    v = v / v.sum()
    s = dict(zip(L12, v))
    bars.append([100 * sum((s[x] for x in xs)) for _, xs, _ in GROUPS])
    labels.append(f'{lab}\nblast {r.blast_pct:.0f}%, RNHI ' + f'{r.RNHI:.1f}'.replace('-', '−'))
bars.append([100 * Hs[xs].sum(axis=1).mean() for _, xs, _ in GROUPS])
labels.append(f'Healthy reference\n(n = {len(hd)} donors)')
axq = fig.add_axes([fx(332), fy(150), fx(100), (150 - 40) / H_PT])
for i, vals in enumerate(bars):
    left = 0.0
    for (name, _, colr), v in zip(GROUPS, vals):
        axq.barh(i, v, left=left, height=0.62, color=colr, lw=0.6, edgecolor='white')
        if v >= 14:
            axq.text(left + v / 2, i, f'{v:.0f}%', ha='center', va='center', fontsize=6.0, color='white')
        left += v
axq.set_ylim(len(bars) - 0.5, -0.5)
axq.set_yticks(range(len(bars)))
axq.set_yticklabels(labels, fontsize=5.9, linespacing=1.15)
axq.set_xlim(0, 100)
axq.set_xticks([0, 25, 50, 75, 100])
axq.set_xticklabels(['0%', '25%', '50%', '75%', '100%'])
axq.tick_params(axis='x', labelsize=6.0, length=2)
axq.tick_params(axis='y', length=0)
axq.set_xlabel('Share of residual normal hematopoiesis', fontsize=6.1, labelpad=1.5)
for s in ['top', 'right']:
    axq.spines[s].set_visible(False)
from matplotlib.patches import Patch
axq.legend(handles=[Patch(color=c, label=n) for n, _, c in GROUPS], ncol=4, loc='lower right', bbox_to_anchor=(1.0, 1.02), fontsize=5.8, handlelength=1.1, handletextpad=0.35, columnspacing=0.7, borderaxespad=0.1, frameon=False)
fig.savefig(FIGS / 'Fig2OPQ_panel.pdf', transparent=True)
fig.savefig(FIGS / 'Fig2OPQ_panel.png', dpi=300)
print(pd.DataFrame(rows).to_string(index=False))
