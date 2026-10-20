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
STATS = _project_io.STATS
FIGS = _PROJECT_ROOT / 'output' / 'figure_2'
FIGS.mkdir(parents=True, exist_ok=True)
_spec = _ilu.spec_from_file_location('figstyle', _PROJECT_ROOT / 'common' / 'figstyle.py')
_fs = _ilu.module_from_spec(_spec)
_sys.modules['figstyle'] = _fs
_spec.loader.exec_module(_fs)
from figstyle import *
import pandas as pd
from scipy.stats import gaussian_kde
d = pd.read_csv(STATS / 'analysis_set_final.csv', usecols=['RNHI', 'pct', 'c_MyeProg', 'c_Mono', 'c_Gran', 'c_Mega', 'c_EarlyEry', 'c_LateEry', 'c_B', 'c_NK', 'c_DC'])
h = pd.read_csv(STATS / 'healthy_donor_RNHI.csv', usecols=['RNHI', 'p_MyeProg', 'p_Mono', 'p_Gran', 'p_Mega', 'p_EarlyEry', 'p_LateEry', 'p_B', 'p_NK', 'p_DC'])
L9 = ['MyeProg', 'Mono', 'Gran', 'Mega', 'EarlyEry', 'LateEry', 'B', 'NK', 'DC']
LAB = {'MyeProg': 'Myeloid prog.', 'Mono': 'Monocyte', 'Gran': 'Granulocytic', 'Mega': 'Megakaryocytic', 'EarlyEry': 'Early erythroid', 'LateEry': 'Late erythroid', 'B': 'B', 'NK': 'NK', 'DC': 'Dendritic'}
W, H = (8.27, 158 / 72)
fig = plt.figure(figsize=(W, H))
axa = fig.add_axes([0.075, 0.24, 0.33, 0.62])
x = np.linspace(-16.5, 0.5, 600)
kp = gaussian_kde(d.RNHI, bw_method=0.12)
kh = gaussian_kde(h.RNHI, bw_method=0.45)
axa.hist(d.RNHI, bins=np.arange(-16.5, 0.6, 0.25), density=True, color='#9ECAE1', alpha=0.55, lw=0)
axa.plot(x, kp(x), color=LOW, lw=1.2, label=f'AML patients (n = {len(d):,})')
axa.plot(x, kh(x), color=HEALTHY, lw=1.2, ls=(0, (4, 2)), label=f'Healthy donors (n = {len(h)})')
med = d.RNHI.median()
hm = h.RNHI.median()
axa.axvline(med, color=LOW, lw=0.7, ls=':')
axa.axvline(hm, color=HEALTHY, lw=0.7, ls=':')
axa.text(med + 0.25, 0.5, 'median ' + f'{med:.2f}'.replace('-', '−'), color=LOW, fontsize=5.8, va='top', ha='left')
axa.text(hm - 0.25, 0.33, 'median ' + f'{hm:.2f}'.replace('-', '−'), color=HEALTHY, fontsize=5.8, va='top', ha='right')
axa.set_xlim(-16.5, 0.5)
axa.set_ylim(0, 0.65)
axa.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:g}'.replace('-', '−')))
axa.set_xlabel('RNHI (negative Aitchison distance to the healthy center)', labelpad=1.5)
axa.set_ylabel('Density')
axa.legend(loc='upper left', fontsize=5.8, handlelength=1.8, borderaxespad=0.1, bbox_to_anchor=(0.3, 1.0))
axa.set_title(f'RNHI in 8 AML cohorts and healthy donors', fontsize=6.6, fontweight='bold', loc='left', pad=2)
NB, NBH, REP = (40, 5, 2)
C = d[['c_' + l for l in L9]].to_numpy(float)
C = C / C.sum(1, keepdims=True)
Hh = h[['p_' + l for l in L9]].to_numpy(float)
hmean = Hh.mean(0)
b = np.minimum((d.pct.to_numpy() * NB).astype(int), NB - 1)
P = np.array([C[b == k].mean(0) for k in range(NB)])[::-1]
hb = np.array_split(np.argsort(-h.RNHI.to_numpy()), NBH)
PH = np.array([Hh[i].mean(0) for i in hb])
lg = lambda M: np.log2(np.maximum(M, 0.0001) / hmean).T
M = np.c_[lg(PH).repeat(REP, axis=1), np.full((9, 2), np.nan), lg(P).repeat(REP, axis=1)]
axb = fig.add_axes([0.53, 0.24, 0.4, 0.62])
cm = plt.get_cmap('RdBu_r').copy()
cm.set_bad('white')
im = axb.imshow(M, aspect='auto', cmap=cm, vmin=-6, vmax=6, interpolation='nearest', rasterized=True)
axb.set_yticks(range(9))
axb.set_yticklabels([LAB[l] for l in L9], fontsize=5.9)
axb.set_xticks([])
axb.tick_params(axis='y', length=0)
for s in ['top', 'right', 'bottom', 'left']:
    axb.spines[s].set_visible(False)
nh = NBH * REP
axb.text(nh / 2 - 0.5, -0.9, 'Healthy\ndonors', ha='center', va='bottom', fontsize=5.6, color=HEALTHY, fontweight='bold')
axb.text(nh + 2 + NB * REP / 2, -0.9, f'AML patients by within-cohort RNHI percentile, high → low (n = {len(d):,})', ha='center', va='bottom', fontsize=5.8, color=LOW, fontweight='bold')
axb.set_xlabel(f'Percentile bins ({NB} AML, {NBH} healthy)', labelpad=1.5)
cax = fig.add_axes([0.94, 0.3, 0.008, 0.5])
cb = fig.colorbar(im, cax=cax)
cb.set_ticks([-6, 0, 6])
cb.set_ticklabels(['≤−6', '0', '≥6'])
cb.ax.tick_params(labelsize=5.4, length=1.5, width=0.5, pad=1)
cb.outline.set_linewidth(0.4)
cb.set_label('log$_2$ ratio to\nhealthy mean', fontsize=5.4, labelpad=1)
fig.savefig(str(FIGS) + '/Fig2AB_panel.pdf')
fig.savefig(str(FIGS) + '/Fig2AB_panel.png', dpi=220)
print('ok')
