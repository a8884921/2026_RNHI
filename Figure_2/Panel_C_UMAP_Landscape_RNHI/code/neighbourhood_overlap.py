from pathlib import Path as _ProjectPath
import sys as _project_sys
_PROJECT_ROOT = next((p for p in _ProjectPath(__file__).resolve().parents if (p / 'common' / 'panel_io.py').is_file()), None)
if _PROJECT_ROOT is None:
    raise RuntimeError('The common analysis module is required.')
_project_sys.path.insert(0, str(_PROJECT_ROOT / 'common'))
import panel_io as _project_io
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.stats import gaussian_kde
R = str(_PROJECT_ROOT / 'output' / 'figure_2' / 'landscape')
_ProjectPath(R).mkdir(parents=True, exist_ok=True)
D = pd.read_csv(_PROJECT_ROOT / 'data' / 'figure_2' / 'landscape' / 'Figure2C_data_frozenUMAP_withRNHI.csv', usecols=['UMAP1', 'UMAP2', 'RNHI'], encoding='utf-8-sig')
K = 30
COLS = {'Low RNHI': '#2166AC', 'Medium RNHI': '#BAB0AC', 'High RNHI': '#FDB42F'}
ORD = ['Low RNHI', 'Medium RNHI', 'High RNHI']
qs = D['RNHI'].quantile([1 / 3, 2 / 3]).to_numpy()
D['ter'] = np.where(D.RNHI <= qs[0], 'Low RNHI', np.where(D.RNHI >= qs[1], 'High RNHI', 'Medium RNHI'))
X = D[['UMAP1', 'UMAP2']].to_numpy(float)
n = len(D)
dm = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1))
np.fill_diagonal(dm, np.inf)
nbr = np.argsort(dm, axis=1)[:, :K]
ter = D['ter'].to_numpy()
is_high = ter == 'High RNHI'
D['frac_high'] = is_high[nbr].mean(1)
Rnh = D['RNHI'].to_numpy()
nbr_rnhi = Rnh[nbr]
D['nbr_min'] = nbr_rnhi.min(1)
D['nbr_max'] = nbr_rnhi.max(1)
D['nbr_span'] = D.nbr_max - D.nbr_min
is_low = ter == 'Low RNHI'
mix = np.mean([is_high[r].any() and is_low[r].any() for r in nbr]) * 100
print('Reproducibility statistics')
print('n = %d | Tercile thresholds %.4f / %.4f' % (n, qs[0], qs[1]))
for g in ORD:
    s = D.loc[D.ter == g, 'frac_high']
    print('  %-13s n=%4d median=%.4f IQR=[%.4f, %.4f]' % (g, len(s), s.median(), s.quantile(0.25), s.quantile(0.75)))
print('  Neighbourhoods containing both High and Low samples = %.1f%%' % mix)
print('  corr(RNHI, UMAP1)=%.3f corr(RNHI, UMAP2)=%.3f' % (np.corrcoef(Rnh, X[:, 0])[0, 1], np.corrcoef(Rnh, X[:, 1])[0, 1]))
gx = np.linspace(0, 1, 2001)
kde = {g: gaussian_kde(D.loc[D.ter == g, 'frac_high'].to_numpy(), bw_method=0.25) for g in ORD}
y = {g: kde[g](gx) for g in ORD}
rows = []
for i, a in enumerate(ORD):
    for b in ORD[i + 1:]:
        ovl = np.trapezoid(np.minimum(y[a], y[b]), gx)
        rows.append(dict(pair=f'{a} vs {b}', OVL=round(float(ovl), 4), median_diff=round(float(D.loc[D.ter == a, 'frac_high'].median() - D.loc[D.ter == b, 'frac_high'].median()), 4)))
ovl3 = float(np.trapezoid(np.minimum(np.minimum(y[ORD[0]], y[ORD[1]]), y[ORD[2]]), gx))
iqr = {g: (float(D.loc[D.ter == g, 'frac_high'].quantile(0.25)), float(D.loc[D.ter == g, 'frac_high'].quantile(0.75))) for g in ORD}
lo3, hi3 = (max((v[0] for v in iqr.values())), min((v[1] for v in iqr.values())))
print('\nOverlap statistics')
for r in rows:
    print('  %-32s overlap coefficient OVL=%.3f' % (r['pair'], r['OVL']))
print('  Three-group overlap area OVL3 = %.3f' % ovl3)
print('  Three-group IQR intersection = [%.3f, %.3f] (width %.3f)' % (lo3, hi3, hi3 - lo3))
print('  Neighbourhood RNHI span: median %.2f, IQR [%.2f, %.2f]; overall RNHI range %.2f (span %.0f%%)' % (D.nbr_span.median(), D.nbr_span.quantile(0.25), D.nbr_span.quantile(0.75), D.RNHI.max() - D.RNHI.min(), 100 * D.nbr_span.median() / (D.RNHI.max() - D.RNHI.min())))
cov_med = 100 * np.mean((D.nbr_min <= D.RNHI.median()) & (D.nbr_max >= D.RNHI.median()))
cov_iqr = 100 * np.mean((D.nbr_min <= D.RNHI.quantile(0.25)) & (D.nbr_max >= D.RNHI.quantile(0.75)))
cov_thr = 100 * np.mean((D.nbr_min <= qs[1]) & (D.nbr_max >= qs[1]))
print('  Neighbourhood coverage: global median (%.2f) %%=%.1f | global IQR %%=%.1f | High threshold (%.2f) %%=%.1f' % (D.RNHI.median(), cov_med, cov_iqr, qs[1], cov_thr))
pd.DataFrame(rows + [dict(pair='Three-group overlap', OVL=round(ovl3, 4), median_diff=np.nan)]).to_csv(f'{R}/Figure2C_density_overlap_stats.csv', index=False)
rcParams.update({'font.family': 'serif', 'font.serif': ['Times New Roman', 'DejaVu Serif'], 'font.size': 9, 'axes.linewidth': 0.7, 'pdf.fonttype': 42})
fig, axes = plt.subplots(2, 1, figsize=(7.4, 8.0), gridspec_kw=dict(height_ratios=[1.28, 1.0], hspace=0.42))
ax = axes[0]
ax.axvspan(lo3, hi3, color='#7A7A7A', alpha=0.13, zorder=0)
ax.axvline(1 / 3, color='#1F1F1F', ls='--', lw=1.0, zorder=1)
for g in ORD:
    ax.fill_between(gx, y[g], color=COLS[g], alpha=0.3, zorder=2)
    ax.plot(gx, y[g], color=COLS[g], lw=1.8, zorder=3, label='%s (n=%d, median %.3f)' % (g, int((D.ter == g).sum()), D.loc[D.ter == g, 'frac_high'].median()))
ax.text(1 / 3 + 0.008, ax.get_ylim()[1] * 0.97, 'random-mixing\nexpected = 1/3', fontsize=8.5, style='italic', va='top', color='#1F1F1F')
ax.annotate('', xy=(lo3, ax.get_ylim()[1] * 0.52), xytext=(hi3, ax.get_ylim()[1] * 0.52), arrowprops=dict(arrowstyle='<->', color='#4A4A4A', lw=1.0))
ax.text((lo3 + hi3) / 2, ax.get_ylim()[1] * 0.54, 'IQR overlap of all three groups\n%.2f–%.2f' % (lo3, hi3), ha='center', va='bottom', fontsize=8, color='#3A3A3A')
_ovl = {r['pair']: r['OVL'] for r in rows}
ax.text(0.985, 0.4, 'overlap coefficient (OVL)\nLow–Medium %.2f   Medium–High %.2f   Low–High %.2f\nall three %.2f' % (_ovl['Low RNHI vs Medium RNHI'], _ovl['Medium RNHI vs High RNHI'], _ovl['Low RNHI vs High RNHI'], ovl3), transform=ax.transAxes, ha='right', va='center', fontsize=8, color='#2B2B2B', bbox=dict(boxstyle='round,pad=0.4', fc='#F7F7F7', ec='#9A9A9A', lw=0.6))
ax.set_xlim(0, 0.85)
ax.set_xlabel('Fraction of the 30 nearest neighbours that are High-RNHI')
ax.set_ylabel('Density')
ax.set_title("A  Density of local high-RNHI neighbourhood fraction, by the patient's own RNHI tercile", loc='left', fontsize=10)
ax.legend(frameon=False, fontsize=8.2, loc='upper left')
ax.grid(axis='y', color='#EFF1F4', lw=0.5)
ax.set_axisbelow(True)
ax2 = axes[1]
sp = D.nbr_span.to_numpy()
sk = gaussian_kde(sp, bw_method=0.28)
sx = np.linspace(0, sp.max() * 1.05, 600)
ax2.fill_between(sx, sk(sx), color='#5B8DB8', alpha=0.32)
ax2.plot(sx, sk(sx), color='#1F4E79', lw=1.8)
for q, ls, lab in [(sp.mean(), '-', 'mean %.2f' % sp.mean()), (np.median(sp), '--', 'median %.2f' % np.median(sp))]:
    ax2.axvline(q, color='#B2182B', ls=ls, lw=1.0)
ax2.text(np.median(sp) + 0.05, ax2.get_ylim()[1] * 0.92, 'median span %.2f RNHI units\n(%.0f%% of the whole patient range %.2f)' % (np.median(sp), 100 * np.median(sp) / (D.RNHI.max() - D.RNHI.min()), D.RNHI.max() - D.RNHI.min()), fontsize=8.2, color='#B2182B', va='top')
ax2.set_xlabel("RNHI span covered by a patient's own 30-neighbourhood  (max − min)")
ax2.set_ylabel('Density')
ax2.set_title('B  Every local neighbourhood covers a wide RNHI interval — %.1f%% contain both High and Low RNHI patients' % mix, loc='left', fontsize=10)
ax2.grid(axis='y', color='#EFF1F4', lw=0.5)
ax2.set_axisbelow(True)
ax2.text(0.985, 0.62, 'neighbourhood intervals overlap:\n%.0f%% contain the global median RNHI\n%.0f%% span the whole global IQR (%.1f–%.1f)\n%.1f%% contain both High and Low patients' % (cov_med, cov_iqr, D.RNHI.quantile(0.25), D.RNHI.quantile(0.75), mix), transform=ax2.transAxes, ha='right', va='center', fontsize=8, color='#2B2B2B', bbox=dict(boxstyle='round,pad=0.4', fc='#F7F7F7', ec='#9A9A9A', lw=0.6))
ax2.text(0.985, 0.16, 'frozen Figure 2B UMAP, k = 30; no re-embedding; RNHI used for grouping only', transform=ax2.transAxes, ha='right', va='center', fontsize=7.4, color='#5A5A5A')
fig.savefig(f'{R}/Figure2C_density_overlap_neighbourRNHI.pdf', dpi=300, bbox_inches='tight')
fig.savefig(f'{R}/Figure2C_density_overlap_neighbourRNHI.png', dpi=300, bbox_inches='tight')
print('Density overlap figures saved.')
