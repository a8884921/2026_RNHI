from pathlib import Path as _ProjectPath
import sys as _project_sys
_PROJECT_ROOT = next((p for p in _ProjectPath(__file__).resolve().parents if (p / 'common' / 'panel_io.py').is_file()), None)
if _PROJECT_ROOT is None:
    raise RuntimeError('The common analysis module is required.')
_project_sys.path.insert(0, str(_PROJECT_ROOT / 'common'))
import panel_io as _project_io
import csv, os, statistics
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
R = str(_PROJECT_ROOT / 'data' / 'figure_1' / 'burden')
OUT = str(_PROJECT_ROOT / 'output' / 'figure_1' / 'burden')
_ProjectPath(OUT).mkdir(parents=True, exist_ok=True)
ps_f = str(_PROJECT_ROOT / 'data' / 'figure_1' / 'burden' / 'per_sample_metrics.csv')
mb_f = str(_PROJECT_ROOT / 'data' / 'figure_1' / 'burden' / 'music_accuracy_metrics_by_burden.csv')
rows = list(csv.DictReader(open(ps_f, encoding='utf-8-sig')))
mb = list(csv.DictReader(open(mb_f, encoding='utf-8-sig')))
BLAB = ['20%', '40%', '60%', '80%', '90%']
order = ['PB20', 'PB40', 'PB60', 'PB80', 'PB90']
xs = np.arange(len(order))

def per_group(key):
    return {b: sorted((float(r[key]) for r in rows if r['burden'] == b)) for b in order}
burden = np.array([float(r['burden_level']) for r in rows])
MET = {k: np.array([float(r[k]) for r in rows]) for k in ('R', 'RMSD', 'mAD')}
MARGIN = {'R': 0.05, 'RMSD': 0.01, 'mAD': 0.01}
rng = np.random.default_rng(20260910)

def bh(ps):
    ks = sorted(ps, key=lambda k: ps[k])
    m = len(ks)
    out = {}
    prev = 1.0
    for i, k in enumerate(reversed(ks)):
        v = min(prev, ps[k] * m / (m - i))
        out[k] = v
        prev = v
    return out
ST = {}
for name, y in MET.items():
    rho = stats.spearmanr(burden, y)
    kw = stats.kruskal(*[y[burden == l] for l in sorted(set(burden))])
    g20, g90 = (y[burden == 0.2], y[burden == 0.9])
    mw = stats.mannwhitneyu(g20, g90, alternative='two-sided', method='exact')
    d = float(np.median(g90) - np.median(g20))
    bs = np.array([np.median(rng.choice(g90, 6, True)) - np.median(rng.choice(g20, 6, True)) for _ in range(20000)])
    lo, hi = (float(v) for v in np.percentile(bs, [2.5, 97.5]))
    pw = {f'{int(l * 100)}%': float(stats.mannwhitneyu(g20, y[burden == l], alternative='two-sided', method='exact').pvalue) for l in sorted(set(burden))[1:]}
    pw = {k: float(v) for k, v in bh(pw).items()}
    ni = lo > -MARGIN[name] if name == 'R' else hi < MARGIN[name]
    ST[name] = dict(rho=float(rho.statistic), p_rho=float(rho.pvalue), p_kw=float(kw.pvalue), p_mw=float(mw.pvalue), dmed=d, ci=[lo, hi], pw=pw, ni=bool(ni), med={f'{int(l * 100)}%': float(np.median(y[burden == l])) for l in sorted(set(burden))})
with open(os.path.join(OUT, 'Fig1F-1_stats_tests.csv'), 'w', newline='', encoding='utf-8') as stream:
    writer = csv.writer(stream)
    writer.writerow(['metric', 'test', 'estimate', 'p_value', 'ci_low', 'ci_high', 'margin', 'passed'])
    for name in ('R', 'RMSD', 'mAD'):
        result = ST[name]
        writer.writerow([name, 'spearman', result['rho'], result['p_rho'], '', '', '', ''])
        writer.writerow([name, 'kruskal_wallis', '', result['p_kw'], '', '', '', ''])
        writer.writerow([name, 'mann_whitney_exact_90_vs_20', result['dmed'], result['p_mw'], result['ci'][0], result['ci'][1], '', ''])
        for burden_group, p_value in result['pw'].items():
            writer.writerow([name, f'mann_whitney_exact_bh_{burden_group}_vs_20', '', p_value, '', '', '', ''])
        writer.writerow([name, 'bootstrap_noninferiority_90_vs_20', result['dmed'], '', result['ci'][0], result['ci'][1], MARGIN[name], result['ni']])
fig = plt.figure(figsize=(13.0, 9.0))
gs = fig.add_gridspec(1, 2, left=0.075, right=0.975, top=0.8, bottom=0.325, wspace=0.24)
ax = fig.add_subplot(gs[0, 0])
Rg = per_group('R')
med = [statistics.median(Rg[b]) for b in order]
meanR = {float(m['burden_level']): float(m['R']) for m in mb}
ax.plot(xs, [meanR[v] for v in (0.2, 0.4, 0.6, 0.8, 0.9)], '--o', color='0.6', lw=1.3, ms=5, zorder=2, label='Group pooled R (module metric)')
for i, b in enumerate(order):
    ax.scatter([i] * 6, Rg[b], s=30, color='#0072B2', alpha=0.8, edgecolor='white', linewidth=0.5, zorder=3, label='Per-pseudo-bulk R (n=6 per group)' if i == 0 else None)
ax.plot(xs, med, '-o', color='#D55E00', lw=2, ms=7, zorder=4, label='Median of per-bulk R')
ax.set_xticks(xs)
ax.set_xticklabels(BLAB)
ax.set_xlabel('Leukemia burden (target AML cell fraction)', fontsize=10)
ax.set_ylabel('Pearson R (truth vs estimate)', fontsize=10)
ax.set_title('(a) Correlation R across leukemia burden', fontsize=11, loc='left', pad=16)
ax.set_ylim(0.6, 0.95)
ax.grid(alpha=0.3, axis='y')
h, l = ax.get_legend_handles_labels()
ax.legend(h, l, frameon=True, fontsize=8, ncol=1, loc='lower left', bbox_to_anchor=(0.0, 1.02, 1.0, 0.0), mode='expand', borderaxespad=0.0, handlelength=1.8)
txt_a = 'Trend across burden (Spearman, n = 30):  rho = %+.3f, p = %.3f\nKruskal-Wallis across the 5 burden groups:  p = %.3f\n90%% vs 20%%:  dR = %+.3f  (bootstrap 95%% CI %+.3f to %+.3f);  exact MWU p = %.3f\nNon-inferiority of 90%% vs 20%% at a %+.2f R margin:  %s' % (ST['R']['rho'], ST['R']['p_rho'], ST['R']['p_kw'], ST['R']['dmed'], ST['R']['ci'][0], ST['R']['ci'][1], ST['R']['p_mw'], MARGIN['R'], 'PASS' if ST['R']['ni'] else 'FAIL')
ax.text(0.015, 0.985, txt_a, transform=ax.transAxes, ha='left', va='top', fontsize=8.2, color='0.1', linespacing=1.5, zorder=6, bbox=dict(boxstyle='round,pad=0.42', facecolor='white', edgecolor='0.7', alpha=0.95))
ax_a = ax
ax = fig.add_subplot(gs[0, 1])
Vr, Vm = (per_group('RMSD'), per_group('mAD'))
for i, b in enumerate(order):
    ax.scatter([i - 0.15] * 6, Vr[b], s=26, color='#009E73', alpha=0.8, edgecolor='white', linewidth=0.4, zorder=3, label='Per-bulk RMSD (n=6 per group)' if i == 0 else None)
    ax.plot([i - 0.15], [statistics.median(Vr[b])], marker='_', ms=15, color='#009E73', lw=2.4)
    ax.scatter([i + 0.15] * 6, Vm[b], s=26, color='#CC79A7', alpha=0.8, edgecolor='white', linewidth=0.4, zorder=3, label='Per-bulk mAD (n=6 per group)' if i == 0 else None)
    ax.plot([i + 0.15], [statistics.median(Vm[b])], marker='_', ms=15, color='#CC79A7', lw=2.4)
ax.set_xticks(xs)
ax.set_xticklabels(BLAB)
ax.set_xlabel('Leukemia burden (target AML cell fraction)', fontsize=10)
ax.set_ylabel('RMSD / mAD (lower is better)', fontsize=10)
ax.set_title('(b) RMSD and mAD across leukemia burden', fontsize=11, loc='left', pad=16)
ax.set_ylim(0.0, 0.1)
ax.grid(alpha=0.3, axis='y')
h, l = ax.get_legend_handles_labels()
ax.legend(h, l, frameon=True, fontsize=8, ncol=1, loc='lower left', bbox_to_anchor=(0.0, 1.02, 1.0, 0.0), mode='expand', borderaxespad=0.0, handlelength=1.8)
txt_b = 'Trend across burden:  RMSD rho = %+.3f, p = %.3f;   mAD rho = %+.3f, p = %.3f\nKruskal-Wallis:  RMSD p = %.3f;   mAD p = %.3f\n90%% vs 20%%:  dRMSD = %+.4f, dmAD = %+.4f;  non-inferiority at +%.2f:  %s / %s' % (ST['RMSD']['rho'], ST['RMSD']['p_rho'], ST['mAD']['rho'], ST['mAD']['p_rho'], ST['RMSD']['p_kw'], ST['mAD']['p_kw'], ST['RMSD']['dmed'], ST['mAD']['dmed'], MARGIN['RMSD'], 'PASS' if ST['RMSD']['ni'] else 'FAIL', 'PASS' if ST['mAD']['ni'] else 'FAIL')
ax.text(0.015, 0.985, txt_b, transform=ax.transAxes, ha='left', va='top', fontsize=8.2, color='0.1', linespacing=1.5, zorder=6, bbox=dict(boxstyle='round,pad=0.42', facecolor='white', edgecolor='0.7', alpha=0.95))
ax_b = ax
fig.text(0.5, 0.965, 'MuSiC deconvolution accuracy across leukemia burden', ha='center', va='top', fontsize=12.5, fontweight='bold')
fig.savefig(os.path.join(OUT, 'Fig1F-1_accuracy_by_burden.png'), dpi=300, bbox_inches='tight')
fig.savefig(os.path.join(OUT, 'Fig1F-1_accuracy_by_burden.pdf'), bbox_inches='tight')
plt.close(fig)
