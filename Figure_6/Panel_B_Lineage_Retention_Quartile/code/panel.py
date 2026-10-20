import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
CODE = HERE.parents[2]
sys.path.insert(0, str(CODE / 'common'))
import panel_io as io
import pandas as pd
_fs = io.figstyle()
globals().update({k: v for k, v in vars(_fs).items() if not k.startswith('_')})
BBOX = (0.636, 0.492, 0.9, 0.968)
cor = pd.read_csv(io.stat('bulk_lineage_correlation.csv')).sort_values('rho', ascending=True)
ret = pd.read_csv(io.stat('lineage_retention_by_quartile.csv'), index_col=0)

def draw(fig):
    panel_letter(fig, 0.53, 0.955, 'B')
    order = [k for k in cor.sort_values('rho', ascending=False).key if k in ret.columns]
    NAMES = dict(zip(cor.key, cor.lineage))
    M = ret.loc[['Q1', 'Q2', 'Q3', 'Q4'], order].to_numpy(float).T
    axb = fig.add_axes([0.665, 0.6, 0.2, 0.325])
    cm = plt.get_cmap('Blues').copy()
    cm.set_bad('white')
    im = axb.imshow(M, aspect='auto', cmap=cm, vmin=0, vmax=100, interpolation='nearest')
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.isnan(M[i, j]):
                continue
            axb.text(j, i, '%.0f' % M[i, j], ha='center', va='center', fontsize=5.6, color='white' if M[i, j] > 55 else INK)
    axb.set_yticks(range(len(order)))
    axb.set_yticklabels([NAMES[k] for k in order], fontsize=6.2)
    axb.set_xticks([0, 1, 2, 3])
    axb.set_xticklabels(['Q1\nlowest', 'Q2', 'Q3', 'Q4\nhighest'], fontsize=5.8, linespacing=1.1)
    axb.tick_params(length=0)
    for s in ['top', 'right', 'bottom', 'left']:
        axb.spines[s].set_visible(False)
    nl = ret.loc[['Q1', 'Q2', 'Q3', 'Q4'], 'n_lineages_mean'].to_numpy()
    nq = ret.loc[['Q1', 'Q2', 'Q3', 'Q4'], 'n'].astype(int).to_numpy()
    axb.set_xlabel('Within-cohort RNHI quartile (n = %d–%d per quartile)' % (nq.min(), nq.max()), labelpad=3)
    axb.set_title('Patients with a detectable lineage (%)', fontsize=6.7, fontweight='bold', loc='left', pad=4)
    axn = fig.add_axes([0.665, 0.505, 0.2, 0.034])
    axn.bar(range(4), nl, color=[Q_COLORS[k] for k in range(4)], width=0.62)
    for j, v in enumerate(nl):
        axn.text(j, v + 0.15, '%.1f' % v, ha='center', va='bottom', fontsize=5.4, color=INK2)
    axn.set_xlim(-0.5, 3.5)
    axn.set_ylim(0, 7.5)
    axn.set_xticks([])
    axn.set_yticks([])
    for s in ['top', 'right', 'left']:
        axn.spines[s].set_visible(False)
    fig.text(0.665 - 0.008, 0.522, 'Detectable lineages\n(mean per patient)', fontsize=5.6, color=INK2, ha='right', va='center', linespacing=1.15)
    cax = fig.add_axes([0.885, 0.69, 0.007, 0.15])
    cb = fig.colorbar(im, cax=cax)
    cb.set_ticks([0, 50, 100])
    cb.ax.tick_params(labelsize=5.4, length=1.5, width=0.5)
    cb.outline.set_linewidth(0.4)
    return axb
