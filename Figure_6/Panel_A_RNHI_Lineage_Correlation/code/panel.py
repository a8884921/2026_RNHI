import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
CODE = HERE.parents[2]
sys.path.insert(0, str(CODE / 'common'))
import panel_io as io
import pandas as pd
_fs = io.figstyle()
globals().update({k: v for k, v in vars(_fs).items() if not k.startswith('_')})
BBOX = (0.012, 0.5, 0.548, 0.968)
cor = pd.read_csv(io.stat('bulk_lineage_correlation.csv')).sort_values('rho', ascending=True)

def draw(fig):
    panel_letter(fig, 0.02, 0.955, 'A')
    axa = fig.add_axes([0.17, 0.56, 0.25, 0.365])
    y = np.arange(len(cor))
    for yi, r in zip(y, cor.itertuples()):
        sig = r.lo > 0 or r.hi < 0
        c = (HIGH if r.rho > 0 else LOW) if sig and r.in_RNHI else '#B0B0B0'
        axa.plot([r.lo, r.hi], [yi, yi], color=c, lw=1.1, solid_capstyle='butt', zorder=2)
        axa.scatter([r.rho], [yi], s=16, color=c, zorder=3, edgecolor='white', linewidth=0.5)
    axa.axvline(0, color=INK, lw=0.6, zorder=1)
    axa.set_yticks(y)
    axa.set_yticklabels(cor.lineage, fontsize=6.4)
    for tl_, inr in zip(axa.get_yticklabels(), cor.in_RNHI):
        if not inr:
            tl_.set_color(MUTED)
    axa.set_ylim(-0.6, len(cor) - 0.4)
    axa.set_xlim(-0.3, 0.45)
    axa.set_xticks([-0.2, 0, 0.2, 0.4])
    axa.set_xticklabels(['−0.2', '0', '0.2', '0.4'])
    axa.set_xlabel('Within-cohort Spearman ρ (95% CI)', labelpad=2)
    axa.spines['left'].set_visible(False)
    axa.tick_params(axis='y', length=0)
    axa.grid(axis='x', color=GRID, lw=0.5)
    axa.set_axisbelow(True)
    tr = matplotlib.transforms.blended_transform_factory(axa.transAxes, axa.transData)
    for yi, r in zip(y, cor.itertuples()):
        axa.text(1.04, yi, fmt_p(r.P, prefix=''), fontsize=5.4, va='center', ha='left', color=INK2, transform=tr, clip_on=False)
        axa.text(1.42, yi, '%d/%d' % (r.n_cohorts_same_sign, r.n_cohorts_evaluable), fontsize=5.4, va='center', ha='center', color=INK2, transform=tr, clip_on=False)
    axa.text(1.04, len(cor) - 0.45, 'P', fontsize=5.8, ha='left', va='bottom', fontweight='bold', color=INK2, transform=tr, clip_on=False)
    axa.text(1.42, len(cor) - 0.45, 'Same\ndirection', fontsize=5.4, ha='center', va='bottom', fontweight='bold', color=INK2, transform=tr, clip_on=False, linespacing=1.0)
    axa.set_title('Bulk AML transcriptomes, 8 cohorts (n = %s)' % format(int(cor.n.iloc[0]), ','), fontsize=6.7, fontweight='bold', loc='left', pad=14)
    fig.text(0.05, 0.522, 'Gray symbols: lineages outside the nine-axis RNHI space (gray labels) or with a 95% CI that includes 0.\nSame direction: cohorts in which the within-cohort coefficient has the sign of the pooled estimate.', fontsize=5.5, color=INK2, ha='left', va='top', linespacing=1.4)
    return axa
