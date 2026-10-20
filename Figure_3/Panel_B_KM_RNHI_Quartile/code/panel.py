import json
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
CODE = HERE.parents[2]
sys.path.insert(0, str(CODE / 'common'))
import panel_io as io
import pandas as pd
_fs = io.figstyle()
globals().update({k: v for k, v in vars(_fs).items() if not k.startswith('_')})
BBOX = (0.012, 0.474, 0.478, 0.766)
km = pd.read_csv(io.stat('km_input.csv'))
R = json.load(open(io.stat('pooled_results.json')))

def draw(fig):
    panel_letter(fig, 0.02, 0.752, 'B')
    axb = fig.add_axes([0.11, 0.575, 0.34, 0.15])
    grp = [(km[km.qt == q].OS_time.values, km[km.qt == q].OS_event.values) for q in [1, 2, 3, 4]]
    draw_km(axb, grp, Q_COLORS, ['Q1 (lowest)', 'Q2', 'Q3', 'Q4 (highest)'], xmax=120, xticks=[0, 24, 48, 72, 96, 120], lw=1.15)
    axb.set_ylabel('Overall survival probability')
    axb.set_xlabel('Months from diagnosis', labelpad=1.5)
    axb.set_title(f"{R['n_cohorts']} cohorts pooled (n = {R['n_total']:,}) by within-cohort RNHI quartile", fontsize=6.7, fontweight='bold', loc='left', pad=3)
    axb.legend(loc='upper right', fontsize=5.9, handlelength=1.6, borderaxespad=0.1, title='RNHI quartile', title_fontsize=5.9, labelspacing=0.25)
    qt, qs = (R['quartile_trend'], R['quartiles'])
    axb.text(0.98, 0.635, 'HR vs Q1: ' + ', '.join((f"{q['Q']} {q['HR']:.2f}" for q in qs)) + '\nStratified log-rank ' + fmt_p(R['strat_logrank_quartile_P']) + '\nHR per quartile ' + f"{qt['HR']:.2f} ({qt['lo']:.2f}–{qt['hi']:.2f})" + '\nP for trend ' + fmt_p(qt['P'], prefix='= '), transform=axb.transAxes, fontsize=5.8, ha='right', va='top', color=INK, linespacing=1.4)
    risk_table(fig, axb, grp, Q_COLORS, ['Q1', 'Q2', 'Q3', 'Q4'], [0, 24, 48, 72, 96, 120], y0=0.523, dy=0.0108)
    return axb
