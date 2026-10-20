import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
CODE = HERE.parents[2]
sys.path.insert(0, str(CODE / 'common'))
import panel_io as io
import pandas as pd
_fs = io.figstyle()
globals().update({k: v for k, v in vars(_fs).items() if not k.startswith('_')})
GROUP = 'Intermediate'
LETTER = 'F'
X0 = 0.415
FIRST = False
BBOX = (0.327, 0.052, 0.66, 0.307)
kme = pd.read_csv(io.stat('km_eln_input.csv'))
eln = pd.read_csv(io.stat('cohort_1_eln_subgroups.csv'))

def draw(fig):
    panel_letter(fig, X0 - 0.08, 0.3, LETTER)
    ax = fig.add_axes([X0 + 0.005, 0.125, 0.225, 0.135])
    s = kme[kme.ELN == GROUP]
    grp = [(s[s.hi_med == 0].OS_time.values, s[s.hi_med == 0].OS_event.values), (s[s.hi_med == 1].OS_time.values, s[s.hi_med == 1].OS_event.values)]
    draw_km(ax, grp, [LOW, HIGH], ['RNHI-low', 'RNHI-high'], xmax=96, xticks=[0, 24, 48, 72, 96], lw=1.1)
    r = eln[eln.ELN == GROUP].iloc[0]
    ax.set_title('Cohort 1, ELN2017 ' + GROUP.lower() + ' (n = ' + str(int(r.n)) + ')', fontsize=6.5, fontweight='bold', loc='left', pad=3)
    ax.text(0.97, 0.97, 'HR ' + fmt_hr(r.HR, r.lo, r.hi) + '\nLog-rank ' + fmt_p(r.P_logrank), transform=ax.transAxes, fontsize=5.8, ha='right', va='top', linespacing=1.4)
    ax.set_xlabel('Months from diagnosis', labelpad=1.5)
    if FIRST:
        ax.set_ylabel('Overall survival probability')
        ax.legend(loc='lower left', fontsize=5.8, handlelength=1.4, borderaxespad=0.3)
    risk_table(fig, ax, grp, [LOW, HIGH], ['Low', 'High'], [0, 24, 48, 72, 96], y0=0.074, dy=0.0108, fontsize=5.6)
    return ax
