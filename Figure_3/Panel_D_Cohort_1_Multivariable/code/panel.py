import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
CODE = HERE.parents[2]
sys.path.insert(0, str(CODE / 'common'))
import panel_io as io
import pandas as pd
_fs = io.figstyle()
globals().update({k: v for k, v in vars(_fs).items() if not k.startswith('_')})
BBOX = (0.012, 0.312, 0.998, 0.482)
mv = pd.read_csv(io.stat('cohort_1_multivariable.csv'))

def draw(fig):
    panel_letter(fig, 0.02, 0.475, 'D')
    fig.text(0.065, 0.462, 'Cohort 1: RNHI-high vs RNHI-low after adjustment for clinical and genetic risk factors', fontsize=6.7, fontweight='bold', ha='left', va='center')
    adj = ['None', 'ELN2017', 'ELN2017 + blast%', 'ELN2017 + blast% + age', 'Model 4 + $\\it{NPM1}$ + $\\it{FLT3}$-ITD', 'Model 5 + $\\it{TP53}$']
    colsd = [(0.0, 'Model', 'left'), (0.06, 'Adjusted for', 'left'), (0.385, 'Patients', 'right'), (0.455, 'Deaths', 'right'), (0.585, 'HR (95% CI)', 'right'), (0.665, 'P', 'right')]
    drows = [{'cells': [str(i + 1), adj[i], f'{int(r.N)}', f'{int(r.events)}', fmt_hr(r.HR, r.lo, r.hi), fmt_p(r.P, prefix='')], 'hr': r.HR, 'lo': r.lo, 'hi': r.hi} for i, r in enumerate(mv.itertuples())]
    axT, axF = forest_table(fig, [0.065, 0.335, 0.89, 0.115], drows, colsd, (0.7, 1.0), xlim=(0.4, 1.25), ticks=(0.4, 0.6, 0.8, 1.0, 1.2), ticklabels=['0.4', '0.6', '0.8', '1.0', '1.2'], color=HIGH, xlabel='Hazard ratio, RNHI-high vs RNHI-low')
    return (axT, axF)
