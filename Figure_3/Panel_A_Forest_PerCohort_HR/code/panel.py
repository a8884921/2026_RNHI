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
BBOX = (0.01, 0.748, 0.998, 0.975)
per = pd.read_csv(io.stat('per_cohort_survival.csv'))
per['cohort'] = io.anonymous_ids(per['cohort'], 'Cohort')
R = json.load(open(io.stat('pooled_results.json')))

def draw(fig):
    cols = [(0.0, 'Cohort', 'left'), (0.205, 'Patients', 'right'), (0.285, 'Deaths', 'right'), (0.395, 'High / Low', 'right'), (0.585, 'HR (95% CI)', 'right'), (0.665, 'P', 'right')]
    rows = [{'cells': [r.cohort, f'{int(r.n):,}', f'{int(r.events):,}', f'{int(r.n_high)} / {int(r.n_low)}', fmt_hr(r.HR, r.lo, r.hi), fmt_p(r.P_wald, prefix='')], 'hr': r.HR, 'lo': r.lo, 'hi': r.hi} for r in per.itertuples()]
    mm = R['meta_median']
    pooled = {'cells': ['Pooled (random effects)', f"{R['n_total']:,}", f"{R['events_total']:,}", '', fmt_hr(mm['HR'], mm['lo'], mm['hi']), fmt_p(mm['P'], prefix='')], 'hr': mm['HR'], 'lo': mm['lo'], 'hi': mm['hi']}
    axT, axF = forest_table(fig, [0.065, 0.775, 0.89, 0.17], rows, cols, (0.7, 1.0), xlim=(0.25, 2.0), ticks=(0.25, 0.5, 1, 2), ticklabels=['0.25', '0.5', '1', '2'], color=HIGH, pooled=pooled, xlabel='Hazard ratio, RNHI-high vs RNHI-low', favour=('Favors RNHI-high', 'Favors RNHI-low'))
    axT.text(0.0, -1.75, f"HR < 1 in {R['n_HRlt1']}/{R['n_cohorts']} cohorts; heterogeneity I$^2$ = {mm['I2']:.0f}%. RNHI-high/low = above/below each cohort's median RNHI.", fontsize=5.9, ha='left', va='center', color=INK2)
    panel_letter(fig, 0.02, 0.958, 'A')
    return (axT, axF)
