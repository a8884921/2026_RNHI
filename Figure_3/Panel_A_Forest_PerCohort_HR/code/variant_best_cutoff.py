import json
import os
import sys
from pathlib import Path
os.environ.setdefault('RNHI_FONT', 'Arial')
HERE = Path(__file__).resolve().parent
PANEL_DIR = HERE.parent
CODE = PANEL_DIR.parents[1]
sys.path.insert(0, str(CODE / 'common'))
import panel_io as io
import pandas as pd
import matplotlib.pyplot as plt
_fs = io.figstyle()
globals().update({k: v for k, v in vars(_fs).items() if not k.startswith('_')})
COLS = [(0.0, 'Cohort', 'left'), (0.205, 'Patients', 'right'), (0.272, 'Deaths', 'right'), (0.392, 'High / Low', 'right'), (0.48, 'Cutoff', 'right'), (0.632, 'HR (95% CI)', 'right'), (0.705, 'P', 'right')]
RECT, FOREST = ([0.065, 0.775, 0.89, 0.17], (0.74, 1.0))
BBOX = (0.01, 0.748, 0.998, 0.975)

def ordinal(p):
    p = int(round(p))
    suf = 'th' if 10 <= p % 100 <= 20 else {1: 'st', 2: 'nd', 3: 'rd'}.get(p % 10, 'th')
    return '%d%s' % (p, suf)

def draw(fig):
    bc = pd.read_csv(io.stat('best_cutoff_by_cohort.csv')).sort_values('n', ascending=False)
    bc['cohort'] = io.anonymous_ids(bc['cohort'], 'Cohort')
    pool = json.load(open(io.stat('best_cutoff_pooled.json')))
    rows = [{'cells': [r.cohort, format(int(r.n), ','), format(int(r.events), ','), '%d / %d' % (r.n_high, r.n_low), ordinal(r.pct), fmt_hr(r.HR, r.lo, r.hi), '%.2e' % r.P], 'hr': r.HR, 'lo': r.lo, 'hi': r.hi} for r in bc.itertuples()]
    pooled = {'cells': ['Pooled (random effects)', format(int(pool['n_total']), ','), format(int(pool['events_total']), ','), '%s / %s' % (format(int(bc.n_high.sum()), ','), format(int(bc.n_low.sum()), ',')), '—', fmt_hr(pool['HR'], pool['lo'], pool['hi']), '%.2e' % pool['P']], 'hr': pool['HR'], 'lo': pool['lo'], 'hi': pool['hi']}
    axT, axF = forest_table(fig, RECT, rows, COLS, FOREST, xlim=(0.25, 2.0), ticks=(0.25, 0.5, 1, 2), ticklabels=['0.25', '0.5', '1', '2'], color=HIGH, pooled=pooled, xlabel='Hazard ratio, RNHI-high vs RNHI-low', favour=('Favors RNHI-high', 'Favors RNHI-low'))
    axT.text(0.0, -1.62, 'HR < 1 in %d/%d cohorts · I² = %.0f%% · Q = %.2f' % (int(pool['n_HRlt1']), int(pool['n_cohorts']), pool['I2'], pool['Q']), fontsize=5.9, ha='left', va='center', color=INK2)
    panel_letter(fig, 0.02, 0.958, 'A')
    return (axT, axF)
if __name__ == '__main__':
    fig = plt.figure(figsize=_fs.A4)
    draw(fig)
    paths = io.save_crop(fig, PANEL_DIR, 'Figure_3_Panel_A_BestCutoff', BBOX)
    print('saved:', *[p.name for p in paths], sep='\n  ')
