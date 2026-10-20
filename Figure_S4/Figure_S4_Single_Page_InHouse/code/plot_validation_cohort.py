from pathlib import Path as _ProjectPath
import sys as _project_sys
_PROJECT_ROOT = next((p for p in _ProjectPath(__file__).resolve().parents if (p / 'common' / 'panel_io.py').is_file()), None)
if _PROJECT_ROOT is None:
    raise RuntimeError('The common analysis module is required.')
_project_sys.path.insert(0, str(_PROJECT_ROOT / 'common'))
import panel_io as _project_io
import json
import importlib.util as _ilu
import sys as _sys
from pathlib import Path as _Path
import numpy as np
import pandas as pd
HERE = _Path(__file__).resolve().parent
ROOT = _PROJECT_ROOT
STATS = _project_io.STATS
FIGS = _PROJECT_ROOT / 'output' / 'figure_s4'
FIGS.mkdir(parents=True, exist_ok=True)
IN = _PROJECT_ROOT / 'data' / 'figure_s4'
_spec = _ilu.spec_from_file_location('figstyle', _PROJECT_ROOT / 'common' / 'figstyle.py')
_fs = _ilu.module_from_spec(_spec)
_sys.modules['figstyle'] = _fs
_spec.loader.exec_module(_fs)
from figstyle import *
_spec = _ilu.spec_from_file_location('rm', _PROJECT_ROOT / 'common' / 'rnhi_math.py')
rm = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(rm)
EPS = 1e-06
hd = pd.read_csv(_PROJECT_ROOT / 'data' / 'reference' / 'healthy_donor_counts.csv', usecols=['sample_id'] + list(rm.LINE9), dtype={'sample_id': str}).set_index('sample_id')
center = rm.comp_center(hd, method='eq_donor_mean', lineage=rm.LINE9)
comp = pd.read_csv(IN / 'validation_composition.csv', usecols=['sample_id'] + list(rm.LINE9), dtype={'sample_id': str}).set_index('sample_id')
clin = pd.read_csv(IN / 'validation_clinical_rnhi.csv', usecols=['sample_id', 'RNHI', 'nz'], dtype={'sample_id': str}).set_index('sample_id')
for _table in (hd, comp, clin):
    _project_io.anonymous_ids(_table.index, 'Sample')
    if not _table.index.is_unique:
        raise ValueError('Sample identifiers must be unique within each table.')
if set(comp.index) != set(clin.index):
    raise ValueError('Composition and clinical sample identifiers must match.')

def rnhi_const(x, eps=EPS):
    x = np.asarray(x, float)
    x = x / x.sum()
    out = np.where(x == 0, eps, x)
    out = out / out.sum()
    lx = np.log(out)
    return -float(np.sqrt(((lx - lx.mean() - center) ** 2).sum()))
rec = pd.Series({s: rnhi_const(comp.loc[s, rm.LINE9].to_numpy()) for s in comp.index})
dmax = float((rec - clin.loc[rec.index, 'RNHI']).abs().max())
assert dmax < 1e-08, f'RNHI not reproduced (max abs difference {dmax:g})'
assert ((comp[rm.LINE9] > 0).sum(axis=1) == clin.loc[comp.index, 'nz']).all()
print(f'RNHI of the 22 patients reproduced from the compositions (max abs difference {dmax:.1e})')
R = json.load(open(STATS / 'validation_survival.json'))
g = pd.read_csv(STATS / 'validation_groups.csv', usecols=['group', 'OS_time', 'OS_event'])
lo_ = g[g.group == 'RNHI-low']
hi_ = g[g.group == 'RNHI-high']
groups = [(lo_.OS_time.values, lo_.OS_event.values), (hi_.OS_time.values, hi_.OS_event.values)]
fig = plt.figure(figsize=A4)
figure_title(fig, 'Supplementary Figure S4')
fig.text(0.065, 0.94, f"Validation cohort: peripheral blood at diagnosis (n = {R['n']}; {R['deaths']} deaths)", fontsize=6.8, fontweight='bold')
fig.text(0.065, 0.927, 'RNHI-high vs RNHI-low by the outcome-derived RNHI threshold', fontsize=6.3, color=INK2)
ax = fig.add_axes([0.16, 0.7, 0.46, 0.205])
xt = [0, 6, 12, 18, 24, 30]
draw_km(ax, groups, [LOW, HIGH], [f"RNHI-low (n = {R['n_low']}; {R['deaths_low']} deaths)", f"RNHI-high (n = {R['n_high']}; {R['deaths_high']} deaths)"], xmax=30, xticks=xt, lw=1.2)
ax.set_xlabel('Months since diagnosis')
ax.set_ylabel('Overall survival')
ax.legend(loc='lower left', fontsize=6.0, handlelength=1.6, borderaxespad=0.4)
ax.text(0.97, 0.97, 'HR (high vs low) ' + fmt_hr(R['HR'], R['lo'], R['hi']) + '\nCox ' + fmt_p(R['P_cox']) + '\nLog-rank ' + fmt_p(R['P_logrank']), transform=ax.transAxes, fontsize=6.0, ha='right', va='top', linespacing=1.45)
risk_table(fig, ax, groups, [LOW, HIGH], ['Low', 'High'], xt, y0=0.645, dy=0.0115, fontsize=5.9)
fig.savefig(str(FIGS) + '/Supplementary_Figure_S4.pdf')
fig.savefig(str(FIGS) + '/Supplementary_Figure_S4.png', dpi=170)
print('ok')
