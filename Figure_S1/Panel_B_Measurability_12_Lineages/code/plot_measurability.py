from pathlib import Path as _ProjectPath
import sys as _project_sys
_PROJECT_ROOT = next((p for p in _ProjectPath(__file__).resolve().parents if (p / 'common' / 'panel_io.py').is_file()), None)
if _PROJECT_ROOT is None:
    raise RuntimeError('The common analysis module is required.')
_project_sys.path.insert(0, str(_PROJECT_ROOT / 'common'))
import panel_io as _project_io
import importlib.util as _ilu
import sys as _sys
from pathlib import Path as _Path
HERE = _Path(__file__).resolve().parent
STATS = _project_io.STATS
FIGS = _PROJECT_ROOT / 'output' / 'figure_s1'
FIGS.mkdir(parents=True, exist_ok=True)
_spec = _ilu.spec_from_file_location('figstyle', _PROJECT_ROOT / 'common' / 'figstyle.py')
_fs = _ilu.module_from_spec(_spec)
_sys.modules['figstyle'] = _fs
_spec.loader.exec_module(_fs)
from figstyle import *
import pandas as pd
r = pd.read_csv(str(STATS) + '/lineage_measurability.csv')
order = ['HSC_MPP', 'LMPP', 'MyeProg', 'Mono', 'Gran', 'Mega', 'EarlyEry', 'LateEry', 'B', 'T', 'NK', 'DC'][::-1]
r = r.set_index('key').loc[order].reset_index()
fig = plt.figure(figsize=(8.27, 4.2))
ax = fig.add_axes([0.2, 0.17, 0.7, 0.74])
y = np.arange(len(r))
cols = ['#3B5998' if k else '#E4572E' for k in r.retained]
ax.barh(y, r.pct_patients_detected, height=0.6, color=cols, lw=0)
ax.axvline(4, color=INK2, lw=0.7, ls=(0, (3, 2)))
ax.text(4.8, len(r) - 0.35, '4% threshold', fontsize=7, color=INK2, va='bottom')
for yi, v in zip(y, r.pct_patients_detected):
    ax.text(max(v + 1.0, 5.2), yi, f'{v:.1f}', va='center', fontsize=7.5, color=INK)
ax.set_yticks(y)
ax.set_yticklabels(r.lineage, fontsize=8.5)
ax.set_xlim(0, 105)
ax.set_xlabel('Patients with a detectable lineage estimate (%)', fontsize=9)
ax.tick_params(axis='x', labelsize=8)
ax.spines['left'].set_visible(False)
ax.tick_params(axis='y', length=0)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color='#3B5998', label='Retained in the RNHI space (9 lineages)'), Patch(color='#E4572E', label='Below 4% (3 lineages)')], loc='upper right', fontsize=8, frameon=False, bbox_to_anchor=(1.0, 0.98))
ax.set_title(f'Measurability of the 12 normal lineages in bulk AML transcriptomes (8 cohorts, n = {int(r.n_patients.iloc[0]):,})', fontsize=9, fontweight='bold', loc='left', pad=6)
fig.savefig(str(FIGS) + '/S1B_panel.pdf')
print('ok')
