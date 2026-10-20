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
HERE = _Path(__file__).resolve().parent
STATS = _project_io.STATS
FIGS = _PROJECT_ROOT / 'output' / 'figure_s2'
FIGS.mkdir(parents=True, exist_ok=True)
_spec = _ilu.spec_from_file_location('figstyle', _PROJECT_ROOT / 'common' / 'figstyle.py')
_fs = _ilu.module_from_spec(_spec)
_sys.modules['figstyle'] = _fs
_spec.loader.exec_module(_fs)
from figstyle import *
import pandas as pd
t = pd.read_csv(STATS / 'citeseq_patient_table.csv', usecols=['frac_CD235a_all', 'erythroid_RNA'])
R = json.load(open(STATS / 'citeseq_concordance.json'))
ERY = '#A50F15'
fig = plt.figure(figsize=A4)
figure_title(fig, 'Supplementary Figure S2')
fig.text(0.065, 0.94, f"CITE-seq cohort: bone marrow at diagnosis (n = {R['n']})", fontsize=6.8, fontweight='bold')
fig.text(0.065, 0.927, 'Erythroid fraction deconvolved from RNA vs CD235a$^{+}$ fraction measured at the protein level', fontsize=6.3, color=INK2)
ax = fig.add_axes([0.16, 0.655, 0.363, 0.257])
x = (100 * t.frac_CD235a_all).to_numpy(float)
y = (100 * t.erythroid_RNA).to_numpy(float)
b1, b0 = np.polyfit(x, y, 1)
n = len(x)
xs = np.linspace(x.min(), x.max(), 200)
resid = y - (b0 + b1 * x)
s = np.sqrt((resid ** 2).sum() / (n - 2))
se = s * np.sqrt(1 / n + (xs - x.mean()) ** 2 / ((x - x.mean()) ** 2).sum())
from scipy import stats as _st
tq = _st.t.ppf(0.975, n - 2)
ax.fill_between(xs, b0 + b1 * xs - tq * se, b0 + b1 * xs + tq * se, color=ERY, alpha=0.12, linewidth=0, zorder=1)
ax.plot(xs, b0 + b1 * xs, color=ERY, lw=1.0, zorder=2)
ax.scatter(x, y, s=22, color=ERY, edgecolor='white', linewidth=0.5, alpha=0.95, zorder=3)
ax.set_xlim(-2, 102)
ax.set_ylim(-1.6, 81.6)
ax.set_xticks([0, 25, 50, 75, 100])
ax.set_yticks([0, 20, 40, 60, 80])
ax.set_xlabel('CD235a$^{+}$ cells, protein\n(% of all cells)')
ax.set_ylabel('Erythroid cells, RNA deconvolution\n(% of all cells)')
ax.text(0.04, 0.97, f"Spearman ρ = {R['spearman_rho']:.2f}\n" + fmt_p(R['P']), transform=ax.transAxes, ha='left', va='top', fontsize=6.0, linespacing=1.45)
fig.savefig(str(FIGS) + '/Supplementary_Figure_S2.pdf')
fig.savefig(str(FIGS) + '/Supplementary_Figure_S2.png', dpi=170)
print('saved Supplementary_Figure_S2')
