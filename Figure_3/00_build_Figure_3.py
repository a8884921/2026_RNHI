import importlib.util
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
CODE = HERE.parent
sys.path.insert(0, str(CODE / 'common'))
import panel_io
import matplotlib.pyplot as plt
fs = panel_io.figstyle()
fig = plt.figure(figsize=fs.A4)
fs.figure_title(fig, 'Figure 3')
for panel_py in sorted(HERE.glob('Panel_*/code/panel.py')):
    name = panel_py.parents[1].name
    spec = importlib.util.spec_from_file_location(f'panel_{name}', panel_py)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.draw(fig)
    print('drew', name)
paths = panel_io.save_full(fig, HERE / '00_output', 'Figure_3')
print('saved:', *[p.name for p in paths], sep='\n  ')
