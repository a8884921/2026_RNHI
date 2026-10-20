import os
import sys
from pathlib import Path
os.environ.setdefault('RNHI_FONT', 'Arial')
HERE = Path(__file__).resolve().parent
PANEL_DIR = HERE.parent
FIGURE = PANEL_DIR.parent.name
LETTER = PANEL_DIR.name.split('_')[1]
CODE = PANEL_DIR.parents[1]
sys.path.insert(0, str(CODE / 'common'))
sys.path.insert(0, str(HERE))
import panel_io
import panel as P
import matplotlib.pyplot as plt
fs = panel_io.figstyle()
fig = plt.figure(figsize=fs.A4)
P.draw(fig)
paths = panel_io.save_crop(fig, PANEL_DIR, '%s_Panel_%s' % (FIGURE, LETTER), P.BBOX)
print('saved:', *[p.name for p in paths], sep='\n  ')
