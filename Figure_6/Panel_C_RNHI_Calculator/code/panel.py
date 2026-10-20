import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
CODE = HERE.parents[2]
sys.path.insert(0, str(CODE / 'common'))
import panel_io as io
_fs = io.figstyle()
globals().update({k: v for k, v in vars(_fs).items() if not k.startswith('_')})
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
BBOX = (0.012, 0.238, 0.968, 0.462)

def draw(fig):
    panel_letter(fig, 0.02, 0.455, 'C')
    axe = fig.add_axes([0.05, 0.255, 0.9, 0.185])
    axe.set_xlim(-1.5, 101.5)
    axe.set_ylim(0, 40)
    axe.axis('off')
    axe.text(0, 39, 'RNHI calculator (web application): one fixed definition applied across data types', fontsize=6.8, fontweight='bold', ha='left', va='top')

    def box(x, y, w, h, title, lines, fc='#F4F6F8', ec='#9AA5B1'):
        axe.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.3,rounding_size=1.2', fc=fc, ec=ec, lw=0.7))
        axe.text(x + w / 2, y + h - 2.0, title, ha='center', va='top', fontsize=6.4, fontweight='bold', color=INK)
        for i, l in enumerate(lines):
            axe.text(x + w / 2, y + h - 7.0 - i * 3.0, l, ha='center', va='top', fontsize=5.8, color=INK2)

    def arrow(x0, x1, y):
        axe.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle='-|>', mutation_scale=8, lw=0.8, color=INK2))
    box(0, 9, 17, 23, '1  Input', ['Bulk RNA-seq or', 'microarray', 'scRNA-seq', '(pseudobulk)', 'Visium spots'])
    arrow(17.8, 21.2, 20.5)
    box(22, 9, 17, 23, '2  Composition', ['Deconvolution with', 'the joint healthy/AML', 'reference;', 'nine measurable', 'normal lineages'])
    arrow(39.8, 43.2, 20.5)
    box(44, 9, 17, 23, '3  Log-ratio space', ['Closure and', 'zero replacement;', 'centered log-ratio', 'transformation'])
    arrow(61.8, 65.2, 20.5)
    box(66, 9, 15, 23, '4  Score', ['Distance to the', 'healthy-donor center;', 'RNHI = −distance', '(0 = healthy-like)'], fc='#FCEFE6', ec=HIGH)
    arrow(81.8, 85.2, 20.5)
    box(86, 9, 14, 23, '5  Output', ['Per-sample RNHI;', 'lineage profile', 'versus healthy;', 'tables, figures'])
    return axe
