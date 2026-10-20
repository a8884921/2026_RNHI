import json
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
CODE = HERE.parents[2]
sys.path.insert(0, str(CODE / 'common'))
import panel_io as io
_fs = io.figstyle()
globals().update({k: v for k, v in vars(_fs).items() if not k.startswith('_')})
BBOX = (0.503, 0.548, 0.97, 0.766)
R = json.load(open(io.stat('pooled_results.json')))

def draw(fig):
    panel_letter(fig, 0.52, 0.752, 'C')
    fig.text(0.555, 0.738, 'Integrity, not quantity, carries survival information', fontsize=6.7, fontweight='bold', ha='left', va='center')
    fig.text(0.555, 0.725, f"Cohort-stratified Cox models (n = {R['n_total']:,}); high vs low = above vs below cohort median", fontsize=5.6, ha='left', va='center', color=INK2)
    colsc = [(0.0, 'Variable', 'left'), (0.5, 'HR (95% CI)', 'right'), (0.635, 'P', 'right')]

    def crow(lab, v, rnhi):
        return {'cells': [lab, fmt_hr(v['HR'], v['lo'], v['hi']), fmt_p(v['P'], prefix='')], 'hr': v['HR'], 'lo': v['lo'], 'hi': v['hi'], 'text_color': HIGH if rnhi else '#555555', 'weight': 'bold' if rnhi else 'normal'}
    crows = [{'cells': ['Separate models', '', ''], 'hr': None, 'text_color': INK2}, crow('RNHI (integrity)', R['strat_median'], True), crow('Normal quantity', R['nq_uni'], False), {'cells': ['Mutually adjusted', '', ''], 'hr': None, 'text_color': INK2}, crow('RNHI (integrity)', R['joint_rnhi'], True), crow('Normal quantity', R['joint_nq'], False)]
    rc = [INK2, HIGH, '#808080', INK2, HIGH, '#808080']
    axTc, axFc = forest_table(fig, [0.555, 0.575, 0.4, 0.125], crows, colsc, (0.68, 1.0), xlim=(0.6, 1.25), ticks=(0.6, 0.8, 1.0, 1.2), ticklabels=['0.6', '0.8', '1.0', '1.2'], row_colors=rc, xlabel='Hazard ratio, high vs low', stripe=False)
    for t in axTc.texts:
        if t.get_text() in ('Separate models', 'Mutually adjusted'):
            t.set_style('italic')
    return (axTc, axFc)
