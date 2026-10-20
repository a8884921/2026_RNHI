import json
import sys
from pathlib import Path
import pandas as pd
HERE = Path(__file__).resolve().parent
CODE = HERE.parents[2]
sys.path.insert(0, str(CODE / 'common'))
import panel_io as io
_fs = io.figstyle()
globals().update({k: v for k, v in vars(_fs).items() if not k.startswith('_')})
STATS = io.STATS
FIGS = HERE.parent / 'output'
FIGS.mkdir(parents=True, exist_ok=True)
F = str(STATS) + '/'
per = pd.read_csv(F + 'per_cohort_survival.csv')
per['cohort'] = io.anonymous_ids(per['cohort'], 'Cohort')
R = json.load(open(F + 'pooled_results.json'))
grp = pd.read_csv(F + 'robustness_groupings.csv')
dfn = pd.read_csv(F + 'robustness_definitions.csv')
lo = pd.read_csv(F + 'robustness_loco.csv')
lo['left_out'] = io.anonymous_ids(lo['left_out'], 'Cohort')
fig = plt.figure(figsize=A4)
figure_title(fig, 'Supplementary Figure S3')
panel_letter(fig, 0.02, 0.952, 'A')
fig.text(0.065, 0.94, 'Continuous RNHI (per 1-SD increase within each cohort): no grouping threshold', fontsize=6.8, fontweight='bold')
cols = [(0.0, 'Cohort', 'left'), (0.22, 'Patients', 'right'), (0.3, 'Deaths', 'right'), (0.585, 'HR per SD (95% CI)', 'right'), (0.665, 'P', 'right')]
rows = [{'cells': [r.cohort, f'{int(r.n):,}', f'{int(r.events):,}', fmt_hr(r.sd_HR, r.sd_lo, r.sd_hi), fmt_p(r.sd_P, prefix='')], 'hr': r.sd_HR, 'lo': r.sd_lo, 'hi': r.sd_hi} for r in per.itertuples()]
ms = R['meta_sd']
pooled = {'cells': ['Pooled (random effects)', f"{R['n_total']:,}", f"{R['events_total']:,}", fmt_hr(ms['HR'], ms['lo'], ms['hi']), fmt_p(ms['P'], prefix='')], 'hr': ms['HR'], 'lo': ms['lo'], 'hi': ms['hi']}
axT, axF = forest_table(fig, [0.065, 0.735, 0.89, 0.19], rows, cols, (0.7, 1.0), xlim=(0.5, 1.5), ticks=(0.5, 0.75, 1.0, 1.5), ticklabels=['0.5', '0.75', '1.0', '1.5'], color=HIGH, pooled=pooled, xlabel='Hazard ratio per 1-SD higher RNHI')
axT.text(0.0, -1.75, f"HR < 1 in {sum(per.sd_HR < 1)}/{len(per)} cohorts; I$^2$ = {ms['I2']:.0f}%.", fontsize=5.9, color=INK2, va='center')
panel_letter(fig, 0.02, 0.69, 'B')
fig.text(0.065, 0.678, f"Alternative ways of grouping patients (cohort-stratified Cox models; n = {R['n_total']:,} unless stated)", fontsize=6.8, fontweight='bold')
order = ['Above vs below cohort median', 'Highest vs lowest quartile', 'Per quartile increase (trend)', 'Per tertile increase (trend)', 'Per 25-percentile increase (continuous rank)', 'Per 1-SD increase (continuous)']
g = grp.set_index('scheme').loc[order].reset_index()
colsb = [(0.0, 'Contrast', 'left'), (0.4, 'Patients', 'right'), (0.585, 'HR (95% CI)', 'right'), (0.665, 'P', 'right')]
rowsb = [{'cells': [r.scheme, f'{int(r.n):,}', fmt_hr(r.HR, r.lo, r.hi), fmt_p(r.P, prefix='')], 'hr': r.HR, 'lo': r.lo, 'hi': r.hi} for r in g.itertuples()]
sp = json.load(open(F + 'robustness_sparsity.json'))['incl_single']
rowsb.append({'cells': ['Median split, incl. single-lineage compositions', f"{int(sp['n']):,}", fmt_hr(sp['HR'], sp['lo'], sp['hi']), fmt_p(sp['P'], prefix='')], 'hr': sp['HR'], 'lo': sp['lo'], 'hi': sp['hi']})
s3 = json.load(open(F + 'robustness_sparsity.json'))['min3_lineages']
rowsb.append({'cells': ['Median split, samples with ≥3 detectable lineages', f"{int(s3['n']):,}", fmt_hr(s3['HR'], s3['lo'], s3['hi']), fmt_p(s3['P'], prefix='')], 'hr': s3['HR'], 'lo': s3['lo'], 'hi': s3['hi']})
forest_table(fig, [0.065, 0.528, 0.89, 0.14], rowsb, colsb, (0.7, 1.0), xlim=(0.6, 1.2), ticks=(0.6, 0.8, 1.0, 1.2), ticklabels=['0.6', '0.8', '1.0', '1.2'], color=HIGH, xlabel='Hazard ratio (higher vs lower RNHI)')
panel_letter(fig, 0.02, 0.5, 'C')
fig.text(0.065, 0.488, 'Alternative RNHI definitions (median split within cohort; cohort-stratified Cox models)', fontsize=6.8, fontweight='bold')
colsc = [(0.0, 'Definition', 'left'), (0.44, 'Spearman ρ with primary', 'right'), (0.585, 'HR (95% CI)', 'right'), (0.665, 'P', 'right')]
rowsc = []
for r in dfn.itertuples():
    lab = r.variant.replace('centre', 'center').replace('Centre', 'Center').replace('Primary definition (9 lineages, 45-donor reference center)', 'Primary definition (9 lineages; 45-donor reference center)')
    rowsc.append({'cells': [lab, f'{r.rho_with_primary:.2f}', fmt_hr(r.HR, r.lo, r.hi), fmt_p(r.P, prefix='')], 'hr': r.HR, 'lo': r.lo, 'hi': r.hi, 'weight': 'bold' if r.variant_key == 'RNHI' else 'normal'})
forest_table(fig, [0.065, 0.24, 0.89, 0.235], rowsc, colsc, (0.7, 1.0), xlim=(0.6, 1.2), ticks=(0.6, 0.8, 1.0, 1.2), ticklabels=['0.6', '0.8', '1.0', '1.2'], color=HIGH, xlabel='Hazard ratio, RNHI-high vs RNHI-low', fs=6.0)
panel_letter(fig, 0.02, 0.198, 'D')
fig.text(0.065, 0.19, 'Leave-one-cohort-out analysis (median split; cohort-stratified Cox models)', fontsize=6.8, fontweight='bold')
colsd = [(0.0, 'Cohort omitted', 'left'), (0.4, 'Patients', 'right'), (0.585, 'HR (95% CI)', 'right'), (0.665, 'P', 'right')]
rowsd = [{'cells': [r.left_out, f'{int(r.n):,}', fmt_hr(r.HR, r.lo, r.hi), fmt_p(r.P, prefix='')], 'hr': r.HR, 'lo': r.lo, 'hi': r.hi} for r in lo.itertuples()]
forest_table(fig, [0.065, 0.04, 0.89, 0.137], rowsd, colsd, (0.7, 1.0), xlim=(0.6, 1.2), ticks=(0.6, 0.8, 1.0, 1.2), ticklabels=['0.6', '0.8', '1.0', '1.2'], color=HIGH, xlabel='Hazard ratio, RNHI-high vs RNHI-low', fs=6.0)
fig.savefig(str(FIGS) + '/Supplementary_Figure_S3.pdf')
fig.savefig(str(FIGS) + '/Supplementary_Figure_S3.png', dpi=170)
print('ok')
