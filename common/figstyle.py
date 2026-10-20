import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import os
from pathlib import Path
_HERE = Path(__file__).resolve().parent

def resolve_fonts():
    requested = os.environ.get('RNHI_FONT', 'Arial').strip() or 'Arial'
    styles = {'regular': ('normal', 'normal'), 'bold': ('bold', 'normal'), 'italic': ('normal', 'italic'), 'bolditalic': ('bold', 'italic')}
    for family in dict.fromkeys([requested, 'Liberation Sans', 'DejaVu Sans']):
        files = {}
        try:
            for key, (weight, style) in styles.items():
                prop = font_manager.FontProperties(family=[family], weight=weight, style=style)
                files[key] = font_manager.findfont(prop, fallback_to_default=False)
        except ValueError:
            continue
        return (family, files)
    raise RuntimeError('No supported font family is available')
FONT, FONT_FILES = resolve_fonts()
for f in FONT_FILES.values():
    font_manager.fontManager.addfont(f)
plt.rcParams.update({'font.family': FONT, 'font.size': 7, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'mathtext.fontset': 'custom', 'mathtext.rm': FONT, 'mathtext.it': FONT + ':italic', 'mathtext.bf': FONT + ':bold', 'mathtext.default': 'regular', 'axes.linewidth': 0.6, 'axes.edgecolor': '#222222', 'axes.labelcolor': '#111111', 'axes.labelsize': 7, 'axes.titlesize': 7.5, 'axes.titleweight': 'bold', 'xtick.labelsize': 6.3, 'ytick.labelsize': 6.3, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.color': '#222222', 'ytick.color': '#222222', 'legend.fontsize': 6.3, 'legend.frameon': False, 'savefig.dpi': 600, 'axes.spines.top': False, 'axes.spines.right': False})
INK = '#111111'
INK2 = '#444444'
MUTED = '#777777'
GRID = '#E6E6E6'
LOW = '#0072B2'
HIGH = '#D55E00'
HEALTHY = '#009E73'
Q_COLORS = ['#0072B2', '#56B4E9', '#E69F00', '#D55E00']
A4 = (8.27, 11.69)

def fmt_p(p, prefix='P = '):
    if p is None or (isinstance(p, float) and math.isnan(p)):
        return prefix + 'NA'
    if p >= 0.1:
        return prefix + f'{p:.2f}'
    if p >= 0.001:
        return prefix + f'{p:.3f}'
    e = int(math.floor(math.log10(p)))
    m = p / 10 ** e
    m = round(m, 1)
    if m >= 10:
        m /= 10
        e += 1
    return prefix + f'${m:.1f}\\times10^{{{e}}}$'

def fmt_hr(hr, lo, hi):
    return f'{hr:.2f} ({lo:.2f}–{hi:.2f})'

def panel_letter(fig, x, y, letter):
    fig.text(x, y, letter, fontsize=18, fontweight='normal', va='top', ha='left', color=INK)

def figure_title(fig, text):
    fig.text(0.019, 0.995, text, fontsize=24, fontweight='bold', va='top', ha='left', color='#231815')

def km_curve(t, e):
    t = np.asarray(t, float)
    e = np.asarray(e, int)
    order = np.argsort(t)
    t = t[order]
    e = e[order]
    uniq = np.unique(t[e == 1])
    S = 1.0
    times = [0.0]
    surv = [1.0]
    for u in uniq:
        n_at = np.sum(t >= u)
        d = np.sum((t == u) & (e == 1))
        S *= 1 - d / n_at
        times.append(u)
        surv.append(S)
    times = np.array(times)
    surv = np.array(surv)
    ct = t[e == 0]
    cs = np.array([surv[np.searchsorted(times, c, side='right') - 1] for c in ct]) if len(ct) else np.array([])
    return (times, surv, ct, cs)

def at_risk(t, grid):
    t = np.asarray(t, float)
    return [int(np.sum(t >= g)) for g in grid]

def draw_km(ax, groups, colors, labels, xmax, xticks, lw=1.1, censor=True):
    for (t, e), c, lab in zip(groups, colors, labels):
        times, surv, ct, cs = km_curve(t, e)
        keep = times <= xmax
        tt = np.append(times[keep], xmax)
        ss = np.append(surv[keep], surv[keep][-1])
        ax.step(tt, ss, where='post', color=c, lw=lw, label=lab, solid_capstyle='butt')
        if censor and len(ct):
            m = ct <= xmax
            ax.plot(ct[m], cs[m], ls='none', marker='|', ms=2.6, mew=0.5, color=c)
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, 1.02)
    ax.set_xticks(xticks)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(['0', '0.25', '0.50', '0.75', '1.00'])

def risk_table(fig, ax, groups, colors, labels, grid, y0, dy, xlabel_offset=0.0, fontsize=5.8, title='No. at risk'):
    bb = ax.get_position()
    x0, x1 = ax.get_xlim()
    fig.text(bb.x0 - 0.028 + xlabel_offset, y0 + dy * 0.9, title, fontsize=fontsize, ha='right', va='center', color=INK2, style='italic')
    for i, ((t, e), c, lab) in enumerate(zip(groups, colors, labels)):
        yy = y0 - i * dy
        fig.text(bb.x0 - 0.028 + xlabel_offset, yy, lab, fontsize=fontsize, ha='right', va='center', color=c, fontweight='bold')
        for g, n in zip(grid, at_risk(t, grid)):
            xf = bb.x0 + (g - x0) / (x1 - x0) * bb.width
            fig.text(xf, yy, str(n), fontsize=fontsize, ha='center', va='center', color=INK)

def forest_table(fig, rect, rows, columns, forest_frac, xlim, ticks, ticklabels=None, color=None, row_colors=None, pooled=None, xlabel='', fs=6.3, hfs=6.5, stripe=True, ms=4.2, favour=None):
    import matplotlib.patches as mpatches
    n = len(rows)
    ybot = -1.9 if pooled else -0.7
    ylim = (ybot, n + 1.25)
    axT = fig.add_axes(rect)
    axT.set_xlim(0, 1)
    axT.set_ylim(*ylim)
    axT.axis('off')
    fx0, fx1 = forest_frac
    axF = fig.add_axes([rect[0] + fx0 * rect[2], rect[1], (fx1 - fx0) * rect[2], rect[3]], sharey=None)
    axF.set_ylim(*ylim)
    axF.set_xscale('log')
    axF.set_xlim(*xlim)
    axF.minorticks_off()
    axF.set_xticks(ticks)
    axF.set_xticklabels(ticklabels if ticklabels else [f'{t:g}' for t in ticks])
    axF.spines['left'].set_visible(False)
    axF.set_yticks([])
    axF.patch.set_alpha(0)
    axF.axvline(1, color=MUTED, lw=0.6, ls=(0, (3, 2)), zorder=0, ymin=0, ymax=(n + 0.5 - ylim[0]) / (ylim[1] - ylim[0]))
    yh = n + 0.97
    for x, h, ha in columns:
        axT.text(x, yh, h, fontsize=hfs, fontweight='bold', ha=ha, va='center', color=INK)
    axT.plot([0, 1], [n + 0.62, n + 0.62], lw=0.5, color=INK2, clip_on=False)
    for i, r in enumerate(rows):
        y = n - i
        if stripe and i % 2 == 0:
            axT.add_patch(mpatches.Rectangle((0, y - 0.5), 1, 1, color='#F3F4F5', zorder=-1, lw=0))
        for (x, h, ha), v in zip(columns, r['cells']):
            axT.text(x, y, v, fontsize=fs, ha=ha, va='center', color=r.get('text_color', INK), fontweight=r.get('weight', 'normal'))
        c = row_colors[i] if row_colors else color
        if r.get('hr') is not None:
            lo, hi = (max(r['lo'], xlim[0]), min(r['hi'], xlim[1]))
            axF.plot([lo, hi], [y, y], color=c, lw=1.0, solid_capstyle='butt', zorder=2)
            for end, inner in ((xlim[1], xlim[1] / 1.12), (xlim[0], xlim[0] * 1.12)):
                if end == xlim[1] and r['hi'] > xlim[1] or (end == xlim[0] and r['lo'] < xlim[0]):
                    axF.annotate('', xy=(end, y), xytext=(inner, y), annotation_clip=False, arrowprops=dict(arrowstyle='-|>', color=c, lw=1.0, mutation_scale=6, shrinkA=0, shrinkB=0))
            axF.plot(r['hr'], y, marker='s', ms=ms, color=c, mec='white', mew=0.5, zorder=3)
    if pooled is not None:
        yp = -1.0
        axT.plot([0, 1], [-0.45, -0.45], lw=0.4, color=MUTED, clip_on=False)
        for (x, h, ha), v in zip(columns, pooled['cells']):
            axT.text(x, yp, v, fontsize=fs, ha=ha, va='center', color=INK, fontweight='bold')
        axF.fill([pooled['lo'], pooled['hr'], pooled['hi'], pooled['hr']], [yp, yp + 0.38, yp, yp - 0.38], color=INK, zorder=3)
    if favour:
        axF.text(0.02, 1.0, favour[0], transform=axF.transAxes, fontsize=5.5, color=MUTED, ha='left', va='bottom', style='italic')
        axF.text(0.98, 1.0, favour[1], transform=axF.transAxes, fontsize=5.5, color=MUTED, ha='right', va='bottom', style='italic')
    axF.set_xlabel(xlabel, fontsize=6.4, labelpad=2)
    return (axT, axF)
