from __future__ import annotations
import importlib.util
import re
import sys
from pathlib import Path
CODE = Path(__file__).resolve().parent.parent
COMMON = CODE / 'common'
STATS = CODE / 'data' / 'statistics'
FIGS = CODE / 'output'

def figstyle():
    key = 'figstyle'
    if key in sys.modules:
        return sys.modules[key]
    spec = importlib.util.spec_from_file_location(key, COMMON / 'figstyle.py')
    mod = importlib.util.module_from_spec(spec)
    sys.modules[key] = mod
    spec.loader.exec_module(mod)
    return mod

def stat(name: str) -> Path:
    p = STATS / name
    if not p.exists():
        raise FileNotFoundError(f'Missing statistics input: {name}')
    return p

def next_version(folder: Path, stem: str) -> int:
    pat = re.compile('^' + re.escape(stem) + '_v(\\d+)\\.(pdf|svg|png)$')
    vs = []
    for f in Path(folder).iterdir():
        m = pat.match(f.name)
        if m:
            vs.append(int(m.group(1)))
    return max(vs) + 1 if vs else 1

def save_crop(fig, panel_dir: Path, stem: str, bbox, dpi=600, svg=True):
    from matplotlib.transforms import Bbox
    w_in, h_in = fig.get_size_inches()
    x0, y0, x1, y1 = bbox
    box = Bbox.from_extents(x0 * w_in, y0 * h_in, x1 * w_in, y1 * h_in)
    out = Path(panel_dir) / 'output'
    out.mkdir(parents=True, exist_ok=True)
    v = next_version(out, stem)
    paths = []
    for ext in ('pdf', 'svg', 'png') if svg else ('pdf', 'png'):
        p = out / f'{stem}_v{v}.{ext}'
        fig.savefig(p, dpi=dpi, bbox_inches=box)
        paths.append(p)
    return paths

def save_full(fig, outdir: Path, stem: str, dpi=600, svg=True):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    v = next_version(outdir, stem)
    paths = []
    for ext in ('pdf', 'svg', 'png') if svg else ('pdf', 'png'):
        p = outdir / f'{stem}_v{v}.{ext}'
        fig.savefig(p, dpi=dpi)
        paths.append(p)
    return paths

def anonymous_ids(values, label='Identifier'):
    result = []
    for value in values:
        token = str(value).strip()
        if re.fullmatch('[1-9][0-9]*(?:\\.0+)?', token) is None:
            raise ValueError(f'{label} values must be positive integer anonymous IDs')
        result.append(str(int(token.split('.')[0])))
    return result
