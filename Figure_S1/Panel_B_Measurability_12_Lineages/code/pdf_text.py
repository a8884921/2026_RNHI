from pathlib import Path as _ProjectPath
import sys as _project_sys
_PROJECT_ROOT = next((p for p in _ProjectPath(__file__).resolve().parents if (p / 'common' / 'panel_io.py').is_file()), None)
if _PROJECT_ROOT is None:
    raise RuntimeError('The common analysis module is required.')
_project_sys.path.insert(0, str(_PROJECT_ROOT / 'common'))
import panel_io as _project_io
import math
import re
import pymupdf as fitz
import importlib.util as _ilu
import sys as _sys
from pathlib import Path as _Path
_HERE = _Path(__file__).resolve().parent
if 'figstyle' not in _sys.modules:
    _spec = _ilu.spec_from_file_location('figstyle', _PROJECT_ROOT / 'common' / 'figstyle.py')
    _fs = _ilu.module_from_spec(_spec)
    _sys.modules['figstyle'] = _fs
    _spec.loader.exec_module(_fs)
_FILES = _sys.modules['figstyle'].FONT_FILES
REG, BOLD, ITAL = (_FILES['regular'], _FILES['bold'], _FILES['italic'])

def _rgb(c):
    return ((c >> 16 & 255) / 255.0, (c >> 8 & 255) / 255.0, (c & 255) / 255.0)

def spans(page):
    out = []
    d = page.get_text('dict')
    for b in d['blocks']:
        if b.get('type') != 0:
            continue
        for l in b['lines']:
            for s in l['spans']:
                s = dict(s)
                s['dir'] = l['dir']
                out.append(s)
    return out

def find_spans(page, needle, exact=False, region=None):
    res = []
    for s in spans(page):
        t = s['text']
        ok = t.strip() == needle if exact else needle in t
        if ok and (region is None or fitz.Rect(s['bbox']).intersects(fitz.Rect(region))):
            res.append(s)
    return res

def replace_span_text(page, needle, new_text_fn, exact=False, region=None, bold=None, size_scale=1.0, dy=0.0, dx=0.0, color=None, pad=0.6):
    targets = find_spans(page, needle, exact=exact, region=region)
    jobs = []
    for s in targets:
        r = fitz.Rect(s['bbox'])
        r.x0 -= pad
        r.y0 -= pad * 0.3
        r.x1 += pad
        r.y1 += pad * 0.3
        page.add_redact_annot(r, fill=False)
        new = new_text_fn(s['text']) if callable(new_text_fn) else new_text_fn
        isbold = 'Bold' in s['font'] if bold is None else bold
        jobs.append((s, new, isbold))
    if targets:
        page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE, text=fitz.PDF_REDACT_TEXT_REMOVE)
    for s, new, isbold in jobs:
        write(page, s['origin'], new, s['size'] * size_scale, color if color else _rgb(s['color']), isbold, s['dir'], dx=dx, dy=dy)
    return len(targets)

def write(page, origin, text, size, color, bold=False, direction=(1.0, 0.0), dx=0.0, dy=0.0, italic=False, tag='A'):
    ff = BOLD if bold else ITAL if italic else REG
    fn = ('LSB' if bold else 'LSI' if italic else 'LSR') + tag
    p = fitz.Point(origin[0] + dx, origin[1] + dy)
    cosv, sinv = direction
    if abs(sinv) < 1e-06 and cosv > 0:
        page.insert_text(p, text, fontsize=size, fontname=fn, fontfile=ff, color=color)
    else:
        ang = math.degrees(math.atan2(-sinv, cosv))
        mat = fitz.Matrix(1, 1).prerotate(ang)
        page.insert_text(p, text, fontsize=size, fontname=fn, fontfile=ff, color=color, morph=(p, mat))

def erase(page, rect, fill=(1, 1, 1), graphics=False):
    page.add_redact_annot(fitz.Rect(rect), fill=fill)
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED if graphics else fitz.PDF_REDACT_LINE_ART_NONE, text=fitz.PDF_REDACT_TEXT_REMOVE)

class Batch:

    def __init__(self, page, tag='B'):
        self.page = page
        self.writes = []
        self.n = 0
        self.tag = tag
        self.done = set()

    def replace(self, needle, new, exact=False, region=None, bold=None, size_scale=1.0, dx=0.0, dy=0.0, color=None, pad=0.6):
        targets = [s for s in find_spans(self.page, needle, exact=exact, region=region) if tuple((round(v, 2) for v in s['bbox'])) not in self.done]
        for s in targets:
            self.done.add(tuple((round(v, 2) for v in s['bbox'])))
            b = fitz.Rect(s['bbox'])
            if abs(s['dir'][1]) < 0.001:
                ym = (b.y0 + b.y1) / 2
                hh = (b.y1 - b.y0) * 0.16
                r = fitz.Rect(b.x0 + 0.25, ym - hh, b.x1 - 0.25, ym + hh)
            else:
                xm = (b.x0 + b.x1) / 2
                ww = (b.x1 - b.x0) * 0.16
                r = fitz.Rect(xm - ww, b.y0 + 0.25, xm + ww, b.y1 - 0.25)
            self.page.add_redact_annot(r, fill=False)
            txt = new(s['text']) if callable(new) else new
            isbold = 'Bold' in s['font'] if bold is None else bold
            if txt == '':
                continue
            self.writes.append((s['origin'], txt, s['size'] * size_scale, color if color else _rgb(s['color']), isbold, s['dir'], dx, dy))
        self.n += len(targets)
        return len(targets)

    def remove_text(self, rect):
        self.page.add_redact_annot(fitz.Rect(rect), fill=False)

    def add_text(self, origin, text, size, color=(0, 0, 0), bold=False, direction=(1.0, 0.0)):
        self.writes.append((origin, text, size, color, bold, direction, 0.0, 0.0))

    def commit(self, images=False):
        self.page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_REMOVE if images else fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE, text=fitz.PDF_REDACT_TEXT_REMOVE)
        for o, t, sz, c, b, d, dx, dy in self.writes:
            write(self.page, o, t, sz, c, b, d, dx=dx, dy=dy, tag=self.tag)
_NUM = re.compile('^[-\\u2212]?\\d+(\\.\\d+)?$')

def fix_minus_signs(batch, page):
    font = fitz.Font(fontfile=REG)
    num = [s for s in spans(page) if _NUM.match(s['text'].strip()) and abs(s['dir'][1]) < 0.001]
    n = 0
    for s in num:
        t = s['text'].strip()
        if not t.startswith('-'):
            continue
        b = fitz.Rect(s['bbox'])
        others = [o for o in num if o is not s]
        right = any((abs(o['bbox'][2] - b.x1) < 0.3 and abs(o['bbox'][0] - b.x0) > 0.3 for o in others))
        left = any((abs(o['bbox'][0] - b.x0) < 0.3 and abs(o['bbox'][2] - b.x1) > 0.3 for o in others))
        new = '−' + t[1:]
        dw = font.text_length(new, fontsize=s['size']) - b.width
        dx = -dw if right else 0.0 if left else -dw / 2
        cx, cy = ((b.x0 + b.x1) / 2, (b.y0 + b.y1) / 2)
        n += batch.replace(t, new, exact=True, region=fitz.Rect(cx - 0.2, cy - 0.2, cx + 0.2, cy + 0.2), dx=dx)
    return n
