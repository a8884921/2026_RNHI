import sys
from pathlib import Path
import fitz
CODE = Path(__file__).resolve().parent.parent
FIGDIR = Path(sys.argv[1]).resolve()
PDF = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else None
import os
PADX = float(os.environ.get('CROP_PAD', 6.0))
if PDF is None:
    cand = list((CODE / 'data' / 'reference_figures').glob(FIGDIR.name.replace('Figure_', 'Figure_') + '.pdf'))
    cand += list((CODE / 'data' / 'reference_figures').glob('Supplementary_' + FIGDIR.name + '.pdf'))
    PDF = cand[0]

def letters(doc):
    out = []
    for b in doc[0].get_text('dict')['blocks']:
        for l in b.get('lines', []):
            for s in l['spans']:
                t = s['text'].strip()
                if len(t) == 1 and t in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' and (s['size'] > 9):
                    out.append((t, s['bbox']))
    return sorted(out, key=lambda r: (round(r[1][1] / 24), r[1][0]))
doc = fitz.open(str(PDF))
page = doc[0]
H, W = (page.rect.height, page.rect.width)
L = letters(doc)
rows = []
for letter, bb in L:
    if rows and abs(bb[1] - rows[-1][0][1][1]) < 24:
        rows[-1].append((letter, bb))
    else:
        rows.append([(letter, bb)])
boxes = {}
for ri, row in enumerate(rows):
    row = sorted(row, key=lambda r: r[1][0])
    y0 = max(0.0, row[0][1][1] - 6)
    y1 = min((r[1][1] for r in rows[ri + 1])) - 2 if ri + 1 < len(rows) else H
    for i, (letter, bb) in enumerate(row):
        x0 = max(0.0, bb[0] - (8 if i == 0 else PADX))
        x1 = row[i + 1][1][0] - 4 if i + 1 < len(row) else W
        boxes[letter] = fitz.Rect(x0, y0, x1, y1)
print(PDF.name, '->', {k: (round(v.width), round(v.height)) for k, v in boxes.items()})
for letter, box in boxes.items():
    hits = sorted(FIGDIR.glob('Panel_%s_*' % letter))
    if not hits:
        print('  (skip %s: folder missing)' % letter)
        continue
    d = hits[0]
    out_dir = d / 'output'
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = '%s_Panel_%s_v1' % (FIGDIR.name.replace('Supplementary_', ''), letter)
    pdf_path = out_dir / (stem + '.pdf')
    if pdf_path.exists():
        print('  keep existing', pdf_path.name)
        continue
    new = fitz.open()
    pg = new.new_page(width=box.width, height=box.height)
    pg.show_pdf_page(pg.rect, doc, 0, clip=box)
    new.save(str(pdf_path), garbage=1, deflate=True)
    new[0].get_pixmap(dpi=600).save(str(out_dir / (stem + '.png')))
    new.close()
    print('  cropped', letter, round(box.width), 'x', round(box.height), 'pt')
doc.close()
