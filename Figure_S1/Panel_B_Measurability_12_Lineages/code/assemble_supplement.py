from pathlib import Path as _ProjectPath
import sys as _project_sys
_PROJECT_ROOT = next((p for p in _ProjectPath(__file__).resolve().parents if (p / 'common' / 'panel_io.py').is_file()), None)
if _PROJECT_ROOT is None:
    raise RuntimeError('The common analysis module is required.')
_project_sys.path.insert(0, str(_PROJECT_ROOT / 'common'))
import panel_io as _project_io
import importlib.util as _ilu
import sys as _sys
import time
from pathlib import Path as _Path
import pymupdf as fitz
HERE = _Path(__file__).resolve().parent
ROOT = _PROJECT_ROOT
SRC = _PROJECT_ROOT / 'data' / 'reference_figures'
FIGS = _PROJECT_ROOT / 'output' / 'figure_s1'
FIGS.mkdir(parents=True, exist_ok=True)
_spec = _ilu.spec_from_file_location('pdfedit', HERE / 'pdf_text.py')
_pe = _ilu.module_from_spec(_spec)
_sys.modules['pdfedit'] = _pe
_spec.loader.exec_module(_pe)
from pdfedit import *
src = fitz.open(SRC / 'Figure_S01.pdf')
p = src[0]
p.add_redact_annot(fitz.Rect(0, 420, 595.276, 841.89), fill=False)
p.apply_redactions(images=fitz.PDF_REDACT_IMAGE_REMOVE, graphics=fitz.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED, text=fitz.PDF_REDACT_TEXT_REMOVE)
B = Batch(p, tag='S1')
B.replace('1 states', lambda t: t.replace('1 states', '1 state'))
B.replace('Supplementary Figure 1', 'Supplementary Figure S1', exact=True)
for old_t, new_t in [('EarlyErythroid', 'Early erythroid'), ('MyeloidProgenitor', 'Myeloid progenitor'), ('HSC_MPP', 'HSC/MPP'), ('LateErythroid', 'Late erythroid'), ('Other_Stromal (excluded)', 'Other stromal (excluded)')]:
    B.replace(old_t, new_t, exact=True)
B.commit()
new = fitz.open()
np_ = new.new_page(width=595.276, height=841.89)
np_.show_pdf_page(fitz.Rect(0, 0, 595.276, 420), src, 0, clip=fitz.Rect(0, 0, 595.276, 420))
b = fitz.open(str(FIGS) + '/S1B_panel.pdf')
w = 595.276
h = w * b[0].rect.height / b[0].rect.width
np_.show_pdf_page(fitz.Rect(0, 436, w, 436 + h), b, 0)
write(np_, (21, 458), 'B', 18, (0.1, 0.1, 0.1), tag='S1new')
new.save(str(FIGS) + '/Supplementary_Figure_S1.pdf', garbage=4, deflate=True)
pix = np_.get_pixmap(matrix=fitz.Matrix(1.6, 1.6))
pix.save(str(FIGS) + '/Supplementary_Figure_S1.png')
print('ok')
