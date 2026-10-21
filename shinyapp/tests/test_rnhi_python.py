# -*- coding: utf-8 -*-
"""Check that the Python implementation reproduces the reference RNHI values.

Run from the application folder:  python3 tests/test_rnhi_python.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

APP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP))
import rnhi_math as rm  # noqa: E402

tab = pd.read_csv(APP / 'data' / 'bulk_compendium_composition.csv', dtype={'sample_id': str}).set_index('sample_id')
ref = pd.read_csv(APP / 'tests' / 'reference_compendium_rnhi.csv', dtype={'sample_id': str}).set_index('sample_id')
healthy = tab[tab.healthy_reference == 1]
aml = tab[tab.healthy_reference == 0]
assert len(healthy) == 45 and len(aml) == 3310

dl = rm.detection_limits(rm.closure(aml[rm.LINEAGES].to_numpy(float)))

# Complete compendium with the fixed detection limits
full = rm.score_table(aml, healthy, dl=dl)
d_full = float(np.max(np.abs(full.loc[ref.index].to_numpy() - ref.RNHI.to_numpy())))
assert d_full < 1e-12, d_full
print('PASS complete compendium reproduces the reference values (max abs difference %.1e)' % d_full)

# Each sample scored on its own with the fixed detection limits
center = rm.comp_center(healthy)
single = np.array([rm.rnhi_score(rm.closure(aml.loc[[s], rm.LINEAGES].to_numpy(float)), center, dl=dl)[0] for s in aml.index])
d_single = float(np.max(np.abs(single - ref.loc[aml.index, 'RNHI'].to_numpy())))
assert d_single < 1e-12, d_single
print('PASS single-sample scoring reproduces the reference values (max abs difference %.1e)' % d_single)
print('All tests passed.')
