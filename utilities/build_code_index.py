import csv
from pathlib import Path
root = Path(__file__).resolve().parent.parent
out = root / 'output'
out.mkdir(parents=True, exist_ok=True)
with (out / 'code_index.csv').open('w', newline='', encoding='utf-8') as stream:
    writer = csv.writer(stream)
    writer.writerow(['figure', 'panel', 'code'])
    for file in sorted(root.glob('Figure_*/**/code/*')):
        if file.suffix in {'.py', '.R', '.sh'}:
            relative = file.relative_to(root)
            writer.writerow([relative.parts[0], relative.parts[1], relative.as_posix()])
