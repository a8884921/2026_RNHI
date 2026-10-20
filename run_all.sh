set -euo pipefail
cd "$(dirname "$0")"
export RNHI_FONT="${RNHI_FONT:-Arial}"
if [ "$#" -gt 0 ]; then
  targets=("$@")
else
  targets=(Figure_*)
fi
count=0
for fig in "${targets[@]}"; do
  [ -d "$fig" ] || continue
  printf '%s\n' "$fig"
  for panel in "$fig"/Panel_*/; do
    [ -f "$panel/code/run_panel.py" ] || continue
    (cd "$panel" && python3 code/run_panel.py)
    count=$((count + 1))
  done
  for build in "$fig"/00_build_*.py; do
    [ -f "$build" ] || continue
    (cd "$fig" && python3 "$(basename "$build")")
    count=$((count + 1))
  done
done
if [ "$count" -eq 0 ]; then
  printf '%s\n' 'No configured panel runner or figure builder was found.'
else
  printf 'Completed %s configured runners and builders.\n' "$count"
fi
