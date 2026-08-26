#!/bin/bash
# Build and run all four systems sequentially
set -e
BASE="$(cd "$(dirname "$0")/.." && pwd)"
for S in PXR_PFECHS PPARa_PFECHS PXR_F-53B PPARa_F-53B; do
  echo ""
  echo "############################################"
  echo "#  $S"
  echo "############################################"
  bash "$BASE/scripts/build_system.sh" "$S"
  bash "$BASE/scripts/run_md.sh" "$S" 0
done
echo ""
echo "All systems done. Next: bash analyze.sh"
