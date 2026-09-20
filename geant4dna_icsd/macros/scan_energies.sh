#!/usr/bin/env bash
set -euo pipefail

OPTION="${1:-2}"
NEVENTS="${2:-100000}"
OUTDIR="results_opt${OPTION}"
mkdir -p "${OUTDIR}"

energy_diameter() {
  case "$1" in
    20|50|100|300|600|1000) echo 8 ;;
    5000|10000) echo 100 ;;
  esac
}

for E in 20 50 100 300 600 1000 5000 10000; do
  D=$(energy_diameter "$E")
  MAC="macros/_tmp_${E}eV.mac"
  cat > "${MAC}" <<MACEOF
/gun/particle e-
/gun/energy ${E} eV
/run/beamOn ${NEVENTS}
MACEOF
  echo "== Running ${E} eV, option ${OPTION}, target ${D} nm =="
  ./nanoICSD "${MAC}" "${OPTION}" "${D}"
  mv icsd_run0.csv "${OUTDIR}/icsd_${E}eV.csv"
  rm -f "${MAC}"
done

echo "Done. Per-event ionization counts in ${OUTDIR}/icsd_*.csv"
echo "Next: python3 analysis/analyze_icsd.py ${OUTDIR}"
