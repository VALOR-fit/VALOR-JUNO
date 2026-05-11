#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${VALOR_EXPT:-}" ]]; then
  echo "[ERROR] VALOR_EXPT is not set. Source JUNO/VALOR setup first."
  exit 1
fi

TAG="${TAG:-ioPred_noData_void_fineDST_pid_contained}"
RUN="${RUN:-1}"
CFG="${VALOR_EXPT}/config/analysis/Atmo_Juno_4cc_fc_2d_perfectdet_ioPred_noData_fineDST_pid_contained.xml"
LOG="/tmp/${TAG}_${RUN}.log"

if [[ ! -f "${CFG}" ]]; then
  echo "[ERROR] Missing config: ${CFG}"
  exit 1
fi

echo "[INFO] Running IO-prediction vs NO-data chi2 check (fine DST PID contained)"
echo "[INFO] Config: ${CFG}"
echo "[INFO] Log: ${LOG}"

valor-run_fit \
  --config "${CFG}" \
  --physics-hypothesis void \
  --param-values-option prefit-errors \
  --osc-prior-file config/priors/void_priors.xml \
  --load-prior-from-file - \
  --fitter minuit+with_simplex_fallback \
  --mesg-thresholds verbosity.xml \
  --tag "${TAG}" \
  --run "${RUN}" \
  --ntoymc 1 \
  --fit-param-opt fix-all \
  --without-stat-fluct \
  --syst-prior-file - \
  --read-data \
  2>&1 | tee "${LOG}"

echo "[INFO] Extracting key lines"
grep -nE "Data file given so will Read|Fit starts with statistic|Fit ends with statistic|Best-fit statistic" "${LOG}" | tail -n 20 || true

echo "[DONE] Completed IO-prediction vs NO-data void chi2 check"
