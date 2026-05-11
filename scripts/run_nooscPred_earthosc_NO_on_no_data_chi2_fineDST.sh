#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${VALOR_EXPT:-}" ]]; then
  echo "[ERROR] VALOR_EXPT is not set. Source JUNO/VALOR setup first."
  exit 1
fi

TAG="${TAG:-nooscPred_earthosc_NO_on_no_data_fineDST}"
RUN="${RUN:-1}"
CFG="${VALOR_EXPT}/config/analysis/Atmo_Juno_4cc_fc_2d_perfectdet_nooscPred_noData_fineDST.xml"
LOG="/tmp/${TAG}_${RUN}.log"

# NO ordering (normal hierarchy) nominal inputs
# Note: these defaults should lie within the ranges in config/params/flux_3fosc.xml
# (otherwise MINUIT will clamp them back inside bounds before the first evaluation).
INPUT_MIXING12="${INPUT_MIXING12:-0.29}"
INPUT_MIXING13="${INPUT_MIXING13:-0.02}"
INPUT_MIXING23="${INPUT_MIXING23:-0.5}"
INPUT_DM12SQ="${INPUT_DM12SQ:-7.01e-5}"
INPUT_DM23SQ="${INPUT_DM23SQ:-2.4e-3}"

# Set CPV as a * pi by default. You can still override with INPUT_CPV directly.
CPV_PI_FACTOR="${CPV_PI_FACTOR:-0.19}"
INPUT_CPV="${INPUT_CPV:-$(awk -v a="${CPV_PI_FACTOR}" 'BEGIN { printf "%.15f", a * atan2(0,-1) }')}"

if [[ ! -f "${CFG}" ]]; then
  echo "[ERROR] Missing config: ${CFG}"
  exit 1
fi

echo "[INFO] Running NOOSC prediction + fitter Earth oscillation (3f_osc), NO ordering, against NO data"
echo "[INFO] Config: ${CFG}"
echo "[INFO] Log: ${LOG}"
echo "[INFO] Physics hypothesis: 3f_osc (SDK parses this from --physics-hypothesis)"
echo "[INFO] Ordering: normal (--normal)"
echo "[INFO] Fit mode: float mixing12, fix remaining oscillation parameters"

echo "[INFO] Osc inputs: mixing12=${INPUT_MIXING12}, mixing13=${INPUT_MIXING13}, mixing23=${INPUT_MIXING23}, dm12sq=${INPUT_DM12SQ}, dm23sq=${INPUT_DM23SQ}, cpv=${INPUT_CPV}"

set +e
valor-run_fit \
  --config "${CFG}" \
  --physics-hypothesis 3f_osc \
  --param-values-option prefit-errors \
  --osc-prior-file config/priors/3Fosc_priors.xml \
  --load-prior-from-file - \
  --fitter minuit+with_simplex_fallback \
  --mesg-thresholds verbosity.xml \
  --tag "${TAG}" \
  --run "${RUN}" \
  --ntoymc 1 \
  --fit-param-opt NONE \
  --without-stat-fluct \
  --syst-prior-file - \
  --read-data \
  --normal \
  --input-mixing12 "${INPUT_MIXING12}" \
  --input-mixing13 "${INPUT_MIXING13}" \
  --input-mixing23 "${INPUT_MIXING23}" \
  --input-dm12sq "${INPUT_DM12SQ}" \
  --input-dm23sq "${INPUT_DM23SQ}" \
  --input-cpv "${INPUT_CPV}" \
  2>&1 | tee "${LOG}"
FIT_RC=${PIPESTATUS[0]}
set -e

if [[ "${FIT_RC}" -ne 0 ]]; then
  if [[ "${FIT_RC}" -eq 139 ]] && grep -q "Fit ends with statistic" "${LOG}"; then
    echo "[WARN] valor-run_fit exited with code 139 after reporting a completed fit."
    echo "[WARN] Treating this as a ROOT shutdown-time cleanup segfault; keeping results."
  else
    echo "[ERROR] valor-run_fit failed with exit code ${FIT_RC}. See ${LOG}"
    exit "${FIT_RC}"
  fi
fi

echo "[INFO] Extracting key lines"
grep -nE "Data file given so will Read|Fit starts with statistic|Fit ends with statistic|Best-fit statistic" "${LOG}" | tail -n 30 || true

echo "[DONE] Completed NOOSC+Earth-osc (NO) chi2 check on NO data (fine DST)"
