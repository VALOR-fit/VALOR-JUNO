#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${VALOR_EXPT:-}" ]]; then
  echo "[ERROR] VALOR_EXPT is not set. Source JUNO/VALOR setup first."
  exit 1
fi

TAG="${TAG:-nooscPred_earthosc_NO_on_no_data_fineDST_detresp}"
RUN="${RUN:-20260522_test}"

CFG="${CFG:-${VALOR_EXPT}/config/analysis/Atmo_Juno_4cc_fc_2d_detresp_nooscPred_noData_fineDST.xml}"

FC_COEFF_DIR="${FC_COEFF_DIR:-${VALOR_EXPT}/data/input/fc_ratio/Total_FC_V2}"
PID_MATRIX_FILE="${PID_MATRIX_FILE:-${VALOR_EXPT}/data/input/detector_response/pid_matrices_8bin_from_user.csv}"

LOG="/scratchfs/juno/heziou/VALOR/softw/JUNO/vmaster/logs/${TAG}_${RUN}.log"

INPUT_MIXING12="${INPUT_MIXING12:-0.307}"
INPUT_MIXING13="${INPUT_MIXING13:-0.0219}"
INPUT_MIXING23="${INPUT_MIXING23:-0.558}"
INPUT_DM12SQ="${INPUT_DM12SQ:-7.53e-5}"
INPUT_DM23SQ="${INPUT_DM23SQ:-2.455e-3}"

CPV_PI_FACTOR="${CPV_PI_FACTOR:-1.19}"
INPUT_CPV="${INPUT_CPV:-$(awk -v a="${CPV_PI_FACTOR}" 'BEGIN { printf "%.15f", a * atan2(0,-1) }')}"

if [[ ! -f "${CFG}" ]]; then
  echo "[ERROR] Missing config: ${CFG}"
  exit 1
fi

for f in \
  "${FC_COEFF_DIR}/coeffs_numu.txt" \
  "${FC_COEFF_DIR}/coeffs_nue.txt" \
  "${FC_COEFF_DIR}/coeffs_antinumu.txt" \
  "${FC_COEFF_DIR}/coeffs_antinue.txt"
do
  if [[ ! -f "$f" ]]; then
    echo "[ERROR] Missing FC coefficient file: $f"
    exit 1
  fi
done

if [[ ! -f "${PID_MATRIX_FILE}" ]]; then
  echo "[ERROR] Missing PID matrix file: ${PID_MATRIX_FILE}"
  exit 1
fi

echo "[INFO] Running NOOSC fine DST + fitter Earth oscillation + JUNO detector response"
echo "[INFO] Config: ${CFG}"
echo "[INFO] FC coeff dir: ${FC_COEFF_DIR}"
echo "[INFO] PID matrix: ${PID_MATRIX_FILE}"
echo "[INFO] Log: ${LOG}"

export LD_LIBRARY_PATH="${VALOR_EXPT}/lib:${LD_LIBRARY_PATH:-}"

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
  --fix-mixing12 --fix-at-mixing12 "${INPUT_MIXING12}" \
  --fix-dm12sq --fix-at-dm12sq "${INPUT_DM12SQ}" \
  --fix-cpv --fix-at-cpv "${INPUT_CPV}" \
  --fix-mixing13 --fix-at-mixing13 "${INPUT_MIXING13}" \
  --fix-mixing23 --fix-at-mixing23 "${INPUT_MIXING23}" \
  --fix-dm23sq --fix-at-dm23sq "${INPUT_DM23SQ}" \
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
grep -nE "DetectorResponse|JunoDetectorResponse|Data file given so will Read|Fit starts with statistic|Fit ends with statistic|Best-fit statistic" "${LOG}" | tail -n 50 || true

echo "[DONE] Completed NOOSC+Earth-osc+detector-response NO chi2 check"