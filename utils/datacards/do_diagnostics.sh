#!/bin/bash
source "$(dirname "$0")/colors.sh"

# Fit diagnostics script

CHANNEL="$1"
YEARS=("Run3")

mkdir -p diagnostics
pushd diagnostics > /dev/null

# Uncomment the options you want to use
EXTRA_OPTS=()
EXTRA_OPTS+=(--saveShapes)
EXTRA_OPTS+=(--saveWithUncertainties)
EXTRA_OPTS+=(--cminDefaultMinimizerStrategy 0)

for YEAR in "${YEARS[@]}"; do
    TAG="${CHANNEL}_${YEAR}"
    CARD="../cards/card_${TAG}.txt"
    WS="../cards/card_${TAG}.root"
    LOGFILE="log_diag_${YEAR}.txt"
    FITDIAGFILE="fitDiagnostics_${TAG}.root"
    SHAPESFILE="prefit_postfit_shapes_${TAG}.root"
    DIFFROOT="diffnuisances_${TAG}.root"

    cecho blue ">> Running FitDiagnostics for ${TAG}"
    combine -M FitDiagnostics ${WS} -n "_${TAG}" "${EXTRA_OPTS[@]}" | tee ${LOGFILE}

    cecho blue ">> Extracting postfit shapes for ${TAG}"
    PostFitShapesFromWorkspace --workspace ${WS} --datacard ${CARD} --output ${SHAPESFILE} \
        --fitresult "${FITDIAGFILE}:fit_b" --postfit --covariance --print
    
    cecho blue ">> Running diffNuisances for ${FITDIAGFILE}"
    python3 "${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py" ${FITDIAGFILE} -g ${DIFFROOT} --all
done

popd > /dev/null