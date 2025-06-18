#!/bin/bash

# Impacts script

CHANNEL="$1"
YEARS=("Run3")

mkdir -p impacts_nocondor
pushd impacts_nocondor > /dev/null

# Uncomment the options you want to use
EXTRA_OPTS=()
# EXTRA_OPTS+=(--saveToys)
# EXTRA_OPTS+=(--rMax 100)
# EXTRA_OPTS+=(--run blind)

for YEAR in "${YEARS[@]}"; do
    TAG="${CHANNEL}_${YEAR}"
    CARD="../../cards/card_${TAG}.root"
    LOGFILE="log_impacts_${YEAR}"
    METHOD="-M Impacts"
    NCORES="32"

    mkdir ${YEAR}
    pushd ${YEAR} > /dev/null

    echo "Running Impacts for ${TAG}"
    # COMMON_OPTS="-t -1 -m 125 --parallel=4 --rMin=-5 --rMax=5 --robustFit 1 --cminDefaultMinimizerStrategy 0 --autoRange 5 --squareDistPoiStep"
    COMMON_OPTS="-t -1 -m 125 --parallel=${NCORES} --rMin=-5 --rMax=5 --robustFit 1 --cminDefaultMinimizerStrategy 0 --autoRange 5 --squareDistPoiStep"
    combineTool.py ${METHOD} -d ${CARD} --doInitialFit ${COMMON_OPTS} | tee ${LOGFILE}_initialFit.txt
    combineTool.py ${METHOD} -d ${CARD} --doFits ${COMMON_OPTS} | tee ${LOGFILE}_Fits.txt
    combineTool.py ${METHOD} -d ${CARD} -o impacts_${TAG}.json ${COMMON_OPTS} | tee ${LOGFILE}_output.txt

    popd > /dev/null
    echo "Plotting Impacts for ${TAG}"
    plotImpacts.py -i ${YEAR}/impacts_${TAG}.json -o impacts_${TAG} # --blind
done

popd > /dev/null