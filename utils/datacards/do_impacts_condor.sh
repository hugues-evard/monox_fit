#!/bin/bash

# Impacts submission to condor script

CHANNEL="$1"
YEARS=("Run3")

mkdir -p impacts
pushd impacts > /dev/null

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

    mkdir ${YEAR}
    pushd ${YEAR} > /dev/null

    COMMON_OPTS="-t -1 -m 125 --parallel=4 --rMin=-5 --rMax=5 --robustFit 1 --cminDefaultMinimizerStrategy 0 --autoRange 5 --squareDistPoiStep"
    # -t -1 -m 125 --robustFit 1

    # this is the same as for nocondor
    echo "Doing initial fit for ${TAG}"
    combineTool.py ${METHOD} -d ${CARD} --doInitialFit ${COMMON_OPTS} | tee ${LOGFILE}_initialFit.txt

    # This gets submitted to condor
    echo "Submitting Impacts for ${TAG}"
    CONDORTAG=task_${YEAR}_${CHANNEL}_${RANDOM}
    combineTool.py ${METHOD} -d ${CARD} --doFits ${COMMON_OPTS} --job-mode condor --task-name ${CONDORTAG} || exit 1

    # Wait for condor jobs to return
    date
    if [ ! -e ${CONDORTAG}*.log ]; then
        exit 1
    fi
    for i in $(seq 1 1 10); do
        condor_wait -debug -status ${CONDORTAG}*.log && break;
        sleep $((i*120))
    done
    date
    sleep 60

    # This is the same as for nocondor
    combineTool.py ${METHOD} -d ${CARD} -o impacts_${TAG}.json ${COMMON_OPTS} | tee ${LOGFILE}_output.txt

    popd > /dev/null
    echo "Plotting Impacts for ${TAG}"
    plotImpacts.py -i ${YEAR}/impacts_${TAG}.json -o impacts_${TAG} --blind
done

popd > /dev/null