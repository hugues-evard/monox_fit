#!/bin/bash
source "$(dirname "$0")/colors.sh"

# Limit script

CHANNEL="$1"
YEARS=("Run3")

mkdir -p limit
pushd limit > /dev/null

# Uncomment the options you want to use
EXTRA_OPTS=()
# EXTRA_OPTS+=(--saveToys)
# EXTRA_OPTS+=(--rMax 100)
# EXTRA_OPTS+=(--run blind)

for YEAR in "${YEARS[@]}"; do
    TAG="${CHANNEL}_${YEAR}"
    CARD="../cards/card_${TAG}.root"
    LOGFILE="log_limit_${YEAR}.txt"
    METHOD="-M AsymptoticLimits -t -1"

    cecho blue "Running AsymptoticLimits for ${TAG}"
    combine ${METHOD} ${CARD} -n "_${TAG}" "${EXTRA_OPTS[@]}" | tee ${LOGFILE}
done

popd > /dev/null