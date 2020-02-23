#!/bin/bash

### Asimov limit
mkdir -p limit
pushd limit
for file in ../cards/*.root; do
    TAG=$(basename $file | sed 's/card_//g;s/.root//g');
    combine -M AsymptoticLimits $file -t -1 -n $TAG | tee log_$TAG.txt & 
done
popd