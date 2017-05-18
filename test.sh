#!/bin/sh

rm ./testResult -rf
./PSSP analyse \
    "./n20t80.pkl.origin" \
    "./testResult/" \
    --num_iters 10 \
    --concentration 1 \

