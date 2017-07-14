#!/bin/sh

rm ./testResult -rf
./pscnaclonal analyse \
    "./n20t80.pkl" \
    "./testResult/" \
    --num_iters 10 \
    --concentration 1 \

