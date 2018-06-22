#!/bin/sh

rm ./testResult -rf
./pscnaclonal analyse \
    "testData/subsimtree.seg.txt.l1000.bed.inputbase.pSCNAClonal.input.pkl" \
    "./testResult/" \
    --num_iters 10 \
    --concentration 1 \

