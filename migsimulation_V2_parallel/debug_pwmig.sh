#!/bin/bash
set -x

idtext=test_migsim_id.txt
id=`cat ${idtext}`
id=$((id+1))
echo ${id} > ${idtext}
totalview /N/u/xinliu/Quarry/mainentr/seism/src/pwmig/migsimulation/migsimulation -a /N/dc/scratch/xinliu/evstacks test_${id}
