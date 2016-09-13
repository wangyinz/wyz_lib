#!/bin/bash
set -x

idtext=test_migsim_id.txt
id=`cat ${idtext}`
id=$((id+1))
echo ${id} > ${idtext}
/N/u/xinliu/Quarry/mainentr/seism/src/pwmig/migsimulation/migsimulation  /N/dc/scratch/xinliu/evstacks test_${id}
