#!/bin/bash

rank_should_profile=false
if [[ -z $ranks_to_profile ]]; then
   ranks_to_profile="0,12"
fi
for i in $(echo $ranks_to_profile | sed "s/,/ /g"); do
   if [[ $i == $OMPI_COMM_WORLD_RANK ]]; then
      rank_should_profile=true
   fi
done
if [[ $ranks_to_profile == all ]]; then
   rank_should_profile=true
fi

if $rank_should_profile; then
   if [[ -n $NSYS_TIMEOUT ]]; then
      timeout_params="-d $NSYS_TIMEOUT"
   fi
   echo "Louis debug NSYS rank $OMPI_COMM_WORLD_RANK"
   set -x
   nsys profile $timeout_params --kill=none --output="nsys.$SLURM_JOBID.$COMMENT$OMPI_COMM_WORLD_RANK.qdrep" $@
else
   $@
fi
