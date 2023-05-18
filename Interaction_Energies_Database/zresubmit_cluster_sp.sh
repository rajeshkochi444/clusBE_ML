#!/bin/bash

for f in N*; do
        cd  $f
        echo $f;
        tot=$(ls VASP_$f/OUTCARS_List  | wc -l)
        echo $tot
        #sed -i  "s/{0..3000}/{${tot}..2000}/g" traj_bare_cluster_sp.sh
        #sbatch traj_bare_cluster_sp.sh
        cd ..
done

