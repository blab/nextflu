#!/usr/bin/env bash

#python src/make_all.py --bin /ebio/ag-neher/share/programs/bin/python --all --lineage $1 --resolution $2
#echo src/$1_process.py --prefix $1_ --resolution $2_$3 --interval $2 $3 -v $4 -r 1.0 --skip genotype_frequencies --lam_HI  1 --lam_pot 0.3 --lam_avi 2
/ebio/ag-neher/share/programs/bin/python src/$1_process.py --prefix $1_ --resolution $2to$3 --interval $2 $3 -v $4 -r 1.0 --skip genotype_frequencies --lam_HI  1 --lam_pot 0.3 --lam_avi 2
