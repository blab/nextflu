#!/usr/bin/env bash

echo src/HI_map_validation.py --prefix $1_ --flutype $1 --resolution $2 --min_aamuts $3 --reg $4 $5 $6 
/ebio/ag-neher/share/programs/bin/python src/HI_map_validation.py --prefix $1_ --flutype $1 --resolution $2 --min_aamuts $3 --reg $4 $5 $6

