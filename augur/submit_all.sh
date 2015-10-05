#!/usr/bin/env bash

qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh H3N2 3y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh H3N2 6y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh H3N2 12y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh H3N2 20y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh H1N1pdm 3y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh H1N1pdm 6y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh H1N1pdm 12y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh Vic 3y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh Vic 6y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh Vic 12y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh Vic 20y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh Yam 3y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh Yam 6y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh Yam 12y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh Yam 20y

#qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh H3N2 1985 1995 50
#qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh H3N2 1990 2000 50
#qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh H3N2 1995 2005 50
#qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh H3N2 2000 2010 20
#qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh H3N2 2005 2016 15
#qsub -cwd -l h_vmem=8G -l h_rt=20:59:00 submit.sh H3N2 1995 2016 10
#qsub -cwd -l h_vmem=8G -l h_rt=20:59:00 submit.sh H3N2 1985 2016 10

