#!/usr/bin/env bash

qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh H3N2 3y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh H3N2 6y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh H3N2 12y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh H3N2 20y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh H1N1pdm 3y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh H1N1pdm 7y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh Vic 3y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh Vic 6y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh Vic 12y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh Vic 20y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh Yam 3y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh Yam 6y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh Yam 12y
qsub -cwd -l h_vmem=8G -l h_rt=10:59:00 submit.sh Yam 20y

