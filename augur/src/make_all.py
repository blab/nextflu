#!/usr/bin/env python
import subprocess,glob, os

run_pipeline = []
for flu in ["H3N2", "H1N1pdm", "H1N1", "Vic", "Yam"]:
	auspice_files = glob.glob('../auspice/data/'+flu+"*json")
	input_files = glob.glob('data/'+flu+'*fasta')
	if len(input_files):
		if len(auspice_files)==0 or min([os.path.getmtime(fname) for fname in auspice_files])<max([os.path.getmtime(fname) for fname in input_files]):
			run_pipeline.append(flu)

#python_binary = '/ebio/ag-neher/share/programs/bin/python'
python_binary = 'python'

if 'H3N2' in run_pipeline:
	call = map(str, [python_binary, 'src/H3N2_process.py', '-v', 50, '-y', 3, '--skip', 'genotype_frequencies', '-r', 1.0, '--prefix', 'data/H3N2_'])
	print call
	subprocess.call(call)
if 'H1N1' in run_pipeline:
	call = map(str, [python_binary, 'src/H1N1historical_process.py', '-v', 20, '--interval', 1990, 2010, '--skip', 'genotype_frequencies', '-r', 1.0, '--prefix', 'data/H1N1_'])
	print call
	subprocess.call(call)
if 'H1N1pdm' in run_pipeline:
	call = map(str, [python_binary, 'src/H1N1pdm_process.py', '-v', 30, '-y', 6, '--skip', 'genotype_frequencies', '-r', 1.0, '--prefix', 'data/H1N1pdm_'])
	print call
	subprocess.call(call)
if 'Vic' in run_pipeline:
	call = map(str, [python_binary, 'src/Vic_process.py', '-v', 30, '-y', 6, '--skip', 'genotype_frequencies', '-r', 1.0, '--prefix', 'data/Vic_'])
	print call
	subprocess.call(call)
if 'Yam' in run_pipeline:
	call = map(str, [python_binary, 'src/Vic_process.py', '-v', 30, '-y', 6, '--skip', 'genotype_frequencies', '-r', 1.0, '--prefix', 'data/Yam_'])
	print call
	subprocess.call(call)
