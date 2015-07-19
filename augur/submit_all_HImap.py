
import os

res='12y'
for reg in [0.1, 0.3, 1, 3, 10]:
    for flu in ['H3N2', 'H1N1pdm', 'Vic', 'Yam']:
        for minaa in [0]: #[0,1,'epi']:
            for model in ["--map_to_tree",""]:
                for train_strains in ['', '--train_strains']:
                    call = ' '.join(['qsub -cwd -l h_vmem=8G -l h_rt=20:59:00 submit_HIvalidation.sh',
                              flu, res, str(minaa), train_strains, str(reg), model])
                    print call
                    os.system(call)
