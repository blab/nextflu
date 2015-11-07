import os

resolutions = ['1y','3y', '6y', '12y']
lineages = ['H3N2', 'H1N1pdm', 'Vic', 'Yam']

for lineage in lineages:
	for res in resolutions:
		for f in ['tree', 'frequencies','meta', 'sequences']:
			fname = '_'.join([lineage, res, f]) + '.json'
			os.system('wget nextflu.org/data/'+fname +' --directory-prefix data') 
