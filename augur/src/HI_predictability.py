from diagnostic_figures import large_effect_mutations, figheight
from itertools import izip
from H3N2_process import H3N2_process, virus_config
import matplotlib.pyplot as plt
plt.ion()

params = {
    'lam_HI':1,
    'lam_avi':2,
    'lam_pot':0.3,
    'prefix':'H3N2_',

}

resolutions = ['1985to1995','1990to2000','1995to2005','2000to2010','2005to2016']
fig, axs = plt.subplots(1,len(resolutions), sharey=True, figsize=(3*figheight, figheight))
cols={}
for res,ax in izip(resolutions,axs):
    params['pivots_per_year'] = 6.0
    params['resolution']=res
    params['time_interval'] = map(float, res.split('to'))

    if params['time_interval'][1]>2015:
        params['time_interval'][1]=2015.8

    # add all arguments to virus_config (possibly overriding)
    virus_config.update(params)
    # pass all these arguments to the processor: will be passed down as kwargs through all classes
    myH3N2 = H3N2_process(**virus_config)
    myH3N2.run(['HI'],
               lam_HI = virus_config['lam_HI'],
               lam_avi = virus_config['lam_avi'],
               lam_pot = virus_config['lam_pot'],
               )

    cols = large_effect_mutations(myH3N2, ax, cols)

