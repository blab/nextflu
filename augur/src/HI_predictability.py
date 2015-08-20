from diagnostic_figures import large_effect_mutations, figheight
from itertools import izip
from H3N2_process import H3N2_process, virus_config
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from scipy.stats import ks_2samp
from fitness_tolerance import *
plt.ion()


params = {
    'lam_HI':1,
    'lam_avi':2,
    'lam_pot':0.3,
    'prefix':'H3N2_',

}

mut_models = True
if mut_models:
    resolutions = ['1985to1995','1990to2000','1995to2005','2000to2010','2005to2016']
    fig, axs = plt.subplots(1,len(resolutions), sharey=True, figsize=(4*figheight, 1.3*figheight))
    cols={}
    HI_distributions_mutations = []
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
        for mut in myH3N2.mutation_effects:
            HI = myH3N2.mutation_effects[mut]
            mutlabel = mut[0]+':'+mut[1][1:]
            if mutlabel in myH3N2.frequencies["mutations"]["global"]:
                HI_distributions_mutations.append([res, mut, HI, np.array(myH3N2.frequencies["mutations"]["global"][mutlabel])])
            else:
                print("no frequencies for ",mut, 'HI', HI)
                continue

    plt.savefig("trajectories_mutations.pdf")

    ### make cumulative distribution of HI titers that fix or don't
    freq_thres = 0.5
    fixed  = np.array([ HI for res, mut, HI, freq in HI_distributions_mutations
                      if freq[0]<0.1 and freq.max()>freq_thres and HI>0.1])
    failed = np.array([ HI for res, mut, HI, freq in HI_distributions_mutations
                      if freq[0]<0.1 and freq.max()<freq_thres and HI>0.1])

    D, p = ks_2samp(fixed, failed)
    print("KS stat:", D, "p-val:",p)
    plt.figure()
    plt.plot(sorted(fixed), np.linspace(0,1,len(fixed)), label = '>'+str(freq_thres)+' n='+str(len(fixed)))
    plt.plot(sorted(failed), np.linspace(0,1,len(failed)), label = '<'+str(freq_thres)+' n='+str(len(failed)))
    plt.xlabel('HI effect')
    plt.ylabel('cumulative distribution')
    plt.legend(loc=4)
    plt.savefig("cumulative_HI_mutations.pdf")

    ### make plot with fraction successful depending on HI effect
    plt.figure()
    HI_threshold = np.array([0.1, 0.3, 0.8, 1.5, 4])
    HI_binc = 0.5*(HI_threshold[:-1]+HI_threshold[1:])
    HI_max = np.array([[HI, freq.max()] for res, mut, HI, freq in HI_distributions_mutations if freq[0]<0.1])
    for freq_thres in [0.5, 0.75, 0.95]:
        frac_success = []
        for HI_lower, HI_upper in zip(HI_threshold[:-1], HI_threshold[1:]):
            ind = (HI_max[:,0]>=HI_lower)&(HI_max[:,0]<HI_upper)
            print(HI_lower, ind.sum())
            frac_success.append((HI_max[ind,1]>freq_thres).mean())
        plt.plot(HI_binc, frac_success, '-o', label = "max freq >"+str(freq_thres))

    plt.legend(loc=2)
    plt.ylim([0,1])
    plt.xlabel('HI effect bins 0.1-0.3, 0.3-0.8, 0.8-1.5, >1.5')
    plt.ylabel('fraction reaching frequency threshold')
    plt.savefig('fraction_successful.pdf')

    ### make cumulative HI on backbone
    HI_backbone = np.array([HI for res, mut, HI, freq in HI_distributions_mutations if freq[0]<0.1 and freq.max()>0.75])
    HI_backbone.sort()
    plt.figure()
    cumHI =  HI_backbone.cumsum()
    plt.plot(HI_backbone, cumHI/cumHI[-1])
    plt.ylabel('fraction of HI due to effects < cutoff')
    plt.xlabel('effect size')
    plt.savefig('cumulative_HI_effects.pdf')

### repeat for tree model
res = '1985to2016'
params['pivots_per_year'] = 3.0
params['resolution']=res
params['time_interval'] = map(float, res.split('to'))

if params['time_interval'][1]>2015:
    params['time_interval'][1]=2015.8

# add all arguments to virus_config (possibly overriding)
virus_config.update(params)
# pass all these arguments to the processor: will be passed down as kwargs through all classes
myH3N2 = H3N2_process(**virus_config)
myH3N2.load()
# assign dates
for node in myH3N2.tree.postorder_internal_node_iter():
    node.num_date = np.min([c.num_date for c in node.child_nodes()])

assign_fitness(myH3N2.tree)
dates_fitness = np.array([(n.num_date, n.tol) for n in myH3N2.tree.postorder_internal_node_iter()])

pivots = myH3N2.tree.seed_node.pivots
HI_vs_max_freq_tree = []
dt=1.0
for node in myH3N2.tree.postorder_internal_node_iter():
    if node.num_date<1987:
        continue
    if node.freq["global"] is not None and node.freq["global"].max()>0.1:
        p = node.parent_node
        cHI = node.dHI
        while  p is not None and (node.num_date - p.num_date)<dt:
            cHI += p.dHI
            p = p.parent_node

        ind = (dates_fitness[:,0]<=node.num_date)&(dates_fitness[:,0]>node.num_date-dt)
        HI_vs_max_freq_tree.append((cHI, (node.tol-np.mean(dates_fitness[ind,1])), np.array(node.freq["global"])))

freq_clusters = []
globbing_thres=0.2
for cHI, tol, freq in HI_vs_max_freq_tree:
    found = False
    for fi, (cHIs, ctol, cfreqs) in enumerate(freq_clusters):
        if np.max(np.abs(freq - np.mean(cfreqs, axis=0)))<globbing_thres:
            freq_clusters[fi][2].append(freq)
            freq_clusters[fi][1].append(tol)
            freq_clusters[fi][0].append(cHI)
            found=True
    if not found:
        freq_clusters.append([[cHI], [tol], [freq]])

freq_thres = 0.75
fixed  = np.array([ max(HI) for HI, tol, freqs in freq_clusters
                  if np.max(np.mean(freqs, axis=0))>freq_thres and max(HI)>0.01])
failed  = np.array([ max(HI) for HI, tol, freqs in freq_clusters
                  if np.max(np.mean(freqs, axis=0))<freq_thres and max(HI)>0.01])

D, p = ks_2samp(fixed, failed)
print("KS stat:", D, "p-val:",p)
plt.figure()
plt.plot(sorted(fixed), np.linspace(0,1,len(fixed)), label = '>'+str(freq_thres)+' n='+str(len(fixed)))
plt.plot(sorted(failed), np.linspace(0,1,len(failed)), label = '<'+str(freq_thres)+' n='+str(len(failed)))
plt.xlabel('HI effect')
plt.ylabel('cumulative distribution')
plt.legend(loc=4)
plt.savefig("cumulative_HI_tree.pdf")

fixed  = np.array([ np.mean(tol) for HI, tol, freqs in freq_clusters
                  if np.max(np.mean(freqs, axis=0))>freq_thres and max(HI)>0.2])
failed  = np.array([ np.mean(tol) for HI, tol, freqs in freq_clusters
                  if np.max(np.mean(freqs, axis=0))<freq_thres and max(HI)>0.2])


D, p = ks_2samp(fixed, failed)
print("KS stat:", D, "p-val:",p)
plt.figure()
plt.plot(sorted(fixed), np.linspace(0,1,len(fixed)), label = '>'+str(freq_thres)+' n='+str(len(fixed)))
plt.plot(sorted(failed), np.linspace(0,1,len(failed)), label = '<'+str(freq_thres)+' n='+str(len(failed)))
plt.xlabel('tolerance')
plt.ylabel('cumulative distribution')
plt.legend(loc=4)
plt.savefig("cumulative_tol_tree.pdf")

################################################################
#### plot tree frequencies
################################################################
fs=14
HI_cutoff=0.3
mycmap = cm.cool
plt.figure(figsize=(3*figheight, figheight))
ax1 = plt.subplot(1,1,1)
for cHIs, tol, cfreqs in freq_clusters:
    if max(cHIs)>HI_cutoff:
        ax1.plot(pivots,np.mean(cfreqs, axis=0), c=mycmap(np.sqrt(np.max(cHIs))/2))

sm = plt.cm.ScalarMappable(cmap=mycmap, norm=plt.Normalize(vmin=0, vmax=2))
# fake up the array of the scalar mappable. Urgh...
sm._A = []
cb = plt.colorbar(sm)
cb.set_ticks(np.sqrt([0, 0.3, 1, 2,4]))
cb.set_ticklabels(map(str, [0, 0.3, 1, 2,4]))
cb.set_label('HI effect', fontsize=fs)
ax1.set_ylabel('frequency', fontsize=fs)
ax1.set_xlabel('year', fontsize=fs)
ax1.tick_params(labelsize=fs)
plt.tight_layout()
plt.savefig("trajectories_tree.pdf")


################################################################
##### add fraction successful
################################################################
plt.figure(figsize=(2.4*figheight, figheight))
ax2 = plt.subplot2grid((1,2),(0,0))
plt.title("tree model", fontsize=fs)
#HI_threshold = np.array([0.1, 0.3, 0.8, 1.5, 4])
HI_threshold = np.array([0.1, 0.3, 1, 4])
HI_binc = 0.5*(HI_threshold[:-1]+HI_threshold[1:])
HI_max = np.array([[np.max(HI), np.max(np.mean(freqs, axis=0))] for HI, tol, freqs in freq_clusters])
for freq_thres in [0.5, 0.75, 0.95]:
    frac_success = []
    for HI_lower, HI_upper in zip(HI_threshold[:-1], HI_threshold[1:]):
        ind = (HI_max[:,0]>=HI_lower)&(HI_max[:,0]<HI_upper)
        print(HI_lower, ind.sum())
        frac_success.append((HI_max[ind,1]>freq_thres).mean())
    ax2.plot(np.arange(len(frac_success))+0.5, frac_success, 'o-', label = "max freq >"+str(freq_thres))

ax2.set_xlabel('HI effect', fontsize=fs)
ax2.set_ylabel('fraction reaching frequency threshold', fontsize=fs)
ax2.tick_params(labelsize=fs)
ax2.set_xticks(np.arange(len(HI_binc))+0.5)
ax2.set_xticklabels([str(lower)+'-'+str(upper) for lower, upper in zip(HI_threshold[:-1], HI_threshold[1:])])
plt.legend(loc=4, fontsize=fs)
plt.ylim([0,1])
plt.xlim([0,len(HI_binc)])

ax3 = plt.subplot2grid((1,2),(0,1))
plt.title("mutation model", fontsize=fs)
HI_max = np.array([[HI, freq.max()] for res, mut, HI, freq in HI_distributions_mutations if freq[0]<0.1])
for freq_thres in [0.5, 0.75, 0.95]:
    frac_success = []
    for HI_lower, HI_upper in zip(HI_threshold[:-1], HI_threshold[1:]):
        ind = (HI_max[:,0]>=HI_lower)&(HI_max[:,0]<HI_upper)
        print(HI_lower, ind.sum())
        frac_success.append((HI_max[ind,1]>freq_thres).mean())
    ax3.plot(np.arange(len(frac_success))+0.5, frac_success, 'o-', label = "max freq >"+str(freq_thres))

ax3.set_xlabel('HI effect', fontsize=fs)
ax3.set_ylabel('fraction reaching frequency threshold', fontsize=fs)
ax3.tick_params(labelsize=fs)
ax3.set_xticks(np.arange(len(HI_binc))+0.5)
ax3.set_xticklabels([str(lower)+'-'+str(upper) for lower, upper in zip(HI_threshold[:-1], HI_threshold[1:])])
plt.legend(loc=4, fontsize=fs)
plt.ylim([0,1])
plt.xlim([0,len(HI_binc)])

plt.tight_layout()
plt.savefig('combined_HI_dynamics.pdf')
