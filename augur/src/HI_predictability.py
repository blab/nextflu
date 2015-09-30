######
# script that explores the predictive power of inferred antigenic change
# It tree and mutation models inferred in intervals of 10years for H3N2
#
######

from collections import defaultdict
from diagnostic_figures import large_effect_mutations, figheight
from itertools import izip
from H3N2_process import H3N2_process, virus_config
from diagnostic_figures import get_slope
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from scipy.stats import ks_2samp
from fitness_tolerance import *
import seaborn as sns
plt.ion()

fig_fontsize=14
fs =fig_fontsize
params = {
    'lam_HI':1.0,
    'lam_avi':2.0,
    'lam_pot':0.3,
    'prefix':'H3N2_',
    'serum_Kc':0.003,
}

def add_panel_label(ax,label, x_offset=-0.1):
    '''add one letter labels to the upper left corner of a figure A, B, C etc '''
    ax.text(x_offset, 0.95, label, transform=ax.transAxes, fontsize=fig_fontsize*1.5)

def select_nodes_in_season(tree, interval):
    '''mark all nodes in a time interval specified by decimalized years, e.g. 2012.34 '''
    for node in tree.leaf_iter(): # mark leafs
        if node.num_date>=interval[0] and node.num_date<interval[1]:
            node.alive=True
            node.n_alive = 1
        else:
            node.alive=False
            node.n_alive = 0
    for node in tree.postorder_internal_node_iter(): # go over all internal nodes: alive iff at least one child alive
        node.alive = any([n.alive for n in node.child_nodes()])
        node.n_alive = np.sum([n.n_alive for n in node.child_nodes()])

def calc_LBI(tree, LBI_tau = 0.0005, attr = 'lb'):
    '''
    traverses the tree in postorder and preorder to calculate the
    up and downstream tree length exponentially weighted by distance.
    then adds them as LBI
    tree -- dendropy tree for whose node the LBI is being computed
    attr     -- the attribute name used to store the result
    '''
    min_bl = 0.00005
    # traverse the tree in postorder (children first) to calculate msg to parents
    for node in tree.postorder_node_iter():
        node.down_polarizer = 0
        node.up_polarizer = 0
        for child in node.child_nodes():
            node.up_polarizer += child.up_polarizer
        bl = max(min_bl, node.edge_length)/LBI_tau
        node.up_polarizer *= np.exp(-bl)
        if node.alive: node.up_polarizer += LBI_tau*(1-np.exp(-bl))

    # traverse the tree in preorder (parents first) to calculate msg to children
    for node in tree.preorder_internal_node_iter():
        for child1 in node.child_nodes():
            child1.down_polarizer = node.down_polarizer
            for child2 in node.child_nodes():
                if child1!=child2:
                    child1.down_polarizer += child2.up_polarizer

            bl =  max(min_bl, child1.edge_length)/LBI_tau
            child1.down_polarizer *= np.exp(-bl)
            if child1.alive: child1.down_polarizer += LBI_tau*(1-np.exp(-bl))

    # go over all nodes and calculate the LBI (can be done in any order)
    max_LBI = 0
    for node in tree.postorder_node_iter():
        tmp_LBI = node.down_polarizer
        for child in node.child_nodes():
            tmp_LBI += child.up_polarizer
        node.__setattr__(attr, tmp_LBI)
        if tmp_LBI>max_LBI:
            max_LBI=tmp_LBI
    return max_LBI

''' mutation model
goes over different intervals and fits the HI model
'''
mut_models = True
save_figs = True
if mut_models:
    resolutions = ['1985to1995','1990to2000','1995to2005','2000to2010','2005to2016']
    fig, axs = plt.subplots(1,len(resolutions), sharey=True, figsize=(4*figheight, 1.3*figheight))
    cols={}
    HI_distributions_mutations = []
    #### make a plot of trajectories colored by HI effect
    for res,ax in izip(resolutions,axs):
        params['pivots_per_year'] = 6.0
        params['resolution']=res
        params['time_interval'] = map(float, res.split('to'))
        #params['time_interval'] = [2015.8-int(res[:-1]), 2015.8]

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

        cols = large_effect_mutations(myH3N2, ax, cols) # plot the mutation trajectories into a multi panel figure
        for mut in myH3N2.mutation_effects: # for each mutation, make a list of mutation, effect and frequency trajectory
            HI = myH3N2.mutation_effects[mut]
            mutlabel = mut[0]+':'+mut[1][1:]
            if mutlabel in myH3N2.frequencies["mutations"]["global"]:
                HI_distributions_mutations.append([res, mut, HI, np.array(myH3N2.frequencies["mutations"]["global"][mutlabel])])
            else:
                print("no frequencies for ",mut, 'HI', HI)
                continue
        print(len(HI_distributions_mutations))
    if save_figs:
        plt.savefig("prediction_figures/"+"trajectories_mutations.pdf")

    ### make cumulative distribution of HI titers that fix or don't
    freq_thres = 0.5    # minimal freq
    HI_threshold =0.1   # minimal HI effect
    fixed  = np.array([ HI for res, mut, HI, freq in HI_distributions_mutations
                      if freq[0]<0.1 and freq.max()>freq_thres and HI>HI_threshold]) # condition on initially rare
    failed = np.array([ HI for res, mut, HI, freq in HI_distributions_mutations
                      if freq[0]<0.1 and freq.max()<freq_thres and HI>HI_threshold])

    D, p = ks_2samp(fixed, failed)
    print("HI distribution of fixed and failed, KS stat:", D, "p-val:",p)
    # plt cumulative distributions
    plt.figure()
    plt.plot(sorted(fixed), np.linspace(0,1,len(fixed)), label = '>'+str(freq_thres)+' n='+str(len(fixed)))
    plt.plot(sorted(failed), np.linspace(0,1,len(failed)), label = '<'+str(freq_thres)+' n='+str(len(failed)))
    plt.xlabel('HI effect')
    plt.ylabel('cumulative distribution')
    plt.legend(loc=4)
    if save_figs:
        plt.savefig("prediction_figures/"+"cumulative_HI_mutations.pdf")

    ################################################################
    ##### fraction successful
    ################################################################

    plt.figure(figsize = (1.6*figheight, figheight))
    ax3 = plt.subplot(1,1,1)
    HI_max = np.array([[HI, freq.max()] for res, mut, HI, freq in HI_distributions_mutations if freq[0]<0.1 and freq.max()>0.1])
    nreps=100
    HI_threshold = np.array([0.0, 0.3, 0.7, 1.2, 4]) #defining HI categories
    HI_binc = 0.5*(HI_threshold[:-1] + HI_threshold[1:])
    for fi,freq_thres in enumerate([0.25, 0.5, 0.75, 0.95]):
        frac_success = []
        stddev_success = []
        for HI_lower, HI_upper in zip(HI_threshold[:-1], HI_threshold[1:]):
            ind = (HI_max[:,0]>=HI_lower)&(HI_max[:,0]<HI_upper)
            vals = HI_max[ind,1]
            tmp = []
            for rep in xrange(nreps):
                tmp_vals = vals[np.random.randint(len(vals), size=len(vals)/2)]
                tmp.append((tmp_vals>freq_thres).mean())
            stddev_success.append(np.std(tmp))
            print(HI_lower, ind.sum())
            frac_success.append((HI_max[ind,1]>freq_thres).mean())
        ax3.errorbar(np.arange(len(frac_success))+0.5+0.03*fi, frac_success,stddev_success, label = "max freq >"+str(freq_thres), lw=2)

    ax3.set_xlabel('HI effect', fontsize=fs)
    ax3.set_ylabel('fraction reaching frequency threshold', fontsize=fs)
    ax3.tick_params(labelsize=fs)
    ax3.set_xticks(np.arange(len(HI_binc))+0.5)
    ax3.set_xticklabels([str(lower)+'-'+str(upper) for lower, upper in zip(HI_threshold[:-1], HI_threshold[1:])])
    plt.legend(loc=8, fontsize=fs)
    plt.ylim([0,1])
    plt.xlim([0,len(HI_binc)])

    plt.tight_layout()
    if save_figs:
        plt.savefig("prediction_figures/"+'fraction_successful.pdf')

    ### make cumulative HI on backbone
    HI_backbone = np.array([HI for res, mut, HI, freq in HI_distributions_mutations if freq[0]<0.1 and freq.max()>0.75])
    HI_backbone.sort()
    plt.figure(figsize = (1.2*figheight, figheight))
    cumHI =  HI_backbone.cumsum()
    plt.plot(HI_backbone, cumHI/cumHI[-1])
    plt.ylabel('fraction of cHI due to effects < cutoff', fontsize = fs)
    plt.xlabel('effect size', fontsize = fs)
    plt.tick_params(labelsize=fs)
    plt.tight_layout()
    if save_figs:
        plt.savefig("prediction_figures/"+'cumulative_HI_effects.pdf')

''' analyze the tree model
This uses the 30 year span and investigates whether antigenic advance (as measured by cHI)
is predictive of clade success.
'''
res = '1985to2016'
params['pivots_per_year'] = 3.0
params['resolution']=res
params['time_interval'] = map(float, res.split('to'))

if params['time_interval'][1]>2015:
    params['time_interval'][1]=2015.8 # set the upper time limit to past the last sequence

# add all arguments to virus_config (possibly overriding)
virus_config.update(params)
# pass all these arguments to the processor: will be passed down as kwargs through all classes
myH3N2 = H3N2_process(**virus_config)
myH3N2.load()
# assign dates to internal nodes as minimum for children. this is needed to calculate the recent
for node in myH3N2.tree.postorder_internal_node_iter():
    node.num_date = np.min([c.num_date for c in node.child_nodes()])

assign_fitness(myH3N2.tree)
dates_fitness = np.array([(n.num_date, n.tol) for n in myH3N2.tree.postorder_internal_node_iter()])
dates_fitness_term = np.array([(n.num_date, n.tol) for n in myH3N2.tree.leaf_iter()])

pivots = myH3N2.tree.seed_node.pivots
HI_vs_max_freq_tree = []
dt=0.5
for node in myH3N2.tree.postorder_internal_node_iter():
    if node.num_date<1987:
        continue
    if node.freq["global"] is not None and node.freq["global"].max()>0.05:
        p = node.parent_node
        cHI = node.dHI
        while  p is not None and (node.num_date - p.num_date)<dt:
            cHI += p.dHI
            p = p.parent_node

        ind = (dates_fitness_term[:,0]<=node.num_date)&(dates_fitness_term[:,0]>node.num_date-dt)
        HI_vs_max_freq_tree.append((cHI, np.array(node.freq["global"])))

freq_clusters = []
globbing_thres=0.25
print("combining trajectories")
for cHI, freq in HI_vs_max_freq_tree:
    found = False
    for fi, (cHIs, cfreqs) in enumerate(freq_clusters):
        max_freq = np.mean(cfreqs, axis=0).max()+0.1
        if np.max(np.abs(freq - np.mean(cfreqs, axis=0)))/max_freq<globbing_thres:
            freq_clusters[fi][1].append(freq)
            freq_clusters[fi][0].append(cHI)
            found=True
    if not found:
        print ("adding new cluster")
        freq_clusters.append([[cHI], [freq]])
print("generated",len(freq_clusters), " trajectory clusters")

freq_thres = 0.75
fixed  = np.array([ max(HI) for HI, freqs in freq_clusters
                  if np.max(np.mean(freqs, axis=0))>freq_thres and max(HI)>0.01])
failed  = np.array([ max(HI) for HI, freqs in freq_clusters
                  if np.max(np.mean(freqs, axis=0))<freq_thres and max(HI)>0.01])

print("split into lists of failed and successful trajectories")

D, p = ks_2samp(fixed, failed)
print("KS stat:", D, "p-val:",p)
plt.figure()
plt.plot(sorted(fixed), np.linspace(0,1,len(fixed)), label = '>'+str(freq_thres)+' n='+str(len(fixed)))
plt.plot(sorted(failed), np.linspace(0,1,len(failed)), label = '<'+str(freq_thres)+' n='+str(len(failed)))
plt.xlabel('HI effect')
plt.ylabel('cumulative distribution')
plt.legend(loc=4)
if save_figs:
    plt.savefig("prediction_figures/"+"cumulative_HI_tree.pdf")

################################################################
#### plot tree frequencies
################################################################
fs=14
HI_cutoff=0.1
mycmap = cm.cool
plt.figure(figsize=(3*figheight, figheight))
ax1 = plt.subplot(1,1,1)
for cHIs, cfreqs in freq_clusters:
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
ax1.set_xlim([1985,2016])
ax1.tick_params(labelsize=fs)
plt.tight_layout()
add_panel_label(ax1, "A", x_offset=-0.07)
if save_figs:
    plt.savefig("prediction_figures/"+"trajectories_tree.pdf")


'''
the following is obsolete code that was used to plot the fraction of successful clades
vs their antigenic effect
'''
################################################################
##### add fraction successful
################################################################
#plt.figure(figsize=(2.4*figheight, figheight))
#ax2 = plt.subplot2grid((1,2),(0,0))
#plt.title("tree model", fontsize=fs)
##HI_threshold = np.array([0.1, 0.3, 0.8, 1.5, 4])
#HI_threshold = np.array([0.0, 0.3, 0.8, 1.5, 4])
#HI_binc = 0.5*(HI_threshold[:-1]+HI_threshold[1:])
#HI_max = np.array([[np.max(HI), np.max(np.mean(freqs, axis=0))] for HI, freqs in freq_clusters])
#for freq_thres in [0.5, 0.75, 0.95]:
#    frac_success = []
#    for HI_lower, HI_upper in zip(HI_threshold[:-1], HI_threshold[1:]):
#        ind = (HI_max[:,0]>=HI_lower)&(HI_max[:,0]<HI_upper)
#        print(HI_lower, ind.sum())
#        frac_success.append((HI_max[ind,1]>freq_thres).mean())
#    ax2.plot(np.arange(len(frac_success))+0.5, frac_success, 'o-', label = "max freq >"+str(freq_thres))
#
#ax2.set_xlabel('HI effect', fontsize=fs)
#ax2.set_ylabel('fraction reaching frequency threshold', fontsize=fs)
#ax2.tick_params(labelsize=fs)
#ax2.set_xticks(np.arange(len(HI_binc))+0.5)
#ax2.set_xticklabels([str(lower)+'-'+str(upper) for lower, upper in zip(HI_threshold[:-1], HI_threshold[1:])])
#plt.legend(loc=4, fontsize=fs)
#plt.ylim([0,1])
#plt.xlim([0,len(HI_binc)])
#
#ax3 = plt.subplot2grid((1,2),(0,1))
#plt.title("mutation model", fontsize=fs)
#HI_max = np.array([[HI, freq.max()] for res, mut, HI, freq in HI_distributions_mutations if freq[0]<0.1 and freq.max()>0.01])
#for freq_thres in [0.25, 0.5, 0.75, 0.95]:
#    frac_success = []
#    for HI_lower, HI_upper in zip(HI_threshold[:-1], HI_threshold[1:]):
#        ind = (HI_max[:,0]>=HI_lower)&(HI_max[:,0]<HI_upper)
#        print(HI_lower, ind.sum())
#        frac_success.append((HI_max[ind,1]>freq_thres).mean())
#    ax3.plot(np.arange(len(frac_success))+0.5, frac_success, 'o-', label = "max freq >"+str(freq_thres))
#
#ax3.set_xlabel('HI effect', fontsize=fs)
#ax3.set_ylabel('fraction reaching frequency threshold', fontsize=fs)
#ax3.tick_params(labelsize=fs)
#ax3.set_xticks(np.arange(len(HI_binc))+0.5)
#ax3.set_xticklabels([str(lower)+'-'+str(upper) for lower, upper in zip(HI_threshold[:-1], HI_threshold[1:])])
#plt.legend(loc=4, fontsize=fs)
#plt.ylim([0,1])
#plt.xlim([0,len(HI_binc)])
#
#plt.tight_layout()
#if save_figs:
#    plt.savefig("prediction_figures/"+'combined_HI_dynamics.pdf')
#
#

#########################################################################
##### the following implements bona fide prediction by estimating HI models
##### and LBI for viruses sampled up to a time cutoff and examining the next season
#########################################################################

gof_by_year = []

alpha = 'ACGT'
def allele_freqs(seqs):
    ''' alignment -> nucleotide frequencies '''
    tmp_seqs = np.array([np.fromstring(seq, 'S1') for seq in seqs])
    af = np.zeros((4,tmp_seqs.shape[1]))
    for ni,nuc in enumerate(alpha):
        af[ni,:] = (tmp_seqs==nuc).mean(axis=0)
    return af

def af_dist(af1, af2):
    ''' average distance between two population characterized by nucleotide frequencies af1 and af2'''
    return 1-(af1*af2).sum(axis=0).mean(axis=0) # distance == 1-probability of being the same

def seq_dist(seq, af):
    ''' average distance betweeen a populiation characterized by nucleotide frequencies and a sequence'''
    ind = np.array([alpha.index(nuc) for nuc in seq]) #calculate the indices in the frequency array that correspond to the sequence state
    return 1.0-np.mean(af[ind, np.arange(len(seq))])  #1 - probability of being the same

# frequency cutoffs for a clade to be included
cutoffs = [0.01, 0.03, 0.05, 0.1]
# set up dictionaries to remember the highest scoring clades for sets defined by each of the cutoffs
LBI_HI_by_date = {cutoff:[] for cutoff in cutoffs}
best_scores = {cutoff:[] for cutoff in cutoffs}
best_LBI = {cutoff:[] for cutoff in cutoffs}
best_HI = {cutoff:[] for cutoff in cutoffs}
best_HI_vs_HI_of_best = {cutoff:[] for cutoff in cutoffs}

# loop over all years we want to include
for year in range(1990, 2015):
    print("#############################################")
    print("### YEAR:",year)
    print("#############################################")
    # train the HI model and remember some basic figures about the fit
    myH3N2.map_HI(training_fraction = 1.0, method = 'nnl1reg', lam_HI=params['lam_HI'], map_to_tree = True,
            lam_pot = params['lam_pot'], lam_avi = params['lam_avi'], cutoff_date = year+2./12.0, subset_strains = False, force_redo = True)

    gof_by_year.append((year, myH3N2.fit_error, myH3N2.tree_graph.shape[0]))

    # take and allele frequency snapshot of future season Sept until June
    select_nodes_in_season(myH3N2.tree, (year+9.0/12, year+18.0/12))
    future_seqs = [node.seq for node in myH3N2.tree.leaf_iter() if node.alive]
    future_af = allele_freqs(future_seqs)

    #current season May until Feb (all previously selected nodes will be erased)
    select_nodes_in_season(myH3N2.tree, (year-7.0/12, year+2.0/12))
    af = allele_freqs([node.seq for node in myH3N2.tree.leaf_iter() if node.alive])
    avg_dist = af_dist(af, future_af)
    max_LBI = calc_LBI(myH3N2.tree, LBI_tau = 0.001)
    total_alive = 1.0*myH3N2.tree.seed_node.n_alive
    # loop over the different frequency cut-offs
    for cutoff in cutoffs:
        # make a list of nodes (clades) that are used for prediction. frequency needs to be >cutoff and <0.95
        nodes = [node for node in myH3N2.tree.postorder_node_iter()
                if node.alive and node.n_alive/total_alive>cutoff and node.n_alive/total_alive<0.95]
        # determine the minimal distance to future, the average cumulative antigenic advance, and the best possible node
        all_distance_to_future = np.array([seq_dist(node.seq, future_af) for node in nodes])
        best = np.argmin(all_distance_to_future)
        min_dist = all_distance_to_future[best]
        current_cHI = np.mean([n.cHI for n in nodes])
        # determine the nodes with the highest LBI and the highest HI
        all_LBI = np.array([n.lb for n in nodes])
        best_LBI_node = nodes[np.argmax(all_LBI)]
        best_HI_node = nodes[np.argmax([n.cHI for n in nodes])]

        # remember the LBI and HI of the best node
        best_scores[cutoff].append([year, nodes[best].lb/max_LBI, nodes[best].cHI-current_cHI])
        # remember the LBI and HI, normalized (d/avg_d), standardized (d-min_d)/(avg_d-min_d), and min_d, avg_d for node with highest LBI
        best_LBI[cutoff].append((year, best_LBI_node.lb/max_LBI, best_LBI_node.cHI-current_cHI,
                        (seq_dist(best_LBI_node.seq, future_af)-min_dist)/(avg_dist-min_dist),
                        seq_dist(best_LBI_node.seq, future_af)/avg_dist, min_dist, avg_dist))
        # remember the LBI and HI, normalized (d/avg_d), standardized (d-min_d)/(avg_d-min_d), and min_d, avg_d for node with highest HI
        best_HI[cutoff].append((year, best_HI_node.lb/max_LBI, best_HI_node.cHI-current_cHI,
                        (seq_dist(best_HI_node.seq, future_af)-min_dist)/(avg_dist-min_dist),
                        seq_dist(best_HI_node.seq, future_af)/avg_dist, min_dist, avg_dist))
        # remember HI of the node closest to the future and the HI of the node iwth the highest HI
        best_HI_vs_HI_of_best[cutoff].append((year, best_HI_node.cHI - current_cHI, nodes[best].cHI - current_cHI))
        print(year, "avg_dist", avg_dist)
        # make a list of the LBI, cHI, and distances for every node in the set belonging to cutoffs. (used for scattering LBI vs HI)
        for node in nodes:
            node_freq = node.n_alive/total_alive
            if node.freq["global"] is not None:
                tmp_freq = np.array(node.freq["global"])
                ii = pivots.searchsorted(year+0.2)
                nii= pivots.searchsorted(year+1.0)
                LBI_HI_by_date[cutoff].append((node, year, node.lb/max_LBI, node.cHI-current_cHI,
                                            (seq_dist(node.seq, future_af)-min_dist)/(avg_dist-min_dist),node_freq,
                                            tmp_freq[ii], tmp_freq[nii], tmp_freq))
            else:
                #print("missing frequency", year, node.n_alive)
                LBI_HI_by_date[cutoff].append((node, year, node.lb/max_LBI, node.cHI-current_cHI,
                                            (seq_dist(node.seq, future_af)-min_dist)/(avg_dist-min_dist), node_freq,
                                            0,0,0))

# make an array out of all values for slicing and plotting
clades = {}
for cutoff in cutoffs:
    best_scores[cutoff] = np.array(best_scores[cutoff])
    best_LBI[cutoff] = np.array(best_LBI[cutoff])
    best_HI[cutoff] = np.array(best_HI[cutoff])
    best_HI_vs_HI_of_best[cutoff] = np.array(best_HI_vs_HI_of_best[cutoff])
    tmp = []
    for cl in LBI_HI_by_date[cutoff]:
        tmp.append(cl[1:-1]) # exclude node (entry 0) and frequencies (entry -1) since they aren't numbers
    clades[cutoff] = np.array(tmp)

############
# make figure that shows the distance to future season for all years
# comparing LBI and HI at different cutoffs
############
cutoff = 0.01
plt.figure(figsize=(2.4*figheight,figheight))
ax=plt.subplot(121)
# -2 == mindist, -1 == avg_dits -3 == dist/avg_dist
ax.plot(best_HI[cutoff][:,0],best_HI[cutoff][:,-2]/best_HI[cutoff][:,-1], label='best',lw=2, c='k')
ax.plot(best_LBI[cutoff][:,0],best_LBI[cutoff][:,-3], label='LBI',lw=2)
ax.plot(best_HI[0.05][:,0],best_HI[0.05][:,-3], label='cHI >0.05',lw=2)
ax.plot(best_HI[0.01][:,0],best_HI[0.01][:,-3], label='cHI >0.01',lw=2)
ax.plot([1990, 2015], [1.0, 1.0], lw=3, c='k', ls='--')
ax.tick_params(labelsize=fs)
add_panel_label(ax, "B", x_offset=-0.12)
ax.set_xlabel('year', fontsize=fs)
ax.set_ylabel('distance to season year/year+1', fontsize=fs)
ax.set_yticks([0,0.5, 1.0, 1.5])
plt.legend(loc=2, fontsize=fs)

############
# second panel showing the distribytion of HI values of the best nodes
############
cols = sns.color_palette(n_colors=5)
symbs = ['o','d','v','^','<']
ax = plt.subplot(122)
ax.hist(best_scores[0.05][:,-1])
add_panel_label(ax, "C")
ax.set_yticks([0,2,4,6, 8])
ax.set_ylim([0,10])
ax.set_ylabel("#years", fontsize=fs)
ax.set_xlabel(r'$cHI-\langle cHI\rangle_{year}$', fontsize=fs)
ax.tick_params(labelsize=fs)

plt.tight_layout()
if save_figs:
    plt.savefig("prediction_figures/"+'LBI_and_HI_vs_distance.pdf')

#fig, axs = plt.subplots(1,3, figsize=(3.0*figheight, figheight))
#lbi_cutoff = 0.2
#for ax, qi in izip(axs,[1,2]):
#    for yi,year in enumerate(range(1990,2015)):
#        ind = (clades[:,0]==year)&((clades[:,-3]>lbi_cutoff)|(clades[:,2]>.5)) #restrict to clades larger than cutoff
#        if ind.sum()==0:
#            continue
#        lstr = str(year) if (year<1998 and qi==1) or (year>=1998 and year<2006 and qi==2) else None
#        ax.scatter(clades[ind,qi], clades[ind,3], c=cols[yi%5], marker=symbs[yi//5], s=50, label=lstr)
#        print(cols[yi%5])
#    x_str = r'$cHI-\langle cHI\rangle_{year}$' if qi==2 else r'$LBI/\max(LBI)$'
#    ax.set_xlabel(x_str, fontsize=fs)
#    ax.tick_params(labelsize=fs)
#    ax.set_xlim((0.2,1.4) if qi==1 else (-3,3))
#    ax.set_xticks((0.25,0.5, 0.75, 1.0) if qi==1 else (-2,-1,0,1,2))
#    ax.set_ylim((-0.2,2.5))
#    if qi<3:
#        ax.set_yticks([0, 0.5,1.0, 1.5, 2.0])
#        ax.set_ylabel(r'distance to season year/year+1', fontsize=fs)
#    ax.legend()
#    add_panel_label(ax, "C" if qi==2 else "B", x_offset=-0.12)
#
#ax = axs[2]
#for yi, (year, lbi, cHI) in enumerate(best_scores):
#    lstr = str(year) if (year>=2006) else None
#    ax.scatter([lbi], [cHI], c=cols[yi%5], marker=symbs[yi//5], s=50, label=lstr)
#ax.set_xlabel(r'$LBI/\max(LBI)$', fontsize=fs)
#ax.set_ylabel(r'$cHI-\langle cHI\rangle_{year}$', fontsize=fs)
#ax.set_xlim([0, 1.1])
#ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
#ax.set_yticks([-0.5, 0, 0.5, 1, 1.5])
#ax.legend()

#################################################################
### plot best HI vs HI of best
################################################################
plt.figure(figsize = (1.2*figheight, figheight))
ax=plt.subplot(111)
for col, cutoff in zip(['b', 'g'], [0.01, 0.05]):
    plt.scatter(best_HI_vs_HI_of_best[cutoff][:,1],
        best_HI_vs_HI_of_best[cutoff][:,2], label = '>'+str(cutoff), s=50, c=col) #, s=50*best_HI[:,-3])
plt.plot([0,3],[0,3])
plt.tick_params(labelsize=fs)
plt.xlabel(r'maximal $cHI-\langle cHI\rangle_{year}$', fontsize=fs)
plt.ylabel(r'successful $cHI-\langle cHI\rangle_{year}$', fontsize=fs)
plt.xticks([0,1,2,3,4])
plt.yticks([-1, 0,1,2,3])
plt.legend(loc=2)
plt.tight_layout()
if save_figs:
    plt.savefig("prediction_figures/"+'best_HI_vs_HI_of_best.pdf')


################################################################
### scatter plot of LBI vs HI
################################################################
plt.figure(figsize=(2.4*figheight, figheight))
mycmap=cm.Set1
cutoff = 0.01
for li,lbi_cutoff in enumerate([0.2, 0.1]):
    ax = plt.subplot(1,2,li+1)
    ind = clades[cutoff][:,-3]>lbi_cutoff #restrict to clades larger than cutoff
    if ind.sum()==0:
        continue
    ax.set_title('clades >'+str(lbi_cutoff)+' frequency')
    ax.scatter(clades[cutoff][ind,1], clades[cutoff][ind,2], c=mycmap((clades[cutoff][ind,0]-1990)/25.0), s=80*(1-clades[cutoff][ind,3])**2)
    ax.set_ylabel(r'$cHI-\langle cHI\rangle_{year}$', fontsize=fs)
    ax.set_xlabel(r'$LBI/\max(LBI)$', fontsize=fs)
    ax.tick_params(labelsize=fs)
    ax.set_yticks([-3,-2,-1,0,1,2])
    add_panel_label(ax, "C" if li else "B", x_offset=-0.15)
    if li:
        sm = plt.cm.ScalarMappable(cmap=mycmap, norm=plt.Normalize(vmin=1990, vmax=2015))
        sm._A = []
        cb = plt.colorbar(sm)
        cb.set_ticks([1990,1995,2000, 2005, 2010, 2015])
        cb.set_label('year', fontsize=fs)

plt.tight_layout()
if save_figs:
    plt.savefig("prediction_figures/"+'LBI_HI.pdf')

