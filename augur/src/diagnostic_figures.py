from glob import glob
from itertools import izip
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('darkgrid')

colors = sns.color_palette(['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
                         '#e31a1c','#fdbf6f','#ff7f00','#cab2d6'])

fs = 12
plt.locator_params(nbins=4)
fmts = ['.pdf','.svg','.png']
figheight = 4
####################################
##### merge accession number files
####################################
def make_combined_accession_number_lists():
    for flu in ['H3N2', 'H1N1pdm', 'Vic', 'Yam']:
        flist = glob.glob('../auspice/data/'+flu+'*accession_numbers.tsv')
        all_accessions=set()
        for fname in flist:
            with open(fname) as infile:
                all_accessions.update([tuple(line.split('\t')[:2]) for line in infile])

        with open('../auspice/data/'+flu+'_all_accession_numbers.tsv', 'w') as ofile:
            for strain, acc in all_accessions:
                ofile.write(strain+'\t'+acc+'\n')




######################################
#### make list of mutation effects across different periods
######################################
def mutation_list(flu = 'H3N2', format='tsv', nconstraints=False):
    import pandas as pd
    flist = glob('../auspice/'+flu+'/*to*/HI_mutation_effects.tsv')
    mutation_effects = {}
    mutation_counts = {}
    for fname in flist:
        interval = fname.split('/')[-2]
        years = map(int, interval.split('to'))
        if years[1]-years[0] not in [10,11]: # skip all non-10year interval
            continue
        tmp_effects, tmp_counts = {},{}
        with open(fname) as infile:
            for line in infile:
                mut, val, count = line.strip().split('\t')
                tmp_effects[mut] = float(val)
                tmp_counts[mut] =  int(count)
        mutation_effects[interval]=tmp_effects
        mutation_counts[interval]=tmp_counts
    M = pd.DataFrame(mutation_effects)
    C = pd.DataFrame(mutation_counts)
    mutation_statistics = []
    for r,c in izip(M.iterrows(), C.iterrows()):
        vals = np.array([x for x in r[1] if not np.isnan(x)])
        mutation_statistics.append((r, c, vals.max(), vals.mean(), vals.std()))
    mutation_statistics.sort(key=lambda x:x[2], reverse=True)
    N = pd.DataFrame.from_items([(x[0][0],x[0][1]) for x in mutation_statistics]).transpose()
    D = pd.DataFrame.from_items([(x[1][0],x[1][1]) for x in mutation_statistics]).transpose()

    if format=='tex':
        sep = ' & '
        le = '\\\\ \n'
    else:
        sep = '\t'
        le = '\n'
    ndigits = 2
    with open('HI_mutation_effects_all.tsv', 'w') as ofile:
        ofile.write("mutation\t"+sep.join(N.columns)+le)
        for r,c in izip(N.iterrows(), D.iterrows()):
            if nconstraints:
                tmp = sep.join([r[0]]+map(str,[(round(x,ndigits),y) for x,y in izip(r[1], c[1])])).replace('nan','---')
            else:
                tmp = sep.join([r[0]]+map(str,[round(x,ndigits) for x,y in izip(r[1], c[1])])).replace('nan','---')
            print tmp
            ofile.write(tmp+le)
    return N, D

def mutation_list_by_position(flu = 'H3N2', format='tsv', nconstraints=False):
    from glob import glob
    import pandas as pd
    from collections import defaultdict
    flist = glob('../auspice/'+flu+'/*to*/HI_mutation_effects.tsv')

    mutation_effects = defaultdict(list)
    intervals = []
    for fname in flist:
        interval = fname.split('/')[-2]
        years = map(int, interval.split('to'))
        if years[1]-years[0] not in [10,11]: # skip all non-10year interval
            continue
        intervals.append(interval)
        with open(fname) as infile:
            for line in infile:
                mut, val, count = line.strip().split('\t')
                positions = map(lambda x:int(x[1:-1]), mut.split('/'))
                for pos in positions:
                    mutation_effects[pos].append((interval, mut, float(val), int(count)))

    positions_by_length = sorted(mutation_effects.items(), key=lambda x:len(x[1]), reverse=True)
    if format=='tex':
        sep = ' & '
        le = '\\\\ \n'
    else:
        sep = '\t'
        le = '\n'
    ndigits = 2
    intervals.sort()
    with open("effects_by_position.tsv", 'w') as ofile:
        for pos, muts in positions_by_length:
            avg_all = np.mean([x[2] for x in muts])
            avg_alone = np.mean([x[2] for x in muts if len(x[1].split('/'))==1])
            ofile.write("#######################\nPosition "+str(pos)+' avg: ' +str(round(avg_all,2))+' avg_alone: ' +str(round(avg_alone,2))+'\n')
            tmp = {(x[1],x[0]):x for x in muts}
            unique_muts = set([x[1] for x in muts])
            for mut in unique_muts:
                ofile.write(mut+sep)
                for interval in intervals:
                    if (mut, interval) in tmp:
                        ofile.write(str((tmp[(mut,interval)][2],tmp[(mut,interval)][3]))+sep)
                    else:
                        ofile.write('---'+sep)
                ofile.write('\n')

######################################
### make big mutation trajectory figure
######################################
def trajectory_figure(flu = 'H3N2', res = '1985to2016'):
    import cPickle as pickle
    from collections import defaultdict
    freqs = pickle.load(open('data/'+'_'.join([flu, res, 'frequencies.pkl'])))
    freqs = freqs['mutations']['global']
    freqs_by_position = defaultdict(list)
    for mut, freq in freqs.iteritems():
        if mut=='pivots':
            pivots = freq
        else:
            freqs_by_position[(mut.split(':')[0], int(mut.split(':')[1][:-1]))].append((mut,freq))
    all_muts = sorted(freqs_by_position.items(), reverse=True,
                      key = lambda x:sum([np.sum(np.abs(np.diff(y[1]))) for y in x[1]]))

    ncols=5
    nrows=6
    fig, axs = plt.subplots(nrows, ncols, sharey=True, sharex=True, figsize = (22,30))
    for mi, (pos, freqs) in enumerate(sorted(all_muts[:nrows*ncols])):
        ax = axs[mi//ncols, mi%ncols]
        for mut, freq in freqs:
            ax.plot(pivots, freq, label = mut[-1])
        ax.set_xlim([1985,2022])
        ax.legend(loc=4, title = pos[0]+' '+str(pos[1]), fontsize=10)
        if mi%ncols==0:
            ax.set_ylabel('frequency')
            ax.set_yticks([0,0.2, 0.4, 0.6, 0.8])
            ax.tick_params(labelsize=14)
        if mi//ncols==nrows-1:
            ax.set_xlabel('year')
            ax.tick_params(labelsize=14)
    plt.tight_layout(w_pad=0.01, h_pad=0.01)

######################################
#### calculate effect size distribution on the trunk of the tree
#### calculate titer variance explained using different cutoffs of trunk effects
######################################
def cumulative_antigenic(myflu, n=10):
    from random import sample
    leaf_sample = sample([leaf for leaf in myflu.tree.leaf_iter()
                          if leaf.num_date>2014], n)


    trunk_effects = []
    trunk_muts = []
    trunk_mut_effects = []
    for node in leaf_sample:
        tmp_muts = []
        tmp_effects = []
        while node.parent_node is not None:
            tmp_effects.append(node.dHI)
            tmp_muts.extend(node.aa_muts.items())
            node = node.parent_node

        tmp_muts = [x for x in tmp_muts if x[1]!='']
        tmp_mut_effects = {}
        for muts in tmp_muts:
            gene = muts[0]
            for pos in muts[1].split(','):
                tmp_mut = (gene, pos)
                if tmp_mut in myflu.mutation_effects:
                    tmp_mut_effects[tmp_mut] = myflu.mutation_effects[tmp_mut]
                else:
                    tmp_mut_effects[tmp_mut] = 0

        trunk_effects.append(tmp_effects)
        trunk_muts.append(tmp_muts)
        trunk_mut_effects.append(tmp_mut_effects)



    plt.figure()
    for eff in trunk_effects[:1]:
        print "sum of effects on trunk", np.sum(eff)
        tmp_eff = np.array(eff)
        plt.hist(tmp_eff[tmp_eff>1e-2],  bins=np.linspace(0,4,21),
                 label='fraction non-neg tree: '+str(np.round(np.mean(tmp_eff>1e-4), 2)), alpha=0.5)

    for eff in trunk_mut_effects[:1]:
        print "sum of mutation effects on trunk:", np.sum(eff.values())
        tmp_eff = np.array(eff.values())
        plt.hist(tmp_eff[tmp_eff>1e-4],  bins=np.linspace(0,4,21),
                 label='fractio non-neg mutations: '+str(np.round(np.mean(tmp_eff>1e-4), 2)), alpha=0.5)

    plt.xlabel('effect size')
    plt.ylabel('number of counts')
    plt.legend()
    plt.savefig(myflu.htmlpath()+ "trunk_effectsize_histogram.png")

    unexplained_variance = []
    cvals = np.linspace(0,2,40)
    for cutoff in cvals:
        myflu.validate(plot=False, cutoff=cutoff)
        unexplained_variance.append([cutoff,myflu.rms_error**2, np.var(myflu.validation.values())])
        print "effect cutoff:", cutoff, unexplained_variance[-1]
    unexplained_variance=np.array(unexplained_variance)
    plt.figure()
    plt.plot(unexplained_variance[:,0], unexplained_variance[:,1]/unexplained_variance[:,2])
    plt.ylabel('fraction of unexplained variance')
    plt.xlabel('effect size cutoff')




######################################
#### make a figure that shows histogram of distance asymmetries and deviations from quartet tests
######################################
def tree_additivity_symmetry(myflu, mtype='tree'):
    reciprocal_measurements = []
    reciprocal_measurements_titers = []
    for (testvir, serum) in myflu.HI_normalized:
        tmp_recip = [v for v in myflu.HI_normalized if serum[0]==v[0] and testvir==v[1][0]]
        for v in tmp_recip:
            val_fwd = myflu.HI_normalized[(testvir,serum)]
            val_bwd = myflu.HI_normalized[v]
            diff_uncorrected = val_fwd - val_bwd
            diff_corrected = (val_fwd - myflu.serum_potency[mtype][serum] - myflu.virus_effect[mtype][testvir])\
                            -(val_bwd - myflu.serum_potency[mtype][v[1]] - myflu.virus_effect[mtype][serum[0]])
            val_bwd = myflu.HI_normalized[v]
            reciprocal_measurements.append([testvir, serum, diff_uncorrected, diff_corrected])
            reciprocal_measurements_titers.append([testvir, serum, val_fwd, val_bwd,
                                                  (val_fwd - myflu.serum_potency[mtype][serum] - myflu.virus_effect[mtype][testvir]),
                                                  (val_bwd - myflu.serum_potency[mtype][v[1]] - myflu.virus_effect[mtype][serum[0]]),
                                                  ])
    plt.figure(figsize=(1.6*figheight, figheight))
    ax = plt.subplot(121)
    plt.text(0.05, 0.93,  ('tree model' if mtype=='tree' else 'mutation model'),
             weight='bold', fontsize=fs, transform=plt.gca().transAxes)
    plt.hist([x[2] for x in reciprocal_measurements],alpha=0.7, label="raw", normed=True)
    plt.hist([x[3] for x in reciprocal_measurements],alpha=0.7, label="tree", normed=True)
    plt.xlabel('distance asymmetry', fontsize=fs)
    ax.tick_params(axis='both', labelsize=fs)
    plt.legend(fontsize=fs)
    plt.tight_layout()

    ####  Analyze all cliques #######################################################
    all_reciprocal = list(set([v[1] for v in reciprocal_measurements_titers]))

    import networkx as nx
    from random import sample
    G = nx.Graph()
    G.add_nodes_from(all_reciprocal)
    for vi,v in enumerate(all_reciprocal):
        for w in all_reciprocal[:vi]:
            if ((v[0], w) in myflu.HI_normalized) and ((w[0], v) in myflu.HI_normalized):
                G.add_edge(v,w)
    print "generated graph of all cliques"
    C = nx.find_cliques(G)
    print "found cliques"
    def symm_distance(v,w):
        res =  myflu.HI_normalized[(v[0], w)] - myflu.virus_effect[mtype][v[0]] - myflu.serum_potency[mtype][w]
        res += myflu.HI_normalized[(w[0], v)] - myflu.virus_effect[mtype][w[0]] - myflu.serum_potency[mtype][v]
        return res*0.5

    additivity_test = {'test':[], 'control':[]}
    n_quartets = 1000
    for clique in C:
        if len(clique)>8:
            for i in xrange(n_quartets):
                Q = sample(clique, 4)
                dists = []
                for (a,b) in [((0,1), (2,3)),((0,2), (1,3)), ((0,3), (1,2))]:
                    dists.append(symm_distance(Q[a[0]], Q[a[1]])+symm_distance(Q[b[0]], Q[b[1]]))
                dists.sort(reverse=True)
                additivity_test['test'].append(dists[0]-dists[1])

                dists = []
                for di in range(3):
                    a,b,c,d = sample(clique,4)
                    dists.append(symm_distance(a, b)+symm_distance(c,d))
                dists.sort(reverse=True)
                additivity_test['control'].append(dists[0]-dists[1])

    ax=plt.subplot(122)
    plt.hist(additivity_test['control'], alpha=0.7,normed=True, bins = np.linspace(0,3,18),
             label = 'control, mean='+str(np.round(np.mean(additivity_test['control']),2)))
    plt.hist(additivity_test['test'], alpha=0.7,normed=True, bins = np.linspace(0,3,18),
             label = 'quartet, mean='+str(np.round(np.mean(additivity_test['test']),2)))
    ax.tick_params(axis='both', labelsize=fs)
    plt.xlabel(r'$\Delta$ top two distance sums', fontsize = fs)
    plt.legend(fontsize=fs)
    plt.tight_layout()

######################################
#### plot large effect mutations
######################################
def large_effect_mutations(myflu, ax=None, cols = None):
    from matplotlib import cm
    if ax is None:
        plt.figure(figsize=(figheight, figheight))
        ax = plt.subplot('111')
    if cols is None:
        cols = {}
    cutoff_freq = 0.25
    HI_cutoff = 0.3
    pivots = myflu.tree.seed_node.pivots
    dt = pivots[1]-pivots[0]
    y3 = int(3.0/dt)
    y1 = int(1.0/dt)
    color_cycle=0
    for mut in myflu.mutation_effects:
        HI = myflu.mutation_effects[mut]
        if HI>HI_cutoff:
            mutlabel = mut[0]+':'+mut[1][1:]
            if mutlabel in myflu.frequencies["mutations"]["global"]:
                mut_freq = np.array(myflu.frequencies["mutations"]["global"][mutlabel])
            else:
                print("no frequencies for ",mut, 'HI', HI)
                continue

            if mut_freq[0]<cutoff_freq:
                print("Plotting ",mut, 'max: ',mut_freq.max(), 'HI', HI)
                if mut not in cols:
                    cols[mut] = colors[color_cycle%len(colors)]
                    color_cycle+=1

                c = cols[mut]
                ax.plot(pivots, mut_freq, lw=2, ls = '--' if mut_freq.max()>0.9 else '-',c=cm.cool(min(np.sqrt(HI-HI_cutoff*0.8),1.5)/1.5))

    ax.set_xlabel('time', fontsize=fs)
    ax.set_ylabel('frequency', fontsize=fs)
    #plt.xlim([-1,2])
    ax.set_ylim([-0.01,1.1])
    ax.tick_params(axis='both', labelsize=fs)
    return cols

######################################
#### analyze correlations between titer distances and sequence/model distances
######################################
def titer_vs_distances(myflu, mtype='tree'):
    from scipy.stats import pearsonr
    from matplotlib import cm
    dists = []
    for (test, serum), val in myflu.HI_normalized.iteritems():
        muts = myflu.get_mutations(serum[0], test)
        ref_node = myflu.node_lookup[serum[0]]
        test_node = myflu.node_lookup[test]
        ref_aaseq = myflu.get_total_peptide(ref_node)
        test_aaseq = myflu.get_total_peptide(test_node)

        hamming_dist = np.sum(np.fromstring(ref_node.seq, 'S1')!=np.fromstring(test_node.seq, 'S1'))
        epidist = myflu.epitope_distance(ref_aaseq, test_aaseq)
        nonepidist = myflu.nonepitope_distance(ref_aaseq, test_aaseq)
        rbsdist = myflu.receptor_binding_distance(ref_aaseq, test_aaseq)
        correction = myflu.virus_effect[mtype][test] + myflu.serum_potency[mtype][serum]
        if mtype=='tree':
            dHI = myflu.predict_HI_tree(test, serum) - correction
        else:
            dHI = myflu.predict_HI_mutations(test, serum) - correction
        val_corrected = val - correction

        dists.append([hamming_dist, len(muts), epidist, rbsdist, dHI, val_corrected, val])

    dists = np.array(dists)
    fig, axs = plt.subplots(1,3,figsize=(3*figheight, figheight))

    for ai, (ax, ci, col_name, gridsize) in enumerate(izip(axs, [2, 3, 4],
                                            ['epitope distance',
                                            'RBS distance', mtype+' model'], [20,6,20])):

        C = pearsonr(dists[:,ci], dists[:,-2])
        ax.hexbin(dists[:,ci], dists[:,-2], gridsize=gridsize, cmap = cm.YlOrRd_r, bins='log') #,s=1, label = u'$r='+str(np.round(C[0],2))+'$')
        ax.tick_params(axis='both', labelsize=fs)
        ax.set_xlabel(col_name, fontsize=fs)
        ax.set_title(u'$r='+str(np.round(C[0],2))+'$', fontsize=fs)
        #ax.legend(loc=4, fontsize=fs)
    axs[0].set_ylabel('HI genetic component', fontsize=fs)
    plt.tight_layout(pad=0.3, w_pad=0.5)
    return dists


