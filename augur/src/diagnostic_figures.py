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
######################################
#### make list of mutation effects across different periods
######################################
def mutation_list(flu = 'H3N2', format='tsv', nconstraints=False):
    from glob import glob
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

######################################
#### make a figure that plots the cumulative antigenic change vs time
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
#### make a figure that compares
#### trunk frequency increase to dHI
######################################
def slope_vs_dHI(myflu):
    plt.figure(figsize=(1.6*figheight, figheight))
    ax = plt.subplot('121')
    slopes = []
    cutoff_freq = 0.25
    tmp_pivots = myflu.tree.seed_node.pivots
    dt = tmp_pivots[1]-tmp_pivots[0]
    y1 = int(1.0/dt)
    y3 = int(3.0/dt)
    for node in myflu.tree.postorder_internal_node_iter():
        tmp_freq = node.freq['global']
        if tmp_freq is not None and tmp_freq[0]<cutoff_freq and np.max(tmp_freq)>cutoff_freq:
            ii, offset, slope = get_slope(tmp_freq, tmp_pivots, cutoff_freq)
            dfreq_sum = tmp_freq[ii:ii+y3].sum()
            if ii+y3<len(tmp_freq):
                dfreq_y1 = tmp_freq[ii+y1]-tmp_freq[ii]
                dfreq_y3 = tmp_freq[ii+y3]-tmp_freq[ii]
            else:
                dfreq_y1 = 0
                dfreq_y3 = 0
            HI=node.dHI
            slopes.append([HI, slope, np.max(tmp_freq), dfreq_y1, dfreq_y3, dfreq_sum])
            plt.plot(tmp_pivots-offset, tmp_freq, lw=2, c=cm.jet(np.sqrt(HI)/2.0), alpha=min(1.0, max(HI, 0.3)))

    plt.xlabel('time', fontsize=fs)
    plt.ylabel('frequency', fontsize=fs)
    plt.xlim([-1,2])
    plt.ylim([-0.01,1.1])
    ax.tick_params(axis='both', labelsize=fs)
    slopes = np.array(slopes)

    ax = plt.subplot('122')
    plt.scatter(slopes[:,1], slopes[:,0])
    plt.xlabel('frequency slope', fontsize=fs)
    plt.ylabel('titer drop', fontsize=fs)
    ax.tick_params(axis='both', labelsize=fs)
    plt.tight_layout()
    from scipy.stats import spearmanr
    print("spearman correlation HI/slope:",spearmanr(slopes[:,:2]))
    print("spearman correlation HI/max:",spearmanr(slopes[:,0],slopes[:,2]))
    print("spearman correlation HI/delta1Y:",spearmanr(slopes[:,0],slopes[:,3]))
    print("spearman correlation HI/delta3Y:",spearmanr(slopes[:,0],slopes[:,4]))
    print("spearman correlation HI/freqsum:",spearmanr(slopes[:,0],slopes[:,5]))
    return slopes

def get_slope(freq, pivots, threshold):
    ii = np.argmax(freq>threshold)
    slope = (freq[ii]-freq[ii-1])/(pivots[ii]-pivots[ii-1])
    offset = pivots[ii-1] + (threshold-freq[ii-1])/slope
    return ii, offset, slope

######################################
#### make a figure that compares
#### trunk frequency increase to dHI
######################################
def slope_vs_mutation(myflu):
    plt.figure(figsize=(1.6*figheight, figheight))
    ax = plt.subplot('121')
    slopes = []
    cutoff_freq = 0.25
    pivots = myflu.tree.seed_node.pivots
    dt = pivots[1]-pivots[0]
    y3 = int(3.0/dt)
    y1 = int(1.0/dt)
    for mut in myflu.mutation_effects:
        HI = myflu.mutation_effects[mut]
        mutlabel = mut[0]+':'+mut[1][1:]
        prevmutlabel = mut[0]+':'+mut[1][1:-1]+mut[1][0]
        if mutlabel in myflu.frequencies["mutations"]["global"] and\
           prevmutlabel in myflu.frequencies["mutations"]["global"]:
            mut_freq = np.array(myflu.frequencies["mutations"]["global"][mutlabel])
            prevmut_freq = np.array(myflu.frequencies["mutations"]["global"][prevmutlabel])
        else:
            continue
        # find sweeps
        cutoff = 0.25 # minimal change over 3 years
        dfreq = mut_freq[y3:]-mut_freq[:-y3]
        #tmp_sweeps = [ii for ii,dx  in enumerate(dfreq[:-1]) if dx>cutoff and dfreq[ii+1]<cutoff]
        tmp_sweeps = [ii for ii,dx  in enumerate(mut_freq[:-y3]) if dx<cutoff and mut_freq[ii+1]>cutoff]
        sweeps=[]
        for ii, si in enumerate(tmp_sweeps):
            # check that one mutation goes up, the prev down
            if np.min((mut_freq+prevmut_freq)[si:si+y3])>0.8:
                #if ii==0 or np.min(mut_freq[tmp_sweeps[ii-1]:si])<0.2:
                sweeps.append(si)
        for si in sweeps:
            tmp_freq = mut_freq[max(0,si-y3):(si+2*y3)]
            tmp_pivots = pivots[max(0,si-y3):(si+2*y3)]
            ii, offset, slope = get_slope(tmp_freq, tmp_pivots, cutoff_freq)
            dfreq_sum = tmp_freq[ii:ii+y3].sum()
            if ii+y3<len(tmp_freq):
                dfreq_y1 = tmp_freq[ii+y1]-tmp_freq[ii]
                dfreq_y3 = tmp_freq[ii+y3]-tmp_freq[ii]
            else:
                dfreq_y1 = 0
                dfreq_y3 = 0
            slopes.append([HI, slope, np.max(tmp_freq), dfreq_y1, dfreq_y3, dfreq_sum])
            plt.plot(tmp_pivots-offset, tmp_freq, lw=2, c=cm.jet(np.sqrt(HI)/2.0), alpha=min(1.0, max(HI, 0.3)))

    plt.xlabel('time', fontsize=fs)
    plt.ylabel('frequency', fontsize=fs)
    plt.xlim([-1,2])
    plt.ylim([-0.01,1.1])
    ax.tick_params(axis='both', labelsize=fs)
    slopes = np.array(slopes)

    ax = plt.subplot('122')
    plt.scatter(slopes[:,1], slopes[:,0])
    plt.xlabel('frequency slope', fontsize=fs)
    plt.ylabel('titer drop', fontsize=fs)
    ax.tick_params(axis='both', labelsize=fs)
    plt.tight_layout()
    from scipy.stats import spearmanr
    print("spearman correlation HI/slope:",spearmanr(slopes[:,:2]))
    print("spearman correlation HI/max:",spearmanr(slopes[:,0],slopes[:,2]))
    print("spearman correlation HI/delta1Y:",spearmanr(slopes[:,0],slopes[:,3]))
    print("spearman correlation HI/delta3Y:",spearmanr(slopes[:,0],slopes[:,4]))
    print("spearman correlation HI/freqsum:",spearmanr(slopes[:,0],slopes[:,5]))
    return slopes


######################################
#### make a figure that shows histogram of distance asymmetries and deviations from quartett tests
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
             label = 'quartett, mean='+str(np.round(np.mean(additivity_test['test']),2)))
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


