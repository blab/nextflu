import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from datetime import datetime
import json
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter

sns.set_style('darkgrid')
plt.ion()

virus = 'H3N2'
#virus = 'H1N1pdm'
#virus = 'Vic'
#virus = 'Yam'


if virus=='H3N2': ########## H3N2
    freqs = json.load(open('../data/H3N2_2y_frequencies.json'))
    clades = ['3c2.a', '3c3.a', '3c3.b']
    mutations = ['HA1:114T', 'HA1:142K', 'HA1:168V', 'HA1:171K']
    #mutations = ['HA1:94H', 'HA1:114T', 'HA1:142K', 'HA1:171K']
    clade_legend = {'panel':0, 'loc':3}
    mut_legend = {'panel':0, 'loc':3}
elif virus=='H1N1pdm': ########## H1N1pdm
    freqs = json.load(open('../data/H1N1pdm_2y_frequencies.json'))
    clades = ['6b.1', '6b.2']
    mutations = ['HA1:84N','HA1:162N','HA1:152T'] #these don't add up to one in Asia, probably due to sketchy sampling.
    clade_legend = {'panel':0, 'loc':3}
    mut_legend = {'panel':0, 'loc':3}
elif virus=='Vic':
    freqs = json.load(open('../data/Vic_2y_frequencies.json'))
    clades = []
    mutations = ['HA1:129D', 'HA1:117V'] # HA1:56K would be good, but it currently isn't computed -> need to lower the threshold.
    clade_legend = {'panel':0, 'loc':3}
    mut_legend = {'panel':0, 'loc':3}
elif virus=='Yam':
    freqs = json.load(open('../data/Yam_2y_frequencies.json'))
    clades = ['2', '3']
    mutations = ['HA1:172Q', 'HA1:251V']
    clade_legend = {'panel':0, 'loc':3}
    mut_legend = {'panel':0, 'loc':3}


offset = datetime(2000,1,1).toordinal()
pivots = [offset+(x-2000)*365.25 for x in  freqs['clades']['global']['pivots']]
regions = ['global', 'NA', 'AS', 'EU', 'OC']
region_label = {'global': 'Global', 'NA': 'N America', 'AS': 'Asia', 'EU': 'Europe', 'OC': 'Oceania'}
cols = sns.color_palette(n_colors=len(regions))
fs=12
months = MonthLocator(range(1, 13), bymonthday=1, interval=3)
monthsFmt = DateFormatter("%b %y")

if len(clades):
    fig, axs = plt.subplots(len(clades), 1, sharex=True, figsize=(8, len(clades)*2))
    for clade, ax in zip(clades, axs):
        for c,region in zip(cols, regions):
            try:
                tmp_freq = freqs['clades'][region][clade]
                if tmp_freq is not None:
                    ax.plot_date(pivots, tmp_freq,'-o', label = region_label[region], c=c, lw=3 if region=='global' else 1)
            except:
                print "skipping", clade, region
        ax.set_xlim([pivots[-1]-700,pivots[-1]+30])
        ax.set_ylim(0,1)
        ax.text(pivots[-1]-700, 0.9, clade)
        ax.tick_params(labelsize=fs)
        ax.set_yticklabels(['{:3.0f}%'.format(x*100) for x in [0, 0.2, 0.4, 0.6, 0.8, 1.0]])
        ax.xaxis.set_major_locator(months)
        ax.xaxis.set_major_formatter(monthsFmt)
    fig.autofmt_xdate()
    fax = fig.add_axes( [0., 0., 1, 1] )
    fax.set_axis_off()
    fax.set_xlim(0, 1)
    fax.set_ylim(0, 1)
    fax.text(0.02, 0.54, "Frequency", rotation='vertical', horizontalalignment='center', verticalalignment='center')
    axs[clade_legend['panel']].legend(loc=clade_legend['loc'], ncol=1, bbox_to_anchor=(1.0, 0.2))
    #plt.tight_layout(h_pad=0.01)
    bottom_margin = 0.2 - 0.03*len(clades)
    plt.subplots_adjust(left=0.12, right=0.84, top=0.96, bottom=bottom_margin)
    plt.savefig('figures/feb-2016/'+virus+'_clades.png')


fig, axs = plt.subplots(len(mutations), 1, sharex=True, figsize=(8, len(mutations)*2))
for mutation, ax in zip(mutations, axs):
    for c,region in zip(cols, regions):
        try:
            tmp_freq = freqs['mutations'][region][mutation]
            if tmp_freq is not None:
                ax.plot_date(pivots, tmp_freq, '-o', label = region_label[region], c=c, lw=3 if region=='global' else 1)
        except:
            print "skipping", mutation, region
    ax.set_xlim([pivots[-1]-700,pivots[-1]+30])
    ax.set_ylim(0,1)
    ax.text(pivots[-1]-700, 0.9, mutation)
    ax.set_yticklabels(['{:3.0f}%'.format(x*100) for x in [0, 0.2, 0.4, 0.6, 0.8, 1.0]]) 
    ax.tick_params(labelsize=fs)
    ax.xaxis.set_major_locator(months)
    ax.xaxis.set_major_formatter(monthsFmt)
fig.autofmt_xdate()
fax = fig.add_axes( [0., 0., 1, 1] )
fax.set_axis_off()
fax.set_xlim(0, 1)
fax.set_ylim(0, 1)
fax.text(0.02, 0.54, "Frequency", rotation='vertical', horizontalalignment='center', verticalalignment='center')
axs[mut_legend['panel']].legend(loc=mut_legend['loc'], ncol=1, bbox_to_anchor=(1.0, 0.2))
#plt.tight_layout(h_pad=0.01)
bottom_margin = 0.2 - 0.03*len(mutations)
plt.subplots_adjust(left=0.12, right=0.84, top=0.96, bottom=bottom_margin)
plt.savefig('figures/feb-2016/'+virus+'_mutations.png')
