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

freqs = json.load(open('../data/'+virus+'_2y_frequencies.json'))
counts = json.load(open('../data/'+virus+'_2y_meta.json'))['virus_stats']
region_names = json.load(open('../data/'+virus+'_2y_meta.json'))['regions']
region_codes = {'EU':['Europe'], 'AS':['China', 'JapanKorea', 'SouthAsia', 'SoutheastAsia', 'WestAsia'],
                'NA':["NorthAmerica"], 'OC':["Oceania"]}


if virus=='H3N2': ########## H3N2
    clades = ['3c2.a', '3c3.a', '3c3.b']
    mutations = ['HA1:114T', 'HA1:142K', 'HA1:168V', 'HA1:171K']
    #mutations = ['HA1:94H', 'HA1:114T', 'HA1:142K', 'HA1:171K']
    clade_legend = {'panel':0, 'loc':3}
    mut_legend = {'panel':0, 'loc':3}
elif virus=='H1N1pdm': ########## H1N1pdm
    clades = ['6b.1', '6b.2']
    mutations = ['HA1:84N','HA1:162N','HA1:152T'] #these don't add up to one in Asia, probably due to sketchy sampling.
    clade_legend = {'panel':0, 'loc':3}
    mut_legend = {'panel':0, 'loc':3}
elif virus=='Vic':
    clades = []
    mutations = ['HA1:129D', 'HA1:117V'] # HA1:56K would be good, but it currently isn't computed -> need to lower the threshold.
    clade_legend = {'panel':0, 'loc':3}
    mut_legend = {'panel':0, 'loc':3}
elif virus=='Yam':
    clades = ['2', '3']
    mutations = ['HA1:172Q', 'HA1:251V']
    clade_legend = {'panel':0, 'loc':3}
    mut_legend = {'panel':0, 'loc':3}


offset = datetime(2000,1,1).toordinal()
pivots = [offset+(x-2000)*365.25 for x in  freqs['clades']['global']['pivots']]
pivots.pop()
regions = ['global', 'NA', 'AS', 'EU', 'OC']
region_label = {'global': 'Global', 'NA': 'N America', 'AS': 'Asia', 'EU': 'Europe', 'OC': 'Oceania'}
cols = sns.color_palette(n_colors=len(regions))
fs=12

n=4
l = len(freqs['clades']['global']['pivots'])
count_array = np.array([c[1:] for c in counts])[-l:,:].T
smoothed_count_array = np.array([np.convolve(np.ones(n, dtype=float)/n, c, mode='same')
                        for c in count_array])
smoothed_count_by_region = {region: np.sum([smoothed_count_array[region_names.index(r)] for r in region_codes[region]], axis=0)
                            for region in regions if region!='global'}
smoothed_count_by_region['global'] = smoothed_count_array.sum(axis=0)




years    = YearLocator()
months = MonthLocator(range(1, 13), bymonthday=1, interval=2)
yearsFmt = DateFormatter('%Y')
monthsFmt = DateFormatter("%b")

if len(clades):
    fig, axs = plt.subplots(len(clades), 1, sharex=True, figsize=(8, len(clades)*2))
    for clade, ax in zip(clades, axs):
        for c,region in zip(cols, regions):
            try:
                tmp_freq = freqs['clades'][region][clade]
                if tmp_freq is not None:
                    tmp_freq = np.array(tmp_freq[:-1])
                    std_dev = np.sqrt(tmp_freq*(1-tmp_freq)/(smoothed_count_by_region[region][:-1]+2))
                    ax.plot_date(pivots, tmp_freq,'-o', label = region_label[region], c=c, lw=3 if region=='global' else 1, clip_on=False)
                    ax.plot_date(pivots, tmp_freq+std_dev,':', c=c, lw=1)
                    ax.plot_date(pivots, tmp_freq-std_dev,':', c=c, lw=1)
            except:
                print "skipping", clade, region
        ax.set_xlim([pivots[0], pivots[-1]])
        ax.set_ylim(0,1)
        ax.text(pivots[0]+5, 0.88, clade)
        ax.set_yticklabels(['{:3.0f}%'.format(x*100) for x in [0, 0.2, 0.4, 0.6, 0.8, 1.0]])
        ax.tick_params(which = 'major', labelsize=fs, pad=20)
        ax.tick_params(which = 'minor', pad=7)
        ax.xaxis.set_major_locator(years)
        ax.xaxis.set_major_formatter(yearsFmt)
        ax.xaxis.set_minor_locator(months)
        ax.xaxis.set_minor_formatter(monthsFmt)
    fig.autofmt_xdate(bottom=0.25, rotation=0, ha='center')
    fax = fig.add_axes( [0., 0., 1, 1] )
    fax.set_axis_off()
    fax.set_xlim(0, 1)
    fax.set_ylim(0, 1)
    fax.text(0.02, 0.54, "Frequency", rotation='vertical', horizontalalignment='center', verticalalignment='center')
    axs[clade_legend['panel']].legend(loc=clade_legend['loc'], ncol=1, bbox_to_anchor=(1.02, 0.2))
    bottom_margin = 0.18 - 0.03*len(clades)
    plt.subplots_adjust(left=0.12, right=0.82, top=0.97, bottom=bottom_margin)
    plt.savefig('figures/feb-2016/'+virus+'_clades.png')


fig, axs = plt.subplots(len(mutations), 1, sharex=True, figsize=(8, len(mutations)*2))
for mutation, ax in zip(mutations, axs):
    for c,region in zip(cols, regions):
        try:
            tmp_freq = freqs['mutations'][region][mutation]
            if tmp_freq is not None:
                tmp_freq = np.array(tmp_freq[:-1])
                std_dev = np.sqrt(tmp_freq*(1-tmp_freq)/(smoothed_count_by_region[region][:-1]+2))
                ax.plot_date(pivots, tmp_freq, '-o', label = region_label[region], c=c, lw=3 if region=='global' else 1, clip_on=False)
                ax.plot_date(pivots, tmp_freq+std_dev,':', c=c, lw=1)
                ax.plot_date(pivots, tmp_freq-std_dev,':', c=c, lw=1)
        except:
            print "skipping", mutation, region
    ax.set_xlim([pivots[0], pivots[-1]])
    ax.set_ylim(0,1)
    ax.text(pivots[0]+5, 0.88, mutation)
    ax.set_yticklabels(['{:3.0f}%'.format(x*100) for x in [0, 0.2, 0.4, 0.6, 0.8, 1.0]])
    ax.tick_params(labelsize=fs, pad=18)
    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_major_formatter(yearsFmt)
    ax.xaxis.set_minor_locator(months)
    ax.xaxis.set_minor_formatter(monthsFmt)
fig.autofmt_xdate(bottom=0.25, rotation=0, ha='center')
fax = fig.add_axes( [0., 0., 1, 1] )
fax.set_axis_off()
fax.set_xlim(0, 1)
fax.set_ylim(0, 1)
fax.text(0.02, 0.54, "Frequency", rotation='vertical', horizontalalignment='center', verticalalignment='center')
axs[mut_legend['panel']].legend(loc=mut_legend['loc'], ncol=1, bbox_to_anchor=(1.02, 0.2))
bottom_margin = 0.18 - 0.03*len(mutations)
plt.subplots_adjust(left=0.12, right=0.82, top=0.97, bottom=bottom_margin)
plt.savefig('figures/feb-2016/'+virus+'_mutations.png')
