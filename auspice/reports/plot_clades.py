import argparse
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from datetime import datetime, timedelta
import json
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter

#################################
# Settings to configure
parser = argparse.ArgumentParser()
parser.add_argument('--lineage', default='h3n2')
args = parser.parse_args()

virus = args.lineage

report = 'sep-2017'
drop = 1
resolution = '2y'

if virus=='h3n2':
    clades = ['3c2.a', '3c3.a', '171K']
    # mutations = ['HA1:121K', 'HA1:92R', 'HA1:131K', 'HA1:31S','HA1:198P', 'HA1:193S']
    # Large subclades
    # 1. 31S, 53N, 144R, 171K, 192T, 197H, characteristic: 197H
    # 2. 121K, 144K, characteristic: 144K
    # 3. 131K, 142K, 261Q, characteristic: 131K
    # 4. 142R, HA2:150E, characteristic: HA2:150E ---> additional 135K
    # 5. 92R, 311Q, characteristic: 92R ---> additional 135K
    mutations = ['HA1:197H', 'HA1:144K', 'HA1:131K', 'HA2:150E', 'HA1:92R']
    clade_legend = {'panel':0, 'loc':3}
    mut_legend = {'panel':0, 'loc':3}
    ymax = 1000
elif virus=='h1n1pdm':
    clades = []
    mutations = ['HA1:74R', 'HA1:295V', 'HA1:164T']
    #mutations = ['HA1:163Q', 'HA1:162N', 'HA1:74R']
    clade_legend = {'panel':0, 'loc':3}
    mut_legend = {'panel':0, 'loc':3}
    ymax = 200
elif virus=='vic':
    clades = []
    mutations = ['HA1:162-', 'HA1:175V', 'HA1:209N']
    clade_legend = {'panel':0, 'loc':3}
    mut_legend = {'panel':0, 'loc':3}
    ymax = 250
elif virus=='yam':
    clades = []
    mutations = ['HA1:251V', 'HA1:211R', 'HA1:76I']
    clade_legend = {'panel':0, 'loc':3}
    mut_legend = {'panel':0, 'loc':3}
    ymax = 250

file_addendum = '_cell_hi'
#file_addendum = ''
file_prefix = '../data/flu_'+virus+'_ha_'+resolution+file_addendum

#################################
# Plotting

print "Plotting dataset " + virus + "/" + resolution + " for " + report + " report"

freqs = json.load(open(file_prefix+'_frequencies.json'))
counts = freqs['counts']
offset = datetime(2000,1,1).toordinal()
pivots = [offset+(x-2000)*365.25 for x in freqs['pivots']]
meta = json.load(open(file_prefix+'_meta.json'))
region_names = meta['geo']['region'].keys()
region_codes = {'EU':['europe'], 'AS':['china', 'south_asia', 'japan_korea','southeast_asia'],
                'NA':["north_america"], 'OC':["oceania"]}

offset = datetime(2000,1,1).toordinal()
regions = ['NA', 'AS', 'EU', 'OC', 'global']
region_label = {'NA': 'N America', 'AS': 'Asia', 'EU': 'Europe', 'OC': 'Oceania', 'global': 'Global'}
cols = sns.color_palette(n_colors=len(regions))
cols.pop(0)
cols.append(sns.color_palette(["#3E547F"], n_colors=1)[0])
fs=12

years = YearLocator()
months = MonthLocator(range(1, 13), bymonthday=1, interval=2)
yearsFmt = DateFormatter('%Y')
monthsFmt = DateFormatter("%b")

sns.set_style('ticks')
plt.ion()
show_errorbars = True

n=2
n_std_dev=2
l = len(pivots)
bins = np.array([c[0] for c in counts])[-l:]
date_bins = []
date_bins=pivots
count_array = np.array([counts[c] for c in counts])[-l:,:].T
count_by_region = {region: sum(x) for region,x in counts.iteritems() if region!='global'}
smoothed_count_array = np.array([np.convolve(np.ones(n, dtype=float)/n, c, mode='same')
                        for c in count_array])
smoothed_count_by_region = {region: np.convolve(np.ones(n, dtype=float)/n, c, mode='same')
                            for region, c in counts.iteritems()}

smoothed_count_by_region.update({r1:smoothed_count_by_region[r2] for r1,r2 in
                            [['north_america', 'NA'], ['europe', 'EU'], ['asia', 'EAS'], ['oceania', 'OC']]})

print "Plotting sample counts"
fig, ax = plt.subplots(figsize=(8, 3))
tmpcounts = np.zeros(len(date_bins[drop:]))
plt.bar(date_bins[drop:], counts['global'][drop:], width=20, linewidth=0, label="Other", color="#bbbbbb", clip_on=False)
for c,region in zip(cols, regions):
    if region!='global':
        plt.bar(date_bins[drop:], counts[region][drop:], bottom=tmpcounts, width=20, linewidth=0,
                label=region_label[region], color=c, clip_on=False)
        tmpcounts += np.array(counts[region][drop:])
ax.set_xlim([pivots[drop-1], pivots[-1]])
ax.set_ylim(0,ymax)
ax.tick_params(axis='x', which='major', labelsize=fs, pad=20)
ax.tick_params(axis='x', which='minor', pad=7)
ax.xaxis.set_major_locator(years)
ax.xaxis.set_major_formatter(yearsFmt)
ax.xaxis.set_minor_locator(months)
ax.xaxis.set_minor_formatter(monthsFmt)
ax.set_ylabel('Sample count', fontsize=fs*1.1)
ax.legend(loc=3, ncol=1, bbox_to_anchor=(1.02, 0.53))
plt.subplots_adjust(left=0.1, right=0.82, top=0.94, bottom=0.22)
sns.despine()
plt.savefig('figures/' + report + '/'+virus+'_counts.png')

if len(clades):
    print "Plotting clade frequencies"
    fig, axs = plt.subplots(len(clades), 1, sharex=True, figsize=(8, len(clades)*2))
    for clade, ax in zip(clades, axs):
        for c,(region, r1) in zip(cols, [('north_america', 'NA'), ('china', 'AS'), ('europe','EU'), ('oceania','OC'), ('global', 'global')]):
            if True: #try:
                tmp_freq = np.array(freqs['%s_%s'%(region, clade)])
                if tmp_freq is not None:
                    std_dev = np.sqrt(tmp_freq*(1-tmp_freq)/(smoothed_count_by_region[r1]+1))
                    ax.plot(pivots, tmp_freq,'-o', label = region_label[r1], c=c, lw=3 if region=='global' else 1)
                    if show_errorbars:
                        ax.fill_between(pivots, tmp_freq-n_std_dev*std_dev, tmp_freq+n_std_dev*std_dev, facecolor=c, linewidth=0, alpha=0.1)
            # except:
            #     print "skipping", clade, region
        ax.set_xlim([pivots[0], pivots[-1]])
        ax.set_ylim(0,1)
        ax.text(pivots[0]+5, 0.88, clade)
        ax.set_yticklabels(['{:3.0f}%'.format(x*100) for x in [0, 0.2, 0.4, 0.6, 0.8, 1.0]])
        ax.tick_params(axis='x', which='major', labelsize=fs, pad=20)
        ax.tick_params(axis='x', which='minor', pad=7)
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
    bottom_margin = 0.22 - 0.03*len(clades)
    plt.subplots_adjust(left=0.12, right=0.82, top=0.97, bottom=bottom_margin)
    sns.despine()
    plt.savefig('figures/' + report + '/'+virus+'_clades.png')

print "Plotting mutation frequencies"
fig, axs = plt.subplots(len(mutations), 1, sharex=True, figsize=(8, len(mutations)*2))
for mutation, ax in zip(mutations, axs):
    for c,region in zip(cols, regions):
        if '%s_%s'%(region, mutation) in freqs:
            tmp_freq = np.array(freqs['%s_%s'%(region, mutation)])
            std_dev = np.sqrt(tmp_freq*(1-tmp_freq)/(smoothed_count_by_region[region]+1))
            ax.plot(pivots[drop:], tmp_freq[drop:], '-o', label = region_label[region], c=c, lw=3 if region=='global' else 1)
            if show_errorbars:
                ax.fill_between(pivots[drop:], (tmp_freq-n_std_dev*std_dev)[drop:], (tmp_freq+n_std_dev*std_dev)[drop:], facecolor=c, linewidth=0, alpha=0.1)
        else:
            print "no data for", mutation, region, "padding with zeros"
            tmp_freq = np.zeros(len(pivots))
            std_dev = np.sqrt(tmp_freq*(1-tmp_freq)/(smoothed_count_by_region[region]+1))
            ax.plot(pivots[drop:], tmp_freq[drop:], '-o', label = region_label[region], c=c, lw=3 if region=='global' else 1)
            if show_errorbars:
                ax.fill_between(pivots[drop:], (tmp_freq-n_std_dev*std_dev)[drop:], (tmp_freq+n_std_dev*std_dev)[drop:], facecolor=c, linewidth=0, alpha=0.1)
    ax.set_xlim([pivots[drop-1], pivots[-1]])
    ax.set_ylim(0,1)
    ax.text(pivots[drop-1]+10, 0.88, mutation, fontsize=fs)
    ax.set_yticklabels(['{:3.0f}%'.format(x*100) for x in [0, 0.2, 0.4, 0.6, 0.8, 1.0]])
    ax.tick_params(axis='x', which='major', labelsize=fs, pad=20)
    ax.tick_params(axis='x', which='minor', pad=7)
    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_major_formatter(yearsFmt)
    ax.xaxis.set_minor_locator(months)
    ax.xaxis.set_minor_formatter(monthsFmt)
fig.autofmt_xdate(bottom=0.25, rotation=0, ha='center')
fax = fig.add_axes( [0., 0., 1, 1] )
fax.set_axis_off()
fax.set_xlim(0, 1)
fax.set_ylim(0, 1)
fax.text(0.02, 0.54, "Frequency", rotation='vertical', horizontalalignment='center', verticalalignment='center', fontsize=fs*1.1)
axs[mut_legend['panel']].legend(loc=mut_legend['loc'], ncol=1, bbox_to_anchor=(1.02, 0.2))
bottom_margin = 0.22 - 0.03*len(mutations)
plt.subplots_adjust(left=0.12, right=0.82, top=0.97, bottom=bottom_margin)
sns.despine()
plt.savefig('figures/'+report+'/'+virus+'_mutations.png')
