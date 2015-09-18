import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from datetime import datetime
import json

sns.set_style('darkgrid')
plt.ion()

virus = 'B/Yam'


if virus=='H3N2': ########## H3N2
    freqs = json.load(open('../data/H3N2_1y_frequencies.json'))
    clades = ['3c2.a', '3c3.a', '3c3.b']
    mutations = ['HA1:159Y', 'HA1:159S', 'HA1:159F']
elif virus=='H1N1pdm': ########## H1N1pdm
    freqs = json.load(open('../data/H1N1pdm_1y_frequencies.json'))
    clades = ['6c', '6b']
    mutations = ['HA1:84N', 'HA2:164G'] #these don't add up to one in Asia, probably due to sketchy sampling.
elif virus=='B/Vic':
    freqs = json.load(open('../data/Vic_1y_frequencies.json'))
    clades = []
    mutations = ['HA1:129D', 'HA1:117V'] # HA1:56K would be good, but it currently isn't computed -> need to lower the threshold.
elif virus=='B/Yam':
    freqs = json.load(open('../data/Yam_1y_frequencies.json'))
    clades = ['3', '2']
    mutations = ['HA1:251V', 'HA1:172Q'] # HA1:56K would be good, but it currently isn't computed -> need to lower the threshold.


offset = datetime(2000,1,1).toordinal()
pivots = [offset+(x-2000)*365.25 for x in  freqs['clades']['global']['pivots']]
regions = ['global', 'NA', 'AS', 'EU', 'OC']
cols = sns.color_palette(n_colors=len(regions))

if len(clades):
    fig, axs = plt.subplots(len(clades), 1, sharex=True)
    for clade, ax in zip(clades, axs):
        for c,region in zip(cols, regions):
            try:
                tmp_freq = freqs['clades'][region][clade]
                if tmp_freq is not None:
                    ax.plot_date(pivots, tmp_freq,'-o', label = region, c=c)
            except:
                print "skipping", clade, region
        ax.legend(loc=2, title=clade)
        ax.set_xlim([pivots[0]-100,pivots[-1]+30])
        ax.set_ylim(0,1)


fig, axs = plt.subplots(len(mutations), 1, sharex=True)
for mutation, ax in zip(mutations, axs):
    for c,region in zip(cols, regions):
        try:
            tmp_freq = freqs['mutations'][region][mutation]
            if tmp_freq is not None:
                ax.plot_date(pivots, tmp_freq, '-o',label = region, c=c)
        except:
            print "skipping", mutation, region
    ax.legend(loc=2, title=mutation)
    ax.set_xlim([pivots[0]-100,pivots[-1]+30])
    ax.set_ylim(0,1)
ax.set_xlabel('time')
