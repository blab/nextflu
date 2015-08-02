import matplotlib.pyplot as plt
import numpy as np
import glob
import seaborn as sns
import cPickle as pickle
from scipy.stats import linregress
from collections import defaultdict
sns.set_style('darkgrid')

fs=14
fmts = ['.pdf', '.svg', '.png']
res = '20y'
figheight=4

#############################################
### plot the cumulative antigenic change
#############################################
plt.figure(figsize=(1.3*figheight, figheight))
ax=plt.subplot(111)
flist = glob.glob('data/*'+res+'_cHI.txt')
for fname in flist:
    flu = fname.split('/')[-1].split('_')[0]
    cHI_trunk = np.loadtxt(fname)
    # regress against time
    R = linregress(cHI_trunk)
    plt.plot(cHI_trunk[:,0], cHI_trunk[:,1], '-o', 
             label= flu+r', $\mu='+str(np.round(R[0],2))+'\pm'+str(np.round(R[-1],2))+'$')

ax.tick_params(axis='both', labelsize=fs)
plt.ylabel(r'cumulative antigenic change $[\log_2]$', fontsize=fs)
plt.xlabel(r'year', fontsize=fs)
plt.legend(loc=2, fontsize = fs-2)
plt.tight_layout()
plt.savefig('cHI_trunk.png')


############################################
### alternative figure tracing cHI along lineages
############################################
flist = glob.glob('data/*'+res+'_cHI_path.pkl')
flist.extend(glob.glob('data/*'+'7y'+'_cHI_path.pkl'))

plt.figure(figsize=(2*figheight, figheight))
ax=plt.subplot(121)
cols = {flu:col for flu, col in zip(['H3N2', 'H1N1pdm', 'Yam', 'Vic'], sns.color_palette(n_colors=4))}
all_path = defaultdict(list)
for fname in flist:
    flu = fname.split('/')[-1].split('_')[0]
    with open(fname) as ifile:
        cHI_path = pickle.load(ifile)
    # add multiple path from leaves to the root
    if flu=='H3N2':
        ax = plt.subplot(121)
    else:
        ax = plt.subplot(122)        
    for p in cHI_path:
        all_path[flu].extend([p[-1,:]])
        ax.plot(p[-1,0], p[-1,1], 'o', ls='none', c=cols[flu])
        ax.plot(p[:,0], p[:,1], ls='-', c=cols[flu], alpha=0.2)

# add linear regression for each flu type and add labels
for flu in all_path:
    if flu=='H3N2':
        ax = plt.subplot(121)
    else:
        ax = plt.subplot(122)        
    R = linregress(np.array(all_path[flu]))
    t = np.array([2007 if flu=='H1N1pdm' else 1993,2015])
    ax.plot(t,R[1]+R[0]*t, c = cols[flu], lw=2, ls='--',
             label= flu+r', $\mu='+str(np.round(R[0],2))+'$')
#plt.ylim([-0.5, 19])
for ai,ax in enumerate([plt.subplot(121), plt.subplot(122)]):
    ax.tick_params(axis='both', labelsize=fs)
    ax.set_xlim([1995,2017])
    ax.set_ylim([-0.5, 4 if ai else 19])
    ax.set_yticks(range(4) if ai else [0,5,10,15])
    if ai==0: ax.set_ylabel(r'cumulative antigenic change $[\log_2]$', fontsize=fs)
    ax.set_xlabel(r'year', fontsize=fs)
    ax.legend(loc=2, fontsize = fs-2)
plt.tight_layout()

for fmt in fmts: plt.savefig('cHI_path'+fmt)

