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
flist = glob.glob('data/*'+res+'_cHI.txt')
figheight=4

plt.figure(figsize=(1.3*figheight, figheight))
ax=plt.subplot(111)
for fname in flist:
    flu = fname.split('/')[-1].split('_')[0]
    cHI_trunk = np.loadtxt(fname)
    R = linregress(cHI_trunk)
    plt.plot(cHI_trunk[:,0], cHI_trunk[:,1], '-o', 
             label= flu+r', $\mu='+str(np.round(R[0],2))+'\pm'+str(np.round(R[-1],2))+'$')

ax.tick_params(axis='both', labelsize=fs)
plt.ylabel(r'cumulative antigenic change $[\log_2]$', fontsize=fs)
plt.xlabel(r'year', fontsize=fs)
plt.legend(loc=2, fontsize = fs-2)
plt.tight_layout()
plt.savefig('cHI_trunk.png')


flist = glob.glob('data/*'+res+'_cHI_path.pkl')
flist.extend(glob.glob('data/*'+'7y'+'_cHI_path.pkl'))

plt.figure(figsize=(1.4*figheight, figheight))
ax=plt.subplot(111)
cols = {flu:col for flu, col in zip(['H3N2', 'H1N1pdm', 'Yam', 'Vic'], sns.color_palette(n_colors=4))}
all_path = defaultdict(list)
for fname in flist:
    flu = fname.split('/')[-1].split('_')[0]
    with open(fname) as ifile:
        cHI_path = pickle.load(ifile)
    for p in cHI_path:
        all_path[flu].extend(p)
        plt.plot(p[:,0], p[:,1], 'o', ls='none', c=cols[flu])

for flu in all_path:
    R = linregress(np.array(all_path[flu]))
    t = np.array([1990,2015])
    plt.plot(t,R[1]+R[0]*t, c = cols[flu], lw=2, ls='-',
             label= flu+r', $\mu='+str(np.round(R[0],2))+'$')
plt.xlim([1995,2017])
plt.ylim([-0.5, 19])
ax.tick_params(axis='both', labelsize=fs)
plt.ylabel(r'cumulative antigenic change $[\log_2]$', fontsize=fs)
plt.xlabel(r'year', fontsize=fs)
plt.legend(loc=2, fontsize = fs-2)
plt.tight_layout()

for fmt in fmts: plt.savefig('cHI_path'+fmt)
