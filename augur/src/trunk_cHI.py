import matplotlib.pyplot as plt
import numpy as np
import glob
import seaborn as sns
from scipy.stats import linregress
sns.set_style('darkgrid')

fs=18
res = '20y'
flist = glob.glob('data/*'+res+'_cHI.txt')

plt.figure()
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
plt.legend(loc=2, fontsize = fs)
plt.tight_layout()
plt.savefig('cHI_trunk.png')
