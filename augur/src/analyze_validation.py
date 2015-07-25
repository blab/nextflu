import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import cPickle
plt.ion()

val_data = []
res='12y'
grid = [0.1, 0.3, 1.0,3.0, 10.0]
for flu in ['H3N2', 'H1N1pdm', 'Vic', 'Yam']:
    for minaa in [0]: #,1,'epi']:
        for hi, lam_HI in enumerate(grid):
            for training in ['measurements', 'virus']:
                fname = 'validation/'+'_'.join([flu, training, 'minaa', str(minaa),'lHI',str(lam_HI)])+'_tree.pkl'
                try:
                    with open(fname) as infile:
                        params,tmp_grid, acc = cPickle.load(infile)
                except:
                    print fname, 'not found'
                else:
                    for pi, lam_pot in enumerate(grid):
                        for ai, lam_avi in enumerate(grid):
                            val_data.append([flu, str(minaa), training, lam_HI, lam_pot, lam_avi] +list(acc[pi,ai]))

val_data = pd.DataFrame(val_data, columns=['lineage', 'minaa', 'training', 'HI_reg', 'pot_reg', 
                                            'avi_reg','rms', 'abs','slope','intercept', 'sym'])



cols = sns.color_palette(n_colors=len(grid))
ls = ['-', '--', '-.']
#training = 'virus'
training = 'measurements'
for flu in ['H3N2', 'H1N1pdm', 'Vic', 'Yam']:
    fig, axs = plt.subplots(2,len(grid), sharey='row')
    plt.suptitle(str([flu, training, minaa]))
    for training in ['virus', 'measurements']:
        for mi, minaa in enumerate([0]): #,1,'epi']):
            for hi, lam_HI in enumerate(grid):
                for pi, lam_pot in enumerate(grid):
                    try:
                        ind = (val_data['lineage']==flu)&(val_data['training']==training)&(val_data['minaa']==str(minaa))&\
                                (val_data['HI_reg']==lam_HI)&(val_data['pot_reg']==lam_pot)
                        print flu, training, minaa, lam_HI, lam_pot, val_data.loc[ind,:]
                        ax =axs[0][hi]
                        ax.plot(grid, val_data.loc[ind,:]['abs'], c=cols[pi], ls = ls[mi], 
                                label='aa='+str(minaa)+', pot='+str(lam_pot),)
                        if hi==len(grid)-1 and training=='virus':
                            ax.legend(ncol=len(ls), loc=4)
                        ax.set_ylim([0,1.5])

                        ax =axs[1][hi]
                        ax.plot(grid, val_data.loc[ind,:]['slope'], c=cols[pi], ls = ls[mi])
                        ax.set_ylim([0,1])
                    except:
                        print 'cant plot'