# %%
%matplotlib inline
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from skbio.stats import distance
from scipy.spatial.distance import pdist
plt.style.use("../lance.txt")
# %%
meta=pd.read_csv('../meta/meta52pc.csv',index_col=0)
table=pd.read_csv('../feature_tables/52N_genus_t70.csv',index_col=0)
# %%
meta['bacto']=0
meta.loc[(.6*meta['PC1']+.7*meta['PC2'] < -.03),'bacto']=1
meta.groupby('bacto').describe()['mbal'] #,'Death','Age','APACHE','Zosyn','ARDS']
# %%
fig,ax=plt.subplots()
markers={'Minimal dysbiosis': 'o', 'Metabolic dysbiosis': 's', 'Pathogen dysbiosis': 'v'}
sns.scatterplot(x=meta['PC1'],y=meta['PC2'],ax=ax,hue=meta['bacto'],style=meta['Enterotype'],markers=markers)
ax.plot([-.4,.3],[.3,-.3],lw=.5)
plt.savefig('../plots/pca_enterotype.png',bbox_inches='tight')
# %%
d = distance.DistanceMatrix( pdist(meta[['PC1','PC2']]),meta.index)

# %%
distance.permanova(d,meta,column='bacto')
'''
method name               PERMANOVA
test statistic name        pseudo-F
sample size                      52
number of groups                  2
test statistic            34.686216
p-value                       0.001
number of permutations          999
'''
# %%
distance.permanova(d,meta,column='month')
'''
method name               PERMANOVA
test statistic name        pseudo-F
sample size                      52
number of groups                  2
test statistic             2.775417
p-value                       0.063
number of permutations          999
'''
# %%
distance.permanova(d,meta,column='Zosyn')
'''
method name               PERMANOVA
test statistic name        pseudo-F
sample size                      52
number of groups                  2
test statistic             3.820463
p-value                       0.028
number of permutations          999
'''
# %%
with open('../negative_permanova.txt','w') as o:
#with open('../positive_permanova.txt','w') as o:
    for i in ['Batch','Gender','ARDS','Shock','diabetes','malignancy','immunosuppressed','Vancomycin','Cipro','year','inhosp']: #,'Enterotype','min','patho','metab',]:
    #for i in ['bacto','month','Zosyn','Enterotype','min','patho','metab',]:
        x=distance.permanova(d,meta,column=i)[3:6]
        print('{}, {} groups'.format(i,x[0]),file=o)
        print('pseudo-F statistic: {:.2f}'.format(x[1]),file=o)
        print('P-value: {:.3f}'.format(x[2]),file=o)
        print('',file=o)


# %%
x = (meta['bacto']==1) & (meta['patho']==1)
x.sum()
# %%
meta['patho'].sum()
# %%
