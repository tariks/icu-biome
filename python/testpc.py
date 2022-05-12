# %%
%matplotlib inline
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from skbio.stats import distance
from skbio.stats.composition import multiplicative_replacement,clr
from scipy.spatial.distance import pdist
from scipy import stats
import cmasher as cmr
plt.style.use("../lance.txt")
# %%
meta=pd.read_csv('../meta/meta52_current.csv',index_col=0)
table=pd.read_csv('../feature_tables/52N_genus_t80.csv',index_col=0)
# %%
def bal(df=table):
    x=df.copy()
    delta=1
    #x=x.loc[:,(x>0).sum()/52 > .25] 
    y = x.divide(x.sum(axis=1),axis=0)
    #delta=y[y>0].min().min()
    #x=y
    print(delta)
    x[x==0]+=np.random.random(x[x==0].shape)*.9*delta +.1*delta
    x[:] = multiplicative_replacement(x,)
    #x[:] = clr(x)
    x[:] = np.log10(x)
    N = x['Enterobacteriaceae'] + x['Anaerococcus']
    D = x['Parasutterella'] + x['Campylobacter']
    b = (N-D)/2
    return b
# %%
meta['mbal2'] = bal()
meta.to_csv('../meta/meta52_current.csv')
sns.histplot(data=meta,x='mbal2',hue='month',binwidth=.2)
meta['mbal2'].describe()
# %%
meta['bacto']=0
meta.loc[(.7*meta['PC2']+.6*meta['PC1'] < -.03),'bacto']=1
meta.groupby('bacto').describe()['mbal'] #,'Death','Age','APACHE','Zosyn','ARDS']
norm=plt.cm.colors.Normalize(vmin=.7,vmax=6,clip=True)
# %%
v='M'
v='D'
v='bacto'
v2='PC1'
v2='mbal'
v2='Death'
v2='logDeath'
meta['D']=0
meta.loc[meta['Death']<24,'D']=1
meta['logDeath'] = np.log2(meta['Death'])
a,b = meta['PC1'], meta['mbal']
a,b = meta.loc[meta[v]==0,v2],meta.loc[meta[v]==1,v2]
#sns.histplot(data=meta,x='mbal',hue='bacto')
print(stats.ranksums(a,b))
print(stats.ttest_ind(a,b))
#stats.pearsonr(a,b)
v2dict={'Death': 'Death-free days', 'mbal': 'Biomarker', 'PC1': 'PC1'}
# %%
h = [plt.scatter([0,1],[0,1],c='#121212',marker='o'),
     plt.scatter([0,1],[0,1],c='#121212',marker='x'),
]
l=['0','1']
def plotstuff(ax,v1,v2,hue,style=None,thres=None):
    order = [0,1] if style else None
    if v2=='mbal':
        cmap=cmr.dusk
    else:
        cmap=cmr.dusk_r
    if v1=='M' or style=='M':
        meta['M']=0
        meta.loc[meta['mbal']>thres,'M']=1
    vdict = {
    'M': 'Biomarker < {}'.format(thres),
    'D': 'Death < 24 days',
    'bacto': 'PCA clusters'
    }
    sns.scatterplot(data=meta,
            x='PC1',y='PC2',ax=ax,
            hue=hue,alpha=.9,
            style=style,style_order=order,palette=cmap,s=30
    #        markers=markers
    )
    if v1=='bacto':
        ax.plot([-.4,.3],[.3,-.3],lw=.8,ls=(5,(5,10)),color='#121212')
    m1 =meta[hue].min()
    m2 = meta[hue].max()
    rang = (m2-m1)*.1
    norm=plt.cm.colors.Normalize(vmin=m1+rang,vmax=m2-rang,clip=True)
    cbar=plt.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap), ax=ax,
        fraction=.03,location='right',aspect=25,shrink=.6,pad=.01
    )
    cbar.outline.set(lw=0)
    ax.get_legend().remove()
    if style:
        ax.legend(h,l,loc='upper right',
        bbox_to_anchor=[.99,.98],
        title=vdict.get(style),title_fontsize=6,
        borderpad=0,borderaxespad=0,markerscale=.75,
        scatteryoffsets=[.5])

    a,b = meta.loc[meta[v1]==0],meta.loc[meta[v1]==1]
    p = stats.ranksums(a[v2],b[v2])
    s = ['Rank sums test:', 
         'Variable: {}'.format(v2dict.get(v2)),
         'Groups: {}'.format(vdict.get(v1)),
         'n1={}, n2={}'.format(a.shape[0],b.shape[0]),
         'Medians: {:.1f}, {:.1f}'.format(a[v2].median(),b[v2].median()),
         'Statistic: {:.2f}'.format(p[0]),
         'P value: {:.3f}'.format(p[1])
    ]
    ax.text(0.01,0.015,'\n'.join(s),fontsize=5,ha='left',va='bottom',
            transform=ax.transAxes)
    return ax
#plt.savefig('../plots/pca_enterotype.png',bbox_inches='tight')
# %%
fig,axs=plt.subplots(5,2)
fig.set(figheight=12,figwidth=6)
plotstuff(axs[0,0],'bacto','mbal','mbal','D')
plotstuff(axs[0,1],'bacto','Death','logDeath','M',thres=1.9)
plotstuff(axs[1,0],'bacto','Death','logDeath','M',thres=2)
plotstuff(axs[1,1],'M','Death','logDeath','M',thres=1.9)
plotstuff(axs[2,0],'M','Death','logDeath','M',thres=2)
plotstuff(axs[2,1],'M','PC1','PC1','M',thres=1.9)
plotstuff(axs[3,0],'M','PC1','PC1','M',thres=2)
plotstuff(axs[3,1],'D','mbal','mbal','D')
plotstuff(axs[4,0],'D','PC1','PC1','D')
for i in axs[:,1]:
    i.set_ylabel('')
axs[4,1].set_axis_off()
plt.savefig('../plots/pca_many.pdf',dpi=300,transparent=True,bbox_inches='tight')
# %%
#distance.permanova(d,meta,column='bacto')

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
meta['M'] = 0
meta.loc[meta['mbal']>2.5,'M'] = 1
print(meta['M'])
distance.permanova(d,meta,column='M')
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
