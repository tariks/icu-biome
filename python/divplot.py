# %%
%matplotlib inline
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

plt.style.use("../lance.txt")
# %%
meta = pd.read_csv('../meta/meta52_current.csv',index_col=0)

# %%
div_metrics = [
# 'ACE (R)',
 'Chao1 (R)',
 'Fisher (D)',
 'Shannon (D)',
 'Pielou (E)',
 'Dominance',
 'Simpson (D)',
 'Simpson (E)',
]



# %%
v='M'
a = meta.loc[meta[v]==0]
b = meta.loc[meta[v]==1]
pvals = {i: stats.ranksums(a[i],b[i],)[1] for i in div_metrics}
pvals = {i: stats.ttest_ind(a[i],b[i],alternative='two-sided')[1] for i in div_metrics}
for k,v in pvals.items():
    print('{}\t{:.3f}'.format(k,v))

# %%
param = dict(boxprops=dict(
                    fill=False,
                    snap=True,
                    capstyle='round',
                    alpha=.8,
                    aa=True,
            )
)
# %%
fig,axs=plt.subplots(2,4)
fig.subplots_adjust(wspace=.3,hspace=.3)
fig.set(figheight=2.75,figwidth=4.1)
for i,v in enumerate(div_metrics):
    ax=axs.flat[i]
    ax.margins(y=.075,x=.2)
    sns.boxplot(data=meta,y=v,x='M',hue='M',ax=ax,
                linewidth=.5,
                fliersize=0,
                dodge=False,
                **param
                )

    sns.swarmplot(data=meta,y=v,x='M',hue='M',ax=ax,
                  size=2.55,alpha=.85,marker='^',
    )    
    ax.set(ylabel='',xlabel='',xticks=[])
    ax.set_title(v,pad=1)
    ax.get_legend().remove()
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    xlim=list(ax.get_xlim())
    xlim[0]-=.1
    xlim[1]+=.1
    ax.set_xlim(xlim)
    p = pvals.get(v)
    print(p)
    if p < .01:
        prefix = '** '
    elif p < .05:
        prefix = '* '
    else:
        prefix = ''
    ax.set_xlabel(prefix+'P= {:.3f}'.format(p),fontsize=6)
    ax.tick_params(labelsize=6)
empty=axs.flat[-1]
empty.set_axis_off()
empty.set(xticks=[],yticks=[])
h,l = axs[0,0].get_legend_handles_labels()
l = ['< 0','> 0']
empty.legend(h,l,loc='upper left',title='MMI',
            borderaxespad=0,borderpad=0
)
plt.savefig('../plots/alphabox_ttest.pdf',dpi=300,transparent=True,bbox_inches='tight')
# %%
