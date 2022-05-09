# %%
%matplotlib inline
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

plt.style.use("../lance.txt")
# %%
meta = pd.read_csv('../meta/meta52pc.csv',index_col=0)

# %%
div_metrics= [
    "ace",
    "chao1",
    "dominance",
    "fisher_alpha",
    "pielou_e",
    "shannon",
    "simpson",
]

# %%
meta['M']=0
meta.loc[meta['mbal']>1.9,'M']=1
#meta.loc[meta['mbal']>4,'M']=
# %%
param = dict(boxprops=dict(
                    fill=False,
                    snap=True,
                    capstyle='round',
            )
)
# %%
fig,axs=plt.subplots(2,4)
fig.subplots_adjust(wspace=.55,hspace=.2)
for i,v in enumerate(div_metrics):
    ax=axs.flat[i]
    ax.margins(y=.2)
    sns.boxplot(data=meta,y=v,x='M',hue='M',ax=ax,
                linewidth=.5,
                fliersize=0,
                dodge=False,
                **param
                )

    sns.swarmplot(data=meta,y=v,x='M',hue='M',ax=ax,
                  size=2,alpha=.8,marker='d'
    )    
    ax.set(ylabel='',xlabel='',xticks=[])
    ax.set_title(v,pad=1)
    ax.get_legend().remove()
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    xlim=list(ax.get_xlim())
    xlim[0]-=.1
    xlim[1]+=.1
    ax.set_xlim(xlim)
    ax.spines
empty=axs[1,3]
empty.set_axis_off()
empty.set(xticks=[],yticks=[])
h,l = axs[0,0].get_legend_handles_labels()
l = ['< 2', '< 4', '> 4']
l = ['< 1.9','> 1.9']
empty.legend(h,l,loc='upper left',title='Biomarker',
            borderaxespad=0,borderpad=0
)
plt.savefig('../plots/alphabox2.pdf',dpi=300,transparent=True,bbox_inches='tight')
# %%
meta['M']
# %%
