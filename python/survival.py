# %%
%matplotlib inline
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import cmasher as cmr
import lifelines
plt.style.use("../lance.txt")
meta=pd.read_csv('../meta/meta52_current.csv',index_col=0)
meta['Sepsis3'] = meta['Sepsis-3'].copy()
del meta['Sepsis-3']
# %%
import warnings
warnings.simplefilter("ignore")

# %%
def cox(v,df=meta,pen=.001,l1=.5,event='month',robust=False):
    f = '+'.join(v) if type(v)==list else v
    modparams=dict(penalizer=pen,l1_ratio=l1)
    mod = lifelines.CoxPHFitter(**modparams)
    params=dict(
                robust=robust,
                formula=f,
    )
    mod.fit(df,duration_col='Death',event_col=event,**params)
    summary = mod.summary
    cols=    ['exp(coef)','exp(coef) lower 95%','exp(coef) upper 95%','p']

    print('\t'.join(['var','exp','low','hi','p']))
    for i in summary.index:
        vals=summary.loc[i,cols]
        print(i+'\t'+'\t'.join(['{:.3f}'.format(vals[j]) for j in cols]))

    aic=mod.AIC_partial_
    log=mod.log_likelihood_
    con=mod.concordance_index_
    print('\nK-fold:(LL, CI)')
    for i in ['log_likelihood','concordance_index']:
        val = lifelines.utils.k_fold_cross_validation(
                    lifelines.CoxPHFitter(**modparams),
                    df,
                    duration_col='Death',
                    event_col=event,
                    scoring_method=i,
                    seed=0,
                    fitter_kwargs=params
            )
        val = ['{: .2f}'.format(j) for j in val]
        print('\t'.join(val))
    print('')
    print('LL: {:.2f}\t CI: {:.2f}\tAIC: {:.2f}'.format(log,con,aic))
    print('LL ratio test: ')
    a=mod.log_likelihood_ratio_test().summary
    print('p: {:.2f}\t\t-log(p): {:.2f}'.format(a.loc[0,'p'],a.loc[0,'-log2(p)']))

    mod.check_assumptions(df)
    print('')
    return summary[cols]

# %%
clin = ['Age','Gender','APACHE','ARDS','Sepsis3','diabetes','malignancy','immunosuppressed']
abx=['Zosyn','Vancomycin','Cipro','Cefepime','Ceftriaxone','Zithromax']
micro = ['mbal','M','PC1','PC1bin']

# %%
meta['M']=0
meta['PC1bin']=0
meta.loc[meta['mbal']>2,'M']=1
meta.loc[meta['PC1']<0,'PC1bin']=1
# %%
meta['D24']=0
meta.loc[meta['Death']<24,'D24']=1
# %%
m = pd.DataFrame(columns=m.columns)
m.index.name='covariate'
for i in micro:
    x=cox(i,pen=.01,event='D24')
    m.loc[i]=0
    for c in x.columns:
        m.loc[i,c]+=x.loc[i,c]
# %%
#m.to_csv('../survival/univariate_clinical_d24.csv')
#m.to_csv('../survival/univariate_abx_d24.csv')
m.to_csv('../survival/univariate_micro_24.csv')

# %%
m

# %%
