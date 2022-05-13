# %%
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use("../lance.txt")
meta = pd.read_csv('../meta/meta52_current.csv',index_col=0)

v = [
 'Gender',
 'Age',
 'ARDS',
 'Shock',
 'APACHE',
 'inhosp',
 'month',
 'year',
 'Death',
 'ICU',
 'Vent',
 'PF',
 'malignancy',
 'diabetes',
 'immunosuppressed',
 'Ceftriaxone',
 'Zosyn',
 'Vancomycin',
 'Cipro',
 'Zithromax',
 'M',
 'ACE (R)',
 'Chao1 (R)',
 'Fisher (D)',
 'Pielou (E)',
 'Shannon (D)',
 'Dominance',
 'Simpson (D)',
 'Simpson (E)',
 'Sepsis3',
 'D24',
 'mmi']
vdict = {i: i for i in v}
vdict['mmi']='MMI'
vdict['M']='MMI > 0'
vdict['D24']='Death < 24 days'
vdict['month']='Death < 28 days'
vdict['year']='Death < 1 year'
vdict['inhosp']='In-hospital death'
vdict['Cipro']='Ciprofloxacin'
vdict['Zithromax']='Azithromycin'
vdict['Zosyn']='Piperacillin-tazobactam'
vdict['immunosuppressed']='Immunosuppressed'
vdict['diabetes']='Diabetes'
vdict['malignancy']='Malignancy'
vdict['PF']='PF ratio'
vdict['Vent']='Ventilator-free days'
vdict['ICU']='ICU-free days'
vdict['Death']='Death-free days'
vdict['APACHE']='APACHE II score'
vdict['Sepsis3']='Sepsis-3'
vdict['Gender']='Gender (male)'
vdict['Shock']='Septic shock'
vdict = {v: k for k,v in vdict.items()}
# %%
def parserank(csv='../bivariate/ranksums.csv',var='Death < 28 days'):
    table = pd.read_csv(csv,index_col=None)
    rows = [i for i in table.index if table.iloc[i,0]==var or table.iloc[i,1]==var]
    table = table.loc[rows]
    cols = ['Test variable',
            '{} (n={}), Median[IQR]'.format(var,table.iloc[0,2]),
            '{} (n={}), Median[IQR]'.format(var.replace('>','<'),table.iloc[0,4]),
            'U statistic',
            'P value',
    ]
    out = pd.DataFrame(columns=cols)
    for i in table.index:
        row = [table.loc[i,'Test variable'],
               table.loc[i,'Median[IQR] (+)'],
               table.loc[i,'Median[IQR] (-)'],
               '{: .2g}'.format(table.loc[i,'U statistic']),
               '{: .2g}'.format(table.loc[i,'P value']),
        ]
        out.loc[len(out.index),cols] = row
    return out
#x=parserank()
#x.to_csv('../bivariate/pretty_ranksums_month.csv',index=None)
x=parserank(var='MMI > 0')
x.to_csv('../bivariate/pretty_ranksums_mmibinary.csv',index=None)
#x.to_csv('../bivariate/pretty_ranksums_24days.csv',index=None)
#parsebiv(var='MMI > 0').to_csv('../bivariate/pretty_ranksums_mmibinary.csv',index=None)
x


# %%
def parserankc(csv='../bivariate/ranksums.csv',var='Death-free days'):
    table = pd.read_csv(csv,index_col=None)
    rows = [i for i in table.index if table.loc[i,'Test variable']==var]
    table = table.loc[rows]
    cols = ['Group variable',
            'n1 (+)',
            'n2 (-)',
            '{} (+), Median[IQR]'.format(var),
            '{} (-), Median[IQR]'.format(var),
            'U statistic',
            'P value',
        ]
    out = pd.DataFrame(columns=cols)
    for i in table.index:
        row = [table.loc[i,'Group variable'],
               table.loc[i,'n1 (+)'],
               table.loc[i,'n2 (-)'],
               table.loc[i,'Median[IQR] (+)'],
               table.loc[i,'Median[IQR] (-)'],
               '{: .2g}'.format(table.loc[i,'U statistic']),
               '{: .2g}'.format(table.loc[i,'P value']),
            ]
        out.loc[len(out.index),cols] = row
    return out
x=parserankc()
x
# %%
x.to_csv('../bivariate/pretty_ranksums_deathfree.csv',index=None)
#x.to_csv('../bivariate/pretty_ranksums_mmi.csv',index=None)
# %%
def parsefish(csv='../bivariate/fisher.csv',var='Death < 28 days'):
    table = pd.read_csv(csv,index_col=None)
    rows = [i for i in table.index if var in table.loc[i,['Variable 1','Variable 2']].to_list()]
    table = table.loc[rows]
    if table.iloc[0,0]==var:
        prev = table.iloc[0,2]
    else:
        prev = table.iloc[0,3]
    cols = ['Variable',
            'Prevalence (%)',
            'With {} (Prevalence={}), %'.format(var,prev),
            'Neither, %',
            'Odds ratio',
            'P value',
        ]
    out = pd.DataFrame(columns=cols)
    for i in table.index:
        other = table.loc[i,'Variable 1'] if var==table.loc[i,'Variable 2'] else table.loc[i,'Variable 2']
        prev = table.loc[i,'N1 (%)'] if var==table.loc[i,'Variable 2'] else table.loc[i,'N2 (%)']
        both = (meta[vdict.get(other)]+meta[vdict.get(var)]==2).sum()/52
        neither = (meta[vdict.get(other)]+meta[vdict.get(var)]==0).sum()/52
        row = [other,
               prev[:-1],
               '{:.0%}'.format(both)[:-1],
               '{:.0%}'.format(neither)[:-1],
               '{: .2g}'.format(table.loc[i,'Odds ratio']),
               '{: .2g}'.format(table.loc[i,'P value']),
            ]
        out.loc[len(out.index),cols] = row
    return out
x=parsefish(var='MMI > 0')
x
# %%
x.to_csv('../bivariate/pretty_fish_mmibin.csv',index=None)

# %%
def parsespear(csv='../bivariate/spearman.csv',var='Death-free days'):
    table = pd.read_csv(csv,index_col=None)
    rows = [i for i in table.index if table.iloc[i,0]==var or table.iloc[i,1]==var]
    table = table.loc[rows]
    cols = ['Variable',
            'Spearman corr. ({})'.format(var),
            'P value',
    ]
    out = pd.DataFrame(columns=cols)
    for i in table.index:
        other = table.loc[i,'Variable 1'] if var==table.loc[i,'Variable 2'] else table.loc[i,'Variable 2']
        row = [other,
               '{: .2f}'.format(table.loc[i,'Spearman corr.']),
               '{:.2g}'.format(table.loc[i,'P value']),
        ]
        out.loc[len(out.index),cols] = row
    return out
x=parsespear(var='MMI')
x
# %%
x.to_csv('../bivariate/pretty_spear_mmi.csv',index=None)
# %%
