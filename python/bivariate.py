# %%
import pandas as pd
from scipy import stats

meta = pd.read_csv("../meta/meta52_current.csv", index_col=0)

# %%
v = [
    "Gender",
    "Age",
    "ARDS",
    "Shock",
    "APACHE",
    "inhosp",
    "month",
    "year",
    "Death",
    "ICU",
    "Vent",
    "PF",
    "malignancy",
    "diabetes",
    "immunosuppressed",
    "Ceftriaxone",
    "Zosyn",
    "Vancomycin",
    "Cipro",
    "Zithromax",
    "M",
    "ACE (R)",
    "Chao1 (R)",
    "Fisher (D)",
    "Pielou (E)",
    "Shannon (D)",
    "Dominance",
    "Simpson (D)",
    "Simpson (E)",
    "Sepsis3",
    "D24",
    "mmi",
]
vdict = {i: i for i in v}
vdict["mmi"] = "MMI"
vdict["M"] = "MMI > 0"
vdict["D24"] = "Death < 24 days"
vdict["month"] = "Death < 28 days"
vdict["year"] = "Death < 1 year"
vdict["inhosp"] = "In-hospital death"
vdict["Cipro"] = "Ciprofloxacin"
vdict["Zithromax"] = "Azithromycin"
vdict["Zosyn"] = "Piperacillin-tazobactam"
vdict["immunosuppressed"] = "Immunosuppressed"
vdict["diabetes"] = "Diabetes"
vdict["malignancy"] = "Malignancy"
vdict["PF"] = "PF ratio"
vdict["Vent"] = "Ventilator-free days"
vdict["ICU"] = "ICU-free days"
vdict["Death"] = "Death-free days"
vdict["APACHE"] = "APACHE II score"
vdict["Sepsis3"] = "Sepsis-3"
vdict["Gender"] = "Gender (male)"
vdict["Shock"] = "Septic shock"
# %%
contin = [
    "Age",
    "APACHE",
    "Death",
    "ICU",
    "Vent",
    "PF",
    "ACE (R)",
    "Chao1 (R)",
    "Fisher (D)",
    "Pielou (E)",
    "Shannon (D)",
    "Dominance",
    "Simpson (D)",
    "Simpson (E)",
    "mmi",
]
outcomes = ["Death", "ICU", "Vent", "month", "inhosp", "year", "D24"]
mmis = ["mmi", "M"]
cats = [i for i in v if i not in contin]
# %%
# ranksums
cols = [
    "Test variable",
    "Group variable",
    "n1 (+)",
    "Median[IQR] (+)",
    "n2 (-)",
    "Median[IQR] (-)",
    "U statistic",
    "P value",
]
out = pd.DataFrame(columns=cols)
for t in contin:
    for g in cats:
        if t in outcomes and g in outcomes:
            pass
        elif t in mmis and g in mmis:
            pass
        else:
            row = [vdict.get(t), vdict.get(g)]
            a, b = meta.loc[meta[g] == 1, t], meta.loc[meta[g] == 0, t]
            for i in [a, b]:
                row += [
                    i.size,
                    "{: ,.2g} [{: ,.2g}-{: ,.2g}]".format(
                        i.median(), i.quantile(0.25), i.quantile(0.75)
                    ),
                ]
            teststat, pval = stats.ranksums(a, b, nan_policy="omit")
            row += [teststat, pval]
            out.loc[len(out.index)] = row
out
# %%
out = out.loc[out["P value"].sort_values().index]
out.to_csv("../bivariate/ranksums.csv", index=None)
# %%
# corr
cols = [
    "Variable 1",
    "Variable 2",
    "Median[IQR] (V1)",
    "Median[IQR] (V2)",
    "Spearman corr.",
    "P value",
]
out = pd.DataFrame(columns=cols)
for t in contin:
    for g in contin:
        if t in outcomes and g in outcomes:
            pass
        elif t == g:
            pass
        else:
            row = [vdict.get(t), vdict.get(g)]
            a, b = meta[t], meta[g]
            for i in [a, b]:
                row += [
                    "{: ,.2g} [{: ,.2g}-{: ,.2g}]".format(
                        i.median(), i.quantile(0.25), i.quantile(0.75)
                    )
                ]
            teststat, pval = stats.spearmanr(a, b, nan_policy="omit")
            row += [teststat, pval]
            out.loc[len(out.index)] = row
out
# %%
out = out.loc[out["P value"].sort_values().index]
out.to_csv("../bivariate/spearman.csv", index=None)
out


# %%
# fish
cols = [
    "Variable 1",
    "Variable 2",
    "N1 (%)",
    "N2 (%)",
    "Odds ratio",
    "P value",
]
out = pd.DataFrame(columns=cols)
for t in cats:
    for g in cats:
        if t in outcomes and g in outcomes:
            pass
        elif t == g:
            pass
        else:
            row = [vdict.get(t), vdict.get(g)]
            a, b = meta[t], meta[g]
            for i in [a, b]:
                row += ["{:.0%}".format(i.sum() / i.size)]
            _, ctab = stats.contingency.crosstab(
                a,
                b,
            )
            teststat, pval = stats.fisher_exact(ctab)
            row += [teststat, pval]
            out.loc[len(out.index)] = row
out
# %%
out = out.loc[out["P value"].sort_values().index]
out.to_csv("../bivariate/fisher.csv", index=None)
out
# %%
def getqs(v):
    vv = v.describe()
    qs = ['50%','25%','75%']
    vv = ['{: .1f}'.format(vv[i]) for i in qs]
    return vv[0] + ' ['+vv[1]+','+vv[2]+']'
# %%

# fish 2
meta['MM'] = 1-meta['M']
vdict['MM'] = 'MMI < 0'

cols = [
    'Variable',
    'Median[IQR] or n',
    'Death < 28 days',
    'Death > 28 days',
    'P value'
]
groups = [j for j in cats if j not in outcomes]
groups = ['Age','Gender','APACHE','ARDS','malignancy','immunosuppressed','diabetes','Sepsis3','Vancomycin','Ceftriaxone','Cipro','Zosyn',
        'mmi','M','Shannon (D)','Chao1 (R)']
out = pd.DataFrame(index=groups,columns=cols)
for i in groups:
    if i in cats:
        a = meta['month']
        b = meta[i]
        asum = b[a==1].sum().astype(int),
        bsum = b[a==0].sum().astype(int),
        tsum = b.sum().astype(int),
        row = [
            vdict.get(i),
            '{} ({:.0%})'.format(tsum[0],tsum[0]/52),
            '{} ({:.0%})'.format(asum[0],asum[0]/52),
            '{} ({:.0%})'.format(bsum[0],bsum[0]/52),
        ]
        ctab = stats.contingency.crosstab(a,b)
        teststat, pval = stats.fisher_exact(ctab[1])
        row.append('{:.2g}'.format(pval))
        out.loc[i] = row
    else:
        a = meta.loc[meta['month']==1,i]
        b = meta.loc[meta['month']==0,i]
        row = [
            vdict.get(i),
            getqs(meta[i]),
            getqs(a),
            getqs(b),
        ]
        teststat, pval = stats.ranksums(a,b,nan_policy='omit')
        row.append('{:.2g}'.format(pval))
        out.loc[i] = row

#out.index=out['Variable'].copy()
#del out['Variable']
out

# %%
out.to_csv("../bivariate/fisher3.csv", index=None)

# %%
