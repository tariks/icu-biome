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

# fish 2

cols = [
    'Variable',
    '(+) & Death < 28 days',
    '(-) & Death < 28 days',
    'Odds ratio',
    'P value'
]
groups = [j for j in cats if j not in outcomes]
out = pd.DataFrame(index=groups,columns=cols)
for i in groups:
    a = meta['month']
    b = meta[i]
    row = [
           vdict.get(i),
           (a+b == 2).astype(int).sum().astype(str),
           a[b==0].sum().astype(int).astype(str),
    ]
    ctab = stats.contingency.crosstab(a,b)
    teststat, pval = stats.fisher_exact(ctab[1])
    row.append('{:.2g}'.format(teststat))
    row.append('{:.2g}'.format(pval))
    out.loc[i] = row
out.index=out['Variable'].copy()
del out['Variable']
out

# %%
out.to_csv("../bivariate/fisher2.csv", index=None)

# %%
