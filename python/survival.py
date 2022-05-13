# %%
import pandas as pd
import matplotlib.pyplot as plt
import lifelines

plt.style.use("../lance.txt")
meta = pd.read_csv("../meta/meta52_current.csv", index_col=0)
# %%
import warnings

warnings.simplefilter("ignore")

# %%
def cox(v, df=meta, pen=0.01, l1=0.5, event="month", robust=True):
    f = "+".join(v) if type(v) == list else v
    modparams = dict(penalizer=pen, l1_ratio=l1)
    mod = lifelines.CoxPHFitter(**modparams)
    params = dict(
        robust=robust,
        formula=f,
    )
    mod.fit(df, duration_col="Death", event_col=event, **params)
    summary = mod.summary
    cols = ["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]

    print("\t".join(["var", "exp", "low", "hi", "p"]))
    for i in summary.index:
        vals = summary.loc[i, cols]
        print(i + "\t" + "\t".join(["{:.3f}".format(vals[j]) for j in cols]))

    aic = mod.AIC_partial_
    log = mod.log_likelihood_
    con = mod.concordance_index_
    print("\nK-fold:(LL, CI)")
    #    for i in ['log_likelihood','concordance_index']:
    #        val = lifelines.utils.k_fold_cross_validation(
    #                    lifelines.CoxPHFitter(**modparams),
    #                    df,
    #                    duration_col='Death',
    #                    event_col=event,
    #                    scoring_method=i,
    #                    seed=0,
    #                    fitter_kwargs=params
    #            )
    #        val = ['{: .2f}'.format(j) for j in val]
    #        print('\t'.join(val))
    print("")
    print("LL: {:.2f}\t CI: {:.2f}\tAIC: {:.2f}".format(log, con, aic))
    print("LL ratio test: ")
    a = mod.log_likelihood_ratio_test().summary
    print("p: {:.2f}\t\t-log(p): {:.2f}".format(a.loc[0, "p"], a.loc[0, "-log2(p)"]))

    #    mod.check_assumptions(df)
    print("")
    summary["model LL"] = log
    summary["model AIC"] = aic
    summary["Concordance"] = con
    return summary[cols + ["model LL", "model AIC", "Concordance"]]


# %%
clin = [
    "Age",
    "Gender",
    "APACHE",
    "ARDS",
    "Sepsis3",
    "diabetes",
    "malignancy",
    "immunosuppressed",
]
abx = [
    "Zosyn",
    "Vancomycin",
    "Cipro",
    "Ceftriaxone",
]
micro = [
    "mmi",
    "M",
]
v = clin + abx + micro

# %%
f = ["APACHE"]
m = cox(["M"] + f, pen=0.1, event="month")
m = cox(["mmi"] + f, pen=0.1, event="month")
m
# %%
m = pd.DataFrame(columns=m.columns)
for i in v:
    x = cox(
        i,
        pen=0.1,
        robust=True,
    )  # event='D24')
    m.loc[i] = 0
    for c in x.columns:
        m.loc[i, c] += x.loc[i, c]
m = m.loc[m["p"].sort_values().index]
# %%
m
# %%
m.to_csv("../survival/univariate.csv")

# %%
from itertools import combinations

cols = [
    "model",
    "variable",
    "exp(coef)",
    "exp(coef) lower 95%",
    "exp(coef) upper 95%",
    "model LL",
    "model AIC",
    "Concordance",
    "p",
]
# %%
m = pd.DataFrame(columns=cols)
for combo in combinations(v, 2):
    for op in ["+", "*"]:
        f = op.join(list(combo))
        x = cox(f, pen=0.1, robust=True, event="month")
        if x["model LL"][0] > -100:
            x["model"] = f
            x["variable"] = x.index.copy().to_list()
            x = x[cols]
            idx = x.index.to_list()
            if len(idx) == 3 and x.iloc[-1, -1] > 0.05:
                pass
            else:
                for row in x.index:
                    m.loc[len(m.index)] = x.loc[row].copy()
# %%
m = m.loc[m["model LL"].sort_values(ascending=False).index]
# m = m.loc[m['model AIC']<200]
m = m.loc[m["model LL"] > -99]
m
# %%
m.to_csv("../survival/bivariate.csv", index=False)

# %%
m = pd.DataFrame(columns=cols)
for combo in combinations(v, 4):
    if "mmi" in combo and "M" in combo:
        pass
    else:
        f = "+".join(list(combo))
        x = cox(f, pen=0.1, robust=True, event="month")
        if x.iloc[0, -3] > -98 and x.iloc[0, -2] < 200:
            x["model"] = f
            x["variable"] = x.index.copy().to_list()
            x = x[cols]
            idx = x.index.to_list()
            for row in x.index:
                m.loc[len(m.index)] = x.loc[row].copy()
# %%
m = m.loc[m["model LL"].sort_values(ascending=False).index]
m
# %%
# m.iloc[:100, :].to_csv("../survival/trivariate.csv", index=False)
# m.iloc[:100, :].to_csv("../survival/quadvariate.csv", index=False)

m.loc[m["model LL"] > -95].to_csv("../survival/quadvariate.csv", index=False)

# %%
