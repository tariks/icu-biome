# %%
import pandas as pd

meta = pd.read_csv("../meta/meta52_current.csv", index_col=0)

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
# vdict = {v: k for k, v in vdict.items()}
# %%
def parserank(csv="../bivariate/ranksums.csv", var="Death < 28 days"):
    table = pd.read_csv(csv, index_col=None)
    rows = [
        i for i in table.index if table.iloc[i, 0] == var or table.iloc[i, 1] == var
    ]
    table = table.loc[rows]
    cols = [
        "Test variable",
        "{} (n={}), Median[IQR]".format(var, table.iloc[0, 2]),
        "{} (n={}), Median[IQR]".format(var.replace("<", ">"), table.iloc[0, 4]),
        "U statistic",
        "P value",
    ]
    out = pd.DataFrame(columns=cols)
    for i in table.index:
        row = [
            table.loc[i, "Test variable"],
            table.loc[i, "Median[IQR] (+)"],
            table.loc[i, "Median[IQR] (-)"],
            "{: .2g}".format(table.loc[i, "U statistic"]),
            "{:.2g}".format(table.loc[i, "P value"]),
        ]
        out.loc[len(out.index), cols] = row
    return out


x = parserank()
x.to_csv("../bivariate/pretty_ranksums_month.csv", index=None)
# x=parserank(var='Death < 24 days')
# x.to_csv('../bivariate/pretty_ranksums_mmibinary.csv',index=None)
# x.to_csv('../bivariate/pretty_ranksums_24days.csv',index=None)
x


# %%
def parserankc(csv="../bivariate/ranksums.csv", var="Death-free days"):
    table = pd.read_csv(csv, index_col=None)
    rows = [i for i in table.index if table.loc[i, "Test variable"] == var]
    table = table.loc[rows]
    cols = [
        "Group variable",
        "n1 (+)",
        "n2 (-)",
        "{} (+), Median[IQR]".format(var),
        "{} (-), Median[IQR]".format(var),
        "U statistic",
        "P value",
    ]
    out = pd.DataFrame(columns=cols)
    for i in table.index:
        row = [
            table.loc[i, "Group variable"],
            table.loc[i, "n1 (+)"],
            table.loc[i, "n2 (-)"],
            table.loc[i, "Median[IQR] (+)"],
            table.loc[i, "Median[IQR] (-)"],
            "{: .2g}".format(table.loc[i, "U statistic"]),
            "{:.2g}".format(table.loc[i, "P value"]),
        ]
        out.loc[len(out.index), cols] = row
    return out


x = parserankc(var="MMI")
x
# %%
# x.to_csv('../bivariate/pretty_ranksums_deathfree.csv',index=None)
x.to_csv("../bivariate/pretty_ranksums_mmi.csv", index=None)
# %%
def parsefish(csv="../bivariate/fisher.csv", var="Death < 28 days"):
    table = pd.read_csv(csv, index_col=None)
    rows = [
        i
        for i in table.index[::2]
        if var in table.loc[i, ["Variable 1", "Variable 2"]].to_list()
    ]
    table = table.loc[rows]
    if table.iloc[0, 0] == var:
        prev = table.iloc[0, 2]
    else:
        prev = table.iloc[0, 3]
    cols = [
        "Variable",
        "n (+)",
        "(+) and {} (n={})".format(var, meta[vdict.get(var)].sum().astype(int)),
        "(-) and {} (n={})".format(
            var.replace("<", ">="), (meta[vdict.get(var)] == 0).sum().astype(int)
        ),
        "Odds ratio",
        "P value",
    ]
    out = pd.DataFrame(columns=cols)
    for i in table.index:
        other = (
            table.loc[i, "Variable 1"]
            if var == table.loc[i, "Variable 2"]
            else table.loc[i, "Variable 2"]
        )
        both = (meta[vdict.get(other)] + meta[vdict.get(var)] == 2).sum()
        neither = (meta[vdict.get(other)] + meta[vdict.get(var)] == 0).sum()
        row = [
            other,
            meta[vdict.get(other)].sum().astype(int),
            "{}".format(both),
            "{}".format(neither),
            "{: .2g}".format(table.loc[i, "Odds ratio"]),
            "{:.2g}".format(table.loc[i, "P value"]),
        ]
        out.loc[len(out.index), cols] = row
    return out


x = parsefish(var="Death < 24 days")
x
# %%
# x.to_csv('../bivariate/pretty_fish_mmibin.csv',index=None)
x.to_csv("../bivariate/pretty_fish_24days.csv", index=None)
# x.to_csv("../bivariate/pretty_fish_month.csv", index=None)

# %%
def parsespear(csv="../bivariate/spearman.csv", var="Death-free days"):
    table = pd.read_csv(csv, index_col=None)
    rows = [
        i for i in table.index if table.iloc[i, 0] == var or table.iloc[i, 1] == var
    ]
    table = table.loc[rows]
    cols = [
        "Variable",
        "Spearman corr. ({})".format(var),
        "P value",
    ]
    out = pd.DataFrame(columns=cols)
    for i in table.index:
        other = (
            table.loc[i, "Variable 1"]
            if var == table.loc[i, "Variable 2"]
            else table.loc[i, "Variable 2"]
        )
        row = [
            other,
            "{: .2f}".format(table.loc[i, "Spearman corr."]),
            "{:.2g}".format(table.loc[i, "P value"]),
        ]
        out.loc[len(out.index), cols] = row
    return out


x = parsespear(var="MMI")
x
# %%
x.to_csv("../bivariate/pretty_spear_mmi.csv", index=None)
# %%
def parsecox(
    csv="../survival/bivariate.csv",
):
    table = pd.read_csv(csv, index_col=None)
    cols = [
        "Model",
        "Model LL",
        "Variable",
        "Hazard ratio",
        "HR 95% CI, lower",
        "HR 95% CI, higher",
        "P value",
    ]
    out = pd.DataFrame(columns=cols)
    current_model = ""
    for i in table.index:
        if "*" not in table.loc[i, "model"]:
            model = table.loc[i, "model"].split("+")
            model = [vdict.get(i) for i in model]
            for j in range(len(model)):
                if " " in model[j]:
                    model[j] = "(" + model[j] + ")"
            model = " + ".join(model)
            if model == current_model:
                m = ["..", ".."]
            else:
                m = [model, "{: .1f}".format(table.loc[i, "model LL"])]
                current_model = model
            other = table.loc[i][-5:].to_list()
            other[0] = vdict.get(other[0])
            for j in range(1, 4):
                other[j] = "{:.1f}".format(other[j])
            other[-1] = "{:.2g}".format(other[-1])
            row = m + other
            out.loc[len(out.index), cols] = row
    return out


# %%
out = parsecox("../survival/quadvariate.csv")
# %%
out.iloc[:40, :]

# %%
out = out.iloc[:40, :]
out.to_csv("../survival/pretty_quad.csv", index=None)

# %%
def parsecox(
    csv="../survival/bivariate.csv",
):
    table = pd.read_csv(csv, index_col=None)
    cols = [
        "Model",
        "Variable",
        "Hazard ratio (95% CI)",
        "P value",
    ]
    out = pd.DataFrame(columns=cols)
    current_model = ""
    for i in table.index:
        if "*" not in table.loc[i, "model"]:
            model = table.loc[i, "model"].split("+")
            model = [vdict.get(i) for i in model]
            for j in range(len(model)):
                if " " in model[j]:
                    model[j] = "(" + model[j] + ")"
            model = " + ".join(model)
            if model == current_model:
                m = [".."]
            else:
                m = [model]
                current_model = model
            other = table.loc[i][-5:].to_list()
            other[0] = vdict.get(other[0])
            for j in range(1, 4):
                other[j] = "{:.1f}".format(other[j])
            other[-1] = "{:.2g}".format(other[-1])
            other = (
                [other[0]]
                + ["{} ({}, {})".format(other[1], other[2], other[3])]
                + [other[-1]]
            )
            row = m + other
            out.loc[len(out.index), cols] = row
    return out


# %%

x4 = parsecox("../survival/quadvariate.csv")
x2 = parsecox("../survival/bivariate.csv")
x3 = parsecox("../survival/trivariate.csv")
# x1 = parsecox('../survival/univariate.csv')

x2
# %%
x = pd.read_csv("../survival/univariate.csv", index_col=0)
x

# %%
x["Model"] = [vdict.get(i) for i in x.index]
x["Variable"] = [vdict.get(i) for i in x.index]
x["Hazard ratio (95% CI)"] = [
    "{: .2g} ({:.2g},{: .2g})".format(
        x.loc[i, "exp(coef)"],
        x.loc[i, "exp(coef) lower 95%"],
        x.loc[i, "exp(coef) upper 95%"],
    )
    for i in x.index
]
x["P value"] = ["{:.2g}".format(x.loc[i, "p"]) for i in x.index]
x = x[["Model", "Variable", "Hazard ratio (95% CI)", "P value"]]
x

# %%
x = x.append(x2, ignore_index=True)
x = x.append(x3, ignore_index=True)
x = x.append(x4, ignore_index=True)
x

# %%

x.to_csv("../survival/pretty_cox.csv", index=None)
# %%
