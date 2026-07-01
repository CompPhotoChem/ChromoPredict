#!/usr/bin/env python3
"""
chromopredict_batch.py
Batch application of the ChromoPredict Woodward-Fieser / Fieser / Fieser-Kuhn
rules to a table of molecules (.xlsx / .csv / .tsv).

Pipeline:
  1. read an .xlsx (or .csv) table with a SMILES column
  2. optionally write the table out as .csv (xlsx -> csv conversion)
  3. run chromopredict.predict() per molecule
       -> new column  lambda_max_calc_nm   (predicted absorption maximum)
       -> new column  rule_system_used     (which empirical rule was applied)
  4. if an experimental column is present: compute error columns
  5. optionally repeat the prediction for several forced rule systems
     and write an error plot (signed error per compound + calc-vs-exp scatter)

Examples:
  python chromopredict_batch.py -input structures.xlsx -to-csv structures.csv \
         -out predictions.csv -plot error.png -rules auto woodward woodward_refine
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# rule-category logic
_CATS = {
    "woodward": {"woodward", "woodward_extended", "woodward_refine", "woodward_coumarin"},
    "fieser": {"fieser"},
    "fieser_kuhn": {"fieser_kuhn"},
}
_RULE2CAT = {rule: cat for cat, rules in _CATS.items() for rule in rules}


def effective_rule(auto, requested):
    """Return the rule that predict() will actually use.

    predict() auto-detects `auto`, if the user forces `requested` from a
    different category, the package warns and falls back to `auto`.
    """
    if not requested:
        return auto
    if auto is None:
        return None
    if _RULE2CAT.get(auto) == _RULE2CAT.get(requested):
        return requested
    return auto

def load_chromopredict(src):
    """Import chromopredict, optionally from a local src/ folder."""
    try:
        import chromopredict  # installed (pip install -e .)
        return chromopredict
    except ImportError:
        if src:
            srcp = Path(src).expanduser().resolve()
            if srcp.exists():
                sys.path.insert(0, str(srcp))
                import chromopredict
                return chromopredict
        raise SystemExit(
            "ERROR: could not import 'chromopredict'.\n"
        )

def read_table(path, sheet):
    p = Path(path)
    ext = p.suffix.lower()
    if ext in (".xlsx", ".xlsm", ".xls"):
        return pd.read_excel(p, sheet_name=sheet if sheet is not None else 0)
    if ext == ".tsv" or ext == ".txt":
        return pd.read_csv(p, sep="\t")
    return pd.read_csv(p)


def _norm(s):
    return str(s).strip().lower().replace(" ", "").replace("_", "")


def pick_col(df, wanted, aliases):
    """Find a column by name/alias, case- and space-insensitive."""
    norm = {_norm(c): c for c in df.columns}
    for cand in [wanted, *aliases]:
        if cand is None:
            continue
        if _norm(cand) in norm:
            return norm[_norm(cand)]
    return None

def predict_one(cp, classify_type, general_rules, smiles, solvent, requested):
    """Return (lambda_max_nm, rule_used, note)."""
    smiles = str(smiles).strip()
    try:
        auto, _ = classify_type(smiles, general_rules)
    except Exception:
        auto = None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            pred, _, _ = cp.predict(smiles, solvent=solvent,
                                    verbose=False, chromlib=requested)
        except Exception as exc:
            return np.nan, (auto or "unclassified"), type(exc).__name__
    eff = effective_rule(auto, requested)
    note = "fallback" if (requested and eff == auto and eff != requested) else ""
    return float(pred), (eff or "unclassified"), note

def make_plot(df, labels, exp, rule_specs, outfile):
    """rule_specs: list of (col_prefix, display_name)."""
    n = len(df)
    idx = np.arange(n)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(max(9, 1.1 * n + 5), 5.2))
    colors = plt.cm.tab10(np.linspace(0, 1, max(3, len(rule_specs))))

    m = len(rule_specs)
    width = 0.8 / m
    stats_lines = []
    for k, (pref, disp) in enumerate(rule_specs):
        err = df[f"error_nm[{pref}]"].to_numpy(dtype=float)
        ax1.bar(idx + (k - (m - 1) / 2) * width, err, width,
                label=disp, color=colors[k], edgecolor="black", linewidth=0.4)
        mae = np.nanmean(np.abs(err))
        rmse = np.sqrt(np.nanmean(err ** 2))
        me = np.nanmean(err)
        stats_lines.append(f"{disp:16s}  MAE={mae:5.1f}  RMSE={rmse:5.1f}  ME={me:+5.1f}")

    ax1.axhline(0, color="black", lw=0.8)
    ax1.set_xticks(idx)
    ax1.set_xticklabels(labels, rotation=40, ha="right", fontsize=8)
    ax1.set_ylabel("signed error  (calc - exp) / nm")
    ax1.set_title("Prediction error per compound")
    ax1.legend(fontsize=8, framealpha=0.9)
    ax1.grid(axis="y", ls=":", alpha=0.5)

    allvals = [exp]
    for pref, _ in rule_specs:
        allvals.append(df[f"lambda_max_calc_nm[{pref}]"].to_numpy(dtype=float))
    allvals = np.concatenate([np.asarray(a, float) for a in allvals])
    lo = np.nanmin(allvals) - 15
    hi = np.nanmax(allvals) + 15
    ax2.plot([lo, hi], [lo, hi], color="black", lw=1, label="y = x")
    ax2.fill_between([lo, hi], [lo - 20, hi - 20], [lo + 20, hi + 20],
                     color="grey", alpha=0.12, label="±20 nm")
    for k, (pref, disp) in enumerate(rule_specs):
        calc = df[f"lambda_max_calc_nm[{pref}]"].to_numpy(dtype=float)
        ax2.scatter(exp, calc, s=55, color=colors[k], edgecolor="black",
                    linewidth=0.5, label=disp, zorder=3)
    ax2.set_xlim(lo, hi); ax2.set_ylim(lo, hi)
    ax2.set_xlabel("experimental λmax / nm")
    ax2.set_ylabel("calculated λmax / nm")
    ax2.set_title("Calculated vs. experimental")
    ax2.legend(fontsize=8, framealpha=0.9)
    ax2.grid(ls=":", alpha=0.5)
    ax2.set_aspect("equal", adjustable="box")

    fig.tight_layout()
    fig.savefig(outfile, dpi=160)
    plt.close(fig)
    return stats_lines

def build_parser():
    p = argparse.ArgumentParser(
        description="Batch Woodward-Fieser / Fieser / Fieser-Kuhn lambda_max prediction.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-input", required=True, help="input table (.xlsx / .csv / .tsv)")
    p.add_argument("-out", default="predictions.csv", help="output CSV with predictions")
    p.add_argument("-to-csv", default=None,
                   help="also write the raw input table to this CSV (xlsx -> csv)")
    p.add_argument("-plot", default="error.png", help="output error plot (.png); '' to skip")
    p.add_argument("-sheet", default=None, help="sheet name/index for xlsx input")

    p.add_argument("-smiles-col", default=None, help="SMILES column name (auto-detected)")
    p.add_argument("-name-col", default=None, help="name/label column (auto-detected)")
    p.add_argument("-cas-col", default=None, help="CAS column (optional)")
    p.add_argument("-exp-col", default=None,
                   help="experimental lambda_max column (auto-detected; enables error columns)")

    p.add_argument("-solvent", default=None,
                   help="solvent SMILES for the correction term (e.g. CCO for ethanol)")
    p.add_argument("-rules", nargs="+", default=["auto"],
                   help="rule systems to evaluate. 'auto' = package auto-classification. "
                        "others: woodward woodward_refine woodward_extended "
                        "woodward_coumarin fieser fieser_kuhn")
    p.add_argument("-chromopredict-src", default=None,
                   help="path to the package 'src' folder if not pip-installed")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)

    cp = load_chromopredict(args.chromopredict_src)
    from chromopredict.strucfeatures import classify_type
    from chromopredict.chrombase import general_rules

    df = read_table(args.input, args.sheet)
    if args.to_csv:
        df.to_csv(args.to_csv, index=False)
        print(f"wrote CSV copy of input -> {args.to_csv}")

    smi_col = args.smiles_col or pick_col(df, "SMILES", ["SMILE", "smiles", "canonical_smiles"])
    if smi_col is None:
        raise SystemExit(f"ERROR: no SMILES column found. Columns: {list(df.columns)}\n"
                         f"  -> pass -smiles-col <name>")
    name_col = args.name_col or pick_col(df, "Name", ["compound", "molecule", "title", "label"])
    cas_col = args.cas_col or pick_col(df, "CAS", ["cas_no", "casrn"])
    exp_col = args.exp_col or pick_col(
        df, "lambda_max_experimental",
        ["λmax experimental", "lambdamaxexp", "expnm", "experimental",
         "lambda_max_exp", "lmax_exp", "abs_max_exp"])

    labels = df[name_col].astype(str).tolist() if name_col else \
        (df[cas_col].astype(str).tolist() if cas_col else
         [f"mol_{i}" for i in range(len(df))])

    rule_specs = []
    for r in args.rules:
        if r.lower() == "auto":
            rule_specs.append(("auto", "auto", None))
        else:
            rule_specs.append((r, r, r))

    print(f"{len(df)} molecules | SMILES='{smi_col}' | exp='{exp_col}' | "
          f"rules={[d for _, d, _ in rule_specs]}")

    for pref, disp, requested in rule_specs:
        preds, used, notes = [], [], []
        for smi in df[smi_col]:
            lam, rule, note = predict_one(cp, classify_type, general_rules,
                                          smi, args.solvent, requested)
            preds.append(lam); used.append(rule); notes.append(note)
        df[f"lambda_max_calc_nm[{pref}]"] = preds
        df[f"rule_system_used[{pref}]"] = used
        if exp_col is not None:
            exp = pd.to_numeric(df[exp_col], errors="coerce")
            df[f"error_nm[{pref}]"] = np.asarray(preds, float) - exp.to_numpy(float)
            df[f"abs_error_nm[{pref}]"] = df[f"error_nm[{pref}]"].abs()

    prim = rule_specs[0][0]
    df["lambda_max_calc_nm"] = df[f"lambda_max_calc_nm[{prim}]"]
    df["rule_system_used"] = df[f"rule_system_used[{prim}]"]
    if exp_col is not None:
        df["error_nm"] = df[f"error_nm[{prim}]"]
        df["abs_error_nm"] = df[f"abs_error_nm[{prim}]"]

    df.to_csv(args.out, index=False)
    print(f"wrote predictions -> {args.out}")

    if exp_col is not None and args.plot:
        exp = pd.to_numeric(df[exp_col], errors="coerce").to_numpy(float)
        specs = [(pref, disp) for pref, disp, _ in rule_specs]
        stats = make_plot(df, labels, exp, specs, args.plot)
        print(f"wrote plot -> {args.plot}")
        print("\n=== error statistics (nm) ===")
        for line in stats:
            print("   " + line)
    elif exp_col is None:
        print("no experimental column detected -> skipped error columns/plot "
              "(pass -exp-col <name> to enable)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
