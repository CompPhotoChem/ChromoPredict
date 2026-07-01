#!/usr/bin/env python3
"""
chromopredict_batch.py
Batch application of the ChromoPredict Woodward-Fieser / Fieser / Fieser-Kuhn
rules to a table of molecules (.xlsx / .csv / .tsv).

flags use a single dash (-input, -smiles-col, -hard-force, ...).

  python chromopredict_batch.py -input structures.xlsx -to-csv structures.csv \
         -out predictions.csv -tidy predictions_tidy.csv -plot error.png \
         -rules auto woodward fieser fieser_kuhn -hard-force -solvent CCO
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

_CATS = {
    "woodward": {"woodward", "woodward_extended", "woodward_refine", "woodward_coumarin"},
    "fieser": {"fieser"},
    "fieser_kuhn": {"fieser_kuhn"},
}
_RULE2CAT = {rule: cat for cat, rules in _CATS.items() for rule in rules}


def soft_effective_rule(auto, requested):
    """What predict() uses in normal mode: falls back to auto across categories."""
    if not requested:
        return auto
    if auto is None:
        return None
    if _RULE2CAT.get(auto) == _RULE2CAT.get(requested):
        return requested
    return auto

def load_chromopredict(src):
    try:
        import chromopredict
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
            "  -> pip install -e .   or   pass -chromopredict-src <path-to>/src"
        )


def enable_hard_force(cp):
    """Patch the private rule selector so the requested rule is always honored,
    overriding the auto-classification."""
    calc = cp.calculator
    key = "__select_mol_type"
    if hasattr(calc, key):
        setattr(calc, key, lambda auto, req: (req if req else auto))
        return True
    return False

def read_table(path, sheet):
    p = Path(path)
    ext = p.suffix.lower()
    if ext in (".xlsx", ".xlsm", ".xls"):
        return pd.read_excel(p, sheet_name=sheet if sheet is not None else 0)
    if ext in (".tsv", ".txt"):
        return pd.read_csv(p, sep="\t")
    return pd.read_csv(p)


def drop_clutter(df):
    """Remove stray unnamed columns and fully-empty columns."""
    keep = []
    for c in df.columns:
        name = str(c)
        if name.startswith("Unnamed:"):
            continue
        if df[c].isna().all():
            continue
        if df[c].dtype == object and df[c].fillna("").astype(str).str.strip().eq("").all():
            continue
        keep.append(c)
    return df[keep].copy()


def _norm(s):
    return str(s).strip().lower().replace(" ", "").replace("_", "")


def pick_col(df, wanted, aliases):
    norm = {_norm(c): c for c in df.columns}
    for cand in [wanted, *aliases]:
        if cand is not None and _norm(cand) in norm:
            return norm[_norm(cand)]
    return None

def predict_one(cp, classify_type, general_rules, smiles, solvent, requested, hard):
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
    if hard and requested:
        used = requested
        note = "forced" if requested != auto else ""
    else:
        used = soft_effective_rule(auto, requested)
        note = "fallback" if (requested and used == auto and used != requested) else ""
    return float(pred), (used or "unclassified"), note

def make_plot(tidy, rule_order, outfile):
    labels = list(dict.fromkeys(tidy["label"]))
    n = len(labels)
    idx = np.arange(n)
    colors = plt.cm.tab10(np.linspace(0, 1, max(3, len(rule_order))))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(max(9, 1.2 * n + 5), 5.2))
    m = len(rule_order)
    width = 0.8 / m
    stats = []
    for k, r in enumerate(rule_order):
        sub = tidy[tidy["rule"] == r].set_index("label").reindex(labels)
        err = sub["error_nm"].to_numpy(dtype=float)
        oob = sub["out_of_domain"].to_numpy(bool)
        err_plot = np.where(oob, np.nan, err)
        ax1.bar(idx + (k - (m - 1) / 2) * width, err_plot, width,
                label=r, color=colors[k], edgecolor="black", linewidth=0.4)
        valid = err[~oob]
        mae = np.nanmean(np.abs(valid)) if valid.size else np.nan
        rmse = np.sqrt(np.nanmean(valid ** 2)) if valid.size else np.nan
        me = np.nanmean(valid) if valid.size else np.nan
        stats.append((r, mae, rmse, me, int(oob.sum())))

    ax1.axhline(0, color="black", lw=0.8)
    ax1.set_xticks(idx)
    ax1.set_xticklabels(labels, rotation=40, ha="right", fontsize=8)
    ax1.set_ylabel("signed error  (calc − exp) / nm")
    ax1.set_title("Error per compound (flagged omitted)")
    ax1.legend(fontsize=8, framealpha=0.9)
    ax1.grid(axis="y", ls=":", alpha=0.5)

    allv = tidy.loc[~tidy["out_of_domain"], ["exp_nm", "calc_nm"]].to_numpy(float)
    allv = allv[np.isfinite(allv).all(axis=1)] if allv.size else allv
    lo, hi = (np.min(allv) - 15, np.max(allv) + 15) if allv.size else (200, 400)
    ax2.plot([lo, hi], [lo, hi], color="black", lw=1, label="y = x")
    ax2.fill_between([lo, hi], [lo - 20, hi - 20], [lo + 20, hi + 20],
                     color="grey", alpha=0.12, label="±20 nm")
    for k, r in enumerate(rule_order):
        sub = tidy[tidy["rule"] == r]
        ok = sub[~sub["out_of_domain"]]
        ax2.scatter(ok["exp_nm"], ok["calc_nm"], s=55, color=colors[k],
                    edgecolor="black", linewidth=0.5, label=r, zorder=3)
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
    return stats


def make_plot_auto_colored(tidy, outfile, base_rule="auto"):
    """Simple single-panel scatter: only the `base_rule` rows, each point colored
    by the rule the classifier actually chose for that molecule (rule_used)."""
    sub = tidy[tidy["rule"] == base_rule].copy()
    if sub.empty: 
        base_rule = tidy["rule"].iloc[0]
        sub = tidy[tidy["rule"] == base_rule].copy()
    for c in ("exp_nm", "calc_nm", "error_nm"):
        sub[c] = pd.to_numeric(sub[c], errors="coerce")

    ok = sub[~sub["out_of_domain"]]
    used_rules = sorted(r for r in sub["rule_used"].dropna().unique())
    cmap = plt.cm.tab10(np.linspace(0, 1, max(3, len(used_rules))))
    color_map = {r: cmap[i] for i, r in enumerate(used_rules)}

    fig, ax = plt.subplots(figsize=(6.6, 6.2))
    vals = ok[["exp_nm", "calc_nm"]].to_numpy(float)
    vals = vals[np.isfinite(vals).all(axis=1)] if vals.size else vals
    lo, hi = (np.min(vals) - 15, np.max(vals) + 15) if vals.size else (150, 400)
    ax.plot([lo, hi], [lo, hi], color="black", lw=1, label="y = x")
    ax.fill_between([lo, hi], [lo - 20, hi - 20], [lo + 20, hi + 20],
                    color="grey", alpha=0.12, label="±20 nm")

    for r in used_rules:
        s = ok[ok["rule_used"] == r]
        ax.scatter(s["exp_nm"], s["calc_nm"], s=55, color=color_map[r],
                   edgecolor="black", linewidth=0.5, label=r, zorder=3)
    fl = sub[sub["out_of_domain"]]
    if not fl.empty:
        ax.scatter(fl["exp_nm"], fl["calc_nm"], s=60, facecolors="none",
                   edgecolors="grey", marker="x", linewidth=1.2,
                   label="out_of_domain", zorder=3)

    err = ok["error_nm"].to_numpy(float)
    mae = np.nanmean(np.abs(err)) if err.size else np.nan
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("experimental λmax / nm")
    ax.set_ylabel("calculated λmax / nm")
    ax.set_title(f"{base_rule}: colored by chosen rule  (MAE {mae:.0f} nm)")
    ax.legend(fontsize=8, framealpha=0.9, title="rule_used")
    ax.grid(ls=":", alpha=0.5)
    fig.tight_layout()
    fig.savefig(outfile, dpi=160)
    plt.close(fig)

    # small stats table by chosen rule
    rows = []
    for r in used_rules:
        e = ok.loc[ok["rule_used"] == r, "error_nm"].to_numpy(float)
        rows.append((r, len(e), np.nanmean(np.abs(e)) if e.size else np.nan,
                     np.nanmean(e) if e.size else np.nan))
    return rows


def build_parser():
    p = argparse.ArgumentParser(
        description="Batch Woodward-Fieser / Fieser / Fieser-Kuhn lambda_max prediction.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-input", required=True, help="input table (.xlsx / .csv / .tsv)")
    p.add_argument("-out", default="predictions.csv", help="clean wide output CSV")
    p.add_argument("-tidy", default=None, help="optional long/tidy CSV (one row per molecule x rule)")
    p.add_argument("-to-csv", default=None, help="write the cleaned input table to this CSV (xlsx -> csv)")
    p.add_argument("-plot", default="error.png", help="error plot (.png); '' to skip")
    p.add_argument("-sheet", default=None, help="sheet name/index for xlsx input")

    p.add_argument("-smiles-col", default=None, help="SMILES column (auto-detected)")
    p.add_argument("-name-col", default=None, help="name/label column (auto-detected)")
    p.add_argument("-cas-col", default=None, help="CAS column (auto-detected)")
    p.add_argument("-exp-col", default=None, help="experimental lambda_max column (auto-detected)")

    p.add_argument("-solvent", default=None, help="solvent SMILES for the correction term (e.g. CCO)")
    p.add_argument("-rules", nargs="+", default=["auto"],
                   help="rule systems: auto woodward woodward_refine woodward_extended "
                        "woodward_coumarin fieser fieser_kuhn")
    p.add_argument("-hard-force", action="store_true",
                   help="always apply the requested rule (override auto-classification, no fallback)")
    p.add_argument("-color-by-rule-used", action="store_true",
                   help="simple single-panel scatter of one rule (default 'auto'), "
                        "each point colored by the rule the classifier chose")
    p.add_argument("-color-base-rule", default="auto",
                   help="which rule's rows to plot when -color-by-rule-used is set")

    p.add_argument("-domain-min", type=float, default=150.0, help="min plausible lambda_max (nm)")
    p.add_argument("-domain-max", type=float, default=600.0, help="max plausible lambda_max (nm)")
    p.add_argument("-keep-flagged", action="store_true",
                   help="include out-of-domain predictions in the error statistics")
    p.add_argument("-round", type=int, default=0,
                   help="decimals for calc/error columns (0 = integer nm, no '.0')")

    p.add_argument("-chromopredict-src", default=None, help="path to package 'src' if not installed")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)

    cp = load_chromopredict(args.chromopredict_src)
    if args.hard_force:
        ok = enable_hard_force(cp)
        print(f"hard-force {'ENABLED' if ok else 'unavailable (selector not found)'}")
    from chromopredict.strucfeatures import classify_type
    from chromopredict.chrombase import general_rules

    df = read_table(args.input, args.sheet)
    before = list(df.columns)
    df = drop_clutter(df)
    dropped = [c for c in before if c not in df.columns]
    if dropped:
        print(f"dropped clutter columns: {dropped}")
    if args.to_csv:
        df.to_csv(args.to_csv, index=False)
        print(f"cleaned CSV copy of input -> {args.to_csv}")

    smi_col = args.smiles_col or pick_col(df, "SMILES", ["SMILE", "canonical_smiles"])
    if smi_col is None:
        raise SystemExit(f"ERROR: no SMILES column. Columns: {list(df.columns)} (use -smiles-col)")
    name_col = args.name_col or pick_col(df, "Name", ["compound", "molecule", "title", "label"])
    cas_col = args.cas_col or pick_col(df, "CAS", ["cas_no", "casrn"])
    exp_col = args.exp_col or pick_col(df, "lambda_max_experimental",
        ["λmax experimental", "lambdamaxexp", "experimental", "lmax_exp", "abs_max_exp"])

    labels = (df[name_col].astype(str) if name_col
              else df[cas_col].astype(str) if cas_col
              else pd.Series([f"mol_{i}" for i in range(len(df))])).tolist()
    exp = pd.to_numeric(df[exp_col], errors="coerce").to_numpy(float) if exp_col else np.full(len(df), np.nan)

    rules = [(None if r.lower() == "auto" else r) for r in args.rules]
    rule_names = [r if r else "auto" for r in rules]
    print(f"{len(df)} molecules | SMILES='{smi_col}' | exp='{exp_col}' | "
          f"rules={rule_names} | hard_force={args.hard_force}")

    id_cols = [c for c in [cas_col, name_col, smi_col, exp_col] if c]
    wide = df[id_cols].copy()
    tidy_rows = []
    rnd = args.round

    for r, rname in zip(rules, rule_names):
        calc, used, flag = [], [], []
        for i, smi in enumerate(df[smi_col]):
            lam, rule_used, _ = predict_one(cp, classify_type, general_rules,
                                            smi, args.solvent, r, args.hard_force)
            oob = (not np.isfinite(lam)) or (lam < args.domain_min) or (lam > args.domain_max)
            calc.append(lam); used.append(rule_used); flag.append(bool(oob))
            err = lam - exp[i]
            row = {"label": labels[i]}
            if cas_col:
                row["CAS"] = df[cas_col].iloc[i]
            row.update({
                "SMILES": smi, "rule": rname, "rule_used": rule_used,
                "calc_nm": round(lam, rnd) if np.isfinite(lam) else np.nan,
                "exp_nm": exp[i],
                "error_nm": round(err, rnd) if np.isfinite(err) else np.nan,
                "abs_error_nm": round(abs(err), rnd) if np.isfinite(err) else np.nan,
                "out_of_domain": bool(oob),
            })
            tidy_rows.append(row)
        calc = np.array(calc, float)
        err = calc - exp
        wide[f"rule_used_{rname}"] = used
        wide[f"calc_{rname}_nm"] = np.round(calc, rnd)
        wide[f"err_{rname}_nm"] = np.round(err, rnd)
        wide[f"absErr_{rname}_nm"] = np.round(np.abs(err), rnd)
        wide[f"flag_{rname}"] = np.where(flag, "out_of_domain", "")

    wide.to_csv(args.out, index=False)
    print(f"wrote clean wide predictions -> {args.out}")

    tidy = pd.DataFrame(tidy_rows)

    if rnd == 0:
        for c in wide.columns:
            if c.startswith(("calc_", "err_", "absErr_")):
                wide[c] = pd.to_numeric(wide[c], errors="coerce").round().astype("Int64")
        for c in ["calc_nm", "exp_nm", "error_nm", "abs_error_nm"]:
            if c in tidy:
                tidy[c] = pd.to_numeric(tidy[c], errors="coerce").round().astype("Int64")
        wide.to_csv(args.out, index=False)

    if args.tidy:
        tidy.to_csv(args.tidy, index=False)
        print(f"wrote tidy predictions -> {args.tidy}")

    if args.keep_flagged:
        tidy_stat = tidy.copy()
        tidy_stat["out_of_domain"] = False

    if exp_col and np.isfinite(exp).any():
        src = tidy if not args.keep_flagged else tidy_stat

        if args.color_by_rule_used:
            if args.plot:
                rows = make_plot_auto_colored(src, args.plot, base_rule=args.color_base_rule)
                print(f"wrote plot -> {args.plot}  (single '{args.color_base_rule}' series, "
                      f"colored by chosen rule)")
                print("\nby chosen rule (rule_used):")
                print(f"   {'rule_used':16s} {'n':>4s} {'MAE':>7s} {'ME':>7s}")
                for r, n, mae, me in rows:
                    print(f"   {r:16s} {n:4d} {mae:7.1f} {me:+7.1f}")
            return 0

        if args.plot:
            stats = make_plot(src, rule_names, args.plot)
            print(f"wrote plot -> {args.plot}")
        else:
            stats = []
            for rname in rule_names:
                sub = src[src["rule"] == rname]
                valid = sub.loc[~sub["out_of_domain"], "error_nm"].to_numpy(float)
                stats.append((rname,
                              np.nanmean(np.abs(valid)) if valid.size else np.nan,
                              np.sqrt(np.nanmean(valid ** 2)) if valid.size else np.nan,
                              np.nanmean(valid) if valid.size else np.nan,
                              int(sub["out_of_domain"].sum())))
        print("\nerror statistics (nm" + (", ALL kept" if args.keep_flagged
              else ", out-of-domain excluded") + ")")
        print(f"   {'rule':16s} {'MAE':>7s} {'RMSE':>7s} {'ME':>7s}  flagged")
        for r, mae, rmse, me, nfl in stats:
            print(f"   {r:16s} {mae:7.1f} {rmse:7.1f} {me:+7.1f}  {nfl}")
    else:
        print("no numeric experimental column -> error columns/plot skipped (-exp-col to set)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
