import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import MACCSkeys, rdFingerprintGenerator

from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error

sys.path.append('../../Woodward_Fieser_Rules-main/src/')
import chromopredict as cp

############################################################
# --- Random Forest Helpers
###########################################################

def get_data_split(X, y, test_size=0.2, random_state=42):

    return train_test_split(X, y, test_size=test_size, random_state=random_state)

def train_rf_model(X_train, X_test, y_train, y_test, random_state=42):
    
    # Train model
    model = RandomForestRegressor(n_estimators=100, random_state=random_state)
    model.fit(X_train, y_train)

    # Predict on test data
    y_pred = model.predict(X_test)
    mae = mean_absolute_error(y_test, y_pred)

    return model, mae

############################################################
# --- Topological Torsions FP Helpers
###########################################################

def get_tt_fp(smiles):

    mol = Chem.MolFromSmiles(smiles)
    ttgen = rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=2048)
    fp = ttgen.GetFingerprint(mol)

    return np.array(fp)

def get_tt_fp_data(df):

    fp_array = np.stack(df['smiles'].map(lambda smi: get_tt_fp(smi)))
    X = pd.DataFrame(fp_array, columns=[f'fp_{i}' for i in range(fp_array.shape[1])])
    y = df['nm_b3lyp']

    return X, y

def get_ttfp_pred(smi, X_test, model):

    unseen_fp = get_tt_fp(smi).reshape(1, -1)
    unseen_df = pd.DataFrame(unseen_fp, columns=X_test.columns)
    pred = model.predict(unseen_df)[0]

    return int(pred)

############################################################
# --- Rooted Topological Torsions FP Helpers
###########################################################

def flatten_and_deduplicate(lst):
    if not lst:
        return []
    # Flatten the list (handling both flat and nested)
    flat = [item for sublist in lst for item in (sublist if isinstance(sublist, list) else [sublist])]
    # Remove duplicates while preserving order
    seen = set()
    return [x for x in flat if not (x in seen or seen.add(x))]

def pattern_matches(smiles,
                    smarts=['[#6]C(=O)[#8,#6,#1]',
                            '[#6,#8,#1]C([#6,#8,#1])=C([#6,#8,#1])C(=O)[#8,#6,#1]',
                            '[#6,#8,#1]C([#6,#8,#1])=C([#6,#8,#1])C([#6,#8,#1])=C([#6,#8,#1])C(=O)[#8,#6,#1]']
                   ):

    mol = Chem.MolFromSmiles(smiles)

    matches = []
    for pattern in smarts:
        mol_pattern = Chem.MolFromSmarts(pattern)
        matches.append([x[0] for x in mol.GetSubstructMatches(mol_pattern)])

    matches_flat = flatten_and_deduplicate(matches)

    ttgen = rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=2048)

    ao = rdFingerprintGenerator.AdditionalOutput()
    ao.AllocateBitPaths()

    fp = ttgen.GetFingerprint(mol,fromAtoms=matches_flat,additionalOutput=ao)

    return np.array(fp)

def get_rooted_fp_data(df):

    fp_array = np.stack(df['smiles'].map(lambda smi: pattern_matches(smi)))
    X = pd.DataFrame(fp_array, columns=[f'fp_{i}' for i in range(fp_array.shape[1])])
    y = df['nm_b3lyp']

    return X, y

def get_rfp_pred(smi, X_test, model):

    unseen_fp = pattern_matches(smi).reshape(1, -1)
    unseen_df = pd.DataFrame(unseen_fp, columns=X_test.columns)
    pred = model.predict(unseen_df)[0]

    return int(pred)

############################################################
# --- MACCS Key FP Helpers
###########################################################

def get_maccs_fp(smiles):

    mol = Chem.MolFromSmiles(smiles)
    fp = MACCSkeys.GenMACCSKeys(mol)

    return np.array(fp)

def get_maccs_fp_data(df):

    fp_array = np.stack(df['smiles'].map(lambda smi: get_maccs_fp(smi)))

    X = pd.DataFrame(fp_array, columns=[f'fp_{i}' for i in range(fp_array.shape[1])])
    y = df['nm_b3lyp']

    return X, y

def get_maccs_pred(smi, X_test, model):

    unseen_fp = get_maccs_fp(smi).reshape(1, -1)
    unseen_df = pd.DataFrame(unseen_fp, columns=X_test.columns)
    pred = model.predict(unseen_df)[0]

    return int(pred)

############################################################
# --- Feature Morgan FP Helpers
###########################################################

def featurize(smiles, radius, fp_size):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(fp_size)

    fmgen = rdFingerprintGenerator.GetMorganGenerator(
        radius=radius,
        fpSize=fp_size,
        atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen()
    )

    fp = fmgen.GetFingerprint(mol)

    return np.array(fp)

def get_morgan_fp_data(df, radius=2, fp_size=1024):

    fp_array = np.stack(df['smiles'].map(lambda smi: featurize(smi, radius, fp_size)))
    X = pd.DataFrame(fp_array, columns=[f'fp_{i}' for i in range(fp_array.shape[1])])
    y = df['nm_b3lyp']

    return X, y

def get_fmfp_pred(case, X_test, model, params):

    unseen_fp = featurize(case, radius=params[0], fp_size=params[1]).reshape(1, -1)
    unseen_df = pd.DataFrame(unseen_fp, columns=X_test.columns)
    pred = model.predict(unseen_df)[0]

    return int(pred)

###########################################################
# --- helper function to get woodward features
###########################################################

def process_smiles_df(df, y_str='nm_b3lyp'):
    """
    Given a dataframe with a 'smiles' column, compute Woodward-Fieser descriptors,
    extract numeric contributions, and return a processed dataframe ready for modeling.
    """

    # Helper functions
    def get_first_value(lst):
        if isinstance(lst, list) and lst and isinstance(lst[0], dict):
            return lst[0].get("value", 0)
        return 0

    def get_descr(smiles):
        _, descr, _ = cp.predict(smiles, solvent=None, verbose=True, draw=False, chromlib='woodward')
        return descr

    # Compute Woodward-Fieser descriptors for each SMILES
    df = df.copy()  # avoid modifying original df
    df['descr'] = df['smiles'].apply(get_descr)

    # Expand 'descr' column into separate columns
    df_descr = pd.concat([df.drop(columns='descr'), df['descr'].apply(pd.Series)], axis=1)

    # Add numeric alpha/beta contributions and total substituents
    df_descr = df_descr.assign(
        alpha_num = df_descr['alpha_all'].apply(get_first_value),
        beta_num  = df_descr['beta_all'].apply(get_first_value),
        Nsub      = df_descr.apply(lambda row: len(row['alpha_all']) + len(row['beta_all']), axis=1)
    )

    # Select columns for modeling
    df_processed = df_descr[['smiles', y_str, 'nm_wf', 'nm_wfr', 'base', 'exo', 'alpha_num', 'beta_num', 'Nsub']]

    return df_processed

############################################################
# --- Plot Helper Function
###########################################################

def plot_multiple_scatter(df, x_col, y_cols, labels=None, markers=None, colors=None, savepath=None):
    """
    Overlay multiple scatterplots on the same axes.

    Parameters
    ----------
    df : pandas.DataFrame
        Data containing the x and y values.
    x_col : str
        Column name for the x-axis.
    y_cols : list of str
        List of column names to plot on the y-axis.
    labels : list of str, optional
        Legend labels for each y column (default: use column names).
    markers : list of str, optional
        Matplotlib marker styles (default: cycle through ['o','X','s','^','D']).
    colors : list of str, optional
        Colors for each y column (default: Matplotlib color cycle).
    savepath : str, optional
        If given, the plot will be saved to this path (e.g., 'figure.pdf' or 'figure.png').
    """

    plt.rcParams.update({'font.size': 16})
    plt.figure(figsize=(6, 6))

    # Default styles
    if labels is None:
        labels = y_cols
    if markers is None:
        markers = ['o', 'X', 's', '^', 'D'] * (len(y_cols) // 5 + 1)
    if colors is None:
        colors = plt.cm.tab10.colors * (len(y_cols) // 10 + 1)

    for i, y_col in enumerate(y_cols):
        plt.scatter(
            df[x_col],
            df[y_col],
            marker=markers[i],
            facecolors='none' if markers[i] == 'o' else colors[i],
            edgecolors=colors[i],
            label=labels[i],
            s=200,
            linewidth=2,
        )

    plt.legend()
    
    # Add identity line
    min_val = min(df[x_col].min(), *(df[y].min() for y in y_cols))
    max_val = max(df[x_col].max(), *(df[y].max() for y in y_cols))
    plt.plot([min_val, max_val], [min_val, max_val], 'k--')

    plt.xlabel('$\lambda_{max}^{ref}$ / nm')
    plt.ylabel('$\lambda_{max}^{pred}$ / nm')
    plt.tight_layout()

    if savepath:
        plt.savefig(savepath, dpi=300, bbox_inches="tight")
    plt.show()


