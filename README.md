
![](chromopredict_header.png)

> From Handbooks to High-Throughput: Automated application of Woodward-Fieser, Fieser and Fieser-Kuhn rules to predict low-energy ππ<sup>∗</sup> absorption maxima of α,β-unsaturated carbonyl compounds, dienes and systems with more than four conjugated carbon-carbon double bonds.
> Application of extended Woodward-Fieser rules to predict the absorption maxima of 3,4,6-substituted coumarins.

## Woodward-Type Examples

<table>
  <tr>
    <td align="left">
      <img src="https://github.com/CompPhotoChem/ChromoPredict/blob/main/examples/CC(%3DO)C1%3DC(C)CCCC1.png" width="300px"><br>
      <b>CC(=O)C1=C(C)CCCC1</b><br>
      <br>
      <b>Experiment:</b> 247 nm<br>
      <b>Woodward:</b> 237 nm<br>
      <b>Woodward refined:</b> 246 nm<br>
    </td>
    <td align="left">
      <img src="https://github.com/CompPhotoChem/Woodward_Fieser_Rules/blob/main/examples/CC1%3DC2CCCCC2CCC1%3DO.png" width="300px"><br>
      <b>CC1=C2CCCCC2CCC1=O</b><br>
      <br>
      <b>Experiment:</b> 239 nm<br>
      <b>Woodward:</b> 242 nm<br>
      <b>Woodward refined:</b> 251 nm<br>
    </td>
    <td align="left">
      <img src="https://github.com/CompPhotoChem/Woodward_Fieser_Rules/blob/main/examples/CC(%3DO)C%3DCC1%3DC(C)CCCC1(C)C.png" width="300px"><br>
      <b>CC(=O)C=CC1=C(C)CCCC1(C)C</b><br>
      <br>
      <b>Experiment:</b> 283 nm<br>
      <b>Woodward:</b> 281 nm<br>
      <b>Woodward refined:</b> 281 nm<br>
    </td>
  </tr>
</table>

## Installation

This project is based on RDKit and the tutorials are provided as Jupyter notebooks.
All dependencies are listed in the `requirements.txt` file from which one can install the necessary packages. 
If you're using a requirements.txt file, navigate to its directory and run:

```
pip install requirements.txt
```

Alternatively, you can install the packages directly:

```
pip install numpy pandas rdkit seaborn jupyterlab notebook py3Dmol
```

## Usage

The basic usage of the tool is outlined in the ![tutorial notebook](https://github.com/CompPhotoChem/ChromoPredict/blob/main/examples/01_tutorial_ChromoPredict.ipynb). 
For example, to predict the absorption maximum of the enon shown above, one can run:

```python

import chromopredict as cp

# original woodward-fieser rules
abs_max, description, image = cp.predict(
  smiles='CC(=O)C1=C(C)CCCC1',
  solvent=None,
  verbose=True, # return increments of all structural features
  chromlib='woodward')

#refined rules by us
abs_max, description, image = cp.predict(
  smiles='CC(=O)C1=C(C)CCCC1',
  solvent=None,
  verbose=True,
  chromlib='woodward_refine')

```
