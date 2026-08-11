# Notebooks

Jupyter notebooks show how to use protein-quest through its Python API and can
be run in cloud notebook environments or locally.

## Available notebooks

| Notebook                        | What you will do                                                                      |
| ------------------------------- | ------------------------------------------------------------------------------------- |
| [Search UniProt](uniprot.ipynb) | Find UniProt accessions and map them to PDB, AlphaFold, EMDB, and partner datasets.   |
| [AlphaFold](alphafold.ipynb)    | Download AlphaFold models, filter on confidence, and visualize structures with Mol\*. |
| [PDBe](pdbe.ipynb)              | Download PDBe structures, extract single chains, and visualize structures with Mol\*. |

## Launch in cloud environments

| Notebook       | Google Colab                                                                                                       | notebooks.egi.eu                      | Binder                                                                                                           | nbgitpuller                                                                                                                          |
| -------------- | ------------------------------------------------------------------------------------------------------------------ | ------------------------------------- | ---------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| Search UniProt | [Open](https://colab.research.google.com/github/haddocking/protein-quest/blob/main/docs/notebooks/uniprot.ipynb)   | [Open](https://notebooks.egi.eu/hub/) | [Open](https://mybinder.org/v2/gh/haddocking/protein-quest/HEAD?urlpath=lab/tree/docs/notebooks/uniprot.ipynb)   | [Generate link](https://nbgitpuller.readthedocs.io/en/latest/link.html?repo=https://github.com/haddocking/protein-quest&branch=main) |
| AlphaFold      | [Open](https://colab.research.google.com/github/haddocking/protein-quest/blob/main/docs/notebooks/alphafold.ipynb) | [Open](https://notebooks.egi.eu/hub/) | [Open](https://mybinder.org/v2/gh/haddocking/protein-quest/HEAD?urlpath=lab/tree/docs/notebooks/alphafold.ipynb) | [Generate link](https://nbgitpuller.readthedocs.io/en/latest/link.html?repo=https://github.com/haddocking/protein-quest&branch=main) |
| PDBe           | [Open](https://colab.research.google.com/github/haddocking/protein-quest/blob/main/docs/notebooks/pdbe.ipynb)      | [Open](https://notebooks.egi.eu/hub/) | [Open](https://mybinder.org/v2/gh/haddocking/protein-quest/HEAD?urlpath=lab/tree/docs/notebooks/pdbe.ipynb)      | [Generate link](https://nbgitpuller.readthedocs.io/en/latest/link.html?repo=https://github.com/haddocking/protein-quest&branch=main) |

- notebooks.egi.eu requires sign-in and VO enrollment before use.

## Run notebooks locally

1. Install Jupyter.

   ```bash
   python -m pip install jupyterlab
   ```

2. Install notebook dependencies.

   ```bash
   python -m pip install protein-quest[nb]
   ```

   (The `[nb]` extra installs `molviewspec` for structure visualization in the
   AlphaFold and PDBe notebooks.)

3. Start Jupyter and open a notebook.

   ```bash
   jupyter lab
   ```
