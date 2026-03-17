# Notebooks

Jupyter notebooks show how to use protein-quest through its Python API and can be
run locally or in cloud notebook environments.

## Avalable notebooks

| Notebook | What you will do |
| --- | --- |
| [Search UniProt](uniprot.ipynb) | Find UniProt accessions and map them to PDB, AlphaFold, EMDB, and partner datasets. |
| [AlphaFold](alphafold.ipynb) | Download AlphaFold models, filter on confidence, and visualize structures with Mol*. |
| [PDBe](pdbe.ipynb) | Download PDBe structures, extract single chains, and visualize structures with Mol*. |

## Launch in cloud environments

Use the links below to open each notebook quickly.

| Notebook | Google Colab | notebooks.egi.eu | Binder | nbgitpuller |
| --- | --- | --- | --- | --- |
| Search UniProt | [Open](https://colab.research.google.com/github/haddocking/protein-quest/blob/main/docs/notebooks/uniprot.ipynb) | [Open hub](https://notebooks.egi.eu/hub/) | [Open](https://mybinder.org/v2/gh/haddocking/protein-quest/HEAD?urlpath=lab/tree/docs/notebooks/uniprot.ipynb) | [Generate link](https://nbgitpuller.readthedocs.io/en/latest/link.html?hub=https://notebooks.egi.eu&repo=https://github.com/haddocking/protein-quest&branch=main) |
| AlphaFold | [Open](https://colab.research.google.com/github/haddocking/protein-quest/blob/main/docs/notebooks/alphafold.ipynb) | [Open hub](https://notebooks.egi.eu/hub/) | [Open](https://mybinder.org/v2/gh/haddocking/protein-quest/HEAD?urlpath=lab/tree/docs/notebooks/alphafold.ipynb) | [Generate link](https://nbgitpuller.readthedocs.io/en/latest/link.html?hub=https://notebooks.egi.eu&repo=https://github.com/haddocking/protein-quest&branch=main) |
| PDBe | [Open](https://colab.research.google.com/github/haddocking/protein-quest/blob/main/docs/notebooks/pdbe.ipynb) | [Open hub](https://notebooks.egi.eu/hub/) | [Open](https://mybinder.org/v2/gh/haddocking/protein-quest/HEAD?urlpath=lab/tree/docs/notebooks/pdbe.ipynb) | [Generate link](https://nbgitpuller.readthedocs.io/en/latest/link.html?hub=https://notebooks.egi.eu&repo=https://github.com/haddocking/protein-quest&branch=main) |

<!-- TODO test links -->

## Notes

- Google Colab and Binder open notebooks directly from this repository.
- notebooks.egi.eu requires sign-in and VO enrollment before use.
- nbgitpuller links depend on a JupyterHub where nbgitpuller is installed.

This section explains how to run notebooks locally.

## Local setup

1. Install Jupyter.

```bash
python -m pip install jupyterlab
```

2. Install notebook dependencies.

For the released package:

```bash
python -m pip install protein-quest[nb]
```
(The `[nb]` extra installs `ipymolstar` for structure visualization in the AlphaFold and PDBe notebooks.)

3. Start Jupyter and open a notebook.

```bash
jupyter lab
```

Then open one of the notebooks (*.ipynb files).

## Runtime dependencies

The first code cell in each notebook installs required packages, including
`protein-quest` and `ipymolstar`.

- Google Colab: best for quick exploration.
- Binder: no local install needed, startup can take a few minutes.
- notebooks.egi.eu: requires EGI account and VO enrollment.
- nbgitpuller: requires a JupyterHub deployment with nbgitpuller installed.
