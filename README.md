# protein-quest

[![Documentation](https://img.shields.io/badge/Documentation-bonvinlab.org-blue?style=flat-square&logo=gitbook)](https://www.bonvinlab.org/protein-quest/)
[![CI](https://github.com/haddocking/protein-quest/actions/workflows/ci.yml/badge.svg)](https://github.com/haddocking/protein-quest/actions/workflows/ci.yml)
[![Research Software Directory Badge](https://img.shields.io/badge/rsd-00a3e3.svg)](https://www.research-software.nl/software/protein-quest)
[![PyPI](https://img.shields.io/pypi/v/protein-quest)](https://pypi.org/project/protein-quest/)
<!-- TODO replace with correct zenodo id -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15632658.svg)](https://doi.org/10.5281/zenodo.15632658)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/7a3f3f1fe64640d583a5e50fe7ba828e)](https://app.codacy.com/gh/haddocking/protein-quest/coverage?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_coverage)

Python package to search/retrieve/filter proteins and protein structures.

It uses

- [Uniprot Sparql endpoint](https://sparql.uniprot.org/) to search for proteins and their measured or predicted 3D structures.
- [Uniprot taxonomy](https://www.uniprot.org/taxonomy?query=*) to search for taxonomy.
- [QuickGO](https://www.ebi.ac.uk/QuickGO/api/index.html) to search for Gene Ontology terms.
- [gemmi](https://project-gemmi.github.io/) to work with macromolecular models.
- [dask-distributed](https://docs.dask.org/en/latest/) to compute in parallel.

An example workflow:

```mermaid
graph TB;
    taxonomy[/Search taxon/] -. taxon_ids .-> searchuniprot[/Search UniprotKB/]
    goterm[/Search GO term/] -. go_ids .-> searchuniprot[/Search UniprotKB/]
    searchuniprot --> |uniprot_accessions|searchpdbe[/Search PDBe/]
    searchuniprot --> |uniprot_accessions|searchaf[/Search Alphafold/]
    searchuniprot -. uniprot_accessions .-> searchemdb[/Search EMDB/]
    searchpdbe -->|pdb_ids|fetchpdbe[Retrieve PDBe]
    searchaf --> |uniprot_accessions|fetchad(Retrieve AlphaFold)
    searchemdb -. emdb_ids .->fetchemdb[Retrieve EMDB]
    fetchpdbe -->|mmcif_files_with_uniprot_acc| chainfilter{Filter on chain of uniprot}
    chainfilter --> |mmcif_files| residuefilter{Filter on chain length}
    fetchad -->|pdb_files| confidencefilter{Filter out low confidence}
    classDef dashedBorder stroke-dasharray: 5 5;
    goterm:::dashedBorder
    taxonomy:::dashedBorder
    searchemdb:::dashedBorder
    fetchemdb:::dashedBorder
```

(Dotted nodes and edges are side-quests.)

## Install

```shell
pip install protein-quest
```

Or to use the latest development version:
```
pip install git+https://github.com/haddocking/protein-quest.git
```

## Usage

The main entry point is the `protein-quest` command line tool which has multiple subcommands to perform actions.

To use programmaticly, see [API documentation](https://www.bonvinlab.org/protein-quest/autoapi/summary/).

### Search Uniprot accessions

```shell
protein-quest search uniprot \
    --taxon-id 9606 \
    --reviewed \
    --subcellular-location-uniprot nucleus \
    --subcellular-location-go GO:0005634 \
    --molecular-function-go GO:0003677 \
    --limit 100 \
    uniprot_accs.txt
```
([GO:0005634](https://www.ebi.ac.uk/QuickGO/term/GO:0005634) is "Nucleus" and [GO:0003677](https://www.ebi.ac.uk/QuickGO/term/GO:0003677) is  "DNA binding")

### Search for PDBe structures of uniprot accessions

```shell
protein-quest search pdbe uniprot_accs.txt pdbe.csv
```

`pdbe.csv` file is written containing the the PDB id and chain of each uniprot accession.

### Search for Alphafold structures of uniprot accessions

```shell
protein-quest search alphafold uniprot_accs.txt alphafold.csv
```

### Search for EMDB structures of uniprot accessions

```shell
protein-quest search emdb uniprot_accs.txt emdbs.csv
```

### To retrieve PDB structure files

```shell
protein-quest retrieve pdbe pdbe.csv downloads-pdbe/
```

### To retrieve AlphaFold structure files

```shell
protein-quest retrieve alphafold alphafold.csv downloads-af/
```

For each entry downloads the summary.json and cif file.

### To retrieve EMDB volume files

```shell
protein-quest retrieve emdb emdbs.csv downloads-emdb/
```

### To filter AlphaFold structures on confidence

Filter AlphaFoldDB structures based on confidence (pLDDT).
Keeps entries with requested number of residues which have a confidence score above the threshold.
Also writes pdb files with only those residues.

```shell
protein-quest filter confidence \
    --confidence-threshold 50 \
    --min-residues 100 \
    --max-residues 1000 \
    ./downloads-af ./filtered
```

### To filter PDBe files on chain of uniprot accession

Make PDBe files smaller by only keeping first chain of found uniprot entry and renaming to chain A.

```shell
protein-quest filter chain \
    pdbe.csv \
    ./downloads-pdbe ./filtered-chains
```

### To filter PDBe files on nr of residues

```shell
protein-quest filter residue  \
    --min-residues 100 \
    --max-residues 1000 \
    ./filtered-chains ./filtered
```

### Search Taxonomy

```shell
protein-quest search taxonomy "Homo sapiens" -
```

### Search Gene Ontology (GO)

You might not know what the identifier of a [Gene Ontology](https://geneontology.org/) term is at `protein-quest search uniprot`.
You can use following command to search for a Gene Ontology (GO) term.

```shell
protein-quest search go --limit 5 --aspect cellular_component apoptosome -
```

##  Model Context Protocol (MCP) server

Protein quest can also help LLMs like Claude Sonnet 4 by providing a [set of tools](https://modelcontextprotocol.io/docs/learn/server-concepts#tools-ai-actions) for protein structures.

![Protein Quest MCP workflow](https://github.com/haddocking/protein-quest/raw/main/docs/protein-quest-mcp.png)

To run mcp server you have to install the `mcp` extra with:

```shell
pip install protein-quest[mcp]
```

The server can be started with:

```shell
protein-quest mcp
```

The mcp server contains an prompt template to search/retrieve/filter candidate structures.

## Contributing

For development information and contribution guidelines, please see [CONTRIBUTING.md](CONTRIBUTING.md).
