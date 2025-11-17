# protein-quest

[![Documentation](https://img.shields.io/badge/Documentation-bonvinlab.org-blue?style=flat-square&logo=gitbook)](https://www.bonvinlab.org/protein-quest/)
[![CI](https://github.com/haddocking/protein-quest/actions/workflows/ci.yml/badge.svg)](https://github.com/haddocking/protein-quest/actions/workflows/ci.yml)
[![Research Software Directory Badge](https://img.shields.io/badge/rsd-00a3e3.svg)](https://www.research-software.nl/software/protein-quest)
[![bio.tools](https://img.shields.io/badge/bio.tools-protein--quest-009fdf.svg)](https://bio.tools/protein-quest)
[![PyPI](https://img.shields.io/pypi/v/protein-quest)](https://pypi.org/project/protein-quest/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16941288.svg)](https://doi.org/10.5281/zenodo.16941288)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/7a3f3f1fe64640d583a5e50fe7ba828e)](https://app.codacy.com/gh/haddocking/protein-quest/coverage?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_coverage)
[![FAIR checklist badge](https://fairsoftwarechecklist.net/badge.svg)](https://fairsoftwarechecklist.net/v0.2?f=31&a=32113&i=32121&r=133)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green)](https://fair-software.eu)
[![Copy/paste detector](https://raw.githubusercontent.com/kucherenko/jscpd/refs/tags/v3.5.10/assets/jscpd-badge.svg?sanitize=true)](https://github.com/kucherenko/jscpd/)


Python package to search/retrieve/filter proteins and protein structures.

It uses

- [Uniprot Sparql endpoint](https://sparql.uniprot.org/) to search for proteins and their measured or predicted 3D structures.
- [Uniprot taxonomy](https://www.uniprot.org/taxonomy?query=*) to search for taxonomy.
- [QuickGO](https://www.ebi.ac.uk/QuickGO/api/index.html) to search for Gene Ontology terms.
- [gemmi](https://project-gemmi.github.io/) to work with macromolecular models.
- [dask-distributed](https://docs.dask.org/en/latest/) to compute in parallel.

The package is used by

- [protein-detective](https://github.com/haddocking/protein-detective)

An example workflow:

```mermaid
graph TB;
    taxonomy[/Search taxon/] -. taxon_ids .-> searchuniprot[/Search UniprotKB/]
    goterm[/Search GO term/] -. go_ids .-> searchuniprot[/Search UniprotKB/]
    searchuniprot --> |uniprot_accessions|searchpdbe[/Search PDBe/]
    searchuniprot --> |uniprot_accessions|searchaf[/Search Alphafold/]
    searchuniprot -. uniprot_accessions .-> searchemdb[/Search EMDB/]
    searchuniprot -. uniprot_accessions .-> searchuniprotdetails[/Search UniProt details/]
    searchintactionpartners[/Search interaction partners/] -.-x |uniprot_accessions|searchuniprot
    searchcomplexes[/Search complexes/]
    searchpdbe -->|pdb_ids|fetchpdbe[Retrieve PDBe]
    searchaf --> |uniprot_accessions|fetchad(Retrieve AlphaFold)
    searchemdb -. emdb_ids .->fetchemdb[Retrieve EMDB]
    fetchpdbe -->|mmcif_files| chainfilter{{Filter on chain of uniprot}}
    chainfilter --> |mmcif_files| residuefilter{{Filter on chain length}}
    fetchad -->|mmcif_files| confidencefilter{{Filter out low confidence}}
    confidencefilter --> |mmcif_files| ssfilter{{Filter on secondary structure}}
    residuefilter --> |mmcif_files| ssfilter
    ssfilter -. mmcif_files .-> convert2cif([Convert to cif])
    ssfilter -. mmcif_files .-> convert2uniprot_accessions([Convert to UniProt accessions])
    classDef dashedBorder stroke-dasharray: 5 5;
    goterm:::dashedBorder
    taxonomy:::dashedBorder
    searchemdb:::dashedBorder
    fetchemdb:::dashedBorder
    searchintactionpartners:::dashedBorder
    searchcomplexes:::dashedBorder
    searchuniprotdetails:::dashedBorder
    convert2cif:::dashedBorder
    convert2uniprot_accessions:::dashedBorder
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

To use programmaticly, see the [Jupyter notebooks](https://www.bonvinlab.org/protein-quest/notebooks) and [API documentation](https://www.bonvinlab.org/protein-quest/autoapi/protein_quest/).

While downloading or copying files it uses a global cache (located at `~/.cache/protein-quest`) and hardlinks to save disk space and improve speed.
This behavior can be customized with the `--no-cache`, `--cache-dir`, and `--copy-method` command line arguments.

### Search Uniprot accessions

```shell
protein-quest search uniprot \
    --taxon-id 9606 \
    --reviewed \
    --subcellular-location-uniprot "nucleus" \
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

For each entry downloads the cif file.

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

### To filter on secondary structure

To filter on structure being mostly alpha helices and have no beta sheets. See the following [notebook](https://www.bonvinlab.org/protein-detective/SSE_elements.html) to determine the ratio of secondary structure elements.

```shell
protein-quest filter secondary-structure \
    --ratio-min-helix-residues 0.5 \
    --ratio-max-sheet-residues 0.0 \
    --write-stats filtered-ss/stats.csv \
    ./filtered-chains ./filtered-ss
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

### Search for interaction partners

Use https://www.ebi.ac.uk/complexportal to find interaction partners of given UniProt accession.

```shell
protein-quest search interaction-partners Q05471 interaction-partners-of-Q05471.txt
```

The `interaction-partners-of-Q05471.txt` file contains uniprot accessions (one per line).

### Search for complexes

Given Uniprot accessions search for macromolecular complexes at https://www.ebi.ac.uk/complexportal
and return the complex entries and their members.

```shell
echo Q05471 | protein-quest search complexes - complexes.csv
```

The `complexes.csv` looks like

```csv
query_protein,complex_id,complex_url,complex_title,members
Q05471,CPX-2122,https://www.ebi.ac.uk/complexportal/complex/CPX-2122,Swr1 chromatin remodelling complex,P31376;P35817;P38326;P53201;P53930;P60010;P80428;Q03388;Q03433;Q03940;Q05471;Q06707;Q12464;Q12509
```

### Search for UniProt details

To get details (like protein name, sequence length, organism) for a list of UniProt accessions.

```shell
protein-quest search uniprot-details uniprot_accs.txt uniprot_details.csv
```

The `uniprot_details.csv` looks like:

```csv
uniprot_accession,uniprot_id,sequence_length,reviewed,protein_name,taxon_id,taxon_name
A0A087WUV0,ZN892_HUMAN,522,True,Zinc finger protein 892,9606,Homo sapiens
```

### Convert structure files to .cif format

Some tools (for example [powerfit](https://github.com/haddocking/powerfit)) only work with `.cif` files and not `*.cif.gz` or `*.bcif` files.

```shell
protein-quest convert structures --format cif --output-dir ./filtered-cif ./filtered-ss
```

### Convert structure files to UniProt accessions

After running some filters you might want to know which UniProt accessions are still present in the filtered structures.

```shell
protein-quest convert uniprot ./filtered-ss uniprot_accs.filtered.txt
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

## Shell autocompletion

The `protein-quest` command line tool supports shell autocompletion using [shtab](https://docs.iterative.ai/shtab).

Initialize for bash shell with:

```shell
mkdir -p ~/.local/share/bash-completion/completions
protein-quest --print-completion bash > ~/.local/share/bash-completion/completions/protein-quest
```

Initialize for zsh shell with:

```shell
mkdir -p ~/.local/share/zsh/site-functions
protein-quest --print-completion zsh > ~/.local/share/zsh/site-functions/_protein-quest
fpath=("$HOME/.local/share/zsh/site-functions" $fpath)
autoload -Uz compinit && compinit
```

## Contributing

For development information and contribution guidelines, please see [CONTRIBUTING.md](CONTRIBUTING.md).
