# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "marimo>=0.23.16",
#     "protein-quest>=1.5.5",
# ]
# ///

import marimo

__generated_with = "0.23.16"
app = marimo.App()


@app.cell(hide_code=True)
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Search on uniprot

    You can search on uniprot to find uniprot accessions and structure identifiers.
    """)
    return


@app.cell
def _():
    # Setup some logging
    import logging
    from pprint import pprint

    logging.basicConfig(level=logging.WARNING)
    # Set to WARNING to see only warnings
    # Set to INFO to see sparql queries
    # Set to DEBUG to see raw results
    return (pprint,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Search for uniprot accessions based on a query
    """)
    return


@app.cell
def _():
    from protein_quest.uniprot import Query, search4uniprot

    return Query, search4uniprot


@app.cell
def _(Query):
    query = Query(
        taxon_id=9606,
        reviewed=True,
        subcellular_location_uniprot="nucleus",
        subcellular_location_go={"GO:0005634"},  # Cellular component - Nucleus
        molecular_function_go={"GO:0003677"},  # Molecular function - DNA binding
    )
    return (query,)


@app.cell
def _(query, search4uniprot):
    uniprot_accessions = search4uniprot(query, limit=200)
    return (uniprot_accessions,)


@app.cell
def _(pprint, uniprot_accessions):
    print(f"Number of Uniprot accessions: {len(uniprot_accessions)}")
    print("First 5:")
    pprint(list(uniprot_accessions)[:5])
    print("Last 5:")
    pprint(list(uniprot_accessions)[-5:])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Find Protein Data Bank (PDB) entries for uniprot entries
    """)
    return


@app.cell
def _():
    from protein_quest.uniprot import search4pdb

    return (search4pdb,)


@app.cell
def _(search4pdb, uniprot_accessions):
    pdb_results = search4pdb(uniprot_accessions, limit=200)
    return (pdb_results,)


@app.cell
def _(pdb_results, pprint):
    pprint(f"Number of PDB results: {len(pdb_results)}")
    pprint("First 5 PDBs of first Uniprot entry:")
    _first_uniprot = next(iter(pdb_results.items()))
    pprint(_first_uniprot[0])
    pprint(list(_first_uniprot[1])[:5])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Find AlphaFold models for uniprot entries
    """)
    return


@app.cell
def _():
    from protein_quest.uniprot import search4af

    return (search4af,)


@app.cell
def _(search4af, uniprot_accessions):
    afresults = search4af(uniprot_accessions, limit=200)
    return (afresults,)


@app.cell
def _(afresults, pprint):
    pprint(f"Number of AlphaFold results: {len(afresults)}")
    first_af = next(iter(afresults.items()))
    pprint(first_af)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Find Electron Microscopy Data Bank (EMDB) entries for uniprot entries
    """)
    return


@app.cell
def _():
    from protein_quest.uniprot import search4emdb

    return (search4emdb,)


@app.cell
def _(search4emdb, uniprot_accessions):
    uniprot_accessions_1 = search4emdb(uniprot_accessions, limit=200)
    return (uniprot_accessions_1,)


@app.cell
def _(pprint, uniprot_accessions_1):
    pprint(f"Number of Uniprot accessions with EMDB entries: {len(uniprot_accessions_1)}")
    _first_uniprot = next(iter(uniprot_accessions_1.items()))
    pprint(_first_uniprot)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Find interaction partners for uniprot entries
    """)
    return


@app.cell
def _():
    from protein_quest.uniprot import search4interaction_partners, search4macromolecular_complexes

    return search4interaction_partners, search4macromolecular_complexes


@app.cell
def _():
    # Helicase SWR1 in yeast
    uniprot_accession = "Q05471"
    return (uniprot_accession,)


@app.cell
def _(search4interaction_partners, uniprot_accession):
    partners = search4interaction_partners(uniprot_accession, limit=100)
    partners
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    To get more information about the complex you can search for the complexes themselves with:
    """)
    return


@app.cell
def _(search4macromolecular_complexes, uniprot_accession):
    complexes = search4macromolecular_complexes([uniprot_accession])
    complexes
    return


if __name__ == "__main__":
    app.run()
