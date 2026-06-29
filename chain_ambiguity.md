<!-- TODO remove this document when ambiguity is resolved -->

# Chain Identifier Audit: label_asym_id vs auth_asym_id

## Scope

This document records where chain identifiers are used in protein-quest and what
identifier system each function expects or returns.

## Gemmi chain systems

8rw8 has:
label=A
auth=B

```python
import gemmi
s = gemmi.read_structure('~/.cache/protein-quest-tests/8rw8_updated.cif.gz')
s[0].find_chain('B')
<gemmi.Chain B with 116 res>
{c for m in s for c in m}
# {<gemmi.Chain B with 116 res>}
c = s[0][0]
c.name
# B
c.subchains()
# [<gemmi.ResidueSpan of 104: A [338(MET) 339(GLU) 340(TYR) ... 441(VAL)]>,
#  <gemmi.ResidueSpan of 1: B [501(A1H3P)]>,
#  <gemmi.ResidueSpan of 1: C [502(A1H3Q)]>,
#  <gemmi.ResidueSpan of 10: D [601(HOH) 602(HOH) 603(HOH) ... 610(HOH)]>]
s.entities[0].subchains
#['A']
atom_site = s.make_mmcif_block().get_mmcif_category("_atom_site.")
atom_site.get('label_asym_id',[])[0]
# 'A'
atom_site.get('auth_asym_id',[])[0]
# 'B'
```

## System legend

- AUTH_CHAIN: Gemmi chain name (model chain id; used by chain lookup and rename
  APIs).
- LABEL_ASYM: mmCIF label asym id (for example \_atom_site.label_asym_id,
  \_struct_ref_seq.pdbx_strand_id, 3D-Beacons chain_ids).
- MIXED: function combines values from multiple sources that may not be in the
  same system.
- EXTERNAL_PDBe: chain token from UniProt/PDBe uniprot_chains text field; not
  explicitly tied to label/auth fields in this code.

## Gemmi evidence used for classification

- Gemmi stores Residue.subchain as mmCIF \_atom_site.label_asym_id
  (include/gemmi/model.hpp).
- Gemmi docs in polyheur state subchain corresponds to label_asym_id and that
  auth_asym_id can be split into label_asym_id units
  (include/gemmi/polyheur.hpp).
- Entity.subchains are subchain ids (label_asym system via residues/entities
  linkage) (include/gemmi/metadata.hpp, include/gemmi/model.hpp).

## Function-by-function findings

### src/protein_quest/structure/chains.py

- find_chain_in_model(model, wanted_chain) -> Chain | None
  - Input: AUTH_CHAIN
  - Output: Gemmi Chain object in AUTH_CHAIN system (Chain.name)
- find_chain_in_structure(structure, wanted_chain) -> Chain | None
  - Input: AUTH_CHAIN or LABEL_ASYM (LABEL_ASYM is translated to AUTH_CHAIN
    before model lookup)
  - Output: AUTH_CHAIN
- resolve_chain_id_to_label(structure, chain_id, chain_system, ...) -> str
  - Input: chain_id interpreted using chain_system flag ("auth" or "label")
  - Output: LABEL_ASYM (canonical id used for downstream Gemmi extraction)
- nr_residues_in_chain(file, chain="A") -> int
  - Input: AUTH_CHAIN
- write_single_chain_structure_file(input_file, chain2keep, output_dir,
  out_chain="A", ...)
  - Input chain2keep: AUTH_CHAIN (Gemmi chain lookup id)
  - Output out_chain: AUTH_CHAIN for chain renaming, then also written into
    entity.subchains during normalization
- \_normalize_single_chain_entities(structure, out_chain)
  - Writes entity.subchains=[out_chain], so creates LABEL_ASYM-style subchain
    values from out_chain text
- retrieve_chain_extraction_provenance(structure)
  - Returns stored chain2keep/out_chain from the chain-filter pipeline (same
    semantic ids used above)

### src/protein_quest/structure/uniprot.py

- struct_ref_seqs_columns_to_records(struct_ref_seqs_columns) ->
  Generator[StructRefSeq]
  - Output StructRefSeq.chain_id from \_struct_ref_seq.pdbx_strand_id =>
    LABEL_ASYM
- selected_struct_ref_seqs_by_chain(structure, accessions) -> dict[str,
  StructRefSeq]
  - Keys: LABEL_ASYM
- structure_to_uniprot(structure) -> Pdb2UniprotMapping
  - Returns MIXED
  - Branch 1: entity.subchains + entity.sifts_unp_acc (subchain-based ids)
  - Branch 2: \_struct_ref_seq.pdbx_strand_id + \_struct_ref UNP rows
    (LABEL_ASYM)
  - In multi-entity structures these may not match one-to-one after
    normalization/roundtrip
- \_append_uniprot_to_structure(structure, chain_uniprot_pairs)
  - Input chain ids expected to match structure.entities[*].subchains keys
  - Writes \_struct_ref_seq.pdbx_strand_id with the same chain string
  - Effective system: subchain/LABEL_ASYM assumptions, but roundtrip
    normalization may remap
- add_uniprot_accessions2structure(structure, pdb2uniprot)
  - Input: `pdb2uniprot` + optional `chain_system` flag (default `auth`)
  - Behavior: resolves incoming chain ids to LABEL_ASYM when `chain_system=auth`
  - Output matching/injection is done against LABEL_ASYM-style chain ids
  - Since structure_to_uniprot is MIXED, this boundary is still partially
    ambiguous in multi-entity/roundtrip cases
- \_rename_chain_based_on_provenance(structure, pdb2uniprot)
  - Rewrites mapping using provenance chain ids from single-chain extraction
    flow

### src/protein_quest/structure/convert.py

- convert_to_cif_file(input_file, output_dir, ..., pdb2uniprot=None,
  chain_system="auth") -> Path
  - When `pdb2uniprot` is provided, forwards `chain_system` to
    add_uniprot_accessions2structure
  - Effective behavior: CSV/auth-style chain ids can be resolved to LABEL_ASYM
    before UniProt injection
- convert_to_cif_files(input_files, ..., pdb2uniprot=None, chain_system="auth")
  - Same as above; forwards chain_system to each per-file conversion

### src/protein_quest/cli/convert.py

- structures(input_dir, ..., uniprots=None, chain_system="auth")
  - Accepts chain_system flag for the `--uniprots` CSV chain column
  - Downstream path: structures -> convert_to_cif_files -> convert_to_cif_file
    -> add_uniprot_accessions2structure
  - This mirrors `filter chain` behavior: `auth` input can be translated to
    label ids at the processing boundary

### src/protein_quest/structure/metadata.py

- structure_metadata(structure, path=None) -> StructureMetadata
  - Pulls selected_struct_ref_seqs_by_chain chain_id (LABEL_ASYM)
  - Resolves chain via find_chain_in_structure(structure, chain_id), which uses
    AUTH_CHAIN lookup
  - Works when those ids coincide; fragile when label/auth differ

### src/protein_quest/filters/chain.py

- filter_file_on_chain((input_file, chain_id), output_dir, out_chain="A", ...)
  -> ChainFilterStatistics
  - Input: chain_id + chain_system (default "auth")
  - Behavior: resolves chain_id to LABEL_ASYM via resolve_chain_id_to_label,
    then passes resolved id to write_single_chain_structure_file
  - Effective extraction lookup remains AUTH_CHAIN inside
    write_single_chain_structure_file/find_chain_in_structure
- filter_files_on_chain(file2chains, ...)
  - Same as above; forwards chain_system to each filter_file_on_chain call
    (sequential and Dask modes)

### src/protein_quest/cli/filter.py

- chain(chains_csv, input_dir, output_dir, ...)
  - Accepts chain_system flag ("auth" default, "label" optional)
  - CSV chain column is interpreted according to chain_system and forwarded to
    filter_files_on_chain
  - Downstream path: chain -> filter_files_on_chain -> filter_file_on_chain ->
    resolve_chain_id_to_label -> write_single_chain_structure_file

### src/protein_quest/mcp_server.py

- extract_single_chain_from_structure(input_file, chain2keep, output_dir,
  out_chain="A") -> Path
  - chain2keep/out_chain semantics inherited from
    write_single_chain_structure_file => AUTH_CHAIN input, renamed output chain

### src/protein_quest/pdbe_3dbeacons/model.py

- Template.chain_id
  - Documented as label_asym_id
- AppUniprotSchemaEntity.chain_ids
  - Documented as label_asym identifiers

### src/protein_quest/pdbe_3dbeacons/search.py

- \_find_chain_for_uniprot(uniprot_accession, summary) -> str
  - Returns entity.chain_ids value => LABEL_ASYM
- flatten_structure_summaries(summaries) -> rows with chain
  - chain field is LABEL_ASYM (from \_find_chain_for_uniprot)

### src/protein_quest/pdbe/result.py

- PdbResult.chain (via \_first_chain_from_uniprot_chains)
  - Output system: EXTERNAL_PDBe (parsed from UniProt/PDBe uniprot_chains
    string)
  - In practice this behaves like author-chain naming for PDBe/UniProt chain
    strings (for example 1F66 returned `C/G=1-128`, not PDB-assigned label asym
    ids)

### src/protein_quest/uniprot.py

- search4pdb(...) and \_flatten_results_pdb(...)
  - Produces PdbResult.uniprot_chains and later PdbResult.chain
  - Chain ids here are EXTERNAL_PDBe from SPARQL `up:chain` values grouped into
    `uniprot_chains`
  - Current observed behavior: these are author-chain style IDs (AUTH-like), not
    PDB-assigned label_asym_id
  - Example: for 1F66 query in tests, search output contains `C/G` while
    label_asym assignment in the structure is `E/I`

## Confirmed issue pattern (1F66)

- The current mapping/injection flow compares and merges chain-UniProt pairs
  across sources with different semantics.
- In 1F66, label/auth drift and Gemmi roundtrip normalization can change where
  an injected pair appears.
- Practical effect: a pair injected as chain C in \_struct_ref_seq may be
  reported back under different chain ids after reconstruction.

## Summary of risk areas

- Highest risk: structure_to_uniprot (MIXED output system).
- Secondary risk: structure_metadata bridging LABEL_ASYM-derived chain ids into
  AUTH_CHAIN chain lookup.
- External input risk: chain ids from PDBe/UniProt search pipeline are
  EXTERNAL_PDBe and currently assumed compatible downstream.

## Implementation status checklist

- [x] `filter chain` command has `chain_system` flag (default `auth`, optional
      `label`).
- [x] `chain_system` is forwarded through `filter_files_on_chain` and
      `filter_file_on_chain`.
- [x] Chain ids are resolved at the filter boundary via
      `resolve_chain_id_to_label` before extraction.
- [x] `find_chain_in_structure` supports label input by translating label to
      auth before Gemmi lookup.
- [x] `convert structures` command has `chain_system` flag (default `auth`,
      optional `label`) for `--uniprots` CSV input.
- [x] `chain_system` is forwarded through
      `convert_to_cif_files`/`convert_to_cif_file` to
      `add_uniprot_accessions2structure`.
- [x] `structure_to_uniprot` still returns MIXED chain-id systems across
      branches.
- [x] `structure_metadata` still bridges LABEL_ASYM-derived ids into AUTH_CHAIN
      lookup and can be fragile when ids differ.
- [ ] External PDBe/UniProt chain ids are still assumed compatible downstream
      without an explicit conversion boundary.

## Recommendation

Adopt one canonical chain system at all function boundaries in structure
workflows, and convert at edges only. The most defensible choice for mmCIF
metadata joins is LABEL_ASYM. If AUTH_CHAIN is needed for direct chain
operations, perform explicit conversion with a helper map and keep both ids
visible in diagnostics.
