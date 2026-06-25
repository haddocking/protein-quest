TODO remove this document when ambiguity is resolved

# Chain Identifier Audit: label_asym_id vs auth_asym_id

## Scope
This document records where chain identifiers are used in protein-quest and what identifier namespace each function expects or returns.

## Namespace legend
- AUTH_CHAIN: Gemmi chain name (model chain id; used by chain lookup and rename APIs).
- LABEL_ASYM: mmCIF label asym id (for example _atom_site.label_asym_id, _struct_ref_seq.pdbx_strand_id, 3D-Beacons chain_ids).
- MIXED: function combines values from multiple sources that may not be in the same namespace.
- EXTERNAL_PDBe: chain token from UniProt/PDBe uniprot_chains text field; not explicitly tied to label/auth fields in this code.

## Gemmi evidence used for classification
- Gemmi stores Residue.subchain as mmCIF _atom_site.label_asym_id (include/gemmi/model.hpp).
- Gemmi docs in polyheur state subchain corresponds to label_asym_id and that auth_asym_id can be split into label_asym_id units (include/gemmi/polyheur.hpp).
- Entity.subchains are subchain ids (label_asym namespace via residues/entities linkage) (include/gemmi/metadata.hpp, include/gemmi/model.hpp).

## Function-by-function findings

### src/protein_quest/structure/chains.py
- find_chain_in_model(model, wanted_chain) -> Chain | None
  - Input: AUTH_CHAIN
  - Output: Gemmi Chain object in AUTH_CHAIN namespace (Chain.name)
- find_chain_in_structure(structure, wanted_chain) -> Chain | None
  - Input: AUTH_CHAIN
  - Output: AUTH_CHAIN
- nr_residues_in_chain(file, chain="A") -> int
  - Input: AUTH_CHAIN
- write_single_chain_structure_file(input_file, chain2keep, output_dir, out_chain="A", ...)
  - Input chain2keep: AUTH_CHAIN
  - Output out_chain: AUTH_CHAIN for chain renaming, then also written into entity.subchains during normalization
- _normalize_single_chain_entities(structure, out_chain)
  - Writes entity.subchains=[out_chain], so creates LABEL_ASYM-style subchain values from out_chain text
- retrieve_chain_extraction_provenance(structure)
  - Returns stored chain2keep/out_chain from the chain-filter pipeline (same semantic ids used above)

### src/protein_quest/structure/uniprot.py
- struct_ref_seqs_columns_to_records(struct_ref_seqs_columns) -> Generator[StructRefSeq]
  - Output StructRefSeq.chain_id from _struct_ref_seq.pdbx_strand_id => LABEL_ASYM
- selected_struct_ref_seqs_by_chain(structure, accessions) -> dict[str, StructRefSeq]
  - Keys: LABEL_ASYM
- structure_to_uniprot(structure) -> Pdb2UniprotMapping
  - Returns MIXED
  - Branch 1: entity.subchains + entity.sifts_unp_acc (subchain-based ids)
  - Branch 2: _struct_ref_seq.pdbx_strand_id + _struct_ref UNP rows (LABEL_ASYM)
  - In multi-entity structures these may not match one-to-one after normalization/roundtrip
- _append_uniprot_to_structure(structure, chain_uniprot_pairs)
  - Input chain ids expected to match structure.entities[*].subchains keys
  - Writes _struct_ref_seq.pdbx_strand_id with the same chain string
  - Effective namespace: subchain/LABEL_ASYM assumptions, but roundtrip normalization may remap
- add_uniprot_accessions2structure(structure, pdb2uniprot)
  - Input expected to match structure_to_uniprot output namespace
  - Since structure_to_uniprot is MIXED, this boundary is currently ambiguous
- _rename_chain_based_on_provenance(structure, pdb2uniprot)
  - Rewrites mapping using provenance chain ids from single-chain extraction flow

### src/protein_quest/structure/metadata.py
- structure_metadata(structure, path=None) -> StructureMetadata
  - Pulls selected_struct_ref_seqs_by_chain chain_id (LABEL_ASYM)
  - Resolves chain via find_chain_in_structure(structure, chain_id), which uses AUTH_CHAIN lookup
  - Works when those ids coincide; fragile when label/auth differ

### src/protein_quest/filters/chain.py
- filter_file_on_chain((input_file, chain_id), output_dir, out_chain="A", ...) -> ChainFilterStatistics
  - Input chain_id passed directly to write_single_chain_structure_file => AUTH_CHAIN
- filter_files_on_chain(file2chains, ...)
  - Same as above

### src/protein_quest/cli/filter.py
- chain(chains_csv, input_dir, output_dir, ...)
  - CSV chain column is forwarded to filter_files_on_chain/write_single_chain_structure_file => AUTH_CHAIN expected

### src/protein_quest/mcp_server.py
- extract_single_chain_from_structure(input_file, chain2keep, output_dir, out_chain="A") -> Path
  - chain2keep/out_chain semantics inherited from write_single_chain_structure_file => AUTH_CHAIN input, renamed output chain

### src/protein_quest/pdbe_3dbeacons/model.py
- Template.chain_id
  - Documented as label_asym_id
- AppUniprotSchemaEntity.chain_ids
  - Documented as label_asym identifiers

### src/protein_quest/pdbe_3dbeacons/search.py
- _find_chain_for_uniprot(uniprot_accession, summary) -> str
  - Returns entity.chain_ids value => LABEL_ASYM
- flatten_structure_summaries(summaries) -> rows with chain
  - chain field is LABEL_ASYM (from _find_chain_for_uniprot)

### src/protein_quest/pdbe/result.py
- PdbResult.chain (via _first_chain_from_uniprot_chains)
  - Output namespace: EXTERNAL_PDBe (parsed from UniProt/PDBe uniprot_chains string)
  - In practice this behaves like author-chain naming for PDBe/UniProt chain strings
    (for example 1F66 returned `C/G=1-128`, not PDB-assigned label asym ids)

### src/protein_quest/uniprot.py
- search4pdb(...) and _flatten_results_pdb(...)
  - Produces PdbResult.uniprot_chains and later PdbResult.chain
  - Chain ids here are EXTERNAL_PDBe from SPARQL `up:chain` values grouped into `uniprot_chains`
  - Current observed behavior: these are author-chain style IDs (AUTH-like), not PDB-assigned label_asym_id
  - Example: for 1F66 query in tests, search output contains `C/G` while label_asym assignment in the structure is `E/I`

## Confirmed issue pattern (1F66)
- The current mapping/injection flow compares and merges chain-UniProt pairs across sources with different semantics.
- In 1F66, label/auth drift and Gemmi roundtrip normalization can change where an injected pair appears.
- Practical effect: a pair injected as chain C in _struct_ref_seq may be reported back under different chain ids after reconstruction.

## Summary of risk areas
- Highest risk: structure_to_uniprot (MIXED output namespace).
- Secondary risk: structure_metadata bridging LABEL_ASYM-derived chain ids into AUTH_CHAIN chain lookup.
- External input risk: chain ids from PDBe/UniProt search pipeline are EXTERNAL_PDBe and currently assumed compatible downstream.

## Recommendation
Adopt one canonical chain namespace at all function boundaries in structure workflows, and convert at edges only. The most defensible choice for mmCIF metadata joins is LABEL_ASYM. If AUTH_CHAIN is needed for direct chain operations, perform explicit conversion with a helper map and keep both ids visible in diagnostics.
