# Coherent Chain–UniProt Pair Handling

This document catalogs every distinct way chain–UniProt-accession pairs appear
in the protein-quest codebase, highlights divergences, and recommends
unification steps.

## 1. Current representations (six variants)

| Representation            | Example / Type                                                                                                       | Where used                                                                                            | Notes                                  |
| ------------------------- | -------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------- | -------------------------------------- |
| **A. Raw string**         | `"A/B=1-42,A/B=50-99"`                                                                                               | `PdbResult.uniprot_chains` (pdbe/result.py:32), CSV I/O                                               | Parsed on demand                       |
| **B. Parsed ranges**      | `tuple[UniprotChainRange, ...]` (each has `chain_ids`, `start`, `end`)                                               | `uniprot_chains.py:10-52`, injection pipeline                                                         | Rich structure, carries ranges         |
| **C. Simple tuple pair**  | `tuple[str, str]` = `(chain_id, uniprot_accession)`                                                                  | `set[tuple[str, str]]` in extraction functions (structure/uniprot.py:109,126,149,308), most pervasive | Minimal, no range info                 |
| **D. Structured mapping** | `UniprotChainMapping` (accession + ranges) → `Pdb2RangeMappings` (pdb→set)                                           | `uniprot_chains.py:56-71`, injection CSV reading                                                      | Carries ranges, used for injection     |
| **E. Separated fields**   | `StructRefSeq.uniprot_accession` + `.chain_id`; `StructureMetadata.uniprot_accession` + `.auth_chain`/`.label_chain` | types.py, metadata.py, resolution.py                                                                  | Pair is implicit, no tuple             |
| **F. Protocol-based**     | `ClusterableStructure` has `uniprot_start`/`end` but no accession field                                              | clustering.py:52-66                                                                                   | Accession implicit in grouping context |

## 2. Key divergences

### A. `Pdb2RawPairs` vs `Pdb2RangeMappings` – previously `Pdb2UniprotMapping` and `Pdb2UniprotChainsMapping`

- `Pdb2RawPairs = dict[str, set[tuple[str, str]]]` (types.py:28) – raw pairs
  from extraction
- `Pdb2RangeMappings = dict[str, set[UniprotChainMapping]]`
  (uniprot_chains.py:70) – structured ranges for injection
- Renamed from the confusingly similar `Pdb2UniprotMapping` /
  `Pdb2UniprotChainsMapping`

### B. `PdbResult` splits the pair across dict key and object

- `PdbResults = dict[str, set[PdbResult]]` – the UniProt accession is the **dict
  key**, not a field on `PdbResult`
- `PdbResult` has `uniprot_chains: str` (raw string) but **no
  `uniprot_accession` field**
- Every consumer must reconstruct: `(accession, pdb_result.chain)` —
  inconsistent with every other dataclass that carries the accession as a field

### C. `structure2uniprot_accessions()` discards chain information

- Returns `set[str]` (just accessions), dropping the chain component
  (structure/uniprot.py:183-195)
- Callers like `structure_metadata()` then re-extract chain info via
  `selected_struct_ref_seqs_by_chain()` or
  `selected_struct_ref_seqs_from_sifts_by_chain()`
- This is an unnecessary round-trip: extract pairs → discard chain → re-extract
  chain

### D. `_mapping_pairs()` round-trips through string format/parse

- Calls `all_chain_ids()` → `format_uniprot_chains()` → `parse_chain_ids()` just
  to extract chain IDs (structure/uniprot.py:308-313)
- Chain IDs already exist in `mapping.chain_ranges[].chain_ids` — can iterate
  directly

### E. `_group_missing_ranges()` and `_rename_chain_based_on_provenance()` flagged as overly complex

- Both have TODO comments noting they may be more complicated than needed (lines
  201, 332)

### F. 3D Beacons `_find_chain_for_uniprot()` returns label_asym IDs but claims auth system

- Source model.py:829 says `chain_ids` are "label_asym identifiers"
- Docstring at search.py:241 claims "auth" system — mismatch

### G. `ResolutionFilterStatistics` has no chain field

- Carries `uniprot_accession: str | None` (resolution.py:51) but **no chain ID**
- Downstream clustering consumers must re-read the structure to get chain info

### H. CSV column schemas diverge across CLI commands

- search.py:70:
  `uniprot_accession, pdb_id, method, resolution, uniprot_chains, chain, chain_length`
- convert.py:86 input: `pdb_id, uniprot_accession, uniprot_chains` (3 columns)
- resolution.py:600: `input_file, id, uniprot_accession, resolution, ...`
- clustering_io.py:83:
  `uniprot_accession, cluster_id, rank_in_cluster, structure_id, ...`
- filter.py:531: `input_file, pdb_id, uniprot_accession, ...`
- No shared CSV schema for pair representation

## 3. Recommended unification steps

1. **~~Merge or clearly rename `Pdb2UniprotMapping` and
   `Pdb2UniprotChainsMapping`~~ DONE**
   - Renamed to `Pdb2RawPairs` and `Pdb2RangeMappings` respectively
   - See commit for full details

2. **Add `uniprot_accession` as a field on `PdbResult`**
   - Remove the dict-key-as-accession pattern; make `PdbResult` self-contained
   - `PdbResult` would become
     `(accession, chain_id, method, resolution, uniprot_chains)` — every other
     dataclass already carries accession as a field

3. **Replace `structure2uniprot_accessions()` with a function that returns
   pairs**
   - Return `set[tuple[str, str]]` instead of `set[str]` so callers don't
     discard and re-extract chain info
   - Update `structure_metadata()` and other callers accordingly

4. **Eliminate the format/parse round-trip in `_mapping_pairs()`**
   - Directly iterate `mapping.chain_ranges` and collect `.chain_ids`
   - Remove dependency on `all_chain_ids()` → `format_uniprot_chains()` →
     `parse_chain_ids()`

5. **Simplify `_group_missing_ranges()` and
   `_rename_chain_based_on_provenance()`**
   - Address the TODO comments; both can likely be simplified significantly

6. **Fix chain ID system mismatch in 3D Beacons search**
   - Either convert label_asym IDs to auth system in
     `_find_chain_for_uniprot()`, or correct the docstring
   - Consider storing the auth chain alongside the label_asym chain in
     `FlattenedUniprotSummary`

7. **Add chain field to `ResolutionFilterStatistics`**
   - Store `chain_id: str | None` alongside `uniprot_accession` so downstream
     consumers don't need to re-read structure files

8. **Unify CSV column schemas across CLI commands**
   - Define a shared schema or base writer that all CLI commands use for pair
     rows
   - Standardize on columns: `uniprot_accession, chain_id, pdb_id, ...`

9. **Consider adding a `Pair` named tuple / dataclass to the public API**
   - `ChainUniprotPair = namedtuple('ChainUniprotPair', ['chain_id', 'uniprot_accession'])`
   - Replace all ad-hoc `tuple[str, str]` usage with the named type for clarity
     and type safety
