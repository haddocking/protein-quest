# Order Cluster Members Implementation Plan

## Goal
Add `order_cluster_members` in `src/protein_quest/pdbe/clustering.py` to rank members in this order:
1. Sequence identity descending (highest first)
2. Resolution ascending (best first)
3. Chain length descending (longest first)
4. PDB ID ascending (deterministic tie-break)

## Confirmed Decisions (Grill-Me)
1. Add a computed cached property for sequence identity on `PdbResult`.
2. Missing/invalid sequence identity is deprioritized behind valid values.
3. Scope for this task: function + unit tests only, no integration into broader filtering flow.
4. TDD process must add one test case at a time (no writing many `pytest.param` cases in one go).

## Files in Scope
1. `src/protein_quest/pdbe/clustering.py`
2. `src/protein_quest/pdbe/result.py`
3. `tests/pdbe/test_clustering.py`
4. `tests/pdbe/test_result.py`

## TDD Sequence (One Behavior Slice at a Time)

### Slice 1: PdbResult.sequence_identity
1. RED: Add one failing test in `tests/pdbe/test_result.py` that reuses `COMMON_CHAIN_CASES`.
   - Compute expected sequence identity as `chain_length / (uniprot_end - uniprot_start + 1)` from one case in `COMMON_CHAIN_CASES`.
2. GREEN: Extend `PdbResult` with cached `sequence_identity` property in `src/protein_quest/pdbe/result.py`.
3. REFACTOR: Keep parsing/derivation local to `PdbResult` and preserve invalid-chain behavior.

### Slice 2: Primary order_cluster_members criterion
1. RED: Add one failing test in `tests/pdbe/test_clustering.py`:
   - higher sequence identity ranks first.
2. GREEN: Add minimal `order_cluster_members` implementation in `src/protein_quest/pdbe/clustering.py` to satisfy this only.
3. REFACTOR: Keep sort-key simple and deterministic.

### Slice 3: Resolution tie-break
1. RED: Add one failing test:
   - when sequence identity ties, lower resolution ranks first.
2. GREEN: Extend sort key with resolution ordering.
3. REFACTOR: Keep missing/invalid resolution behavior aligned with current project semantics.

### Slice 4: Chain-length tie-break
1. RED: Add one failing test:
   - when identity and resolution tie, longer chain ranks first.
2. GREEN: Extend sort key with chain-length ordering.
3. REFACTOR: Keep invalid chain-length fallback behavior unchanged.

### Slice 5: Deterministic tie-break
1. RED: Add one failing test:
   - when all rank fields tie, sort by ID ascending.
2. GREEN: Add ID as final sort key element.

### Slice 6: Missing/invalid sequence identity behavior
1. RED: Add one failing test:
   - missing/invalid sequence identity ranks after valid values.
2. GREEN: Add sequence-identity priority bucket in the sort key.

## Verification Plan
1. Run targeted tests after each slice:
   - `pytest tests/pdbe/test_clustering.py -q`
   - `pytest tests/pdbe/test_result.py -q`
2. Run broader package tests at end:
   - `pytest tests/pdbe -q`
3. Confirm no integration changes were made outside this scope.

## Non-Goals
1. Do not wire `order_cluster_members` into `filter_pdb_results_on_resolution` in this task.
2. Do not alter pipeline ordering or clustering algorithm behavior beyond introducing the new helper and its tests.

TODO remove this file once implementation is complete and tested.
