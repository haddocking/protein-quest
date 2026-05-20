# Resolution-Length Clustering Plan

## Goal
Improve PDB result filtering so selected models are not only high quality (resolution) but also provide useful sequence spread.

## Agreed Behavior
1. Use PdbResult.uniprot_start and PdbResult.uniprot_end for spread logic.
2. For split ranges like A=398-459,A=74-386,A=520-584,A=1-53, use uniprot_start=1 and uniprot_end=584.
3. Ignore full UniProt sequence length for this feature.
4. Selection objective is novelty coverage first, then quality tie-breaks.
5. If fewer distinct regions exist than top, continue filling remaining slots up to top.
6. Broad models may outrank nested narrow models when they add better coverage utility.

## Definitions
- Model span: [uniprot_start, uniprot_end] computed on PdbResult.
- Covered union: the union of selected model spans.
- Novelty: number of residues in a candidate span not yet covered by the covered union.

## Selection Algorithm (Per UniProt Accession)
1. Validate top > 0, else raise ValueError.
2. For all pdb chains belonging to same uniprot accession
3. Cluster them based on start/stop ranges
4. For each cluster, order by sequence identity (highest first), resolution ascending (best first), length descending (longest first).
5. Pick top from each cluster. For example given top 100 and 4 clusters take top 25 from each cluster

Use scipy clustering
compare function of 2 pdb result objects by how many residues two pdb chains share.

## Invalid Data Handling
1. Invalid chain ranges (for example A=-) fail uniprot_start/uniprot_end and contribute zero novelty.
2. Existing chain-length fallback behavior remains unchanged.
3. Missing or invalid resolution stays deprioritized as in current logic.

## Implementation Scope
- Keep function signature unchanged: filter_pdb_results_on_resolution(pdb_results, top).
- Keep pipeline order unchanged:
  1. filter_pdb_results_on_chain_length first.
  2. spread-aware filter_pdb_results_on_resolution second.
- Do not add new external dependencies for interval/graph logic.

## Files To Touch
- src/protein_quest/pdbe/result.py
- tests/pdbe/test_result.py

## Test Plan
1. Run tests/pdbe/test_result.py as primary contract suite.
2. Confirm issue-102-inspired domain/overlap tests pass.
3. Add explicit test for split-range span using uniprot_start=1 and uniprot_end=584.
4. Add tie determinism tests where novelty and resolution match.
5. Run broader tests under tests/pdbe to catch regressions.

## Notes
- This plan intentionally avoids normalization by full UniProt length.
- Novelty-driven iterative selection resolves the overlap-cluster conflict where one-component clustering would collapse useful broad models.

TODO remove this file once implementation is complete and tested.
