# Resolution-Length Clustering Plan

## Goal
Improve PDB result filtering so selected models are not only high quality (resolution) but also provide useful sequence spread.

## Agreed Behavior
1. Use sequence envelopes from uniprot_chains for spread logic.
2. For split ranges like A=398-459,A=74-386,A=520-584,A=1-53, use envelope [1,584].
3. Ignore full UniProt sequence length for this feature.
4. Selection objective is novelty coverage first, then quality tie-breaks.
5. If fewer distinct regions exist than top, continue filling remaining slots up to top.
6. Broad models may outrank nested narrow models when they add better coverage utility.

## Definitions
- Envelope: [min_start, max_end] parsed from all valid ranges in one uniprot_chains value.
- Covered union: the union of envelopes from already selected models.
- Novelty: number of residues in a candidate envelope not yet covered by the covered union.

## Selection Algorithm (Per UniProt Accession)
1. Validate top > 0, else raise ValueError.
2. For each PdbResult, parse envelope from uniprot_chains.
3. Repeatedly pick one model until selected count reaches top or candidates are exhausted.
4. Rank candidates by:
   1. Highest novelty first.
   2. Better resolution next (valid numeric lower is better; missing/invalid is worst).
   3. Larger chain length next (existing behavior).
   4. Stable deterministic PDB ID tie-break.
5. Add selected model to output and update covered union with its envelope.
6. Remove selected model from candidate pool and continue.

## Invalid Data Handling
1. Invalid chain ranges (for example A=-) produce no envelope and contribute zero novelty.
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
3. Add explicit test for split-range envelope [1,584].
4. Add tie determinism tests where novelty and resolution match.
5. Run broader tests under tests/pdbe to catch regressions.

## Notes
- This plan intentionally avoids normalization by full UniProt length.
- Novelty-driven iterative selection resolves the overlap-cluster conflict where one-component clustering would collapse useful broad models.

TODO remove this file once implementation is complete and tested.
