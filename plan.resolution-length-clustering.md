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
7. Sequence identity is computed as chain_length / (uniprot_end - uniprot_start + 1). Higher is better.
8. Hard filters are strict by default for predictability and auditability.
9. Optional fill-to-top fallback can be enabled to relax only spread constraints in a fixed order.

## Definitions
- Model span: [uniprot_start, uniprot_end] computed on PdbResult.
- Covered union: the union of selected model spans.
- Novelty: number of residues in a candidate span not yet covered by the covered union.
- Sequence identity: chain_length / (uniprot_end - uniprot_start + 1), a float in (0, 1]. A value < 1 means the PDB chain has gaps relative to the UniProt span it covers. Returns 0.0 when chain_length or span cannot be determined.

## Selection Algorithm (Per UniProt Accession)
1. Validate top > 0, else raise ValueError.
2. Use each PdbResult.uniprot_start and PdbResult.uniprot_end as the candidate span.
3. Repeatedly pick one model until selected count reaches top or candidates are exhausted.
4. Before each pick, recompute novelty for every remaining candidate relative to the current covered union. This is always dynamic regardless of priority_profile.
5. Rank candidates by:
   1. Highest novelty first.
   2. Larger chain length next (more coverage is preferred when novelty is tied).
   3. Better resolution next (valid numeric lower is better; missing/invalid is worst).
   4. Higher sequence identity next (higher is better; 0.0 for invalid/error cases is worst).
   5. Stable deterministic PDB ID tie-break.
6. Add selected model to output and update covered union with its span.
7. Remove selected model from candidate pool and continue.
8. Filter policy:
  1. strict (default): apply all configured hard filters as strict gates, even if fewer than top are returned.
  2. fill_to_top (optional): if selected count < top after strict pass, relax in this order:
    1. Lower min_novelty_residues to 0.
    2. Raise max_overlap_ratio to 1.0.
    3. Disable per_region_cap.
    4. Keep max_resolution and min_sequence_identity strict.
      5. When multiple zero-novelty candidates compete, novelty is tied at 0 for all; ranking falls through to the next criterion in the priority profile naturally. No special case is required.
4. Missing or invalid resolution stays deprioritized as in current logic when `max_resolution` is not set.
5. When `max_resolution` is set, missing or invalid resolution is treated as failing the threshold and the model is discarded. Applies in both `strict` and `fill_to_top` modes.

## User Knobs (Steering Passed/Discarded PDBs)
These knobs are optional and default to current behavior when unset.

1. `max_resolution` (Knob #1: resolution ceiling)
  - Meaning: discard models with numeric resolution greater than the threshold.
  - Missing/invalid resolution handling: when `max_resolution` is set, missing or non-numeric resolution is treated as failing the threshold in both `strict` and `fill_to_top` modes. Rationale: `max_resolution` is a quality gate and unknown quality must not pass it silently.
  - Suggested default: `None` (no hard discard).

2. `min_sequence_identity` (Knob #3)
  - Meaning: discard models with `sequence_identity < threshold`.
  - Uses `sequence_identity = chain_length / (uniprot_end - uniprot_start + 1)`.
  - Invalid chain/range yields `0.0` and will be filtered out if threshold > 0.
  - Suggested default: `0.0`.

3. `min_novelty_residues` (Knob #4: novelty threshold)
  - Meaning: during iterative selection, a candidate must add at least this many newly covered residues to be eligible.
  - Applies per pick against the current covered union.
  - Suggested default: `0` (preserves current fill behavior).

4. `max_overlap_ratio` (Knob #5: overlap tolerance)
  - Meaning: optionally cap allowed overlap with already selected coverage.
  - Formula: `overlap_ratio = overlap_len / candidate_span_len` (candidate-relative).
  - Semantics: measures what fraction of the *incoming* candidate is redundant. A large model that partially overlaps can still pass because most of its span is novel. Consistent with agreed behavior #6 (broad models favoured when they extend coverage).
  - Candidate is ineligible when `overlap_ratio > threshold`.
  - Suggested default: `1.0` (no extra overlap restriction).

5. `per_region_cap` (Knob #8: per-region diversity cap)
  - Meaning: limit how many selected models can come from the same region.
  - Region definition: transitive overlap clustering using connected components. Two candidates are connected when their spans overlap by any amount (overlap_len > 0). A region is the full connected component under this relation. No internal threshold required.
  - Region membership is computed once over all candidates at the start of selection, before any model is picked. It is static: adding models to the output does not change region assignments.
  - Suggested default: `None` (disabled).

6. `priority_profile` (Knob #9: tie-break priority profile)
  - Meaning: user-selectable ranking preference among novelty, resolution, sequence identity, and chain length.
  - Novelty is always computed dynamically relative to the current covered union at each pick step, regardless of profile. The profile controls ranking priority only, not how novelty is computed.
  - Preserve deterministic PDB ID as final tie-break in all profiles.
  - Suggested built-ins:
    - `novelty_first` (current plan default): novelty > chain_length > resolution > identity > id.
    - `quality_first`: resolution > identity > novelty > chain_length > id.
    - `identity_first`: novelty > identity > chain_length > resolution > id.

7. `filter_policy` (strict vs fill behavior)
  - Meaning: controls whether hard filters may be relaxed to reach top.
  - Values:
    - `strict` (default): do not relax filters.
    - `fill_to_top`: relax only spread constraints using the fixed fallback order in Selection Algorithm.
  - Suggested default: `strict`.

### Backward Compatibility
- Existing function `filter_pdb_results_on_resolution(pdb_results, top)` stays valid.
- If no knobs are provided, behavior remains equivalent to the current plan semantics.
- Knobs can be introduced via optional parameters or a config object in a follow-up without breaking current callers.

## API Sketch

```python
from dataclasses import dataclass
from typing import Literal

PriorityProfile = Literal["novelty_first", "quality_first", "identity_first"]


@dataclass(frozen=True)
class ResolutionSelectionOptions:
  max_resolution: float | None = None
  min_sequence_identity: float = 0.0
  min_novelty_residues: int = 0
  max_overlap_ratio: float = 1.0
  per_region_cap: int | None = None
  priority_profile: PriorityProfile = "novelty_first"
  filter_policy: Literal["strict", "fill_to_top"] = "strict"


def filter_pdb_results_on_resolution(
  pdb_results: PdbResults,
  top: int,
  options: ResolutionSelectionOptions | None = None,
) -> PdbResults:
  ...
```

### Compatibility Wrapper Option

```python
def filter_pdb_results_on_resolution(
  pdb_results: PdbResults,
  top: int,
) -> PdbResults:
  return filter_pdb_results_on_resolution_with_options(
    pdb_results,
    top=top,
    options=ResolutionSelectionOptions(),
  )
```

### Example Usage

```python
# Keep current behavior
results = filter_pdb_results_on_resolution(pdb_results, top=3)

# Strict quality and spread
results = filter_pdb_results_on_resolution(
  pdb_results,
  top=5,
  options=ResolutionSelectionOptions(
    max_resolution=3.0,
    min_sequence_identity=0.7,
    min_novelty_residues=25,
    max_overlap_ratio=0.8,
    per_region_cap=2,
    priority_profile="novelty_first",
    filter_policy="strict",
  ),
)

# Try to fill to top by relaxing spread-only constraints
results = filter_pdb_results_on_resolution(
  pdb_results,
  top=5,
  options=ResolutionSelectionOptions(
    max_resolution=3.0,
    min_sequence_identity=0.7,
    min_novelty_residues=25,
    max_overlap_ratio=0.8,
    per_region_cap=2,
    priority_profile="novelty_first",
    filter_policy="fill_to_top",
  ),
)
```

### Validation Rules (Sketch)
- `top > 0`
- if set: `max_resolution > 0`
- `0.0 <= min_sequence_identity <= 1.0`
- `min_novelty_residues >= 0`
- `0.0 <= max_overlap_ratio <= 1.0`
- if set: `per_region_cap > 0`
- `filter_policy in {"strict", "fill_to_top"}`

## Use Cases

These use cases are drawn from https://github.com/haddocking/protein-quest/issues/102 and are
encoded as `test_on_sequence_domains` and `test_on_sequence_overlaps` in `tests/pdbe/test_result.py`.

**Scale:** each `=` represents ~25 residues. Positions are approximate for readability.

---

### Use Case 1: Distinct Domains (test_on_sequence_domains)

**Problem:** Eight models exist but six of the eight best-resolution models cover only two sub-regions.
Pure resolution ranking returns three models all covering 500–1000 (m6, m7, m8) — useless redundancy.

```
Residue: 1        250  400  500                   1000
         |         |    |    |                      |
1AAA  3.6[==========]
2BBB  5.4[==========]
3CCC  2.1[==========]
4DDD  8.1          [========]
5EEE  4.6          [========]
6FFF  1.3                   [=====================]
7GGG  1.4                   [=====================]
8HHH  1.6                   [=====================]
```

**Old behaviour** (pure resolution, top=3): keeps 6FFF, 7GGG, 8HHH — all in the same domain.

**New behaviour** (novelty-first, top=3): 

```
Pick 1 — covered={}
  novelty: 6FFF=7GGG=8HHH=501, 1AAA=2BBB=3CCC=250, 4DDD=5EEE=201
  three-way tie at 501 → resolution: 6FFF(1.3) wins
  covered: {500-1000}

Pick 2 — covered={500-1000}
  novelty: 1AAA=2BBB=3CCC=250, 4DDD=5EEE=201, 7GGG=8HHH=0
  three-way tie at 250 → resolution: 3CCC(2.1) wins
  covered: {1-250, 500-1000}

Pick 3 — covered={1-250, 500-1000}
  novelty: 4DDD=5EEE=150 (251-400 uncovered), rest=0
  tie at 150 → resolution: 5EEE(4.6) beats 4DDD(8.1)
  covered: {1-400, 500-1000}

Selected: 3CCC, 5EEE, 6FFF  ✓
```

```
Residue: 1        250  400  500                   1000
         |         |    |    |                      |
3CCC  2.1[==========]            ← pick 2
5EEE  4.6          [========]   ← pick 3
6FFF  1.3                   [=====================] ← pick 1
                    ^^^^^^^^
                    gap 401-499 (no model available)
```

---

### Use Case 2: Overlapping Broad Models (test_on_sequence_overlaps)

**Problem:** Broad models m9 (1–600) and m10 (1–1000) have poor resolution individually but together
with m6 they give the most useful coverage for docking. Pure resolution would discard them in favour
of small high-resolution models.

```
Residue: 1        250  400  500   600              1000
         |         |    |    |     |                 |
1AAA  3.6[==========]
4DDD  8.1          [========]
6FFF  1.3                   [=======================]
9III  4.2[========================]
10JJJ 1.4[========================================]
```

**New behaviour** (novelty-first, top=3):

```
Pick 1 — covered={}
  novelty: 10JJJ=1000, 9III=600, 6FFF=501, 1AAA=250, 4DDD=201
  10JJJ wins outright
  covered: {1-1000}

Pick 2 — covered={1-1000}
  novelty: all=0 — fall to next tiebreaker
  **⚠ plan order says resolution next: 6FFF(1.3) wins  → picks 6FFF**
  **but test expects 9III here (broader chain_length=600 > 6FFF=501)**

Pick 3 — novelty: all=0
  **⚠ plan order (resolution) picks 1AAA(3.6) over 9III(4.2)**
  **but test expects 9III**

Expected: 9III, 10JJJ, 6FFF
```

Verification with corrected order (novelty → chain_length → resolution → identity → id):

```
Pick 2 — novelty=0, chain_length: 9III=600 > 6FFF=501 > 1AAA=250 > 4DDD=201
  → 9III selected.

Pick 3 — novelty=0, chain_length: 6FFF=501 > 1AAA=250 > 4DDD=201
  → 6FFF selected.

Selected: 10JJJ, 9III, 6FFF  ✓
```

## Implementation Scope
- Immediate implementation keeps function signature unchanged: filter_pdb_results_on_resolution(pdb_results, top).
- API sketch with `options` is a follow-up evolution path, not required for this implementation pass.
- Keep pipeline order unchanged:
  1. filter_pdb_results_on_chain_length first.
  2. spread-aware filter_pdb_results_on_resolution second.
- Do not add new external dependencies for interval/graph logic.

## Files To Touch
- src/protein_quest/pdbe/result.py
  - ~~Add cached_property `sequence_identity` to PdbResult~~ — done.
  - Add helper `_sequence_identity_or_zero(entry: PdbResult) -> float` (mirrors `_chain_length_or_zero`).
  - Add helper `_compute_novelty(covered_intervals: list[tuple[int, int]], span: tuple[int, int]) -> int`:
    returns the number of residues in `span` not covered by any interval in `covered_intervals`.
  - Add helper `_build_region_map(candidates: Iterable[PdbResult]) -> dict[PdbResult, int]`:
    assigns each candidate a region id using transitive connected-component clustering over spans (two candidates connected when their spans overlap by any amount).
  - Rewrite body of `filter_pdb_results_on_resolution` as an iterative pick loop:
    - Each iteration recomputes novelty for all remaining candidates against the current covered union.
    - Builds a per-pick sort key using the priority profile order.
    - Enforces min_novelty_residues, max_overlap_ratio, per_region_cap, max_resolution, min_sequence_identity gates.
    - Applies filter_policy fallback when strict pass underflows top.
  - `_sort_resolution_key` is replaced by a per-pick ranking tuple builder; remove the old function.
- tests/pdbe/test_result.py

## Test Plan
1. Run tests/pdbe/test_result.py as primary contract suite.
2. Confirm issue-102-inspired domain/overlap tests pass.
3. Add explicit test for split-range span using uniprot_start=1 and uniprot_end=584.
4. Add tests for PdbResult.sequence_identity computed correctly (continuous range and split range).
5. Add tie determinism tests where novelty and resolution match but identity differs.
6. Add test: PDB with identity=0.0 (invalid chain range) ranks below any PDB with valid identity.
7. Run broader tests under tests/pdbe to catch regressions.

## Notes
- This plan intentionally avoids normalization by full UniProt length.
- Sequence identity is derived from existing PdbResult fields — no SPARQL or data pipeline changes needed.
- Novelty-driven iterative selection resolves the overlap-cluster conflict where one-component clustering would collapse useful broad models.

TODO remove this file once implementation is complete and tested.
