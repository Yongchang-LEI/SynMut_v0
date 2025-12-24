# Dinucleotide Bounding & max.dinu Integration Proposal

## 1. Objectives
- certify attainable upper/lower bounds for any junction dinucleotide (e.g., CpG) before mutation.
- integrate bounds with `dinu_to()`'s `max.dinu` / `min.dinu` switches so the greedy mutators know their feasible targets.
- expose the results through both user-facing knobs (vignette + API) and internal optimizers (goals modulation / `opt_mode`).

## 2. Current Mechanism Recap
- `R/dinu_to.R` already exposes `max.dinu` / `min.dinu` and two strategies:
  - `dinu_to.not.keep()` rewrites every mutable codon to an "optimal" synonym returned by `get_optimal_codon()` when `keep = FALSE`.
  - `dinu_to.keep()` performs swap-only adjustments on 3'→1' junctions so codon usage is preserved when `keep = TRUE`.
- The methods succeed in **monotonic** maximization/minimization but do not report what range of junction counts is reachable nor whether a user target is infeasible.

## 3. Gap Analysis vs. Requirements Document
| Requirement from bounds memo | Status in repo | Missing pieces |
| --- | --- | --- |
| Fast analytical bounds per junction | ❌ not implemented | Need per-junction capability probes prior to running greedy edits. |
| Exact relaxed bounds (DP-Relax) | ❌ | Need codon-by-codon DP over synonym sets. |
| Multiset-aware bounds | ⚠️ Partial (`keep = TRUE` swaps) | Need MILP or Lagrangian outer bounds + MC inner bounds for swap-only sequences. |
| Integration into goals modulation | ❌ | Need API returning `[L_d, U_d]` and hooking into tolerance mapping. |
| Reporting / validation | ❌ | Need logging & CI checks showing MC inner ≤ DP outer, etc. |

## 4. Proposed Solution Stack
1. **Analytical Junction Probes** (O(n))
   - Precompute, for each boundary i→i+1, whether synonyms exist that can host the desired 3'/1' bases. Sum them to form loose outer bounds (`analytical_upper`, `analytical_lower`).
   - Reuse codon tables already assembled in `dinu_to()` (i.e., `seqinr::ucoweight`).

2. **DP-Relax Engine** (Exact relaxed bounds)
   - Implement `dinu_bounds_dp_relax(seq, aa, region, target_dinu, maximize = TRUE, numcode)` that runs Viterbi-style DP across synonyms, weighting transitions by whether their junction equals the target dinucleotide.
   - Acts as the certified bound when `keep = FALSE` (or when multiset is ignored); also as an outer bound otherwise.
   - Returns both best and worst achievable counts (`max_dp`, `min_dp`).

3. **Swap-Only Enhancements**
   - Wrap the existing `dinu_to.keep()` logic in a Monte Carlo harness:
     - Randomly permute codons within AA groups (respecting conserved regions).
     - Apply `dinu_to.keep()` greedily to polish.
     - Track observed extrema (`mc_inner_min`, `mc_inner_max`).
   - Optionally add a Lagrangian relaxer that penalizes codons for violating multiset counts but still solves via DP; iterate multipliers for tighter outer bounds.

4. **MILP Windows (Optional Certification)**
   - For challenging segments, build a `glpkAPI`/`ompr` model with decision variables `x_{i,c}` and adjacency `z_{i,c,d}` to enforce exact multiset counts.
   - Solve windows of 400–800 codons, sum their maxima as an outer bound, and optionally stitch solutions together.

5. **API & UX Integration**
   - New exported helper `compute_dinu_bounds(regioned_dna, dinu = "CG", keep = FALSE, methods = c("analytical","dp_relax","mc"), ...)` returning:
     ```r
     list(
       lower = ..., upper = ...,
       outer_certified = list(analytical = ..., dp_relax = ..., lagrangian = ..., milp = ...),
       inner_empirical = list(mc = ...),
       diagnostics = ...
     )
     ```
   - `dinu_to()` consults `compute_dinu_bounds()` when a numeric target (e.g., "max CpG to >= k") is supplied, and emits warnings if the target exceeds the feasible `[L, U]`.
   - Goals-modulation uses normalized residual `(X_d - target)/(U - L)` to scale penalties.

6. **Validation / Reporting**
   - Unit tests ensuring `mc_inner_max <= dp_relax_upper + ε` and `mc_inner_min >= dp_relax_lower - ε` for relaxed cases.
   - Vignette updates (this document + tutorials) explaining how to interpret `max.dinu` with bounds.
   - Logging hooks so benchmarking scripts can print `% of theoretical max` for SA/2-Opt runs.

## 5. Implementation Roadmap
1. **Foundations (1 sprint)**
   - Implement analytical probe + DP-Relax helper.
   - Add `compute_dinu_bounds()` API + tests covering CpG/TpA use-cases with/without conserved regions.
2. **Swap-Only Coverage (1–1.5 sprints)**
   - Monte Carlo harness leveraging `dinu_to.keep()`.
   - Optional Lagrangian relaxer for tighter outer bounds when multiset fixed.
3. **Advanced Certification (as needed)**
   - MILP window solver and reporting utilities.
4. **Integration & Docs (0.5 sprint)**
   - Wire bounds into `dinu_to()` & `goals_modulation` workflows.
   - Publish vignette + user guide (linking to https://guhaogao.com/files/SynMut/index.html#synonymous-mutants-mimicking-a-specific-codon-usage-pattern).

## 6. Expected Benefits
- Early infeasibility detection for ambitious `max.dinu` targets (e.g., CpG = 50 ± 1 when only 35 are possible).
- Quantitative benchmarking of heuristic optimizers against proven bounds.
- Improved user confidence via certified ranges and normalized progress metrics.

