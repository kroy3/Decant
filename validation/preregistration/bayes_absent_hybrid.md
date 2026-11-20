# Pre-registration: `bayes_absent` with the hybrid cap, on unused mixtures

Written before any of the data below was downloaded. The commit that adds this
file, together with the frozen code, is the timestamp. Amendments are listed at
the end, each made before the result it could affect, with its reason.

## What changed, and why

`bayes_absent` was judged on three fresh 10x mixtures (BENCHMARK.md 1g) and
failed two pre-stated criteria: mean retain 0.9837 against DecontX 0.9847, and
on PBMC 4k it removed more counts too high to be ambient (likely doublets)
than `bayes_cell`.

Development since then used only data already seen: the hgmm 1k/5k/10k v3
mixtures and PBMC 4k.

- PBMC: the extra loss came entirely from the cap test, not from the absence
  posterior. With the same cap test, `bayes_absent` and `bayes_cell` gave
  identical PBMC numbers.
- The new cap test, `cap_test = "hybrid"`, sets an entry aside if it passes
  the per-cluster test, or the per-entry test with at least 5 UMIs. On the
  v3 mixtures, mean retain/removal: cluster 0.9817 / 0.8620, hybrid
  0.9817 / 0.8615. On PBMC, hybrid matches `bayes_cell` exactly
  (Ig/monocyte/NK: 12.7 / 16.1 / 32.1% of non-ambient counts removed).
- Retain is unchanged. Its loss on hgmm 10k comes mostly from native-species
  genes the posterior calls absent in a cluster. Half of those are at least
  3-fold enriched in another cluster of the same species: that is
  within-species ambient RNA, which this scorer counts as loss. The rest are
  low-count genes whose call follows the prior.
- Lowering the absence prior traded removal for retain along the same
  frontier (0.9832 / 0.8524). That is a tuning knob, not a fix, and was not
  adopted.

**Prediction:** on the v3 mixtures, retain is about 0.2 points below DecontX.
The retain criterion is therefore expected to fail again.

## Data (not yet downloaded)

The 10x human/mouse 1:1 mixtures not used so far:

- `hgmm_1k`: Cell Ranger 2.1.0, v2 chemistry.
- `hgmm_1k_v2`: Cell Ranger 3.0.0, v2 chemistry.
- `hgmm_100`: Cell Ranger 2.1.0, about 100 cells.

Before scoring, the two 1k sets are checked for being the same library: the
share of filtered barcodes in common and the total UMIs. If more than half of
the barcodes are shared, only `hgmm_1k_v2` is scored. This check reads no
corrected counts.

The pipeline is unchanged from `benchmark_removal_heldout.R`:

- cells called by Cell Ranger, cross-species doublets excluded;
- empty droplets with 1 to 100 UMIs;
- Seurat clusters at resolution 0.1;
- the same per-gene retain/removal scorer (Cargnelli et al. 2026);
- SoupX and DecontX with the same wrappers.

## Pre-stated criteria

Primary: the mean over `hgmm_1k` and `hgmm_1k_v2`, or `hgmm_1k_v2` alone if
the two are one library.

1. `bayes_absent` mean removal >= DecontX mean removal.
2. `bayes_absent` mean retain >= DecontX mean retain.
3. The Smart-seq2 negative control: retain >= 0.99.
4. The simulation gate: counts are never added, and mean L1 is no worse
   than `bayes_cell` by more than 0.002.
5. PBMC 4k: counts that cannot be ambient are removed no more than with
   `bayes_cell`. This is not a held-out test, since PBMC was used to choose
   the hybrid cap; it is reported as a regression check.

`hgmm_100` is reported, but not scored: there are too few cells for stable
clusters.

`bayes_absent` becomes the default only if criteria 1 to 4 pass. If criterion
2 fails, it stays opt-in and the README states the trade-off in numbers.
Whatever the outcome, nothing is tuned on these data afterwards.

## Amendments
None. Library check: hgmm_1k and hgmm_1k_v2 share 5 of about 1,000 filtered
barcodes, so both are scored as pre-stated. Results: BENCHMARK.md 1h.

## Decision after the result (not an amendment)

Criterion 2 failed, so under the rule above `bayes_absent` should have stayed
opt-in. After seeing the result the maintainer made it the default anyway.
That overrides the pre-registered rule and is recorded as such, with its
reasons and known costs, in BENCHMARK.md 1h. The result itself is unchanged:
criterion 2 failed.
