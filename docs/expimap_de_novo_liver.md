# Liver De Novo expiMap Programs

## Design

The bounded ARCHS4 liver reference (1,000 leakage-excluded samples, 9,319
shared genes, and 1,140 Reactome terms) was used as the frozen reference. The
231 OSDR liver samples were mapped as a query with ten additional,
unconstrained expiMap programs. These programs are query-specific and do not
alter the official Reactome architecture.

All query-extension runs used raw counts, negative-binomial reconstruction,
150 mapping epochs, and GPU `cuda:0`. The de novo program decoder has L1
sparsity controlled by `gamma_ext`.

## Results

| `gamma_ext` | Nonzero genes per program | Zero-weight programs | Best aggregate program FDR | Best study-aware FDR |
| --- | --- | --- | --- | --- |
| 0.7 | 3,808-5,427 | 0 | 0.888 | 0.852 |
| 1.5 | 643-5,621 | 0 | 0.870 | 0.860 |
| 3.0 | 0-75 | 2 | 0.866 | 0.860 |

The low-penalty run is too diffuse for an interpretable program definition.
The high-penalty run produces compact decoder gene sets, but no de novo
program is associated with spaceflight after multiple-testing correction.
The closest nominal result is `unconstrained_4` in every run; even in the
compact run its Welch FDR is 0.866 and its accession-aware FDR is 0.860.

This is therefore a negative discovery result, not evidence for a new liver
spaceflight program. The compact run's program gene loadings and post-hoc
Reactome association tables are under
`outputs/expimap_archs4_reference_osdr_query_liver/query_denovo10_gamma3_150epoch/de_novo_analysis/`.

## Limitation

The original installed scArches HSIC implementation was numerically unstable for
this 1,140-term latent model: its normalizer overflowed, yielding NaN latent
values. These completed runs used L1 sparsity only. Later tutorial-style runs
use a repo runtime patch with a stable log-gamma bandwidth calculation when
HSIC is enabled. Any future de novo claim should add multi-seed stability
selection and, ideally, a larger ARCHS4 reference before biological
interpretation.
