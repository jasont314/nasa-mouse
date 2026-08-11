# Adaptive Generative Benchmark Decisions

This log records why bounded experiments were continued, stopped, or added. It is
not a substitute for the resolved configurations and per-run manifests.

Entries below preserve the metric language used when each decision was made. The
former scalar composite was retired on 2026-07-16; current selection requires every
paper-aligned metric to pass independently, so no historical composite value can
promote a model.

## 2026-07-17

1. **Use a replicated within-study split for the conditional-generation target.**
   The unseen-accession experiments tested study extrapolation as well as synthetic
   fidelity and failed. The replacement split contains 781 train, 536 validation,
   and 293 locked-test profiles across 75 accessions, stratified within
   accession/tissue/condition. This answers generation for represented studies and
   is explicitly not an unseen-study benchmark.
2. **Promote the factorized ModelDDIM after repeated validation.** Starting from the
   exact ARCHS4 tissue model, the selected adapter uses study and material
   conditioning, rank-512 domain LoRA, and correlation-regularized refinement. Its
   train-only condition-blind mean/covariance calibration passed all broad fidelity
   metrics and pooled FLT/GC recovery across the required fraction of six validation
   generations. Validation muscle recovery was unstable, so muscle remained a
   diagnostic rather than part of broad finalist selection.
3. **Lock the final protocol before opening test.** Commit `7e8dec7` fixed four
   generation seeds, a 75% repeat threshold, independent metric gates, classifier
   utility analysis, and overwrite refusal. The one-time test then passed the broad
   rule: mean Corr/P/R/F1/AA were 0.977/0.998/0.997/0.997/0.458, FD ratio was 0.075,
   and pooled FLT/GC recovery passed three of four generations. All four muscle
   accession diagnostics passed, but exact concordant LOO-stable gene hits remained
   0, 0, 0, and 1.
4. **Do not claim paper parity or augmentation benefit.** All test generations fall
   below the separate 0.98 Corr target despite passing the finite-sample Corr floor.
   Real-plus-synthetic FLT/GC training scored 0.734 balanced accuracy versus 0.754
   for real-only training. The accepted scope is therefore within-study conditional
   simulation, not strict paper parity, unseen-study generation, or improved
   downstream classification.
5. **Reject skeletal-muscle specialist branches.** ARCHS4-initialized, pooled-model
   transfer, frozen-control, and effect-refined specialist runs were less faithful
   and more seed-sensitive than the pooled factorized generator. They are retained
   as negative evidence and are not promoted.

## 2026-07-16

1. **Promote the completed exact ARCHS4 ModelDDIM baseline.** It has common-score
   composite 0.9571, direct precision/recall 0.9655/0.8645, synthetic-to-real tissue
   balanced accuracy 0.8687, memorization fraction 0.0035, and diversity SD ratio
   0.842. It passes both fixed diversity and memorization gates.
2. **Stop the first broad WGAN attempt.** Its model architecture was correct, but
   the named native preprocessing profile inserted CPM and the early-stop checks
   occurred after epochs 5, 10, and so on. The partial epoch-5 gamma was 0.9731. The
   output is retained as an interrupted NASA CPM adaptation, not paper-native evidence.
3. **Rerun strict WGAN preprocessing and source timing.** The replacement applies
   direct `log1p` plus train-fitted gene z-score and checks gamma after epochs 1, 6,
   11, and so on. Gamma rose from 0.9579 to 0.9641 through epoch 11, so the run was
   continued under released-code patience.
4. **Use exact ModelDDIM for the OSDR conditional arm.** The configurable 512-unit
   adapter is a useful smoke-test proxy but cannot answer paper-architecture transfer.
   The exact extension uses 1,610 API profiles, accession grouping, full-transcriptome
   TPM, 974 landmarks, and 48 tissue-by-FLT/GC classes.
5. **Reject the strict broad WGAN as a fidelity finalist.** Released-code timing
   restored epoch 11 and stopped at epoch 61 after 3,533 seconds. The run has mean/SD
   correlations 0.323/0.563, precision/recall 0.586/0.490, and adversarial accuracy
   0.950. It passes diversity and memorization but fails the fixed fidelity gate, so
   additional seeds would not resolve the model-class mismatch.
6. **Reject direct OSDR ModelDDIM.** It reached mean correlation 0.804 but SD
   correlation 0.153, precision/recall 0.196/0.700, composite 0.367, and global
   FLT-minus-GC delta correlation 0.125. Tissue consistency was useful, but fidelity
   and condition-effect gates failed; the locked test and extra seeds remain blocked.
7. **Continue exact ARCHS4 transfer before changing architecture or duration.** The
   reference MaxAbs distribution is compatible with OSDR, all tissue weights are
   explicitly mapped, and new FLT/GC columns are learnable during fine-tuning. The
   first AMP overflows exposed an accounting bug; the clean restart advances EMA,
   scheduler, and global step only after a successful optimizer step and records each
   skip. Validation, not training loss, decides whether transfer warrants replication.
8. **Reject the reference-only DDIM transfer baseline.** It recovers the pooled
    FLT/GC effect better than direct training (delta correlation 0.487, direction
    agreement 0.637) and has mean/SD correlations 0.937/0.841, but nearest-neighbor
    two-sample accuracy is 0.973. Its fidelity composite is 0.607 and synthetic
    classifier balanced accuracy is 0.529 versus 0.537 from real training data. It
    therefore fails the fixed fidelity gate despite passing condition, diversity, and
    memorization gates.
9. **Test and reject function-preserving DDIM transfer.** The new expansion copies
    shared tissue embeddings and condition-input columns and adjusts condition-layer
    biases so every expanded shared-tissue class is numerically identical to the
    source denoiser before fine-tuning. Exact parity is covered by a numerical test.
    Validation at 1,000, 3,000, and 5,000 fine-tune epochs peaked at epoch 3,000 with
    composite 0.623, mean/SD correlations 0.937/0.850, precision/recall F1 0.625,
    and two-sample accuracy 0.960. The condition effect passes, but local real-versus-
    synthetic separation still fails fidelity, so no extra seeds or locked test are
    permitted.
10. **Stop stochastic DDIM sampling as a rescue branch.** At the completed
    function-preserving checkpoint, eta 0.25 and 0.5 reduce the composite to 0.602
    and 0.587 and leave two-sample accuracy at 0.967 and 0.969. Eta 1.0 is not run
    because the monotonic degradation already falsifies the proposed remedy.
11. **Prefer CPM/log/z-score for the bounded WGAN transfer test.** Direct OSDR WGAN
    training with strict preprocessing has composite 0.424; the NASA CPM adaptation
    improves it to 0.468 and raises SD correlation from 0.603 to 0.766. Both pass the
    pooled condition-effect gate but fail fidelity, chiefly because two-sample
    accuracy is 1.000/0.998 and precision/recall F1 is 0.250/0.238. Only the stronger
    CPM branch proceeds to a 100-reference/100-fine-tune-epoch ARCHS4 transfer screen.
12. **Reject ARCHS4-pretrained WGAN transfer and strengthen the condition gate.**
    After 100 ARCHS4 epochs and 66 OSDR fine-tuning epochs, the model has composite
    0.437, mean/SD correlations 0.894/0.784, precision/recall F1 0.071, and nearest-
    neighbor two-sample accuracy 1.000. Its pooled FLT/GC delta correlation of 0.805
    is misleading: the accession-aware meta-effect correlation is -0.022 and zero
    of nine held-out tissues passes the same correlation/direction thresholds.
    Pooled conditional recovery can therefore reflect tissue/accession composition.
    Future augmentation requires both pooled and accession-aware effect gates. This
    post-screen correction does not change the WGAN decision because fidelity had
    already failed; no additional WGAN seeds or locked-test evaluation are run.
13. **Reject all nine matched liver harmonization arms.** The comparison fixes the
    same 119/50/70 train/validation/locked-test profiles and 974 genes for no
    correction, Ilangovan per-study normalization, mentor two-stage scaling, ComBat,
    ComBat-seq, three official MBatch methods, and MOBER. No arm passes every
    independent fidelity gate, and all fail both pooled FLT/GC and accession-aware
    effect recovery. Mentor two-stage is the closest balanced result at Corr 0.348,
    precision/recall/F1 0.440/1.000/0.611, AA 0.690, and FD ratio 0.977. MOBER's Corr
    0.808 is not promotable because precision is 0.260, AA is 0.770, and FD ratio is
    33.311. The OSD-379 test partition remains locked and no arm advances to five
    seeds.
14. **Bound GPU parallelism by throughput rather than allocated memory.** Six exact
    OSDR ModelDDIM jobs fit concurrently on the A100 40 GB at about 6.2 GiB process
    residency each, but saturate GPU compute and slow each job substantially. Run
    distinct pending arms concurrently while they improve aggregate throughput; do
    not spend free memory on repeat seeds before a single-seed arm passes all gates.
