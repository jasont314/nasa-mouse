# Adaptive Generative Benchmark Decisions

This log records why bounded experiments were continued, stopped, or added. It is
not a substitute for the resolved configurations and per-run manifests.

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
5. **Gate GeneJEPA by A100 feasibility.** The queued pilot keeps the exact 768-wide,
   512-latent, 24-block architecture but runs only two profiles. It will not be called
   paper-native duration; a longer representation baseline is scheduled only if the
   exact architecture fits and one-step runtime is tractable.
6. **Reject the strict broad WGAN as a fidelity finalist.** Released-code timing
   restored epoch 11 and stopped at epoch 61 after 3,533 seconds. The run has mean/SD
   correlations 0.323/0.563, precision/recall 0.586/0.490, and adversarial accuracy
   0.950. It passes diversity and memorization but fails the fixed fidelity gate, so
   additional seeds would not resolve the model-class mismatch.
7. **Keep GeneJEPA representation-only and bound its duration.** Exact batch 92 with
   accumulation two fits the A100 40 GB at 29.88 GB and takes 6.57 seconds per 184
   profile exposures. A literal 50-million-exposure run is about 20 single-GPU days.
   Any shorter mouse run is a declared duration adaptation and GeneJEPA guidance is
   attempted only if its held-out tissue representation improves over expression and
   the existing practical screen.
8. **Reject direct OSDR ModelDDIM.** It reached mean correlation 0.804 but SD
   correlation 0.153, precision/recall 0.196/0.700, composite 0.367, and global
   FLT-minus-GC delta correlation 0.125. Tissue consistency was useful, but fidelity
   and condition-effect gates failed; the locked test and extra seeds remain blocked.
9. **Continue exact ARCHS4 transfer before changing architecture or duration.** The
   reference MaxAbs distribution is compatible with OSDR, all tissue weights are
   explicitly mapped, and new FLT/GC columns are learnable during fine-tuning. The
   first AMP overflows exposed an accounting bug; the clean restart advances EMA,
   scheduler, and global step only after a successful optimizer step and records each
   skip. Validation, not training loss, decides whether transfer warrants replication.
10. **Reject the reference-only DDIM transfer baseline.** It recovers the pooled
    FLT/GC effect better than direct training (delta correlation 0.487, direction
    agreement 0.637) and has mean/SD correlations 0.937/0.841, but nearest-neighbor
    two-sample accuracy is 0.973. Its fidelity composite is 0.607 and synthetic
    classifier balanced accuracy is 0.529 versus 0.537 from real training data. It
    therefore fails the fixed fidelity gate despite passing condition, diversity, and
    memorization gates.
11. **Test and reject function-preserving DDIM transfer.** The new expansion copies
    shared tissue embeddings and condition-input columns and adjusts condition-layer
    biases so every expanded shared-tissue class is numerically identical to the
    source denoiser before fine-tuning. Exact parity is covered by a numerical test.
    Validation at 1,000, 3,000, and 5,000 fine-tune epochs peaked at epoch 3,000 with
    composite 0.623, mean/SD correlations 0.937/0.850, precision/recall F1 0.625,
    and two-sample accuracy 0.960. The condition effect passes, but local real-versus-
    synthetic separation still fails fidelity, so no extra seeds or locked test are
    permitted.
12. **Stop stochastic DDIM sampling as a rescue branch.** At the completed
    function-preserving checkpoint, eta 0.25 and 0.5 reduce the composite to 0.602
    and 0.587 and leave two-sample accuracy at 0.967 and 0.969. Eta 1.0 is not run
    because the monotonic degradation already falsifies the proposed remedy.
13. **Prefer CPM/log/z-score for the bounded WGAN transfer test.** Direct OSDR WGAN
    training with strict preprocessing has composite 0.424; the NASA CPM adaptation
    improves it to 0.468 and raises SD correlation from 0.603 to 0.766. Both pass the
    pooled condition-effect gate but fail fidelity, chiefly because two-sample
    accuracy is 1.000/0.998 and precision/recall F1 is 0.250/0.238. Only the stronger
    CPM branch proceeds to a 100-reference/100-fine-tune-epoch ARCHS4 transfer screen.
14. **Gate GeneJEPA guidance with grouped tissue prediction.** The practical mouse
    adaptation reached held-out-series tissue balanced accuracy 0.453, versus 0.839
    from the same preprocessed expression. The exact architecture receives one
    43,744-exposure replacement-sampled epoch as a duration screen. Guidance is
    considered only if its representation exceeds both baselines; otherwise the
    estimated 20-GPU-day paper-duration run and guided-diffusion branch are stopped.
15. **Reject ARCHS4-pretrained WGAN transfer and strengthen the condition gate.**
    After 100 ARCHS4 epochs and 66 OSDR fine-tuning epochs, the model has composite
    0.437, mean/SD correlations 0.894/0.784, precision/recall F1 0.071, and nearest-
    neighbor two-sample accuracy 1.000. Its pooled FLT/GC delta correlation of 0.805
    is misleading: the accession-aware meta-effect correlation is -0.022 and zero
    of nine held-out tissues passes the same correlation/direction thresholds.
    Pooled conditional recovery can therefore reflect tissue/accession composition.
    Future augmentation requires both pooled and accession-aware effect gates. This
    post-screen correction does not change the WGAN decision because fidelity had
    already failed; no additional WGAN seeds or locked-test evaluation are run.
16. **Stop GeneJEPA guidance after the exact-architecture duration screen.** The
    768-wide, 512-latent, 24-block model completed 43,744 ARCHS4 exposures in 1,074
    training seconds on the A100, with 31.68 GB peak allocated memory and finite
    loss. On the balanced held-out-series figure cohort, its tissue probe reached
    0.703 balanced accuracy and 0.701 macro F1, improving on the practical model's
    0.453/0.423 but remaining below expression's 0.839/0.840. Embedding and UMAP
    silhouettes were -0.176 and -0.215. It therefore fails the predeclared guidance
    gate. This one-epoch screen made only 238 optimizer updates, below the paper's
    2,000-step EMA warmup and far below 50 million paper-duration exposures, so the
    result is a bounded-compute stopping decision rather than evidence against a
    fully trained GeneJEPA. GeneJEPA still has no expression decoder; no guided-
    diffusion branch or longer representation run is started.
