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
