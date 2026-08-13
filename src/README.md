# Source packages

The repository keeps project-owned code in six packages. The split reflects the
development history and the different model runtimes; it is not six independent
projects.

| Package | Responsibility |
|---|---|
| [`nasa_mouse_glare`](nasa_mouse_glare/README.md) | NASA OSDR API ingestion, tissue normalization, TMS preparation, GLARE adaptation, MOBER integration, and exploratory validation |
| [`nasa_mouse_expimap`](nasa_mouse_expimap/README.md) | Reactome-mask preparation, ARCHS4 reference training, OSDR query mapping, pathway analysis, and expiMap paper builds |
| [`nasa_mouse_generative`](nasa_mouse_generative/README.md) | Shared configuration, cohort assembly, preprocessing, harmonization adapters, experiment matrices, and common metrics |
| [`nasa_mouse_diffusion`](nasa_mouse_diffusion/README.md) | Conditional diffusion training and the paper-parity DDIM implementation used for the final synthetic analysis |
| [`nasa_mouse_wgan`](nasa_mouse_wgan/README.md) | Conditional WGAN-GP training, generation, calibration, and comparison analyses |
| [`nasa_mouse_internship_report`](nasa_mouse_internship_report/README.md) | Figures and document build for the combined internship report |

`nasa_mouse_generative` owns the common benchmark contract. The diffusion and
WGAN packages contain model-specific implementations because their training and
evaluation code differ substantially.

The complete dependency and data-flow explanation is in
[`docs/codebase_guide.md`](../docs/codebase_guide.md). Each package README lists
its modules and distinguishes final workflow code from supporting or historical
experiments.

Install all six packages from the repository root:

```bash
python -m pip install -e .
```

Use `python -m <package>.<module> --help` for module-specific options. The final
run sequence is documented in `outputs/COMMANDS.md`, while
`configs/generative/README.md` identifies the shared and model-specific config
directories.
