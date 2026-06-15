# NASA Mouse GLARE

This repository vendors GLARE in `src/glare` and adversarial gene expression
code in `src/adversarial-gene-expression`, while project-specific mouse
spaceflight code lives in `src/nasa_mouse_glare`.

Run workflow commands from this directory:

```bash
cd path/to/nasa-mouse
conda activate nasa
export PYTHONPATH=src
```

See [src/nasa_mouse_glare/README.md](src/nasa_mouse_glare/README.md) for the
TMS pretraining and OSDR fine-tuning workflow, including strict
GLARE-compatible `.mtx`/CSV exports.
