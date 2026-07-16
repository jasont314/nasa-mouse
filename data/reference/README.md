# Mouse Reference Annotation

`gencode_vM39_mouse_gene_lengths.tsv` contains versionless mouse Ensembl gene
lengths calculated as the union of all annotated exon intervals in the GENCODE
M39 primary-assembly GTF. The generative benchmark uses these lengths for TPM
normalization in the Lacan diffusion paper-native preprocessing arm.

The adjacent manifest records the official source URL, source SHA-256, generation
time, length definition, and number of genes. Regenerate the table with:

```bash
PYTHONPATH=src python -m nasa_mouse_generative gene-lengths \
  --gtf assets/reference/gencode.vM39.primary_assembly.annotation.gtf.gz \
  --output data/reference/gencode_vM39_mouse_gene_lengths.tsv \
  --source-url https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M39/gencode.vM39.primary_assembly.annotation.gtf.gz
```
