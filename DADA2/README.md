# DADA2 — Aurora PS137

16S rRNA amplicon analysis for cruise PS137 (Aurora). Standard DADA2 pipeline unless noted below.

---

## Non-default parameters

| Parameter | Value | Note |
|---|---|---|
| `TRUNC_LEN_F` | 280 | |
| `TRUNC_LEN_R` | 280 | |
| `MAX_EE_F` | 3.0 | |
| `MAX_EE_R` | 3.0 | |
| `MIN_OVERLAP` | 12 | |
| `MIN_ASV_LENGTH` | 390 | |
| `MAX_ASV_LENGTH` | 450 | |
| `POOLING` | TRUE | Chimera removal also pooled |
| `BINNED_QUALITY_SCORES` | AUTO | For NovaSeq/NextSeq binned scores |
| `MULTIPLE_RUNS` | FALSE | Single run |
| `N_THREADS` | 35 | |

**Taxonomy database:** SILVA 138.1 (`silva_nr99_v138.1_train_set`, species assignment via `silva_species_assignment_v138.1`)

---

## Output files

```
tables/
  seqtab_final_TRUE_SingleRun.csv   # ASV × sample abundance table
  taxonomy_TRUE_SingleRun.csv       # ASV taxonomy (Kingdom → Species)
  tracking_tech_per_run_TRUE.csv    # Read counts per step (input→filter→denoise→merge)
  tracking_final_TRUE_SingleRun.csv # Final chimera removal summary per sample

rds_objects/
  seqtab_final_TRUE_SingleRun.rds   # seqtab as R object
  taxonomy_TRUE_SingleRun.rds       # taxonomy as R object

qc_plots/
  quality_profile_fwd_SingleRun.pdf
  quality_profile_rev_SingleRun.pdf
  error_rates_TRUE_SingleRun_fwd.pdf

logs/
  config.txt            # Full run configuration
  ANALYSIS_REPORT_TRUE.txt
```
