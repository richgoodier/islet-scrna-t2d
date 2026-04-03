# GitHub Repo Reorganization — Task Tracker

**Session:** 2026-04-02
**Goal:** Transform this bioinformatics course project into a professional, fully reproducible GitHub repository suitable for a data science job portfolio.

**Constraints:**

- Never delete or alter existing project file contents
- Move out-of-place files to `_archive/` (gitignored), not delete them
- Use subagents for parallel research tasks where possible
- After each task: review and correct before moving on
- Final step: iterate corrections until none remain

---

## Key Findings (for README)

- **Data:** ENA accession PRJEB15401; 142 beta-cell samples (97 healthy, 45 T2D)
- **Genes analyzed:** 62,710 raw → 16,361 after filtering (present in ≥5 cells)
- **Model:** Logistic regression (L1 penalty, liblinear) on log2-RPKM features; top predictor by coefficient: MICOS10P3 (note: README previously stated IRF2BPL — verify against current notebook output)
- **Key finding:** Significant expression differences in INS, GCG confirmed; GSEA identified T2D mellitus and insulin resistance pathways enriched
- **11123248.json:** Zenodo metadata file for ENA dataset — moved to `_archive/`
- **manipulate_beta_cells.xlsx:** Sample filter spreadsheet (created by teammate); data already extracted to `files/healthy_beta.txt` and `files/t2d_beta.txt`

## Team

- Richard Goodier *(repo author)*
- Group 3: Dinesh Sambhaji Pradhan, Asmitha Nagajothi Purushotam, Vedant Kulkarni, Richard Goodier

---

## Task Breakdown

### Phase 1 — Structure & Organization

- [x] Create `pipeline/` folder; move all `.nf` scripts and `build_star_index.sh` into it
- [x] Move `shell_scripts/` contents under `pipeline/shell_scripts/`
- [x] Rename `paper, figures, tables/` → `references/`
- [x] Rename `Individual Write Up/` → `writeup/`
- [x] Move `misc/bbab563.pdf` → `references/`
- [x] Create `_archive/raw_fastq/` and move both root-level FASTQ files there
- [x] Create `_archive/` and move `manipulate_beta_cells.xlsx`, `11123248.json` there
- [x] Create `_archive/generated/` and move `multiqc.pdf`, `pipeline.pdf` there
- [x] Create `docs/` and move `Project - BINF 6310.pdf/.docx`, `notes.txt`, `set up env.txt` there

### Phase 2 — GitHub Essentials

- [x] Create `.gitignore` (excludes `_archive/`, large matrices, FASTQ, BAM, OS files)
- [x] Create `environment.yml` (conda environment spec with all Python/R dependencies)
- [x] Create `LICENSE` (MIT)
- [x] Update `CLAUDE.md` with session goal

### Phase 3 — Documentation

- [x] Create root `README.md` (project overview, pipeline summary, key results, team credits)
- [x] Create `pipeline/README.md` (HPC pipeline instructions)
- [x] Create `analysis/README.md` (local analysis instructions, notebook descriptions)
- [x] Create `files/README.md` (data download instructions, ENA accession)
- [x] Create `gene_matrix/README.md` (how to regenerate count matrices)
- [x] Create `qc_reports/README.md` (QC summary)

### Phase 4 — Git Init & Final Review

- [x] Initialize git repository (`git init`, initial commit — 118 files)
- [x] Review full repo structure — verify nothing broken, all paths logical
- [x] Check `.gitignore` is correctly excluding large files
- [x] Verify all new markdown files render correctly (no broken links)
- [x] Iterate corrections (fixed CLAUDE.md pipeline/ paths, STAR version notes)

### Phase 5 — Streamlining & Professionalization

<!-- presentation/ -->

- [x] Archive `presentation/speaker_notes.txt` — internal prep notes, not portfolio-appropriate
- [x] Archive `presentation/fig3_1.jpg`, `fig3_1_legend.jpg`, `fig6.jpg`, `fig8.jpg` — copyrighted figures from source paper, already captured in the PDF

<!-- qc_reports/ -->

- [x] Archive entire `qc_reports/multiqc_data/` directory — HTML report is self-contained; raw backing files not used downstream

<!-- rpkm_values/ -->

- [x] Convert `rpkm_values/calculate_rpkm.txt` → `calculate_rpkm.R` — fix extension, replace hardcoded `/scratch/<username>/` path, reformat as proper R script

<!-- analysis/ -->

- [x] Rename `analysis/analysis.ipynb` → `expression_analysis.ipynb`
- [x] Rename `analysis/analysis_1.ipynb` → `ml_classification.ipynb`
- [x] Update `analysis/README.md` to reflect new notebook names

<!-- references/ -->

- [x] Verify `references/data_paper/*.xlsx` files are unused as analysis input (confirmed unused)
- [x] Archive all PDFs and xlsx files in `references/`
- [x] Add References section to root `README.md` with full citations and DOI links
- [x] Remove empty `references/` directory once cleared

<!-- housekeeping -->

- [x] Fix stale Team section in TASKS.md
- [x] Final review of all READMEs for broken links after renames/removals (fixed gene_matrix/README.md, CLAUDE.md)

### Phase 6 — Future Tasks

- [x] Create `analysis/ANALYSIS_REPORT.md` — distills findings with embedded figures
- [x] Comprehensive review and rename of all file names for clarity (20 files renamed; all references updated; documented in file_rename_review.md)
- [x] Review each `.ipynb` notebook and clean up unnecessary code blocks

#### Notebook Cleanup Checklist

**1. `pipeline.ipynb`** ✓

- [x] Cell 0: Fix typo "Pipline" → "Pipeline"
- [x] Cell 4: Delete stray terminal output (`python-3.12.9` / `conda-25.3.0`)
- [x] Cell 12: Remove commented-out `conda install cutadapt` line
- [x] Cell 13: Rewritten — featureCounts/subread install note
- [x] Cell 20: Replace hardcoded course path with placeholder
- [x] Cell 26: Remove inline STAR index command
- [x] Cell 28: Replace email and course paths with placeholders
- [x] Cell 31: Replace `/scratch/<username>/...` paths with `<username>` placeholders
- [x] Cell 32: Delete stray Nextflow version output
- [x] Cell 35: Delete duplicate compute-node access instructions
- [x] Cell 37: Replace `/scratch/<username>/...` (teammate path) with `<username>` placeholder
- [x] Cell 39: Replace absolute featureCounts.nf path
- [x] Cell 41: Replace hardcoded path → `../gene_matrix/counts_combined.txt`
- [x] Cell 42: Delete empty cell

**2. `files/create_download_scripts.ipynb`** ✓

- [x] Cell 11: Replace hardcoded course path and `explorer` hostname with `<username>@<hpc-hostname>` placeholders
- [x] Cell 12: Raw bash commands in code cell → converted to markdown block with placeholders
- [x] Cell 14: `nano` editor call + hardcoded scp paths → converted to markdown
- [x] Cell 19: Bash script in code cell → converted to markdown

**3. `files/donor_samples.ipynb`** ✓

- [x] Cell 15: Removed commented-out dead code line (`#links = chosen_samples["submitted_ftp"].to_numpy()`)
- [x] Cell 19: Deleted empty cell

**4. `files/read_counts.ipynb`** ✓

- [x] Cell 5: Replaced stub `# Plot the reads` comment with a markdown note explaining the visualization context

**5. `gene_matrix/combine_tables.ipynb`** ✓

- [x] Cell 6: Removed dead first `pd.concat(...)` that was immediately overwritten
- [x] Cell 15: Removed orphaned debug cell (displayed only `Length` column)
- [x] Cell 16: Deleted empty cell

**6. `qc_reports/qc_analysis.ipynb`** ✓

- [x] Cell 6: Deleted empty cell
- [x] Cell 7: Renamed variable `bad_files` → `flagged_samples`
- [x] Cell 8: Deleted empty cell
- [x] Cell 9: Deleted empty cell

**7. `analysis/ensembl_to_ids.ipynb`** ✓

- [x] Cell 9: Removed commented-out save operation (`# df_out_with_condition.to_csv(...)`) — superseded by Cell 11

**8. `analysis/differential_expression.ipynb`** ✓

- [x] No critical issues — notebook is clean. Inline comments are documentation, not dead code.

**9. `analysis/expression_analysis.ipynb`** ✓

- [x] Cell 3: Removed `# df.columns.tolist()` — leftover debug line

**10. `analysis/ml_classification.ipynb`** ✓

- [x] Cell 23: Deleted empty cell
- [x] Cell 29: Deleted empty cell

### Phase 7 — Directory Cleanup & Organization

- [x] Create `analysis/figures/` subfolder; move all 10 output PNGs there — keeps notebooks and reports distinct from figure outputs
- [x] Update `analysis/ANALYSIS_REPORT.md` — all 7 image references updated to `figures/` prefix
- [x] Update `analysis/README.md` — figures section updated to `figures/` prefix; table separator spacing fixed (MD060)
- [x] Rename `files/` → `samples/` — more informative name reflecting content (sample metadata + download scripts)
- [x] Update `README.md` — directory tree and figure link updated
- [x] Update `CLAUDE.md` — all `files/` references updated to `samples/`
- [x] Update `pipeline/README.md` — cross-link to `../samples/README.md` updated
- [ ] No subfolders created within `samples/` — 14 files are navigable with README; moving metadata CSVs/TXTs would require updating notebook pd.read_csv paths

### Phase 8 — Personal Information Scrub

Full scan of all committable files. Findings and resolutions:

- [x] `pipeline/build_star_index.sh` — replaced `<your-email>` placeholder (email was `goodier.r@northeastern.edu`), replaced hardcoded HPC log/error/mkdir paths with `./genome_index/logs/`; updated conda env name
- [x] `pipeline/star_align.nf` — replaced 5 instances of `/scratch/<username>/` (was `goodier.r`) with `/scratch/<username>/` placeholder
- [x] `samples/count_reads.py` — replaced full WSL path (`\\wsl.localhost\Ubuntu\home\richgoodier\...`) with `/path/to/fastq/files`
- [x] `samples/create_download_scripts.ipynb` — removed collaborator first name ("Vedant") from cell 0; replaced with neutral "The sample list was filtered..."
- [x] `pipeline.ipynb` — replaced `explorer` HPC hostname in cell 19 with `<hpc-hostname>` placeholder
- [x] `gene_matrix/featureCounts_healthy/counts.txt` — replaced teammate username in featureCounts command header comment
- [x] `gene_matrix/featureCounts_t2d/counts.txt` — same
- [x] `qc_reports/multiqc_report.html` — replaced `/scratch/<username>/` in embedded analysis path
- [x] `TASKS.md` — removed actual usernames from task description text
- [x] Kept intentional: `LICENSE` (copyright holder name), `README.md` Acknowledgments (standard attribution), `TASKS.md` Team section (appropriate for portfolio)
- [x] Noted: `presentation/presentation_slides.pdf` is binary — cannot edit without source file; names likely appear inside

### Last updated: 2026-04-02
