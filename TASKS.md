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
- [ ] Create `environment.yml` (conda environment spec with all Python/R dependencies)
- [ ] Create `LICENSE` (MIT)
- [x] Update `CLAUDE.md` with session goal

### Phase 3 — Documentation

- [ ] Create root `README.md` (project overview, pipeline summary, key results, team credits)
- [ ] Create `pipeline/README.md` (HPC pipeline instructions)
- [ ] Create `analysis/README.md` (local analysis instructions, notebook descriptions)
- [ ] Create `files/README.md` (data download instructions, ENA accession)
- [ ] Create `gene_matrix/README.md` (how to regenerate count matrices)
- [ ] Create `qc_reports/README.md` (QC summary)

### Phase 4 — Git Init & Final Review

- [x] Initialize git repository (`git init`, initial commit — 118 files)
- [x] Review full repo structure — verify nothing broken, all paths logical
- [x] Check `.gitignore` is correctly excluding large files
- [x] Verify all new markdown files render correctly (no broken links)
- [x] Iterate corrections (fixed CLAUDE.md pipeline/ paths, STAR version notes)

---

## Key Findings (for README)

- **Data:** ENA accession PRJEB15401; 142 beta-cell samples (97 healthy, 45 T2D)
- **Genes analyzed:** 62,710 raw → 16,361 after filtering (present in ≥5 cells)
- **Model:** Logistic regression (L1 penalty, liblinear) on log2-RPKM features; top predictor: IRF2BPL (slide 18 of presentation)
- **Key finding:** Significant expression differences in INS, GCG confirmed; GSEA identified T2D mellitus and insulin resistance pathways enriched
- **11123248.json:** Zenodo metadata file for ENA dataset — moved to `_archive/`
- **manipulate_beta_cells.xlsx:** Sample filter spreadsheet (created by teammate); data already extracted to `files/healthy_beta.txt` and `files/t2d_beta.txt`
- **No explicit model accuracy metric** recorded in notebooks; feature importances shown in `analysis/analysis_1.ipynb`

## Team

- Richard Goodier *(repo author)*
- Group 3 contributors *(names pending permission — see `submission/` folder for full list)*

---

*Last updated: 2026-04-02*
