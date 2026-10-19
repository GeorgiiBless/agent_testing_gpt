# Work Summary and Continuation Handoff

## Purpose

This folder contains an in-progress Python workflow for Gaia open-cluster
analysis. The intended scientific project is to identify members of the Pleiades,
Praesepe, and M67, create Gaia colour-magnitude diagrams, and estimate cluster
distances.

This document is a handoff for another model or developer. It describes what
actually exists and has been run, not only what was originally planned.

## Workspace

Project root:

`.`

The folder is not currently a Git repository. Do not assume that commits,
branches, or remote history are available.

Important existing files and directories:

- `Cluster_project.py` - the implemented single-file command-line pipeline.
- `plan.md` - the broader original design and requirements.
- `requirements.txt` - Python dependencies.
- `data_vizier_m67/` - cached M67 query and raw catalogue.
- `outputs_vizier_m67/` - generated M67 tables, plots, metadata, and report.
- `outputs/`, `outputs_gavo_smoke/`, `data_vizier/` - present in the directory
  listing but not inspected in detail during the prior summary; inspect before
  modifying or deleting anything.
- `.venv/` - local Python virtual environment.

## Implemented code

The main script is:

`Cluster_project.py`

It defines configurations for:

- Pleiades: RA 56.75 deg, Dec 24.1167 deg, radius 2 deg.
- Praesepe: RA 130.05 deg, Dec 19.67 deg, radius 2 deg.
- M67: RA 132.825 deg, Dec 11.8 deg, radius 2 deg.

The script currently:

1. Queries or loads cached data from Gaia DR3 table
   `gaiadr3.gaia_source_lite`.
2. Selects a documented set of Gaia source-lite astrometry, photometry, and
   quality fields, with a 2-degree cone and `phot_g_mean_mag <= 18` ADQL
   predicate.
3. Uses the Gaia column names directly as normalized names such as
   `source_id`, `parallax`, `pmra`, `pmdec`, and Gaia photometry names.
4. Applies finite-value and quality filters.
5. Computes `bp_rp`.
6. Applies a fixed parallax correction:
   `parallax_corrected = parallax + 0.017` mas.
7. Scales corrected parallax, `pmra`, and `pmdec` with
   `sklearn.preprocessing.RobustScaler`.
8. Uses DBSCAN for astrometric clustering.
9. Chooses the component closest to the configured expected cluster locus,
   with sky-position and component-size terms in the heuristic cost.
10. Assigns a continuous Gaussian-like `membership_probability` and retains
    only DBSCAN component sources with probability at least 0.5 in
    `member_flag`.
11. Estimates a weighted cluster parallax and distance modulus.
12. Computes simple extinction-corrected HR-diagram quantities:
    `abs_g = phot_g_mean_mag - distance_modulus - 0.86 * A_V`
    and `bp_rp_0 = bp_rp - 0.42 * A_V`.
13. Writes astrometric diagnostic, HR-diagram, and colour--magnitude PNG
    files.
14. Computes a simple reproducible turn-off proxy from the brightest selected
    stars in a colour range.
15. Writes a member CSV, metadata JSON, a summary CSV, and a short Markdown
    report.

The command-line interface can be inspected with:

```bash
python Cluster_project.py --help
```

Important options include:

- `--clusters`
- `--cache-dir`
- `--output-dir`
- `--refresh`
- `--offline`
- `--dbscan-eps`
- `--min-samples`
- `--verbose`

## Actual completed run

Only M67 has been processed in the existing outputs. The query metadata is:

`data_vizier_m67/m67_query.json`

The raw cached catalogue is:

`data_vizier_m67/m67_raw.ecsv`

Query details:

- Source: CDS VizieR.
- Catalogue: `I/345/gaia2`.
- Centre: RA 132.825 deg, Dec 11.8 deg.
- Cone radius: 2 deg.
- Query version: `1`.
- Retrieval timestamp: `2026-10-06T09:37:00.986568+00:00`.
- Rows returned: 71,390.

Generated outputs:

- `outputs_vizier_m67/m67_members.csv`
- `outputs_vizier_m67/m67_diagnostics.png`
- `outputs_vizier_m67/m67_hr_diagram.png`
- `outputs_vizier_m67/m67_metadata.json`
- `outputs_vizier_m67/cluster_summary.csv`
- `outputs_vizier_m67/report.md`

Results recorded in `outputs_vizier_m67/cluster_summary.csv`:

| Quantity | Value |
|---|---:|
| Queried sources | 71,390 |
| Quality-filtered sources | 18,645 |
| Selected members | 12,989 |
| Corrected parallax | 0.991172 mas |
| Parallax error floor/result | 0.015 mas |
| Distance modulus | 10.019254 |
| Extinction prior (`A_V`) | 0.13 |
| Metallicity prior (`[Fe/H]`) | 0.03 |
| Turn-off absolute-G proxy | -1.317354 |
| Turn-off BP-RP proxy | 1.148300 |

The metadata file confirms:

- Membership label selected: `0`.
- Median selected-member `pmra`: `-4.248`.
- Median selected-member `pmdec`: `-2.997`.
- Median selected-member corrected parallax: `0.9718`.
No isochrone fit or age estimate is part of the current pipeline; the outputs
are based on Gaia astrometry, DBSCAN membership, and distance estimates.

## Dependency state

`requirements.txt` contains:

- `astroquery>=0.4.7`
- `astropy>=6.0`
- `matplotlib>=3.7`
- `numpy>=1.24`
- `pandas>=2.0`
- `scikit-learn>=1.3`
- `scipy>=1.10`
- `tabulate>=0.9`

Before running Python commands, prefer the workspace's configured interpreter
or the existing `.venv`. Do not install packages unless a validation command
shows that a required dependency is missing or the dependency manifest is
changed.

## Gaia archive availability check

Verified on 2026-10-13 using the workspace interpreter
`.venv/bin/python`:

- `astroquery` version: `0.4.11`.
- `astroquery.gaia` imports successfully.
- `Gaia.load_table("gaiadr3.gaia_source_lite")` succeeded.
- The table metadata reported 51 columns, including `source_id`, `ra`, `dec`,
  `parallax`, `parallax_error`, `pmra`, `pmdec`, `ruwe`, and Gaia photometry
  fields.
- A live asynchronous query succeeded:
  `SELECT TOP 1 source_id, ra, dec, parallax FROM gaiadr3.gaia_source_lite`.
- The query returned one row with source ID
  `137333365999306624`.

Therefore, the Gaia archive and the requested `gaiadr3.gaia_source_lite` table
were reachable at the time of this check. The archive also reported a planned
deployment outage on 2026-10-14 from 11:00 to 12:00 CEST.

## Continuation update

The query path now uses the verified Gaia DR3 `gaiadr3.gaia_source_lite`
table through `astroquery.gaia`. Every cluster uses a 2-degree cone and the
ADQL query limits results to `phot_g_mean_mag <= 18`. The exact selected
columns and normalized ADQL are recorded in each query metadata JSON.

Existing VizieR/DR2 caches are treated as stale because their catalog,
query-version, and selection metadata do not match the Gaia DR3 query. They
are not silently reused.

## Planned but not implemented

`plan.md` describes a substantially more complete design. The following
planned pieces do not currently exist or are not complete:

- `pyproject.toml`, package layout, README, and separate source modules.
- Configuration file for clusters and scientific thresholds.
- Gaia DR3 table schema validation before querying.
- Explicit ADQL provenance, checksums, archive retries, and query validation.
- Versioned Gaia DR3 parallax zero-point correction.
- Full Gaia quality metadata beyond the fields available in
  `gaia_source_lite`, such as visibility/transit counts.
- HDBSCAN preference and measurement-error-aware Mahalanobis membership.
- External membership validation.
- Rigorous extinction coefficients and uncertainty propagation.
- Isochrone-grid likelihood/posterior output and credible intervals.
- Binary modelling, bootstrap uncertainty, sensitivity analysis, and
  statistical/systematic uncertainty separation.
- Three-cluster end-to-end validation.
- Unit tests and synthetic integration tests.
- Staged CLI commands such as `query`, `members`, `plot`, `fit`, and `report`.

## Important scientific caveats

1. The query now uses Gaia DR3 `gaiadr3.gaia_source_lite`, but the code still
   uses a fixed parallax correction rather than a versioned Gaia DR3
   zero-point implementation.
2. The parallax correction is a fixed constant and should not be treated as a
   current Gaia DR3 zero-point implementation.
3. Membership is heuristic DBSCAN clustering in three astrometric dimensions.
4. The selected M67 component has a median proper motion notably different from
   the configured expected `pmra` (`-4.248` versus `-10.9` mas/yr). This is a
   warning that the current component-selection result needs scientific
   validation before being trusted.
5. The selected sample is very large (12,989 of 18,645 quality-filtered
   sources), which may indicate field contamination or an overly broad
   clustering criterion.
6. The turn-off values are proxies, not a formal turn-off fit.
7. No physical cluster age has been estimated in the saved run.
8. The existing report is only a compact summary table; it is not the full
   report envisioned in `plan.md`.

## Recommended continuation order

### First: reproduce and inspect safely

1. Inspect the unreviewed directories `outputs/`, `outputs_gavo_smoke/`, and
   `data_vizier/`.
2. Run the script in offline mode against the existing M67 cache, writing to a
   new temporary output directory rather than overwriting the saved results.
3. Compare the regenerated summary and member counts with the recorded values.
4. Inspect the M67 diagnostic plot and astrometric distributions before making
   scientific claims.

### Second: fix core correctness issues

1. Investigate why the selected M67 astrometric component does not match the
   configured expected proper motion.
2. Confirm catalogue column semantics and whether the VizieR data really match
   the assumptions in the code.
3. Improve component selection and add explicit diagnostics for all DBSCAN
   components.
4. Preserve raw and derived columns while making quality thresholds explicit.
5. Add tests for cache loading, quality filtering, component selection,
   distance calculations, HR quantities, and turn-off extraction.

### Third: add a controlled age fit

1. Obtain or provide a documented Gaia-passband isochrone CSV with the four
   required columns.
2. Run the optional isochrone fit on a validated member sample.
3. Save the best-fit parameters and model metadata separately from the
   no-fit baseline.
4. Do not report an age until the fit, distance, membership, and model
   assumptions have been checked.

### Fourth: expand scope

1. Process Pleiades and Praesepe only after the M67 pipeline is validated.
2. Add reproducibility metadata, uncertainty estimates, and sensitivity tests.
3. Refactor into the modular structure described by `plan.md` only when the
   single-file behavior is covered by tests.

## Safe working rules

- Do not delete or overwrite existing cached data or outputs without explicit
  user approval.
- Use new output directories for experiments.
- Treat `plan.md` as a specification, not evidence that a feature exists.
- Treat `outputs_vizier_m67/report.md` and
  `outputs_vizier_m67/cluster_summary.csv` as the authoritative record of the
  saved M67 run.
- If changing scientific assumptions, update this handoff and the related
  documentation.
