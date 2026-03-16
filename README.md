# A twin-based analysis of proteomic organ aging

DOI: `TBD`

This folder contains the analysis scripts used for the paper.

## Scripts

- `clock_developer.Rmd`
  - Builds organ-age clocks from proteomic data.
  - Produces fitted models and prediction outputs used in downstream analyses.

- `clock_functions.R`
  - Helper functions used by `clock_developer.Rmd`.
  - Includes shared modeling and utility code.

- `baseline_analyses.Rmd`
  - Baseline cross-sectional analyses of organ-age measures and outcomes.
  - Includes main baseline statistical models and summary figures/tables.

- `incident_analyses.Rmd`
  - Longitudinal and incident-event analyses.
  - Includes time-to-event and follow-up based models.

- `multimorbidity_analyses.Rmd`
  - Analyses of multimorbidity burden in relation to organ aging.
  - Includes count/survival analyses and corresponding figures/tables.

## Notes

- Set input and output paths in each script before running.
- Data files are not included in this repository.
