# Public reproducibility analysis

This package reproduces the principal statistical analyses reported in the
manuscript on complete blood count Mahalanobis distance and survival after
definitive concurrent chemoradiotherapy.

## Folder structure

Place the files as follows:

```text
public_analysis.R
README.md
data/
  analysis_dataset.csv
  selected_CBC_long.csv
  reference_mean.csv
  reference_covariance.csv
  data_dictionary_eng.xlsx
```

## Running the analysis

Open `public_analysis.R` in RStudio and run the complete file, or run:

```r
source("public_analysis.R")
```

The script requires R and the `survival` package. It creates plain CSV outputs
under `results/public_analysis/`. It does not create figures, HTML files, PDFs,
or RDS model objects.

## Analyses reproduced

- Mahalanobis distance calculated from Hb, MCV, RDW, ANC, ALC, AMC, PLT, and MPV
- Reverse Kaplan-Meier median follow-up
- Paired changes from baseline to treatment weeks 2-4
- Continuous adjusted Cox and proportional-hazards analyses
- Baseline time-varying effects using `D x log(t/12)`
- Added value of week-specific D beyond baseline D
- Time-dependent cumulative/dynamic OS ROC analyses at 12-60 months
- Paired bootstrap AUC comparisons and baseline cutoff stability
- Kaplan-Meier log-rank and adjusted cutoff-based models
- Minimum-Hb pretreatment CBC sensitivity analysis

All survival times are measured from treatment completion. Bootstrap analyses
use 1,000 resamples and fixed random seeds to reproduce the reported results.

The legacy variables named `baseline_NLR_gt5` and `minHb_NLR_gt5` in
`analysis_dataset.csv` represent an NLR cutoff of 5. They are not used by the
public script. The reported comparison at NLR >= 3 is recalculated directly
from the continuous baseline NLR.
