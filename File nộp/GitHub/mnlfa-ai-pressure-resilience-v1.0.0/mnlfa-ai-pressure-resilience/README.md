# Replication Archive: AI Integration Pressure, Pedagogical Autonomy, and Faculty Sustainable Performance

## Citation

> Tran Quang Canh (2025). AI Integration Pressure, Pedagogical Autonomy, and Faculty Sustainable Performance: The Role of Digital Resilience in Vietnamese Higher Education. *Computers & Education*.

- **Author**: Tran Quang Canh
- **Email**: canhtq@uef.edu.vn
- **ORCID**: https://orcid.org/0000-0001-6513-9319
- **Affiliation**: Faculty of Business Administration, Ho Chi Minh City University of Economics and Finance (UEF), Vietnam

---

## Overview

This repository contains the anonymized dataset, R analysis scripts, and supplementary output files for the above paper. The study examined whether AI integration pressure (AIP) and pedagogical autonomy (PA) exhibit quadratic effects on faculty sustainable performance (FSP), and whether digital resilience (DR) moderates these relationships as a continuous boundary condition, using a two-stage factor score regression approach within the Moderated Nonlinear Factor Analysis (MNLFA) framework.

**Sample**: N = 522 full-time faculty members across five Vietnamese universities (data collected March–June 2025).

---

## Repository Structure

```
mnlfa-ai-pressure-resilience/
│
├── README.md                        # This file
├── LICENSE                          # MIT License
│
├── data/
│   └── mnlfa_data.csv               # Anonymized survey dataset (N = 522)
│
├── code/
│   ├── MNLFA_lavaan.R               # Main analysis: CFA, DIF, structural models,
│   │                                #   bootstrap, robustness checks
│   └── MNLFA_Supplementary.R        # Supplementary analyses: descriptives,
│                                    #   VIF, Harman CMB, simple slopes,
│                                    #   Fornell-Larcker, robustness comparison
│
└── results/
    └── [Excel output files]         # Generated automatically when scripts are run
```

---

## Data Description

`data/mnlfa_data.csv` contains item-level responses for 20 survey items across four constructs:

| Variable | Items | Description |
|----------|-------|-------------|
| AIP1–AIP5 | 5 items | AI Integration Pressure |
| PA1–PA5 | 5 items | Pedagogical Autonomy |
| DR1–DR5 | 5 items | Digital Resilience |
| FSP1–FSP5 | 5 items | Faculty Sustainable Performance |

All items used a five-point Likert scale (1 = Strongly Disagree, 5 = Strongly Agree). No personally identifiable information is included. Cases with identical response patterns across all 20 items were removed prior to archiving (final N = 522).

---

## Requirements

### R version
R 4.4.1 or higher

### Required R packages
```r
install.packages(c(
  "lavaan",      # CFA and SEM
  "semTools",    # Reliability and validity indices
  "writexl",     # Export to Excel
  "readxl",      # Read Excel files
  "mgcv",        # GAM smoothers
  "robustbase"   # MM-estimator robust regression
))
```

---

## How to Reproduce

### Main analysis
```r
Rscript --vanilla code/MNLFA_lavaan.R \
  data/mnlfa_data.csv \
  results/ \
  2000 42
```

Arguments:
1. Path to input CSV
2. Path to output directory
3. Number of bootstrap replications (recommended: 2000)
4. Random seed (recommended: 42)

### Supplementary analysis
```r
Rscript --vanilla code/MNLFA_Supplementary.R \
  data/mnlfa_data.csv \
  results/
```

### Expected runtime
- Main analysis: approximately 2–5 minutes (bootstrap with 2,000 replications)
- Supplementary analysis: approximately 1 minute

---

## Key Results

| Hypothesis | Description | Result |
|------------|-------------|--------|
| H1 | AIP² → FSP (inverted-U) | Not supported (*b* = −0.006, *p* = .729) |
| H2 | PA² → FSP (concave) | Not supported (*b* = +0.002, *p* = .899) |
| H3 | AIP² × DR moderates AIP–FSP curve | **Supported** (*b* = −0.060, *p* < .001) |
| H4 | PA² × DR moderates PA–FSP curve | **Supported** (*b* = −0.045, *p* < .001) |
| H5 | DR → FSP (direct positive effect) | **Supported** (*b* = +0.521, *p* < .001) |

CFA fit: χ²(164) = 668.762, CFI = 0.930, TLI = 0.919, RMSEA = 0.077 [0.071, 0.083], SRMR = 0.062

Model C: R² = 0.741, Adj-R² = 0.736, N = 522

---

## License

This repository is released under the [MIT License](LICENSE). You are free to use, modify, and distribute the code with attribution.

---

## Acknowledgements

Data collection was conducted at five Vietnamese universities. The author thanks all faculty participants for their time and responses.
