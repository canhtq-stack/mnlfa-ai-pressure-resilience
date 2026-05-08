# Data Folder

Place `mnlfa_data.csv` in this folder before running the analysis scripts.

## File Description

`mnlfa_data.csv` — Anonymized survey dataset (N = 522)

### Variables

| Column | Type | Description |
|--------|------|-------------|
| AIP1 | integer (1–5) | AI Integration Pressure item 1 |
| AIP2 | integer (1–5) | AI Integration Pressure item 2 |
| AIP3 | integer (1–5) | AI Integration Pressure item 3 |
| AIP4 | integer (1–5) | AI Integration Pressure item 4 |
| AIP5 | integer (1–5) | AI Integration Pressure item 5 |
| PA1  | integer (1–5) | Pedagogical Autonomy item 1 |
| PA2  | integer (1–5) | Pedagogical Autonomy item 2 |
| PA3  | integer (1–5) | Pedagogical Autonomy item 3 |
| PA4  | integer (1–5) | Pedagogical Autonomy item 4 |
| PA5  | integer (1–5) | Pedagogical Autonomy item 5 |
| DR1  | integer (1–5) | Digital Resilience item 1 |
| DR2  | integer (1–5) | Digital Resilience item 2 |
| DR3  | integer (1–5) | Digital Resilience item 3 |
| DR4  | integer (1–5) | Digital Resilience item 4 |
| DR5  | integer (1–5) | Digital Resilience item 5 |
| FSP1 | integer (1–5) | Faculty Sustainable Performance item 1 |
| FSP2 | integer (1–5) | Faculty Sustainable Performance item 2 |
| FSP3 | integer (1–5) | Faculty Sustainable Performance item 3 |
| FSP4 | integer (1–5) | Faculty Sustainable Performance item 4 |
| FSP5 | integer (1–5) | Faculty Sustainable Performance item 5 |

### Preprocessing Applied
- Responses with fewer than one semester of teaching experience: excluded
- Responses with more than 20% missing items: excluded
- Responses with straight-lining within any single construct block: excluded
- Responses with identical answers across all 20 items: excluded
- Final analytic sample: N = 522, no missing data
