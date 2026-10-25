# Body Adiposity Index and Dysglycemia in Hausa Women

[![Code license: MIT](https://img.shields.io/badge/Code%20license-MIT-blue.svg)](LICENSE) [![Data license: CC BY 4.0](https://img.shields.io/badge/Data%20license-CC%20BY%204.0-lightgrey.svg)](LICENSE-DATA)

Data and R code for the article:

> Asuku AY, Bashir MA, Muhammed J. **Body Adiposity Index Outperforms Waist Circumference and Waist-to-Height Ratio in Identifying Dysglycemia Among Hausa Women in Northern Nigeria.** *Scientific African* (under review).

The code reproduces every table, figure and number reported in the manuscript and its supplementary material from the de-identified dataset.

## Study summary

A cross-sectional study of 249 Hausa women aged 35 years and older, without known diabetes, recruited at two teaching hospitals in Kano, Nigeria (119 premenopausal, 130 postmenopausal). Seven anthropometric indices were compared for their ability to identify dysglycemia (HbA1c ≥5.7%):

| Abbreviation | Index |
|---|---|
| BAI | Body adiposity index = hip circumference (cm) / height (m)^1.5 − 18 |
| WHtR | Waist-to-height ratio |
| WC | Waist circumference |
| HC | Hip circumference |
| NC | Neck circumference |
| BMI | Body mass index |
| WHR | Waist-to-hip ratio |

**Main findings.** With age and menopausal status in the model, BAI discriminated dysglycemia better than WHtR and WC (AUC 0.882 vs 0.860 and 0.849). Waist-based indices added no information once BAI was in the model. BAI was not superior to WHtR for diabetes-range HbA1c (≥6.5%). Menopausal status was associated with higher HbA1c independently of age and BAI.

## Repository contents

| File | Description |
|---|---|
| `clean_data.csv` | De-identified dataset (one row per participant) |
| `Analysis_script.R` | All analyses; creates the tables, figures and inline values |
| `Article.Rmd` | Manuscript; sources `Analysis_script.R` and knits to Word |
| `Supplementary.Rmd` | Supplementary Tables S1–S3 and Figures S1–S2 |
| `mybibfile.bib` | References cited in the original version of the manuscript |
| `new_references.bib` | Methods references added in the revision |
| `references.docx` | for formatting the article |
| `LICENSE` | MIT License for the code |
| `LICENSE-DATA` | CC BY 4.0 license for the data |
| `scientific-african.csl` | Citation style (Elsevier numeric with titles, the parent of the official Scientific African CSL) |
| `figures/` | Created when the script runs: figures as PDF and 300-dpi TIFF |

## Data dictionary (`clean_data.csv`)

One row per participant (n = 249), no missing values.

| Column | Description | Units |
|---|---|---|
| `id` | Sequential participant number (not an identifier) | — |
| `age_years` | Age | years |
| `menopausal_status` | `Premenopausal` or `Postmenopausal` (≥12 consecutive months of amenorrhea) | — |
| `age_at_menarche_years` | Age at menarche | years |
| `height_cm` | Height | cm |
| `weight_kg` | Weight | kg |
| `waist_circumference_cm` | Waist circumference, midway between lowest rib and iliac crest | cm |
| `hip_circumference_cm` | Hip circumference, widest point over the greater trochanters | cm |
| `neck_circumference_cm` | Neck circumference, just below the thyroid cartilage | cm |
| `bmi_kg_m2` | Body mass index = weight / height² (2 d.p.) | kg/m² |
| `whr` | Waist-to-hip ratio = WC / HC (3 d.p.) | ratio |
| `whtr` | Waist-to-height ratio = WC / height (3 d.p.) | ratio |
| `bai_pct` | Body adiposity index = HC (cm) / height (m)^1.5 − 18 (2 d.p.) | % |
| `hba1c_pct` | Glycated hemoglobin, ion-exchange HPLC | % |

## Analyses

`Analysis_script.R` runs the following, in order:

1. **Descriptives** by menopausal status (Welch's t-test, χ² test) — Table 1
2. **Correlations** with HbA1c, with Meng–Rosenthal–Rubin tests of each index against BAI and Fisher z tests between groups — Supplementary Table S1
3. **Discrimination of dysglycemia**: logistic models (age + status + index), DeLong comparisons of AUC, NRI and IDI, bootstrap optimism correction, index-alone AUCs and Youden cut-offs in premenopausal women — Table 2, Figure 1
4. **Incremental information** for continuous HbA1c: ΔR² over age and status with bootstrap CIs, information added beyond BAI, index × status interactions — Table 3, Supplementary Figure S1
5. **Menopausal status, age and HbA1c**: linear and spline age adjustment, adjustment for adiposity, age-overlap (45–55 years) sensitivity analysis — Table 4, Figure 2
6. **Model comparison**: linear regression, LASSO, random forest and SVR against a three-variable linear model, using 10 × 5-fold repeated cross-validation with Nadeau–Bengio corrected tests — Supplementary Table S2
7. **LASSO bootstrap selection stability** and random forest permutation importance — Supplementary Table S3
8. **Decision-curve analysis** for dysglycemia — Supplementary Figure S2

## How to reproduce

### Requirements

- R ≥ 4.3
- R packages:

```r
install.packages(c("tidyverse", "pROC", "glmnet", "ranger", "kernlab",
                   "patchwork", "flextable", "officer", "rmarkdown", "bookdown"))
```

The results were produced with tidyverse 2.0.0, pROC 1.18.5, glmnet 4.1.8, ranger 0.16.0, kernlab 0.9.32 and patchwork 1.2.0. Other versions should give the same estimates, with possible small differences in bootstrap and cross-validation results.

### Run the analysis only

```r
source("Analysis_script.R")
```

This takes about one minute. It creates the result objects (`table1`–`table5`, `tableS1`, `tableS2`, `auc_dys`, `inc_res`, `status_res` and others) and writes the figures to `figures/`.

### Rebuild the supplementary material

```r
rmarkdown::render("Supplementary.Rmd")
```

### Rebuild the manuscript

`Article.Rmd` reads its references from `mybibfile.bib` and `new_references.bib` and formats them with `scientific-african.csl`, so no reference manager is needed. It uses `jas_template.docx` as the Word reference template; if you don't have that file, delete the `reference_docx:` line from the YAML header and Word's default styles will be used. Then run:

```r
rmarkdown::render("Article.Rmd")
```

### Reproducibility

The random seed is fixed (`set.seed(2026)`) at the top of `Analysis_script.R`. The number of resamples can be changed there:

| Setting | Default | Used for |
|---|---|---|
| `N_BOOT` | 1000 | Bootstrap CIs for ΔR² |
| `N_OPT` | 300 | Optimism correction of AUC |
| `N_LASSO` | 500 | LASSO selection stability |
| `CV_FOLDS` × `CV_REPEATS` | 5 × 10 | Repeated cross-validation |

## Ethics

The study was approved by the research ethics committees of Murtala Muhammad Teaching Hospital (MMSHZ//GEN/0831/I), Muhammad Abdullahi Wase Teaching Hospital (MAWTH/HMB/GEN/456/I) and the Kano State Ministry of Health (NHREC/17/03/2018). All participants gave written informed consent. The dataset contains no names, dates or other direct identifiers.

## Citation

If you use these data or code, please cite the article above. Citation details will be updated on publication.

## License

- **Code** (`Analysis_script.R`, `Article.Rmd`, `Supplementary.Rmd`): [MIT License](LICENSE)
- **Data** (`clean_data.csv`): [Creative Commons Attribution 4.0 International (CC BY 4.0)](LICENSE-DATA)
- **Citation style** (`scientific-african.csl`): third-party file from the [Citation Style Language project](https://citationstyles.org/), under its own CC BY-SA 3.0 license

You may reuse the code and data, including commercially, provided you credit the authors and cite the article.

## Contact

M. A. Bashir (corresponding author), Department of Anatomy, Bayero University Kano, Nigeria — musabashir34@gmail.com
