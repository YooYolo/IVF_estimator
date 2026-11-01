# IVF Estimator & Statistics Dashboard

A publicly accessible interactive **Streamlit** data dashboard for exploring In Vitro Fertilization (IVF) cycle outcomes, ovarian reserve markers, pairwise correlation statistics with multiple-testing adjustment, and personalized Nearest-Neighbor (k-NN) attrition funnel estimations.

---

## Overview

The **IVF Estimator** helps patients, researchers, and clinicians visualize real-world IVF statistics across both peer-reviewed clinical research and community-reported cycle data. Users can explore how baseline parameters (Age, AMH, FSH, stimulation dosages) relate to oocyte yield, fertilization, blastocyst development, and PGT-A euploidy rates, and estimate individual cycle funnels using the closest matching historical cycles.

---

## Core Dashboard Features

### 1. Flexible Dataset Display, Color Coding & Bilingual UI (English / 中文)
- **Bilingual Switch:** Toggle seamlessly between **English** and **Chinese (中文)** in the sidebar, with succinct parameter and stage labels.
- **Dataset Selection:** View the **Study Dataset** alone, the **Community Dataset** alone, or **Both Datasets** overlaid.
- **Visual Encoding:** When displaying both datasets together on scatterplots and distributions:
  - **Study Dataset (Frontiers in Endocrinology, 2025):** Smaller black dots (`#000000`, `size=5`, `alpha=0.3`)
  - **Community Dataset (Live Google Sheet):** Prominent orange dots (`#FF8C00`, `size=8`, `alpha=0.85`)
- **Wide-Format Summary Statistics:** Displays dataset summary metrics (`N`, `Mean`, `Median`, `SD`) in a compact wide table across all key parameters.

### 2. All-Pairs Correlation Analysis & Interactive Scatterplot Explorer
- **Comprehensive Pairwise Correlations:** Separately computes pairwise correlations (Pearson $r$ and Spearman $\rho$) across all numerical parameters within each dataset.
- **Multiple-Testing Adjustment:** Automatically adjusts $p$-values across all tested parameter pairs using **Benjamini–Hochberg False Discovery Rate (FDR $q$-value)**.
- **Interactive Scatterplot:** Choose any pair of $X$ and $Y$ parameters to inspect their bivariate distribution, trendline, and correlation metrics.
- **Selectable Reference Statistics:** Toggle which dataset is used to report headline correlation statistics (defaults to the **Study Dataset**). Includes a full searchable/sortable correlation matrix and pairwise summary table.

### 3. Nearest-Neighbor (k-NN) Attrition Funnel Estimator
- **Personalized Funnel Prediction:** Matches on **Age** and **AMH** (with default $k = 10$, user-adjustable) to estimate expected counts and stage-to-stage conversion rates across the IVF attrition funnel:
  $$\text{Eggs Retrieved} \rightarrow \text{Mature (MII)} \rightarrow \text{Fertilized (2PN)} \rightarrow \text{Blastocysts} \rightarrow \text{Euploid Embryos}$$
- **Configurable Reference Pool:** Choose to find nearest neighbors from:
  - **Study Dataset Only**
  - **Community Dataset Only**
  - **Blend (Both Datasets Combined)**
- **Interactive Neighbor Inspection & Pruning:** Displays the **actual cycle records of the $k$ nearest neighbors** (including source dataset, standardized distance, baseline markers, and cycle outcomes) and allows users to **exclude/deselect specific neighbors** based on their stats (or filter out outliers/missing stages) to dynamically recompute the funnel estimate.

---

## Data Sources

### 1. Clinical Study Dataset (Local)
- **File:** [`data/data Diminished ovarian reserve is associated to euploidy rate.xlsx`](data/data%20Diminished%20ovarian%20reserve%20is%20associated%20to%20euploidy%20rate.xlsx)
- **Source Study:** Carnesi E, Castellano S, Albani E, Busnelli A, Smeraldi A, Bulbul O, Morenghi E, Immediata V, and Levi-Setti PE (2025). *Diminished ovarian reserve is associated to euploidy rate: a single center study.* **Frontiers in Endocrinology**, 15:1535776.
- **Article Link:** [https://www.frontiersin.org/journals/endocrinology/articles/10.3389/fendo.2024.1535776/full](https://www.frontiersin.org/journals/endocrinology/articles/10.3389/fendo.2024.1535776/full)
- **Description:** Single-center retrospective cohort dataset (777 cycles / 773 patients undergoing IVF with PGT-A). AMH is reported in **ng/mL**.

### 2. Community-Contributed Dataset (Live Google Sheet)
- **Live Spreadsheet:** [Google Sheets Link (`gid=920796724`)](https://docs.google.com/spreadsheets/d/1sHdyPNEOtVj8Yi_vtrNq4T2b3LaqdqYU/edit?gid=920796724#gid=920796724)
- **CSV Export Endpoint:** `https://docs.google.com/spreadsheets/d/1sHdyPNEOtVj8Yi_vtrNq4T2b3LaqdqYU/export?format=csv&gid=920796724`
- **Description:** Anonymized (numeric IDs), crowd-sourced IVF cycle outcomes collected from an online support group. AMH is reported in **ng/mL** (US standard).
- **Columns:**
  `Name`, `Egg retrieval round`, `age`, `AMH`, `FSH`, `Menopur`, `Egg retreived`, `Egg matured`, `fertilized`, `Blastocyst`, `Euploid`, `FSH+Menopur`, `Mature rate`, `Fertilize rate`, `Blastocyst rate`, `Eu rate`, `Mosaic`, `Male`, `Female`, `Blastocyst/egg`, `PGT-a`, `Fertilization`, `Trigger`, `IVF reason`, `Note`

---

## Harmonized Parameter Mapping

To enable unified visualizations, correlation analyses, and blended Nearest-Neighbor estimation, parameters from both datasets are mapped to a common schema:

| Standardized Parameter | Unit / Type | Frontiers Study (`data/*.xlsx`) | Community Google Sheet (`gid=920796724`) |
| :--- | :--- | :--- | :--- |
| `record_id` | String | Column 0 (`ICH 1`, `ICH 2`, ...) | `Name` + `Egg retrieval round` |
| `age` | Years | `paz_eta_donna` | `age` |
| `amh` | ng/mL | `esa_last_amh` | `AMH` |
| `fsh_baseline` | mIU/mL | `esa_last_fsh` | *(Study only)* |
| `afc` | Count | `ecg_afc` | *(Study only)* |
| `bmi` | $\text{kg/m}^2$ | `ana_bmi_donna` | *(Study only)* |
| `gonadotropin_total_iu` | IU | `mon_gonad_dosaggio` | `FSH+Menopur` (or `FSH` + `Menopur`) |
| `fsh_dose_iu` | IU | *(Community only)* | `FSH` |
| `menopur_dose_iu` | IU | *(Community only)* | `Menopur` |
| `eggs_retrieved` | Count | `ovo_prelevati` | `Egg retreived` |
| `eggs_mature` | Count | `ovo_post_m2` | `Egg matured` |
| `fertilized` | Count | `ovo_2pn` | `fertilized` |
| `blastocysts` | Count | `emb_biopsati` | `Blastocyst` |
| `euploid` | Count | `emb_biopsia_pgs_euploide` | `Euploid` |
| `mosaic` | Count | *(Community only)* | `Mosaic` |
| `mature_rate` | Proportion ($0\text{–}1$) | `ovo_post_m2 / ovo_prelevati` | `Egg matured / Egg retreived` (or `Mature rate`) |
| `fertilization_rate` | Proportion ($0\text{–}1$) | `ovo_2pn / ovo_post_m2` | `fertilized / Egg matured` (or `Fertilize rate`) |
| `blastocyst_rate` | Proportion ($0\text{–}1$) | `emb_biopsati / ovo_2pn` | `Blastocyst / fertilized` (or `Blastocyst rate`) |
| `euploid_rate` | Proportion ($0\text{–}1$) | `emb_biopsia_pgs_euploide / emb_biopsati` | `Euploid / Blastocyst` (or `Eu rate`) |
| `blast_per_egg_rate` | Proportion ($0\text{–}1$) | `emb_biopsati / ovo_prelevati` | `Blastocyst / Egg retreived` (or `Blastocyst/egg`) |
| `euploid_per_egg_rate` | Proportion ($0\text{–}1$) | `emb_biopsia_pgs_euploide / ovo_prelevati` | `Euploid / Egg retreived` |
| `trigger` | Category | `mon_hcg_des` | `Trigger` |
| `ivf_reason` | Category | `pma_indicazione_des` | `IVF reason` |

---

## Tech Stack & Local Setup

- **Framework:** [Streamlit](https://streamlit.io/)
- **Data & Statistics:** `pandas`, `numpy`, `scipy`, `openpyxl`
- **Interactive Charts:** `plotly`

### Running Locally (with Conda)

```bash
# 1. Create the conda environment from environment.yml
conda env create -f environment.yml

# 2. Activate the conda environment
conda activate ivf_estimator

# 3. Launch the Streamlit dashboard
streamlit run app.py
```

*(Alternatively, to create it manually: `conda create -n ivf_estimator python=3.11 -y && conda activate ivf_estimator && pip install -r requirements.txt`)*

---

## Disclaimer

This dashboard and repository are intended strictly for **informational and educational purposes** and do not constitute medical advice, diagnosis, or treatment recommendations. Always consult a qualified reproductive endocrinologist or healthcare provider regarding fertility treatments.
