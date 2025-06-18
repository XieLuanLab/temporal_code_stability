# Code for "Temporal coding carries more stable cortical visual representations than firing rate over time" 

This repository contains MATLAB code for reproducing key analyses and figures from the manuscript titled:  
**"Temporal coding carries more stable cortical visual representations than firing rate over time"**

The demos provided here implement core steps such as temporal code fitting, representational drift analysis, decoding, unit tracking, and longitudional dimensionality reduction for populational neural code. 

Please kindly cite the mentioned work if these code help you.  

---

## 1. Installation

1. Install **MATLAB 2023b** (tested on Linux OS) following official MATLAB installation instructions
2. Estimated installation time: **~1 hour**

---

## 2. Demos

### 2.1 Demo A – Temporal Code Fitting

- **Run:** `main_TemporalCodeFitting.m`
- **Purpose:** Fits stimulus-coding temporal components from spike trains. Reproduces **Fig. 4a–d**.
- **Input:** Loads the variable `VS_Store` (1×15 cell array). Each cell contains a 50 trials × 16 classes × 520 time bins matrix (1 ms resolution), with 0 meaning no spike and 1 indicating a spike.  
- Also loads source data for the rest of **Fig. 4**.
- **Expected Output:** Complete **Fig. 4**
- **Run Time:** ~3 minutes

---

### 2.2 Demo B – Representational Drift Index

- **Run:** `main_RepresentationalDriftIndex.m`
- **Download Data:** [https://doi.org/10.6084/m9.figshare.28877813](https://doi.org/10.6084/m9.figshare.28877813)
- **Purpose:** Calculates representational drift index across 15 time intervals based on fitted temporal components, following our interpretations from [Marks and Goard, 2021](https://www.nature.com/articles/s41467-021-25436-3).  
  Change the variable `nVS` to 2, 3, or 5 to generate additional figure panels.
- **Expected Output:** **Fig. 4f**
- **Run Time:** ~10 minutes

---

### 2.3 Demo C – UMAP Embedding of Population Activity

- **Run:** `main_Umap.m`
- **Download Data:** [https://doi.org/10.6084/m9.figshare.28877813](https://doi.org/10.6084/m9.figshare.28877813)
- **Purpose:** Generates low-dimensional representations of population activity using firing rates and fitted temporal components. Combines data across individual animals for comparison between rate and temporal codes.
- **Expected Output:** **Fig. 5**, **SFig. 10a–d**
- **Run Time:** ~20 minutes

---

### 2.4 Demo D – Decoding Across Days

- **Run:** `main_Decode_fix_slide_shuffle.m` **or** `main_Decode_rate_components.m`
- **Download Data:** [https://doi.org/10.6084/m9.figshare.28877813](https://doi.org/10.6084/m9.figshare.28877813)
- **Purpose:** Trains stimulus decoders on past days and tests on future days. Compares performance using rate codes vs temporal codes.
- **Expected Output:** **Fig. 1j**, **SFig. 3a–c**, **SFig. 10e–f**
- **Run Time:** 300+ minutes

---

### 2.5 Demo E – GLM Fit Against Behavior

- **Run:** `main_Rate_Component_GLM_FitAgainst_Behavior.m`
- **Download Data:** [https://doi.org/10.6084/m9.figshare.28877813](https://doi.org/10.6084/m9.figshare.28877813)
- **Purpose:** Fits generalized linear models (GLMs) using rate/temporal components to predict behavioral variables, see captions of SFig18d for details.
- **Expected Output:** **SFig. 18d**
- **Run Time:** Overnight

---

### 2.6 Demo F – Manual Unit Tracking Across Days

- **Run:** `main_Tracking.m`
- **Download Data:** [https://doi.org/10.6084/m9.figshare.28877813](https://doi.org/10.6084/m9.figshare.28877813)
- **Purpose:** Tracks units from the same 32-channel (2×16) electrode array (with closed channels represented as blanks) across 4 example days.
- Users are repeatedly prompted to merge (press `1`) or reject (press `0` or `Enter`) unit groups after viewing overlaid waveform plots.
- Example: For groups (1,2,3) and (4,5,6), six panels are shown—(1,2,3), (4,5,6), (1,4), (1,6), (3,4), (3,6)—including intra- and inter-group combinations involving the first and last temporally appearing units.
- Press `2` to view additional data: waveform movies, autocorrelograms (ACGs), and cross-correlograms (CCGs, when applicable: any two units from the same session).
- Advanced inputs:
  - `Nx2` array like `[1 3; 2 5]` to block merges between units 1 vs 3 and 2 vs 5.
  - `[-3 -3]` blocks unit 3 from merging with any other units.
  - Inputs with complex numbers (e.g., `[2 4i]`, `[-1 -1i]`) tentatively reject specified merges but may prompt reconsideration in future iterations.
- **Expected Output:** Procedures described in **SFig. 19** and Supplementary Methods. Waveform overlay plots of tracked units.
- **Run Time:** ~120 minutes, depending on manual input frequency

---

## 3. Instructions for Use

- Follow the demos above by running each `main_*.m` file in MATLAB.
- Resulting figures will be either displayed or saved to the `results/` folder.

---

