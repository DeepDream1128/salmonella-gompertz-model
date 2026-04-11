# Salmonella Gompertz Growth Model — Slide Outline

---

## Slide 1: Title Slide
- **Title**: Predictive Modeling of *Salmonella enteritidis* Growth Under Sinusoidal Temperature
- Team members, course name, date

---

## Slide 2: Introduction & Background
- *Salmonella enteritidis* is a major foodborne pathogen in egg products
- Temperature fluctuations during storage/transport affect microbial growth
- Reference: Gumudavelli et al. (2007), *J. Food Sci.* 72:M254–62
- **Objective**: Use a modified Gompertz model to predict *Salmonella* growth under sinusoidal heating, estimate model parameters via OLS, and evaluate model adequacy

---

## Slide 3: Model Description
- **Primary model** (differential form for changing temperature):
  - $\frac{d\log N}{dt} = \mu \cdot K \cdot C \cdot \exp(-K)$, where $K = \exp(-\mu(t-M))$
- **Secondary model** (growth rate as function of temperature):
  - $\mu = a(T - T_{min})^2 (1 - \exp[b(T - T_{max})])$
- **Dependent variable**: $\log_{10} N(t)$ — log microbial concentration (log cfu/mL)
- **Independent variable**: $t$ — time (hr)
- **7 candidate parameters**:

| Parameter | Symbol | Description | Units | Initial Guess |
|-----------|--------|-------------|-------|---------------|
| A | $A$ | Initial log microbial number | log cfu/mL | log₁₀(400) ≈ 2.60 |
| C | $C$ | Difference between upper and lower asymptote | log cfu/mL | 11 |
| M | $M$ | Inflection point time (max growth rate) | hr | 7.5 |
| a | $a$ | Secondary model regression coefficient | °C⁻² hr⁻¹ | 0.000338 |
| b | $b$ | Secondary model regression coefficient | °C⁻¹ | 0.275 |
| T_min | $T_{min}$ | Minimum growth temperature | °C | 6 |
| T_max | $T_{max}$ | Maximum growth temperature | °C | 46.3 |

---

## Slide 4: Data Overview
- Growth data: log₁₀(CFU/mL) at multiple time points (some with replicates), n = 20
- Temperature data: sinusoidal profile T(t) recorded every minute
- Two independent time vectors: one for growth data, one for temperature data
- 📊 **Fig**: `fig07_dual_axis_logN_temp.png`

---

## Slide 5: Forward Problem — Initial Guess
- All 7 parameter guesses from literature (Gumudavelli et al., 2007)
- ODE solved with ode45; temperature interpolated via interp1
- 📊 **Fig**: `fig01_initial_guess_vs_data.png`

---

## Slide 6: SSC Analysis — All 7 Parameters (Parameter Selection Step 1)
- Scaled sensitivity coefficients (SSC) for all 7 parameters:

| Parameter | A | C | M | a | b | T_min | T_max |
|-----------|-----|------|------|------|------|-------|-------|
| Max\|SSC\| | 2.61 | 5.11 | 5.45 | 2.20 | **0.88** | **0.83** | 9.02 |

- **Decision**: T_min (0.83) and b (0.88) have the lowest SSC → model is least sensitive to these two → **fix at literature values**
- T_min = 6°C, b = 0.275 (from Gumudavelli et al.)
- 📊 **Fig**: `fig02_SSC_7params.png`

---

## Slide 7: Round 1 — Estimate 5 Parameters [A, C, M, a, T_max]
- Fixed: T_min = 6°C, b = 0.275
- Estimated: A, C, M, a, T_max (p = 5)
- 📊 **Fig**: `fig04_SSC_5params_round1.png` — SSC with initial guesses

---

## Slide 8: Round 1 Results & Parameter Selection Step 2
- nlinfit OLS estimation (5 params)
- Report: parameter estimates, RMSE, cond(J), det(J'J), CI, correlation matrix
- **Key finding**: T_max CI is wide, SSC shape coupled with other params → still poorly identifiable due to collinearity
- **Decision**: Fix T_max at Round 1 estimate → reduce to 4 parameters
- 📊 **Fig**: `fig05_SSC_5params_estimated.png`

---

## Slide 9: Parameter Selection Summary

| Parameter | SSC Rank | Identifiable? | Action | Value Used |
|-----------|----------|---------------|--------|------------|
| A | 2.61 | ✅ Yes | **Estimate** | OLS result |
| C | 5.11 | ✅ Yes | **Estimate** | OLS result |
| M | 5.45 | ✅ Yes | **Estimate** | OLS result |
| a | 2.20 | ✅ Yes | **Estimate** | OLS result |
| b | 0.88 | ❌ Low SSC | Fix (literature) | 0.275 °C⁻¹ |
| T_min | 0.83 | ❌ Lowest SSC | Fix (literature) | 6 °C |
| T_max | 9.02 | ⚠️ High SSC but collinear | Fix (Round 1 est.) | from Round 1 |

- **Stepwise elimination**: 7 → 5 (fix T_min, b) → 4 (fix T_max)
- **Final estimable parameters**: [A, C, M, a]
- Rationale: SSC magnitude alone is not sufficient; collinearity (correlation matrix, cond(J)) must also be checked

---

## Slide 10: OLS Inverse Problem — Final Results (4 params)
- **Final parameter estimates**:

| Parameter | Estimate | Std Error | Rel Error | 95% CI |
|-----------|----------|-----------|-----------|--------|
| A | ~2.24 | ~0.20 | ~8.8% | [1.83, 2.66] |
| C | ~13.39 | ~0.08 | ~0.6% | [13.23, 13.55] |
| M | ~3.53 | ~0.19 | ~5.4% | [3.13, 3.93] |
| a | ~0.000360 | ~0.00003 | ~8.3% | [0.0003, 0.0004] |

- Correlation matrix: M and a highly correlated (R ≈ -0.998), but both identifiable (narrow CI, distinct physical meaning)
- RMSE ≈ 0.37, relRMSE ≈ 6.3%

---

## Slide 11: Fitted Curve with Asymptotic CB and PB
- 📊 **Fig**: `fig06_OLS_fit_CB_PB.png`
- 95% simultaneous confidence bands (CB) and prediction bands (PB)

---

## Slide 12: Residual Analysis
- 📊 **Fig**: `fig08_residual_scatter.png` + `fig09_residual_histogram.png`
- **5 standard statistical assumptions**:
  1. ✅ Zero mean errors (mean residual ≈ 0.064)
  2. ⚠️ Constant variance — scatter plot shows slight pattern
  3. ⚠️ Uncorrelated errors — runs = 6 < min required 10.5 → autocorrelation present
  4. ✅ Normal distribution — histogram approximately symmetric
  5. ✅ Correct model — fitted curve tracks data well

---

## Slide 13: Final SSC (4 params, estimated)
- 📊 **Fig**: `fig10_SSC_4params_final.png`
- All 4 parameters have significant SSC → confirms identifiability
- SSC shapes are distinct → parameters are distinguishable

---

## Slide 14: Optimal Experimental Design
- 📊 **Fig**: `fig11_OED_delta_Cii.png`
- Δ criterion increases monotonically → more data always helps (no single optimum)
- C_ii (parameter variance) decreases with more data points for all 4 parameters
- Practical recommendation: ≥ 20 equally-spaced points for adequate estimation

---

## Slide 15: Bootstrapping of Residuals
- 1000 bootstrap iterations (residual resampling), 1000/1000 valid
- **Bootstrap 95% CI**:

| Parameter | Bootstrap 95% CI | Asymptotic 95% CI |
|-----------|-----------------|-------------------|
| A | [2.13, 2.53] | [1.83, 2.66] |
| C | [12.76, 14.40] | [13.23, 13.55] |
| M | [3.41, 3.64] | [3.13, 3.93] |
| a | [0.000340, 0.000381] | [0.0003, 0.0004] |

- 📊 **Fig**: `fig12_bootstrap_CB_PB.png`

---

## Slide 16: Asymptotic vs Bootstrap Comparison

| Metric | Asymptotic | Bootstrap |
|--------|-----------|-----------|
| CB width (avg) | 1.10 | 0.36 |
| PB width (avg) | 3.04 | 1.80 |

- **Bootstrap bands are significantly narrower**
- Asymptotic method overestimates uncertainty due to M-a correlation inflating covariance matrix
- Bootstrap provides more realistic uncertainty for this nonlinear model

---

## Slide 17: Conclusions
- Modified Gompertz model successfully predicts *Salmonella* growth under sinusoidal temperature
- **Parameter selection**: SSC analysis + collinearity check identified 4 of 7 parameters (A, C, M, a) as estimable; b, T_min, T_max fixed at literature values via stepwise elimination
- OLS fit is adequate (relRMSE ≈ 6.3%), though residual autocorrelation (runs test) suggests potential for model improvement
- Bootstrap CI are tighter than asymptotic CI → asymptotic approximation is conservative for this nonlinear problem
