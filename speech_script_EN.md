# Speech Script — English (~15 min)

> Predictive Modeling of *Salmonella enteritidis* Growth Under Sinusoidal Temperature

---

## Slide 1 — Title

Good morning/afternoon, everyone. My presentation today is titled "Predictive Modeling of *Salmonella enteritidis* Growth Under Sinusoidal Temperature." We use a modified Gompertz model combined with Ordinary Least Squares parameter estimation to predict how this pathogen grows when the temperature changes over time.

---

## Slide 2 — Table of Contents

Here is a brief overview of the topics I will cover. We will start with the introduction and background, then describe the mathematical model. Next, I will present the data overview, followed by the forward problem. After that, we will discuss parameter selection using Scaled Sensitivity Coefficient analysis, OLS estimation results, model diagnostics, optimal experimental design, and bootstrapping. Finally, I will conclude with the key takeaways.

---

## Slide 3 — Introduction & Background

*Salmonella enteritidis* is a major foodborne pathogen, especially associated with egg products. Food safety modeling is critical because during real-world storage and transportation, the temperature is never truly constant — refrigeration cycles, ambient conditions, and handling all cause fluctuations. Understanding how bacteria respond to these dynamic temperature environments is essential for predicting shelf life and ensuring food safety.

Most traditional predictive microbiology models assume isothermal conditions. However, the real world is rarely that simple. This is where our work comes in: we model the growth under a sinusoidal temperature profile, which better approximates the periodic heating and cooling that occurs in practice.

Our work is built upon the study by Gumudavelli et al., published in the *Journal of Food Science* in 2007. They provide both the mathematical framework — the modified Gompertz model with a secondary temperature-dependent growth rate equation — and the initial parameter values we use as starting points.

Our objective is: use a modified Gompertz model to predict *Salmonella* growth under sinusoidal heating, estimate model parameters via OLS, and evaluate model adequacy.

---

## Slide 4 — Model Description (Equations)

The mathematical model has two components.

The primary model describes the growth kinetics in differential form. The rate of change of log N is given by mu times K times C times exp(-K), where K equals exp(-mu times (t minus M)). Here, C is the total growth range — the difference between the upper and lower asymptotes; M is the inflection point time where growth rate is maximum; and mu is the specific growth rate.

The secondary model connects mu to the instantaneous temperature T: mu = a times (T minus T_min) squared times (1 minus exp of b times (T minus T_max)). This equation has important biological meaning: growth only occurs between T_min and T_max. Below T_min or above T_max, the growth rate drops to zero.

The dependent variable is log base 10 of N(t), in units of log CFU per mL. The independent variable is time t in hours. Because temperature follows a sinusoidal profile, mu is time-varying, making the ODE non-autonomous and requiring numerical integration.

---

## Slide 5 — 7 Candidate Parameters

The combined model has seven candidate parameters. Let me walk through them:

- A: initial log microbial number, in log CFU/mL, initial guess log₁₀(400) ≈ 2.60
- C: upper–lower asymptote difference, in log CFU/mL, initial guess 11
- M: inflection point time, in hours, initial guess 7.5
- a: secondary model coefficient, in °C⁻² hr⁻¹, initial guess 0.000338
- b: secondary model coefficient, in °C⁻¹, initial guess 0.275
- T_min: minimum growth temperature, initial guess 6°C
- T_max: maximum growth temperature, initial guess 46.3°C

All initial guesses are taken directly from Gumudavelli et al. 2007. A natural question is: can we reliably estimate all seven from our sinusoidal growth data? As we will see, the answer is no — some parameters are not identifiable from this dataset.

---

## Slide 6 — Data Overview

Our experimental dataset consists of 20 growth measurements — log CFU per mL — collected at discrete time points over approximately 22 hours. The temperature was recorded every minute and follows a sinusoidal profile. We maintain two separate time vectors: one for the 20 growth data points, and one for the approximately 1300 temperature measurements. The temperature time vector must span at least as long as the growth data time vector to allow proper interpolation during ODE integration.

Figure 1 shows the dual-axis plot with log₁₀ N on the left axis (both observed and predicted) and temperature T on the right axis, all versus time.

---

## Slide 7 — Forward Problem — Initial Guess

Before parameter estimation, we solve the forward problem: given the 7 literature parameter values, what does the model predict? The ODE is solved using MATLAB's ode45 integrator, and at each time step the current temperature is obtained by interpolating via interp1 from the recorded sinusoidal profile.

As Figure 2 shows, the initial guess already captures the general S-shaped growth trend reasonably well — the curve goes from about 2.6 to around 8 over 17 hours. However, the fit is not precise, which motivates the inverse problem: adjusting parameters to better match the observed data.

---

## Slide 8 — SSC Analysis — All 7 Parameters (Step 1)

Now we arrive at a critical step: determining which parameters can be reliably estimated. We use the Scaled Sensitivity Coefficient, defined as SSC = beta times partial Y over partial beta. This dimensionless quantity measures how much the model output changes when each parameter is perturbed proportionally.

If a parameter has a high Max|SSC|, the data is highly sensitive to it — small changes lead to large changes in prediction, so it can be estimated accurately. Conversely, low SSC means the data cannot constrain that parameter.

Looking at the Max|SSC| values: A is 2.61, C is 5.11, M is 5.45, a is 2.20, b is 0.88, T_min is 0.83, and T_max is 9.02.

T_min has the lowest at 0.83, and b is 0.88. Both are below 1, the general threshold for estimability. Therefore, we fix T_min at 6°C and b at 0.275, both at their literature values.

---

## Slide 9 — Parameter Selection Summary

This table summarizes the complete parameter selection process. We go from 7 to 5 estimable parameters:

- A (Max|SSC| = 2.61): identifiable → estimate via OLS
- C (5.11): identifiable → estimate
- M (5.45): identifiable → estimate
- a (2.20): identifiable → estimate
- b (0.88): low SSC → fix at 0.275 °C⁻¹
- T_min (0.83): lowest SSC → fix at 6°C
- T_max (9.02): identifiable → estimate

The direct elimination takes us from 7 to 5 parameters. The final estimable parameter set is [A, C, M, a, T_max].

---

## Slide 10 — OLS Results

Here are the final OLS estimation results for the five parameters.

A is estimated at 2.358, with standard deviation 0.146, relative error 6.2%, and 95% confidence interval [2.046, 2.670].

C is 11.062, σ = 2.013, relative error 18.2%, CI [6.771, 15.354].

M is 5.541 hours, σ = 0.528, relative error 9.5%, CI [4.414, 6.667].

a is 3.0×10⁻⁴, σ = 1.1×10⁻⁴, relative error 38.3%, CI [1.0×10⁻⁴, 6.0×10⁻⁴].

T_max is 45.72°C, σ = 9.15, relative error 20.0%, CI [26.23, 65.22].

Fixed parameters: T_min = 6°C, b = 0.275.

Key statistics: condition number of J is 3.14×10⁵, RMSE = 0.232, relative RMSE = 3.84%, R² = 0.989.

Notable correlations: R(C, T_max) = 0.98, R(C, a) = -0.97, R(a, T_max) = -0.94, R(C, M) = 0.92. All 95% confidence intervals exclude zero, confirming estimability.

---

## Slide 11 — Fitted Curve with Asymptotic CB and PB

Figure 6 shows the OLS fitted curve together with two types of uncertainty bands. The narrower shaded region is the 95% simultaneous confidence band — representing uncertainty in the mean prediction. The wider band is the 95% prediction band, which accounts for both parameter uncertainty and observation noise.

All 20 data points fall within the prediction band, which is encouraging for model adequacy.

---

## Slide 12 — Residual Analysis

Residual analysis checks the five standard OLS assumptions. Looking at the residual scatter plot and histogram:

1. Zero mean: satisfied — mean residual ≈ 0.000
2. Constant variance: borderline — some heteroscedasticity may be present
3. Uncorrelated residuals: not satisfied — runs test gives 6 runs, below the minimum required 10.5, indicating autocorrelation
4. Normality: satisfied — histogram is approximately bell-shaped
5. Correct model: satisfied — overall trend is well captured

The residual autocorrelation is the main concern and suggests the model structure could potentially be improved with a more complex growth model or different secondary model.

---

## Slide 13 — SSC Convergence — Before & After Estimation

This is an important validation step. We compare the SSC profiles computed before estimation (using initial guesses) and after estimation (using OLS-optimized parameters).

A, C, and M show stable SSC curves — similar shapes before and after. a and T_max show SSC amplitudes reduced by roughly 50% after estimation, because the initial guesses were relatively far from the final OLS values. However, all parameters still have Max|SSC| ≥ 1, confirming identifiability of the five-parameter set. No further parameter reduction is required.

---

## Slide 14 — Optimal Experimental Design

For optimal experimental design, we follow Beck & Arnold (1977), Equation 8.3.5. We define the information matrix C as a function of experiment duration t:

C_ij(t) = (1/t) ∫₀ᵗ X'_i · X'_j dτ

where X' denotes the scaled sensitivity coefficients. Delta = det(C) serves as the D-optimality criterion.

Delta shows a clear maximum at approximately t ≈ 17 hours, indicating the optimal experiment duration. Beyond this point, the growth curve approaches its plateau and the model becomes less sensitive to dynamic parameters.

In the C_ii plot: M has the largest diagonal element, peaking around 10 hours then declining. C grows steadily. a has the smallest C_ii throughout — consistent with it being the hardest parameter to constrain.

Our actual experiment spans about 17 hours, aligning well with this optimal duration.

---

## Slide 15 — Bootstrapping of Residuals

The asymptotic confidence intervals rely on a linear approximation of the model around the parameter estimates. For nonlinear models, this can be inaccurate. Bootstrap provides a more robust alternative.

We performed 600 bootstrap iterations using residual resampling. 522 out of 600 yielded valid estimates.

Bootstrap 95% confidence intervals compared with asymptotic:

| Parameter | Bootstrap 95% CI | Asymptotic 95% CI |
|-----------|-----------------|-------------------|
| A | [2.093, 2.606] | [2.046, 2.670] |
| C | [10.27, 12.94] | [6.771, 15.354] |
| M | [4.652, 5.836] | [4.414, 6.667] |
| a | [2.5e-4, 4.4e-4] | [1.0e-4, 6.0e-4] |
| T_max | [43.26, 152.3] | [26.23, 65.22] |

Figure 11 shows the bootstrap confidence and prediction bands.

---

## Slide 16 — Asymptotic vs Bootstrap Comparison

The quantitative comparison is striking:

| Metric | Asymptotic | Bootstrap |
|--------|-----------|-----------|
| CB width (avg) | 0.867 | 0.365 |
| PB width (avg) | 2.119 | 1.274 |

Bootstrap CB is 2.4 times narrower; bootstrap PB is 1.7 times narrower. This indicates the asymptotic method overestimates uncertainty — the high condition number (3.14×10⁵) inflates the covariance matrix. Bootstrap provides more realistic uncertainty quantification for this nonlinear model.

---

## Slide 17 — Conclusions

To summarize:

First, the modified Gompertz model successfully predicts *Salmonella enteritidis* growth under sinusoidal temperature conditions.

Second, SSC analysis led to parameter selection: fix T_min and b (low SSC); estimate [A, C, M, a, T_max] — a 7-to-5 reduction. R(C, T_max) = 0.980, but all 95% CI exclude zero.

Third, OLS fit quality: condition number 3.14×10⁵, relative RMSE = 3.84%, R² = 0.989. Residual autocorrelation (runs = 6 < 10.5) remains a limitation.

Fourth, bootstrap CI are narrower — 2.4× for CB, 1.7× for PB — indicating asymptotic methods overestimate uncertainty due to high parameter correlation.

Thank you for your attention. I am happy to take any questions.

---

## Slide 18 — References

Our primary reference is Gumudavelli V, Subbuh J, Thippareddi H, Velugoti PR, Froning G. 2007. Growth and Inactivation of *Salmonella* in Liquid Whole Egg. *J. Food Sci.* 72:M254–62.
