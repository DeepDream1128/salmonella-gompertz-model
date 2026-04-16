# Speech Script — English (~15 min)

> Predictive Modeling of *Salmonella enteritidis* Growth Under Sinusoidal Temperature

---

## Slide 1 — Title

Good morning/afternoon, everyone. My presentation today is titled "Predictive Modeling of *Salmonella enteritidis* Growth Under Sinusoidal Temperature." We use a modified Gompertz model combined with Ordinary Least Squares parameter estimation to predict how this pathogen grows when the temperature changes over time.

---

## Slide 2 — Table of Contents

Here is a brief overview of the topics I will cover. We will start with the background and objectives, then describe the mathematical model and the experimental data. Next, I will walk through how we solved the forward problem and selected identifiable parameters using Scaled Sensitivity Coefficient analysis. After that, we will look at the OLS estimation results, model diagnostics, optimal experimental design, and a bootstrap comparison. Finally, I will conclude with the key takeaways.

---

## Slide 3 — Introduction & Background

*Salmonella enteritidis* is one of the most common foodborne pathogens worldwide, especially associated with egg and poultry products. Food safety modeling is critical because during real-world storage and transportation, the temperature is never truly constant — it fluctuates due to refrigeration cycles, ambient conditions, and handling. Understanding how bacteria respond to these dynamic temperature environments is essential for predicting shelf life and ensuring food safety.

Most traditional predictive microbiology models assume isothermal — that is, constant temperature — conditions. However, the real world is rarely that simple. This is where our work comes in: we model the growth under a sinusoidal temperature profile, which better approximates the periodic heating and cooling that occurs in practice.

Our work is built upon the study by Gumudavelli and colleagues, published in the *Journal of Food Science* in 2007. They provide both the mathematical framework — the modified Gompertz model with a secondary temperature-dependent growth rate equation — and the initial parameter values that we use as starting points.

Our three-fold objective is: first, to predict *Salmonella* growth under sinusoidal heating using the modified Gompertz model; second, to estimate model parameters via the OLS method; and third, to evaluate how adequate the model is through residual diagnostics and uncertainty analysis.

---

## Slide 4 — Model Description (Equations)

The mathematical model has two components. The primary model describes the growth kinetics in ODE form. Specifically, the rate of change of log N is given by mu times K times C times exp(-K), where K itself is an exponential function: K = exp(-mu times (t - M)). Here, C is the total growth range (the difference between the upper and lower asymptotes), M is the time at which the growth rate is maximum (the inflection point), and mu is the specific growth rate.

The secondary model connects mu to the instantaneous temperature T. It follows the form: mu = a times (T minus T_min) squared times (1 minus exp of b times (T minus T_max)). This equation has important biological meaning: growth only occurs when the temperature is between T_min and T_max. Below T_min or above T_max, the growth rate drops to zero, reflecting the fundamental biological limits of the organism.

The key point is that because temperature changes over time — following a sinusoidal profile in our case — the growth rate mu is itself time-varying. This makes the ODE non-autonomous and requires numerical integration.

---

## Slide 5 — 7 Candidate Parameters

The combined primary and secondary model has seven candidate parameters. A is the initial log microbial concentration. C is the difference between the upper and lower asymptotes of the growth curve. M is the inflection point time, where growth is fastest. The parameter a is a coefficient that controls how steeply the growth rate increases with temperature, and b controls the suppression near T_max. T_min and T_max define the biological temperature limits for growth.

All initial guesses are taken directly from the Gumudavelli et al. 2007 literature: for example, A is approximately log base 10 of 400 (about 2.60), a is 0.000338, b is 0.275, T_min is 6 degrees Celsius, and T_max is 46.3 degrees Celsius.

A natural question is: can we reliably estimate all seven from our sinusoidal growth data? The answer, as we will see, is that we must fix two of them and estimate five — not all seven are independently well constrained by this particular dataset.

---

## Slide 6 — Data Overview

Our experimental dataset consists of 20 growth measurements — log CFU per mL — collected at discrete time points over approximately 22 hours. The temperature was recorded every minute and follows a sinusoidal profile, starting from around 3 degrees, rising to about 43 degrees, and cycling back.

An important implementation detail: we maintain two separate time vectors. One corresponds to the 20 growth data points; the other corresponds to the approximately 1300 temperature measurements. The temperature time vector must span at least as long as the growth data time vector to allow proper interpolation during ODE integration.

Figure 1 shows the dual-axis plot with log N on the left axis and temperature on the right. You can see the sigmoidal growth curve overlaid with the sinusoidal temperature.

---

## Slide 7 — Forward Problem — Initial Guess

Before doing any parameter estimation, we first solve the forward problem: given the 7 literature parameter values, what does the model predict? The ODE is solved using MATLAB's ode45 integrator, and at each time step the current temperature is obtained by interpolating (interp1) from the recorded sinusoidal profile.

As Figure 2 shows, the initial guess already captures the general S-shaped growth trend reasonably well — the curve goes from about 2.6 to around 8 over 17 hours. However, the fit is not precise, which motivates the inverse problem: adjusting parameters to better match the observed data.

---

## Slide 8 — SSC Analysis — All 7 Parameters

Now we arrive at a critical step: determining which parameters can be reliably estimated from this dataset. We use the Scaled Sensitivity Coefficient, defined as SSC = beta times partial Y over partial beta. This dimensionless quantity measures how much the model output changes when each parameter is perturbed proportionally.

If a parameter has a high Max |SSC|, it means the data is highly sensitive to that parameter — small changes in its value lead to large changes in the model prediction. Such parameters can be estimated accurately. Conversely, a low SSC means the data is insensitive to the parameter — its value barely affects the prediction, so the data cannot constrain it.

Looking at the table: T_min has the lowest Max |SSC| at 0.83, and b is 0.88. Both are well below 1, which is generally considered the threshold for estimability. In practice, this means that if you try to fit T_min or b freely, the optimizer will struggle — the objective function is nearly flat in those directions.

Therefore, we fix T_min at 6°C and b at 0.275, both at their literature values.

You'll also notice that T_max has a high SSC of 9.02. It is also strongly coupled with C in the sensitivity structure. On the next slide I will explain how we still estimate T_max jointly with the other primary parameters rather than fixing it at the literature value.

---

## Slide 9 — Parameter Selection Summary

This table summarizes the complete parameter selection. We go directly from 7 to 5 parameters.

T_min and b are fixed due to their low SSC — the data simply cannot inform us about their values. We do not fix T_max: it is estimated together with A, C, M, and a. The correlation R(C, T_max) equals 0.98, so C and T_max remain strongly coupled; even so, the asymptotic 95% confidence interval for T_max excludes zero — so T_max is estimable and we retain it in the fitted set rather than fixing it at the literature value.

Key correlations in the Jacobian-based linearization include R(C, T_max) = 0.98, R(C, a) = -0.97, and R(a, T_max) = -0.94. These high correlations inflate standard errors for a and T_max, as we will see in the OLS results, but they do not force us to fix T_max.

The condition number of the Jacobian is 3.14 times 10 to the 5th — still below the usual 10 to the 6th concern threshold, indicating a numerically workable fit for the estimated set A, C, M, a, and T_max.

---

## Slide 10 — OLS Results

Here are the final OLS estimation results for the five estimated parameters — A, C, M, a, and T_max.

A is 2.358, with standard error 0.146, 6.2% relative error, and 95% confidence interval [2.046, 2.670] — so the initial bacterial concentration is about 10^2.36 CFU/mL.

C, the total growth range, is 11.062 log units, standard error 2.013, 18.2% relative error, and interval [6.771, 15.354].

M, the inflection point, is 5.541 hours, standard error 0.528, 9.5% relative error, and interval [4.414, 6.667].

a, the secondary-model coefficient, is 3.0 times 10 to the minus 4, standard error 1.1 times 10 to the minus 4, 38.3% relative error, and interval from 1.0 times 10 to the minus 4 to 6.0 times 10 to the minus 4.

T_max is 45.72 degrees Celsius, standard error 9.15, 20.0% relative error, and interval [26.23, 65.22].

a and T_max have the largest relative errors, consistent with strong correlations among a, C, and T_max — in particular R(C, T_max) = 0.98, R(C, a) = -0.97, and R(a, T_max) = -0.94.

The overall fit quality is excellent: RMSE is 0.232, relative RMSE is 3.84%, R-squared is 0.989, and the Jacobian condition number is 3.14 times 10 to the 5th, consistent with the identifiability discussion on the previous slide.

---

## Slide 11 — Fitted Curve with Asymptotic CB and PB

Figure 6 shows the OLS fitted curve together with two types of uncertainty bands. The narrower shaded region is the 95% simultaneous confidence band, which represents uncertainty in the mean prediction — if we repeated the experiment many times, the true mean curve would fall within this band 95% of the time. The wider band is the 95% prediction band, which accounts for both parameter uncertainty and observation noise — individual future observations would be expected to fall within this band.

All 20 data points fall within the prediction band, which is encouraging for model adequacy.

---

## Slide 12 — Residual Analysis

Residual analysis checks the five standard OLS assumptions. Looking at the scatter plot and histogram:

Zero mean: satisfied — the mean residual is approximately 0.000, essentially zero relative to the data range. Constant variance: borderline — there may be slight heteroscedasticity, with residuals appearing somewhat larger at certain time points. Uncorrelated residuals: not satisfied — the runs test gives only 6 runs versus a minimum required of 10.5. This means consecutive residuals tend to have the same sign, indicating autocorrelation. This is a known limitation and suggests the model structure could potentially be improved — perhaps with a more complex growth model or a different secondary model. Normality: satisfied — the histogram is approximately bell-shaped. Correct model: generally satisfied — the overall trend is well captured.

The autocorrelation is the main concern and is worth noting as a direction for future improvement.

---

## Slide 13 — SSC Convergence

This is an important validation step. We compare the SSC profiles computed before estimation — using the initial guess from the first round — and after estimation — using the OLS-optimized parameters.

For A, C, and M, the shapes of the SSC curves are similar before and after — the model remains most sensitive to these structural parameters in much the same way. For a and T_max, the SSC amplitudes are reduced by roughly fifty percent after estimation, because the initial guesses for those two were relatively far from the final OLS values; once the parameters move toward the data-supported region, local sensitivities shrink.

Despite that reduction, every parameter still has Max absolute SSC greater than or equal to 1, which confirms identifiability of the five-parameter set. No further parameter reduction is required.

---

## Slide 14 — Optimal Experimental Design

For optimal experimental design we follow Beck and Arnold (1977), Equation 8.3.5. We compute the information matrix C as a function of experiment duration t: C_ij(t) = (1/t) times the integral from 0 to t of X prime sub i times X prime sub j d tau, where X prime denotes the scaled sensitivity coefficients. We set Delta equal to the determinant of C. The integral is evaluated numerically using cumtrapz on a fine time grid spanning the full temperature data range.

Delta shows a clear maximum at approximately t equals 17 hours, indicating the optimal experiment duration for this model and parameter set. Beyond this point, additional measurement time yields diminishing information content per unit time — the growth curve approaches its plateau, and the model becomes less sensitive to the dynamic parameters.

In the C_ii plot, M has the largest diagonal element, peaking around 10 hours then declining — consistent with M governing the inflection timing, which is most informative during the active growth phase. C increases steadily as the model output approaches the upper asymptote. A remains roughly constant. T_max grows gradually, and a has the smallest C_ii throughout — consistent with it being the hardest parameter to constrain.

The fact that our actual 20-measurement experiment spans about 17 hours aligns well with the optimal duration identified by this analysis.

---

## Slide 15 — Bootstrapping

The asymptotic confidence intervals we computed earlier rely on a linear approximation of the model around the parameter estimates. For nonlinear models, this approximation can be inaccurate. Bootstrap provides a more robust alternative.

We performed 600 bootstrap iterations using residual resampling: we take the estimated residuals, resample them with replacement, add them back to the fitted values to create pseudo-data, and re-estimate the parameters. 522 runs yielded valid estimates; 1 failed to converge, and 77 were discarded as outliers — so 522 of 600 usable draws.

Bootstrap 95% confidence intervals are: A from 2.093 to 2.606; C from 10.27 to 12.94; M from 4.652 to 5.836; a from 2.5 times 10 to the minus 4 to 4.4 times 10 to the minus 4; and T_max from 43.26 to 152.3 degrees Celsius. The interval for T_max is very wide, reflecting its strong correlation with C — the same coupling seen in the Jacobian correlations.

The table compares these bootstrap intervals to the asymptotic ones from OLS. Figure 11 shows the resulting bootstrap confidence and prediction bands.

---

## Slide 16 — Asymptotic vs Bootstrap Comparison

The quantitative comparison is striking. The bootstrap confidence band has an average width of 0.365, compared to 0.867 for the asymptotic band — about 2.4 times narrower. The prediction band is 1.274 versus 2.119 — about 1.7 times narrower.

This pattern indicates that the asymptotic linear approximation overstates band widths for this nonlinear, strongly coupled parameterization. The bootstrap, by empirically resampling residuals, yields tighter uncertainty envelopes for both mean and predictive intervals in this application.

---

## Slide 17 — Conclusions

To summarize:

First, the modified Gompertz model with the temperature-dependent secondary model successfully predicts *Salmonella enteritidis* growth under realistic sinusoidal temperature conditions.

Second, through SSC analysis and collinearity diagnosis, we performed a direct 7-to-5 parameter reduction. T_min and b remain fixed at literature values due to low sensitivity; T_max is estimated together with A, C, M, and a, with R(C, T_max) = 0.98 indicating strong coupling.

Third, the OLS fit is strong: R-squared of 0.989, relative RMSE of 3.84%, and Jacobian condition number 3.14 times 10 to the 5th. a and T_max show the largest relative errors, consistent with high inter-parameter correlations. The runs test result of 6 versus 10.5 still suggests residual autocorrelation — a direction for future model refinement.

Fourth, bootstrap analysis on 522 valid replicates shows asymptotic bands are wider: about 2.4 times for confidence bands and 1.7 times for prediction bands compared with bootstrap — underscoring the value of resampling for this nonlinear fit.

Thank you for your attention. I am happy to take any questions.

---

## Slide 18 — References

Our primary reference is Gumudavelli et al. 2007, "Growth and Inactivation of *Salmonella* in Liquid Whole Egg," published in the *Journal of Food Science*, volume 72, pages M254-262.
