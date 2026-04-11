Gompertz Growth Model

Therefore, the primary model I recommend you use is:

$\frac{d\log N(t)}{dt} = \mu KC\exp( - K)$

The secondary model for growth rate as a function of temperature is Eq.
2 in the attached paper:

$\mu = a(T - T_{\min})^2(1- \exp[b(T - T_{\max})]$

In summary, the parameters to be considered for estimation are:

*A, C, M, a, b, T~min~*, and *T~max~*, where *a* and *b* are regression
coefficients, and *T~min~* and *T~max~* are the minimum and maximum
temperatures where the microorganism will grow.

You may use the following guesses for plotting the SSC:

*A* = log(400 cfu/ml); *C* = 11 log cfu/ml; *M* = 7.5 hr; *a* = 0.000338
(ºC^2^); *b* = 0.275 ºC^-1^; *T~min~* = 6 ºC; *T~max~* = 46.3 ºC.
Examine the SSCs to determine whether all the parameters can be
estimated. If there are parameters that cannot be estimated, you can set
them at a constant value, and then renumber your remaining parameters
that will be estimated.

I recommend using inv_soln5.mlx as a templates for your ordinary least
squares code. For 7.d. in your Final Project Instructions at the end of
the Syllabus, plot log*N*(t) observed (left axis), log*N*(t) predicted
(left axis), and *T*(t) (right axis) all on the same plot using MATLAB's
yyaxis left and yyaxis right commands.

The log*N*(t) data are in the attached "Salmonella sin growth.xlsx." The
temperature data *T*(t) are in "Salmonella sin growth Temps.xlsx." These
data (log*N*(t) and *T*(t) ) are shown in Figure 5c, left side
(Sinusoidal heating) in Gumudavelli V, Subbuh J, Thippareddi H, Velugoti
PR, Froning G. 2007. Dynamic predictive model for growth of *Salmonella
enteritidis* in egg yolk. *J. Food Sci.* 72:M254--62. Note that you must
keep vectors for two times: 1) time connected to the log*N*(t) data; and
2) time connected the temperature data. I suggest using different names
for these two vectors. The largest value of the 2) vector must be \>=
the largest value of the 1) vector. All your plots should have time in
hours on the x-axis. Use "readmatrix" to import data from Excel.

Some N(t) (CFU/mL) data are missing for 0, 0.5, 1, and 2 hr. You must
read in all the data. You can make one column with duplicate times (hr)
and one column with ALL the CFU/mL, and then read in the data from these
two new columns. You will have only one value of CFU/mL for times 0,
0.5, 1, and 2 hr, and two CFU/mL values for the other times.

Use MATLAB's built-in function interp1(tobs, Tobs, time) to interpolate
temperatures from the temperature data, where the time-temperature data
are (tobs, Tobs), and "time" is the desired time where ode45 needs the
temperature. If you are using a script with function files at the bottom
(that's what I do) you can use a global statement to allow the script
and the function(s) to share the time-temperature data. Download
global_example.m file at D2L\>\>MATLAB files\>\>Inverse codes.

Syllabus：

Final presentations will be given orally, with all team members
presenting, and uploaded

electronically by April 21, 2026.

Each student group must submit final project powerpoint presentation
electronically, with

supplemental MATLAB codes "published." Each group must give a \~15-min
presentation of a

final project. Any student not attending ALL the presentations without
legitimate excuse will

automatically lose 25% of the final project grade. Project must be done
in MATLAB. If the

group does not have its own project from its research or from papers, I
will assign a project by

the second week of the course.

Here is the outline for your .ppt presentation:

Introduction: Give the background of your problem as to what the problem
is, and an overview

of the experimental methods, even though your group did not do the
experiment. State your

objectives.

1\. Choose a nonlinear model (one or more ODEs or PDEs, or equation(s)
in explicit form)

with at least 2 parameters.

a\. Good: Take model from a paper;

b\. Better: Take model from paper related to your research;

c\. Best: Take model from your own research.

d\. Alternative: Ask Dr. Dolan for model and data.

2\. Make reasonable guesses of parameter values.

Forward Problem:

3\. Make reasonable guesses of the parameters. Use the parameters to
predict Y. Plot Y and

the data to make sure the parameters are reasonably close to the true
parameters.4. Plot the scaled sensitivity coefficients for each
parameter and Ypredicted on the same

plot.

a\. Which parameters can be estimated, and why?

b\. State which parameters will be estimated most accurately, in order.

Ordinary Least Squares (OLS) Inverse Problem:

5\. Obtain data for your model.

a\. Good: Generate the data using appropriate assumptions of random
error.

b\. Better: Take data from a journal article. Can use MATLAB's GRABIT,
or some

other plot digitizer

c\. Best: Collect data from your own research.

6\. Estimate parameters using OLS and obtain all statistical results.
Report all statistical

measures, and interpret the results as to what insights we can gain:

a\. Parameter standard errors, relative errors, confidence intervals,
and correlation

matrix;

b\. Root mean square error of the dependent variable;

c\. Fitted curve and the observed data on the same plot. Show the
asymptotic CB and

PBs.

d\. Residual scatter plot and residual histogram. List whether the first
five standard

statistical assumptions were met or not.

e\. Final scaled sensitivity coefficients.

Optimal Experimental Design (Forward Problem):

7\. Plot the delta criterion and plot all C 11 , C22 , ...Cpp (p = \# of
parameters) curves, assuming

fixed parameter values, to determine an optimal experimental design. If
there is no

optimum, show how you know that.

Bootstrapping:

8\. Using either bootstrapping of data or bootstrapping of residuals,
report the 95%

confidence interval for each parameter individually, and plot the
bootstrap confidence

bands (CBs) and prediction bands (PBs). Comment on whether the bootstrap
CBs and

PBs were larger or smaller than the asymptotic CBs and PBs
