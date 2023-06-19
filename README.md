# bivpoisson_ate
A Stata Command for Average Treatment Effects Estimation with Correlated Count-Valued Outcomes

When we encounter correlated count-valued outcomes y1 in {0,1,...,M}, y2 in {0,1,...,M}, the identification and estimation of average treatment effects (ATEs) need to take into account the correlation structure of the data generating process. As illustrated by Fisher, Terza, Zhang (2022), Stata command "bivpoisson" estimates the deep parameters in count-valued seemingly unrelated regression (count SUR) model. Our model affords greater precision and accuracy in terms of deep parameter estimations in comparison to single-equation Poisson model (by Stata "poisson" command). Post-Estimation Command bivpoisson_ate supports the estimation of ATEs in our count SUR model. We provide formulas for the conditional means and the ATEs of outcomes as functions of deep parameter estimates. We show, by MC simulations, that our bivpoisson_ate affords greater precision and accuracy in terms of ATEs in comparison to the ATEs estimated using "poisson" estimated parameters. We allow the treatment variable to be binary, and we plan to extend it to allow count-valued treatment. An example is provided to estimate the ATEs of private insurance status on the numbers of physician office visits, and non-physician health professional office visits within 2-week. User will specify: (1) outcome y1, (2) outcome y2, (3) policy variable, (4) a vector of control variables.![image](https://github.com/zhangyl334/bivpoisson_ate/assets/47189909/c9975e2a-b179-44d5-8847-b4a7737108f4)
