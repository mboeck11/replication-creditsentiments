# replication-creditsentiments

<<<<<<< HEAD
Replication files for Böck, M. and T. O. Zörner (2021) *The Impact of Credit Market Sentiment Shocks*
=======
Replication files for Böck, M. and T. O. Zörner (2022) *The Impact of Credit Market Sentiment Shocks*.
>>>>>>> 0b64740 (revision dec22)

In order to reproduce all the graphics visible in the paper, just run **main.r** in the *scripts* folder. Estimations are based on 25.000 MCMC draws where the first 15.000 are discarded. Hence, the script takes a considerable amount of time to run through. The number of saved draws (*draws*) and discared draws (*burnins*) can be adapted by the user. The script then reproduces the following

- Figures:
  + Figure 1. Baa Bond - Treasury Credit Spread and Its Diagnostic Expectations.
  + Figure 2. Impulse Response to a Credit Market Sentiment Shock.
  + Figure 3. Credit Market Sentiment Regimes.
  + Figure 4. Impulse Response to a Credit Market Sentiment Shock (Diagnostic Expectations).
  + Figure 5. Impulse Response to a Credit Market Sentiment Shock (Cholesky Decomposition).
  + Figure 6. Distribution of Shock on Impact for Different Expectation Formation Mechanisms.
  + Figure 7. Comparison to a Factor-Augmented TVAR.
  + Figure F1. Alternative Ordering of the VAR - Linear Model.
  + Figure F2. Alternative Ordering of the VAR - Threshold Model.
  + Figure F3. Enriching the Information Set with More Factors.

- Tables:
  + Table E1. Convergence Statistics.

**Abstract** This paper investigates the role of credit market sentiments and investor beliefs on credit cycle dynamics and their propagation to business cycle fluctuations. Using US data from 1968 to 2014, we show that credit market sentiments are indeed able to detect asymmetries in a small-scale macroeconomic model. By exploiting recent developments in behavioral finance on expectation formation in financial markets, we are able to identify an unexpected credit market news shock exhibiting different impacts in an optimistic and pessimistic credit market environment. While an unexpected movement in the optimistic regime leads to a rather low to muted impact on output and credit, we find a significant and persistent negative impact on those variables in the pessimistic regime. Therefore, this article departs from the current literature on the role of financial frictions for explaining business cycle behavior in macroeconomics and argues in line with recent theoretical contributions on the relevance of expectation formation mechanisms as a source of macroeconomic instability.

**Links** [WU Working Paper](https://epub.wu.ac.at/7087/)
