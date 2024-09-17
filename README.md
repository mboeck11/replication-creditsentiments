# replication-creditsentiments


Replication files for Boeck, M. and Zörner, T. O. (2023)
=======
Replication files for Boeck, M. and Zörner, T. O. (2023) *The Impact of Credit Market Sentiment Shocks*, Journal of Money, Credit and Banking, forthcoming.

In order to reproduce all the graphics visible in the paper, just run **main.r** in the *scripts* folder. Estimations are based on 25.000 MCMC draws where the first 15.000 are discarded. Hence, the script takes a considerable amount of time to run through. The number of saved draws (*draws*) and discarded draws (*burnins*) can be adapted by the user. The script then reproduces the following

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
  + Figure F4. Internal Instrument Approach - Linear Model.
  + Figure F5. Internal Instrument Approach - Threshold Model.

- Tables:
  + Table E1. Convergence Statistics.

**Abstract** This paper investigates the role of credit market sentiments and investor beliefs on credit cycle dynamics and their propagation to business cycle fluctuations. Using US data from 1968 to 2014, we find that credit market sentiments are indeed able to detect asymmetries in a small-scale macroeconomic model. An unexpected credit market news shock exhibits different impacts in an optimistic and pessimistic credit market environment. While an unexpected movement in the optimistic regime leads to a rather low to muted impact on output and credit, we find a significant negative impact on those variables in the pessimistic regime. The findings highlight the relevance of expectation formation mechanisms as a source of macroeconomic instability.

**Links** [(Latest Version Dec 2022)](https://mboeck11.github.io/papers/BZ2023JMCB.pdf) [(WU Working Paper Jul 2019)](https://epub.wu.ac.at/7087/)
