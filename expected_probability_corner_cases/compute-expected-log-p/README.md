# Lower bound

We have the following lower bound using Jensen's inequality:

`(1-s)^(L-Nmut) <= expected probability`

If we try to compute expectation of log probability, this will simply give us the lower bound when comparing to some threshold. Some sample comparisons are as follows:

```
(10000, 21, 0.2, 0.1, 6.020276123141471e-05,  0.001)
(10000, 31, 0.2, 0.1, 0.3522408842623004,     0.4915)
(10000, 21, 0.3, 0.1, 0.5551667263098173,     0.6275)
(10000, 31, 0.3, 0.1, 0.9835141089434607,     0.988)
```
