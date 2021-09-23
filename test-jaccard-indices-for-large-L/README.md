### Test the statistical significance of the jaccard confidence interval

We noticed that for small L, the % of times true mutation rate lies within our calculated confidence interval is very off (not close to 96%).

```
(k, mut.rate, L, s, %success, num_success, num_failure, above_point_est, below_point_est)
(21, 0.01, 100000, 0.1, 0.717, 717, 283, 484, 516)
(31, 0.01, 100000, 0.1, 0.744, 744, 256, 471, 529)
(21, 0.01, 10000, 0.4, 0.655, 655, 345, 444, 556)
(31, 0.01, 10000, 0.4, 0.706, 706, 294, 458, 542)
```

This occurs when mutation rate is small (p = 0.01). Therefore, we decided to make L really large, and see if that improves results. We generated the scaled jaccard indices for L = 10^9 using a mutation model written in C++ (to reduce the time to run the experiments). You can see the directory by going to `faster-mutation-model`. Then, we collected those scaled jaccard indices, and calculated the confidence interval for mutation rate for each of these indices.

#### Dealing with very large L
Our implementation of the confidence interval for mutation rate (calculated from a scaled jaccard index) involves calculating the third moment of N_mut. To do so, we need to use the Hypergeometric function. There are four instances where we had to use it:

1. Hypergeometric (1, 2 + k - L, k - L, 1)
1. Hypergeometric (1, 2 + k - L, k - L, 1 - p)
1. Hypergeometric (1, 1 + 2k - L, -1 + 2k - L, 1)
1. Hypergeometric (1, 1 + 2k - L, -1 + 2k - L, 1 - p)

When L is really large, these do not converge to a finitite value. We checked this in python (returns `nan`) and Mathematica (keeps running the loop forever).

#### Quick fix

To make a quick fix, we used this property of Hypergeometric:

```
Hypergeometric(a, b, b, z) = 1/(1-z)^a for all b
```

When L is really large, we can assume that the second and third arguments are the same. Using this, Hypergeometric (1, 2 + k - L, k - L, 1) will be 1/0. With L relly large, we simply return infinity whenever this happens.


#### Results

Unfortunately, even with this fix and really large L, the accuracy of the confidence interval does not really improve! The accuracy is only 35%.

```
(k, mut.rate, L, s, %success, num_success, num_failure, above_point_est, below_point_est)
(21, 0.01, 1000000000, 0.1, 0.3515704154002026, 347, 640, 476, 511)
```

