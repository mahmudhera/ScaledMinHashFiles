## Overview
The codes in this directory are used to calculate the expected probability of the cases where there is nothing common between the sketches.

To get nothing common in two sketches, there are $`L-N_{mut}`$ k-mers that are not mutated by the mutation process, and all of these need to be thrown out while taking a sketch. The probability of that happening is $`(1-s)^{L - N_{mut}}`$. Next, using mathematics, we can calculate the expectation of this probability.

The python file implements this expected probability. The pdf files are sample plots.

### plot_expected_probability_of_pathological_cases.python

This is the code used to plot the expected probability. The other files are used to make the whole thing run.

### check_expected_probability_against_simulations.py

In this script, the code to run the following experiment resides: since we do not know if the formula (for expected probability of having nothing common in the sketches) is correct, we calculate the probability of finding nothing in common in the two sketches (before and after the mutation process) from multiple simulations, and compare the value with what our formula spits out.

