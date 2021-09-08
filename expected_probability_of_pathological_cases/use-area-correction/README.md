Sep 8, 2021

The problem we have is: we have an expression for the expected probability of nothing common in sketches. The formula has some issues. At scale_factor = 0, the probability rises, then falls. Conceptually, it should just be 1.0, as we are making sure that there is nothing common in the sketches.

In this nb file, we try area correction to see if the problem is solved.

### What is area correction?

When we approximate a discrete distribution (for example, N_mut) with a continuous distribution (such as the normal distribution), we often have to do area correction. By using that, probability(N_mut = n) = intergration of continuous PDF from n-0.5 to n+0.5.

### What are we doing in this file?

In this file, we calculate the precise probability that N_mut = 0, 1 and 2. Then, we calculate the same probability from normal distribution. Finally, we try to understand if working with area correction can make things better.