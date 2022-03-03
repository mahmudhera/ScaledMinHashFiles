In this file, we calculate and return the PDF of Nmut given L, k, and p.

## DP
This is a dynamic programming algorithm. To understand the algorithm, we notice that if we need a base to be mutated, the probability is `p`. If we require the base to be non-mutated, the probability we need is `1-p`. If we do not care if a base is mutated or not, we can simply ignore it (the probability will simply be `p + 1 - p = 1`).

First, let us express the mutated sequence as a bit string. If original and mutated sequences match at position `i`, then we have a `0` at that position in our bit string. Otherwise, we would have a `1` indication a mutation. Given such a bit string, it is easy to count the number of k-mers that has been mutated.

### Recurrence for number of ways mutations can happen

```
c[l, x, z] = # of ways there can be x mutations in a bit-string of length l, where there are z trailing zeros
```

Strings of length `k` can be:

```
0000...0
1000...0
x100...0
xx10...0
xxx1...0
.
.
xxxx...1
```

Thus, we have the following base conditions:

```
c[l=k, x=0, z=k  ] = 1

c[l=k, x=1, z=k-1] = 1
c[l=k, x=1, z=k-2] = 2
c[l=k, x=1, z=k-3] = 4
.
.
.
c[l=k, x=1, z=k-k] = 2^(k-1)
```

Now, we will use solutions for `l=l` to compute solutions for `l=l+1` as follows:
`c[l, x, z]` has `z` trailing zeros. We can add a 0 or a 1 to it.
If we add 0, then if `z < k-1`, only then we have a new mutated kmer.
Otherwise, we have a new non-mutated kmer. In both scenarios, the number of trailing zeros will increment by 1.

On the other hand, if we add 1, we will end up with a new mutated kmer, and the resulting
bit string will have no trailing zeros (it will end in a 1).

```
start with c[l+1, x, z] = 0
for every x and z, do the following:
    if z < k-1:
        c[l+1, x+1, z+1] += c[l, x, z]
    else:
        c[l+1, x, z+1] += c[l, x, z]
    c[l+1, x, 0] += c[l, x, z]
```

### Translating to probablities
From the number of ways we can have `x` mutations, we can easily make transition to the probability of having `x` mutations as follows:

```
p[l, x, z] = probability that there can be x mutations in a bit-string of length l, where there are z trailing zeros

p[l=k, x=0, z=k  ] = (1-p)^k

p[l=k, x=1, z=k-1] = p (1-p)^(k-1)
p[l=k, x=1, z=k-2] = p (1-p)^(k-2)
p[l=k, x=1, z=k-3] = p (1-p)^(k-3)
.
.
.
p[l=k, x=1, z=k-k] = p

start with p[l+1, x, z] = 0
for every x and z, do the following:
    if z < k-1:
        p[l+1, x+1, z+1] += p[l, x, z] * (1-p)
    else:
        p[l+1, x, z+1] += p[l, x, z] * (1-p)
    p[l+1, x, 0] += p[l, x, z] * p
```

## Correctness
Verified using formulas for three moments, and this pdf.
For L = 50, k = 21, p = 0.1, we get:

Expectation from formula:
44.5290505434
Expectation from pdf:
44.5290505434

2nd moment from formula:
2042.78475811
2nd moment from pdf:
2042.78475811

3rd moment from formula:
95435.2403145
3rd moment from pdf:
95435.2403145

## Running time
Runs in order of `L^3`
