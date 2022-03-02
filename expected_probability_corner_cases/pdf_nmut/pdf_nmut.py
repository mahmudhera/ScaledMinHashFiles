import numpy as np

def perform_exhaustive_counting(k, p, L):
    arr = np.zeros((2, L+1))
    # arr[i,j] indicates c[k,i,j]
    arr[0, k] = 1 # (1-p)**k
    for i in range(1,k+1):
        arr[1, k-i] = 2**(i-1) # p * (1-p)**(k-i)
    return arr

k = 21
p = 0.1
L = 1000
arr = perform_exhaustive_counting(k, p, L)

curr = np.zeros((L+1, L+1))
next = np.zeros((L+1, L+1))
curr[0] = arr[0]
curr[1] = arr[1]

print(curr)

for l in range(k+1, L+1):
    next[next>0] = 0
    for x in range(L):
        for z in range(L):
            if z < k-1:
                next[x+1, z+1] += curr[x, z]
            else:
                next[x, z+1] += curr[x, z]
            next[x+1, 0] += curr[x, z]
    next, curr = curr, next
print(curr)
print(int(sum(sum(curr))))
print(2**L)
'''
print(sum(sum(arr)))
print(2**k)

print(arr)
L = k+1
arr_L = np.zeros((k+1, k+2))
for x in range(L-k+1):
    for z in range(k+1):
        if z < k-1:
            arr_L[x+1, z+1] += arr[x,z]
        else:
            arr_L[x, z+1] += arr[x,z]
        arr_L[x+1, 0] += arr[x,z]
        print(x, z)
print(arr_L)
print(sum(sum(arr_L)))
'''
'''
[[  0.   0.   0.   0.   0.   0.   1.]
 [  1.   0.   0.   0.   0.   1.   0.]
 [ 31.  16.   8.   4.   2.   0.   0.]
 [  0.   0.   0.   0.   0.   0.   0.]
 [  0.   0.   0.   0.   0.   0.   0.]
 [  0.   0.   0.   0.   0.   0.   0.]]
'''
