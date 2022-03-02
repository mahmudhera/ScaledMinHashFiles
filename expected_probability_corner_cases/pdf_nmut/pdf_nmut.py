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
L = 50
arr = perform_exhaustive_counting(k, p, L)

c = np.zeros((L+1, L+1, L+1))
c[k][0] = arr[0]
c[k][1] = arr[1]

for l in range(k+1, L+1):
    for x in range(L):
        for z in range(L):
            if z < k-1:
                c[l, x+1, z+1] += c[l-1, x, z]
            else:
                c[l, x, z+1] += c[l-1, x, z]
            c[l, x+1, 0] += c[l-1, x, z]
print(c[L])
print(int(sum(sum(c[L]))))
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
