# %%

from collections import Counter

def solution(A):

    keys = Counter(A).keys() # equals to list(set(words))
    values = Counter(A).values() # counts the elements' frequency

    for i in range(len(values)):

        n_times = list(values)[i]

        if n_times%2==0:
            pass
        elif n_times%2 != 0:
            return list(keys)[i]

A = [500111222, 3, 9, 3, 9, 7, 7] 

solution(A)

#arr = np.array(list)


# %%
A = [9,3,9,3,9,7,9]
A[0] = 9  
A[1] = 3  
A[2] = 9
A[3] = 3  
A[4] = 9  
A[5] = 7
A[6] = 9

def solution(A):
    A.sort()
    N = len(A)
    i = 0
    while i < N:
         try:
           if A[i] == A[i+1]:
               i+=2
           else:
               return A[i]
         except:
            return A[N-1]
    pass
solution(A)
# %%
