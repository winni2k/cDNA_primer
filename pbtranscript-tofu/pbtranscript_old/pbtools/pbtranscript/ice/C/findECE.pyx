from libc.stdlib cimport malloc, free

cdef inline double double_min(double a, double b): return a if a <= b else b
cdef inline double double_max(double a, double b): return a if a >= b else b


def alignment_has_large_nonmatch(list ece_arr, int penalty, int min_len):
    """
    penalty of (-)1: 50%
    penalty of (-)2: 66%
    penalty of (-)4: 80%
    penalty of (-)9: 90%

    Return True when alignment has large non-matches not explained by low base QVs
    (in other words, "reject" as an isoform hit and don't put in the same cluster)
    """
    cdef int i
    cdef int L = len(ece_arr) + 1
    cdef double *arr2 = <double *>malloc(L * sizeof(double))

    # replacing the following:
    # ece_arr = ece_arr * (penalty + 1)
    # s = [0] + list(ece_arr - penalty)
    # calling: findECE(s, len(s), min_len, True)
    arr2[0] = 0
    for i in xrange(1, L):
        arr2[i] = (<double>ece_arr[i-1] * (penalty + 1)) - penalty
        
    answer = len(findECE2(arr2, L, min_len, True)) > 0
    free(arr2)

    return answer # fix this later to something faster & better

cdef findECE2(double *s, int L, int min_ece_length, give_up_as_soon_one_found=False):
    """
    s           --- the vector containing the -1 and +4s
    L           --- just the length of s
    min_ece_length --- the minimum ece length

    NOTE: s[0] should be 0, the extra 0 we need for the ECE alg
          so the real data should start from s[1], s[2]....

    Returns all valid ECEs (>= min_ece_length) in s, which could be overlapping
    Call remove_overlapping_bests afterwards to process overlapping ECEs
    """
    bests = []
    cdef int i, j
    cdef double *r = <double *>malloc(L * sizeof(double)) # r[i] = sum(s[0], s[1]...s[i])
    cdef double *X = <double *>malloc(L * sizeof(double)) # X[i] = min(r[0], r[1]...r[i])
    cdef double *Y = <double *>malloc(L * sizeof(double)) # Y[i] = max(r[i], r[i+1]...)

    r[0] = 0
    for i in range(1, L):
        r[i] = r[i-1] + s[i]

    X[0] = 0
    for i in range(1, L):
        X[i] = double_min(X[i-1], r[i])

    Y[L-1] = r[L-1]
    for i in range(L-2, -1, -1):
        Y[i] = double_max(Y[i+1], r[i])

    i = j = 0
    while j < L:
        if j == L-1 or Y[j+1] < X[i]:
            if j - i >= min_ece_length:
                if give_up_as_soon_one_found:
                    free(r)
                    free(X)
                    free(Y)
                    return [(i+1, j)]
                else: bests.append((i+1, j))
            j += 1
            while j < L and i < L and Y[j] < X[i]:
                i += 1
        else:
            j += 1

    free(r)
    free(X)
    free(Y)

    return bests


def findECE(list s, int L, int min_ece_length, give_up_as_soon_one_found=False):
    """
    s           --- the vector containing the -1 and +4s
    L           --- just the length of s
    min_ece_length --- the minimum ece length

    NOTE: s[0] should be 0, the extra 0 we need for the ECE alg
          so the real data should start from s[1], s[2]....

    Returns all valid ECEs (>= min_ece_length) in s, which could be overlapping
    Call remove_overlapping_bests afterwards to process overlapping ECEs
    """
    assert len(s) == L and s[0] == 0

    bests = []
    cdef int i, j
    cdef double *r = <double *>malloc(L * sizeof(double)) # r[i] = sum(s[0], s[1]...s[i])
    cdef double *X = <double *>malloc(L * sizeof(double)) # X[i] = min(r[0], r[1]...r[i])
    cdef double *Y = <double *>malloc(L * sizeof(double)) # Y[i] = max(r[i], r[i+1]...)

    r[0] = 0
    for i in range(1, L):
        r[i] = r[i-1] + s[i]

    X[0] = 0
    for i in range(1, L):
        X[i] = double_min(X[i-1], r[i])
    
    Y[L-1] = r[L-1]
    for i in range(L-2, -1, -1):
        Y[i] = double_max(Y[i+1], r[i])

    i = j = 0
    while j < L:
        if j == L-1 or Y[j+1] < X[i]:
            if j - i >= min_ece_length: 
                if give_up_as_soon_one_found:
                    free(r)
                    free(X)
                    free(Y)
                    return [(i+1, j)]
                else: bests.append((i+1, j))
            j += 1
            while j < L and i < L and Y[j] < X[i]:
                i += 1
        else:
            j += 1

    free(r)
    free(X)
    free(Y)

    return bests

def test():
    s = [+1 if int(x)==1 else -4 for x in list('0111011111011011101101')]
    if findECE([0]+s, len(s)+1, 10) == [(1, 10), (2, 17), (6, 20)]:
        print("Test passed.")
    else:
        print("Test not passed.")

