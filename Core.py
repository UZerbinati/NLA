from __init__ import *

# Redefine norms using np.array allowing to process
# as well array not only defined using numpy.
def norm2(A,HighPrec=False):
    if A.NP == None:
        if A.MP.rows == 1:
            return mpmath.norm(A.MP,2)
        elif A.MP.cols == 1:
            return mpmath.norm(A.MP,2)
        else:
            U, s, V = mpmath.svd_c(A.MP)
            return s[0]
    elif HighPrec == False:
        return np.linalg.norm(np.array(A.NP),2)
    else:
        if A.MP.rows == 1:
            return mpmath.norm(A.MP,2)
        elif A.MP.cols == 1:
            return mpmath.norm(A.MP,2)
        else:
            U, s, V = mpmath.svd_c(A.MP)
            return s[0]

def frobeniusNorm(A):
    if A.NP == None:
        return mpmath.mnorm(A.MP,'F')
    elif HighPrec == False:
        return np.linalg.norm(np.array(A.NP),'fro')
    else:
        return mpmath.mnorm(A.MP,'F')
def Memory_Occupied(M,Verbose=False):
    Size = sys.getsizeof(M)
    if Verbose:
        if Size < 1000:
            print(Size,'Bytes')
        elif Size > 1000:
            print(Size/1000,'KBytes')
        elif Size > 1000000:
            print(Size/1000000,'MBytes')
    return Size
