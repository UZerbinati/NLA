from __init__ import *

#Singular Value Decompostion
def SVD(A,l=None,RealOnly=False):
    if A.onlyNP == True:
        U = Matrix([0],onlyNP=True)
        S = Matrix([0],onlyNP=True)
        V = Matrix([0],onlyNP=True)
        m = A.NP.shape[0]
        n = A.NP.shape[1]
        if (l == None):
            l = min([m,n])
        if l > m or l > n:
            raise NLAMathError('l can\'t be greater then min{m,n}')
        U.NP, s, V.NP = np.linalg.svd(A.NP, full_matrices=True)
        U.NP = U.NP[0:m,0:l]
        S.NP = np.diag(s[0:l])
        V.NP = V.NP[0:l,0:n]
        #Not trasposing V, alredy done by np.linalg.svd
        return U,S,V
    elif A.onlyMP == True:
        U = Matrix([0],onlyMP=True)
        S = Matrix([0],onlyMP=True)
        V = Matrix([0],onlyMP=True)
        if RealOnly == False:
            m = A.MP.rows
            n = A.MP.cols
            if (l == None):
                l = min([m,n])
            if l > m or l > n:
                raise NLAMathError('Error you should cut more element then the size of the matrix')
            else:
                U.MP,s,V.MP = mpmath.svd_c(A.MP)
                U.MP = U.MP[0:m,0:l]
                S.MP = mpmath.diag(s[0:l])
                V.MP = V.MP[0:n,0:l]
                return U,S,V
        elif RealOnly == True:
            m = A.MP.rows
            n = A.MP.cols
            if (l == None):
                l = min([m,n])
            if l > m or l > n:
                raise NLAMathError('Error you should cut more element then the size of the matrix')
            else:
                U.MP,s,V.MP = mpmath.svd_c(A.MP)
                U.MP = U.MP[0:m,0:l]
                S.MP = mpmath.diag(s[0:l])
                V.MP = V.MP[0:n,0:l]
                return U,S,V
    else:
        U = Matrix([0])
        S = Matrix([0])
        V = Matrix([0])
        m = A.NP.shape[0]
        n = A.NP.shape[1]
        if (l == None):
            l = min([m,n])
        if l > m or l > n:
            raise NLAMathError('Error you should cut more element then the size of the matrix')
        U.NP, s, V.NP = np.linalg.svd(A.NP, full_matrices=True)
        U.NP = U.NP[0:m,0:l]
        S.NP = np.diag(s[0:l])
        V.NP = V.NP[0:l,0:n]
        if RealOnly == False:
            m = A.MP.rows
            n = A.MP.cols
            if (l == None):
                l = min([m,n])
            if l > m or l > n:
                raise NLAMathError('Error you should cut more element then the size of the matrix')
            else:
                U.MP,s,V.MP = mpmath.svd_c(A.MP)
                U.MP = U.MP[0:m,0:l]
                S.MP = mpmath.diag(s[0:l])
                V.MP = V.MP[0:n,0:l]
        elif RealOnly == True:
            m = A.MP.rows
            n = A.MP.cols
            if (l == None):
                l = min([m,n])
            if l > m or l > n:
                raise NLAMathError('Error you should cut more element then the size of the matrix')
            else:
                U.MP,s,V.MP = mpmath.svd_c(A.MP)
                U.MP = U.MP[0:m,0:l]
                S.MP = mpmath.diag(s[0:l])
                V.MP = V.MP[0:n,0:l]
        return U,S,V
