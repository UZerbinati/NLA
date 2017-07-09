from __init__ import *
import Core
from copy import *
#Singular Value Decompostion
def SVD(A,l=None,RealOnly=False):
    if A.onlyNP == True:
        U = Matrix([[0]],onlyNP=True)
        S = Matrix([[0]],onlyNP=True)
        V = Matrix([[0]],onlyNP=True)
        m = A.rows
        n = A.cols
        if (l == None):
            l = min([m,n])
        if l > m or l > n:
            raise NLAMathError('l can\'t be greater then min{m,n}')
        U.NP, s, V.NP = np.linalg.svd(A.NP, full_matrices=True)
        U.NP = U.NP[0:m,0:l]
        U.Array =U.NP.tolist()
        U.rows = U.NP.shape[0]
        U.cols = U.NP.shape[1]
        S.NP = np.diag(s[0:l])
        S.rows = S.NP.shape[0]
        S.cols = S.NP.shape[1]
        S.Array = S.NP.tolist()
        V.NP = V.NP[0:l,0:n]
        V.rows = V.NP.shape[0]
        V.cols = V.NP.shape[1]
        V.Array =V.NP.tolist()
        #Not trasposing V, alredy done by np.linalg.svd
        return U,S,V
    elif A.onlyMP == True:
        U = Matrix([[0]],onlyMP=True)
        S = Matrix([[0]],onlyMP=True)
        V = Matrix([[0]],onlyMP=True)
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
                U.Array =U.MP.tolist()
                S.MP = mpmath.diag(s[0:l])
                S.Array =S.MP.tolist()
                V.MP = V.MP[0:n,0:l]
                V.Array =V.MP.tolist()
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
                U.Array =U.MP.tolist()
                S.MP = mpmath.diag(s[0:l])
                S.Array =S.MP.tolist()
                V.MP = V.MP[0:n,0:l]
                V.Array =V.MP.tolist()
                return U,S,V
    else:
        U = Matrix([[0]],onlyNP=False)
        S = Matrix([[0]],onlyNP=False)
        V = Matrix([[0]],onlyNP=False)
        m = A.NP.shape[[0]]
        n = A.NP.shape[1]
        if (l == None):
            l = min([m,n])
        if l > m or l > n:
            raise NLAMathError('Error you should cut more element then the size of the matrix')
        U.NP, s, V.NP = np.linalg.svd(A.NP, full_matrices=True)
        U.NP = U.NP[0:m,0:l]
        U.Array =U.NP.tolist()
        S.NP = np.diag(s[0:l])
        S.Array =S.NP.tolist()
        V.NP = V.NP[0:l,0:n]
        V.Array =V.NP.tolist()
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
                U.Array =U.MP.tolist()
                S.MP = mpmath.diag(s[0:l])
                S.Array =S.MP.tolist()
                V.MP = V.MP[0:n,0:l]
                V.Array =V.MP.tolist()
        elif RealOnly == True:
            m = A.MP.rows
            n = A.MP.cols
            if (l == None):
                l = min([m,n])
            if l > m or l > n:
                raise NLAMathError('Error you should cut more element then the size of the matrix')
            else:
                U.MP,s,V.MP = mpmath.svd_r(A.MP)
                U.MP = U.MP[0:m,0:l]
                U.Array =U.MP.tolist()
                S.MP = mpmath.diag(s[0:l])
                S.Array =S.MP.tolist()
                V.MP = V.MP[0:n,0:l]
                V.Array =V.MP.tolist()
        return U,S,V
def QR(A,Reduced=False):
    if A.onlyNP == True:
        m = A.NP.shape[0]
        n = A.NP.shape[1]
        Q = Matrix([[0] * n for _ in range(m)])
        R = Matrix([[0] * n for _ in range(m)])
        Q.NP,R.NP = np.linalg.qr(A.NP)
        if Reduced==True:
            while (R.NP[-1,:] == np.zeros(n)).all():
                Q.NP=Q.NP[0:-1,:]
                R.NP=R.NP[0:-1,:]
        Q.Array=Q.NP.tolist()
        R.Array=R.NP.tolist()
        return Q,R
    elif A.onlyMP == True:
        m = A.MP.rows-1
        n = A.MP.cols
        Q = Matrix([[0] * n for _ in range(m)],onlyMP=True)
        R = Matrix([[0] * n for _ in range(m)],onlyMP=True)
        Q.MP,R.MP = mpmath.qr(A.MP)
        if Reduced==True:
            while (mphveceq(R.MP[m,:],mpmath.zeros(1,n))):
                Q.MP=Q.MP[0:-1,:]
                R.MP=R.MP[0:-1,:]
                m = m-1
        Q.Array=Q.MP.tolist()
        R.Array=R.MP.tolist()
        return Q,R
    else:
        Q = Matrix([[0] * n for _ in range(m)],onlyNP=False)
        R = Matrix([[0] * n for _ in range(m)],onlyNP=False)
        Q.NP,R.NP = np.linalg.qr(A.NP)
        Q.MP,R.MP = mpmath.qr(A.MP)
        if Reduced==True:
            while (R.NP[-1,:] == np.zeros(n)).all():
                Q.NP=Q.NP[0:-1,:]
                R.NP=R.NP[0:-1,:]
            while (mphveceq(R.MP[m,:],mpmath.zeros(1,n))):
                Q.MP=Q.MP[0:-1,:]
                R.MP=R.MP[0:-1,:]
                m = m-1
        Q.Array=Q.NP.tolist()
        R.Array=R.NP.tolist()
        Q.Array=Q.MP.tolist()
        R.Array=R.MP.tolist()
        return Q,R

def LU(A,Pivoting=False):
    U = A
    m = A.rows
    L = Core.MatrixGallery("I",A.rows,A.cols)
    P = Core.MatrixGallery("I",A.rows,A.cols)
    if A.rows != A.cols:
        print("LU Factorization only with square matrix (I'm lazy not yet implemented not square)!")
    if A.Sparse == True:
        print("LU Factorization only works with dense matrix please convert it first")
    if A.onlyMP == True:
        print("LU Factorization not yet implemented for high precision")
    if Pivoting:
        if A.onlyNP == True:
            U = copy(U.NP)
            A = copy(A.NP)
            L = copy(L.NP)
            P = copy(P.NP)
            for k in range(m-1):
                M = [np.absolute(U[i,k]) for i in range(k,m)]
                i = np.argmax(M)+k
                #print(M)
                #print(U[i,k])

                tmp1 = copy(U[k,k:m])
                tmp2 = copy(U[i,k:m])
                U[k,k:m] = copy(tmp2)
                U[i,k:m] = copy(tmp1)
                del tmp1, tmp2

                tmp1 = copy(L[k,0:k])
                tmp2 = copy(L[i,0:k])
                L[k,0:k] = copy(tmp2)
                L[i,0:k] = copy(tmp1)
                del tmp1, tmp2

                tmp1 = copy(P[k,:])
                tmp2 = copy(P[i,:])
                P[k,:] = copy(tmp2)
                P[i,:] = copy(tmp1)
                del tmp1, tmp2

                for j in range(k+1,m):
                    L[j,k]=copy(U[j,k]/U[k,k])
                    U[j,k:m]=U[j,k:m]-L[j,k]*U[k,k:m]
            U = Matrix(U.tolist(),onlyNP=True)
            L = Matrix(L.tolist(),onlyNP=True)
            P = Matrix(P.tolist(),onlyNP=True)
        return [P.T(),L,U]
    else:
        for i in range(A.cols-1):
            for j in range(i+1,A.cols):
                L[j,i]=U[j,i]/U[i,i]
                if A.onlyNP == True:
                    #print(U.NP[j,i:A.cols])
                    U.NP[j,i:A.cols]=U.NP[j,i:A.cols]-U.NP[i,i:A.cols]*L[j,i]
                    U = Matrix(U.NP.tolist(),onlyNP=True)
        return [L,U]
