from __init__ import *
import Factorization
import random
import scipy
import scipy.sparse as sparse
# Redefine norms using np.array allowing to process
# as well array not only defined using numpy.
def norm2(A,HighPrec=False):
    if A.onlyMP == True:
        if A.MP.rows == 1:
            return mpmath.norm(A.MP,2)
        elif A.MP.cols == 1:
            return mpmath.norm(A.MP,2)
        else:
            U, s, V = mpmath.svd_c(A.MP)
            return s[0]
    elif A.Sparse == True:
        return np.linalg.norm(A.S.todense())
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
def MConditioning(A):
    if A.onlyMP == True:
        print "Not yet implemente for high precision matrix"
    elif A.Sparse == True:
        return np.linalg.cond(A.S.todense())
    else:
        return np.linalg.cond(A.NP)
def MatrixRank(A):
    if A.onlyNP == True:
        return np.linalg.matrix_rank(A)
    elif A.onlyMP == True:
        #TO BE IMPLEMENTED
        pass
def frobeniusNorm(A):
    if A.onlyMP == True:
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
def MatrixFromColumns(q):
    A = []
    for i in range(len(q)):
        A.append(q[i].tolist())
    return Matrix(A,onlyNP=q[0].onlyNP,onlyMP=q[0].onlyMP,Sparse=q[0].Sparse).T()
def MatrixGallery(option,n,m=0,OnlyNP=True,onlyMP=False,StartingPoint=0,SparseDensity=0.01,Sigma=0):
    """
    _____________________________________________
    |           DENSE MATRIX GALLERY            |
    _____________________________________________
    """
    ONP = OnlyNP
    OMP = onlyMP
    A = []
    if option == "Zeros" or option=="zeros":
        if m == 0:
            A = Matrix([[0.0]*n for _ in range(n)],onlyNP=ONP,onlyMP=OMP)
        else:
            A = Matrix([[0.0]*m for _ in range(n)],onlyNP=ONP,onlyMP=OMP)
    elif option == "Random Diag" or option == "random diag" or option=="random_diag" or option == "rdiag":
        if m == 0:
            A = Matrix([[0.0]*n for _ in range(n)],onlyNP=ONP,onlyMP=OMP)
            m = n
        else:
            A = Matrix([[0.0]*m for _ in range(n)],onlyNP=ONP,onlyMP=OMP)
        for i in range(n):
            for j in range(m):
                if j == (i+StartingPoint):
                    A[i,j]=random.random()
    elif option == "Random Tridiagonal" or option == "r3diag":
            if m == 0:
                A = Matrix([[0.0]*n for _ in range(n)],onlyNP=ONP,onlyMP=OMP)
                m = n
            else:
                A = Matrix([[0.0]*m for _ in range(n)],onlyNP=ONP,onlyMP=OMP)
            X = [-1,0,1]
            for x in X:
                for i in range(n):
                    for j in range(m):
                        if j == (i+x):
                            A[i,j]=random.random()
    elif option == "Square Random" or option == "qrand":
                    if m == 0:
                        R = np.random.rand(n,n)
                        A = Matrix(R.tolist(),onlyNP=ONP,onlyMP=OMP)
                    else:
                        R = np.random.rand(n,m)
                        A = Matrix(R.tolist(),onlyNP=ONP,onlyMP=OMP)
    elif option == "qrandsigma" or option == "Square Random Sigma":
        mu, sigma = 0, Sigma # mean and standard deviation
        s = np.random.normal(mu, sigma, (n,n))
        A = Matrix(s,onlyNP=ONP, onlyMP=OMP)
        #print(s)
    elif option == "Eye" or option == "eye" or option=="I":
        if m == 0:
            A = Matrix([[0.0]*n for _ in range(n)],onlyNP=ONP,onlyMP=OMP)
            m = n
        else:
            A = Matrix([[0.0]*m for _ in range(n)],onlyNP=ONP,onlyMP=OMP)
        for i in range(n):
            for j in range(m):
                if j == (i+StartingPoint):
                    A[i,j]=1.0
    elif option == "Ones" or option == "ones" or option=="1v":
        if m == 0:
            A = Matrix([[1]*n for _ in range(n)],onlyNP=ONP,onlyMP=OMP)
            m = n
        else:
            A = Matrix([[1]*m for _ in range(n)],onlyNP=ONP,onlyMP=OMP)
    elif option == "Random Vector" or option=="rvector":
        A = Matrix([[0] for _ in range(n)],onlyNP=ONP,onlyMP=OMP)
        for i in range(n):
            A[i,0] = random.random()
    elif option == "DenseWathen" or option=="dense_wathen":
        A,_ = wathen_csr(n,m)
        A = A.S.todense()
        return Matrix(A.tolist(),onlyNP=ONP,onlyMP=OMP)
    elif option == "Square Sparse Random Diag" or option == "square sparse random diag" or option=="square_sparse_random_diag" or option == "qsrdiag":
        if m == 0:
            diag = np.random.rand(n-np.absolute(StartingPoint))
            diag = np.diag(diag,StartingPoint)
            #print(diag)
            A = Matrix(diag,Sparse=True)
    elif option == "Square Sparse Random Tridiag" or option == "square sparse random tridiag" or option=="square_sparse_random_tridiag" or option == "qsr3diag":
        if m == 0:
            diag = np.random.rand(n)
            diag = np.diag(diag,0)
            #print(diag)
            A = Matrix(diag,Sparse=True)
            diag = np.random.rand(n-1)
            diag = np.diag(diag,1)
            #print(diag)
            B = Matrix(diag,Sparse=True)
            diag = np.random.rand(n-1)
            diag = np.diag(diag,-1)
            #print(diag)
            C = Matrix(diag,Sparse=True)
            A=A+B+C
    elif option == "Wathen" or option=="wathen":
        A,_ = wathen_csr(n,m)
    elif option == "sdones":
        A=Matrix(np.ones(n),Sparse=True)
    elif option == "Sparse Random Vector" or option == "sparse random vector" or option=="sparse_random_vector" or option == "srvector":
        #print(sparse.rand(n,1,density=SparseDensity,format="csr").__class__.__name__)
        A=Matrix(sparse.rand(n,1,density=SparseDensity,format="csr"),Sparse=True)
    else:
        print("This matrix isn't in the gallery.")
        return None
    if m==0:
        A.cols=n
        A.rows=n
    else:
        if option == "Wathen" or option=="wathen":
            A.cols = 3*n*m+2*n+2*m+1
            A.rows = 3*n*m+2*n+2*m+1
        else:
            A.cols=m
            A.rows=n
    return A
