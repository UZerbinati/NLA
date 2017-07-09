from __init__ import *
import Core
import Factorization

def GS(A,QR=False):
    V = Matrix([[0]*A.cols]*A.rows)
    Q = Matrix([[0]*A.cols]*A.rows)
    R = Matrix([[0]*A.cols]*A.rows)
    for j in Slice(0,A.cols):
        V[0:4,j] = A[0:4,j]
        for i in Slice(0,j-1,True):
            R[i,j]=Q[0:4,i] * A[0:4,j]
            D = Q[0:4,i] * R[i,j]
            V[0:4,j]=V[0:4,j] - D
        R[j,j]=Core.norm2(V[0:4,j])
        Q[0:4,j]= V[0:4,j] / R[j,j]
    if QR == True:
        return Q,R
    else:
        return Q
def MGS(A,QR=False):
    print(A)
    V = Matrix([[0]*A.cols]*A.rows)
    Q = Matrix([[0]*A.cols]*A.rows)
    R = Matrix([[0]*A.cols]*A.rows)
    for i in Slice(0,A.cols):
        V[0:4,i]=A[0:4,i]
    for i in Slice(0,A.cols):
        R[i,i]=Core.norm2(V[0:4,i])
        #print(R[i,i])
        Q[0:4,i]=V[0:4,i]/R[i,i]
        #print(Q[0:4,i])
        for j in Slice(i+1,A.cols):
            R[i,j]=Q[0:4,i].T()*V[0:4,j]
            V[0:4,j]=V[0:4,j]-Q[0:4,i]*R[i,j]
    if QR == True:
        return Q,R
    else:
        return Q
def isUT(A):
    UT = True
    #print(A.rows,A.cols)
    for j in Slice(0,A.cols):
        for i in Slice(0,A.rows):
            if i < j:
                #print(i,j)
                if (A[j,i]==0):
                    UT = UT and True
                else:
                    UT = UT and False
    return UT
def isD(A):
    D = True
    #print(A.rows,A.cols)
    for j in Slice(0,A.cols):
        for i in Slice(0,A.rows):
            if i != j:
                #print(i,j)
                if (A[j,i]==0):
                    D = D and True
                else:
                    D = D and False
    return D
def LinSys(A,b):
    x = Matrix([[0] for _ in range(A.rows)],onlyNP=A.onlyNP,onlyMP=A.onlyMP)
    #print(x)
    #print(A.cols)
    #print(A.rows)
    #print(b.cols)
    #print(b.rows)
    if A.cols == A.rows:
        if isD(A):
            for i in Slice(0,A.cols):
                for j in Slice(0,A.rows):
                    if i == j:
                        #print(b[2,0])
                        x[i,0] = b[i,0]/A[i,j]
                        #print(x[0,i])
            return x
        else:
            Q,R = Factorization.QR(A)
            x = Matrix([[0] for _ in range(A.rows)],onlyNP=A.onlyNP,onlyMP=A.onlyMP)
            #print("----------|Matrix Q|-------")
            #print(Q)
            #print("----------|Matrix R|-------")
            #print(R)
            """
            Ax=b -> QRx=b -> Rx = bQ^T
            """
            b = Q.T()*b
            #print("----------|Vector bQ^T|----------")
            #print(b)
            if isUT(R):
                #print("[*]Proceding with a backword substitution ...")
                for i in Slice(0,b.rows):
                    d=0
                    i=b.rows-i-1
                    """
                    x_{i}=\frac{b_i-\sum_{i=1}^{n}a_{i,j}x_{j}}{a_{i,i}}
                    """
                    for j in Slice(0,R.cols):
                        d=d+R[i,j]*x[j,0]
                    #print(j,i)
                    #print(R[1,1])
                    x[i,0]=(b[i,0]-d)/R[i,i]
                #print("----------|Solution|-------")
                return(x)

            else:
                print("[Error]Something went wrong during QR-Factorization")
    else:
        Q,R = Factorization.QR(A)
        x = Matrix([[0] for _ in range(A.rows)],onlyNP=A.onlyNP,onlyMP=A.onlyMP)
        #print("----------|Matrix Q|-------")
        #print(Q)
        #print("----------|Matrix R|-------")
        #print(R)
        """
        Ax=b -> QRx=b -> Rx = bQ^T
        """
        b = Q.T()*b
        #print("----------|Vector bQ^T|----------")
        #print(b)
        if isUT(R):
            #print("[*]Proceding with a backword substitution ...")
            for i in Slice(0,b.rows):
                d=0
                i=b.rows-i-1
                """
                x_{i}=\frac{b_i-\sum_{i=1}^{n}a_{i,j}x_{j}}{a_{i,i}}
                """
                for j in Slice(0,R.cols):
                    d=d+R[i,j]*x[j,0]
                #print(j,i)
                #print(R[1,1])
                x[i,0]=(b[i,0]-d)/R[i,i]
            #print("----------|Solution|-------")
            return(x)

        else:
            print("[Error]Something went wrong during QR-Factorization")
def Vander(x,m,onlyMP=False,onlyNP=True):
    A = Matrix([[0] * m for _ in Slice(0,len(x))],onlyMP=onlyMP,onlyNP=onlyNP)
    for i in range(0,len(x)):
        for j in range(0,m):
            if type(x) == np.ndarray:
                #print(type(x))
                A[i,j]=np.power(x[i],j)
            if type(x) == list:
                A[i,j]=x[i]**j
    return(A)
def flipLR(M):
    Array = M.Array
    #print(M)
    A = Matrix([[0]*M.cols for _ in range(M.rows)],onlyNP=M.onlyNP,onlyMP=M.onlyMP)
    for i in Slice(0,A.cols):
        A[:,i]=M[:,M.cols-1-i]
    return (A)
def flipUD(M):
    #print(M)
    if M.cols > 1:
        A = Matrix([[0]*M.cols for _ in range(M.rows)],onlyNP=M.onlyNP,onlyMP=M.onlyMP)
        B = Matrix([[0]*M.cols for _ in range(M.rows)],onlyNP=M.onlyNP,onlyMP=M.onlyMP)
        for i in Slice(0,A.rows):
            A[i,:]=M[M.cols-1-i,:]
        #print(A)
        for i in range(0,A.rows):
            B[i,:]=A[i,0]
        return (B)
    elif M.cols == 1:
        A = Matrix([[0]*M.cols for _ in range(M.rows)],onlyNP=M.onlyNP,onlyMP=M.onlyMP)
        #print(M[2,0])
        for i in Slice(0,M.rows):
            #print(i,M.rows-i)
            A[i,0] = M[M.rows-i-1,0]
        #print(A)
        return(A)
def FilterVector(v,k):
    """
    The vector must be horizontal to be filtered !
    """
    Array = []
    for element in v.Array:
        if element[0] != k:
            Array.append(element[0])
    #print("...",Array)
def ToList(A):
    Array = [];
    if A.onlyNP == True:
        Array= A.NP.tolist()
    elif A.onlyMP == True:
        Array= A.MP.tolist()
    else:
        Array.append(A.NP.tolist())
        Array.append(A.MP.tolist())
    return Array
