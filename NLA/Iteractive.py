from __init__ import *
from tqdm import *
from copy import *
import Core

def CG(A,b,maxit,tol,sequence=False):
    """
    _________
    ! WIRED !
    _________
    If the sequence option is used the array of matrix that the function reutrn,
    is made by matrix of the same order of A having as rows the vector we are
    interested in (only with dense Matrix).
    """
    x = []
    r = []
    p = []
    s = []
    x.append(Matrix([[0] for _ in range(A.rows)],onlyNP=A.onlyNP,onlyMP=A.onlyMP,Sparse=A.Sparse))
    r.append(InRow(b,Sparse=A.Sparse))
    p.append(r[0])
    i = 1
    while(True):
        #print(p[i-1].S.shape)
        #print(float((p[i-1]*A*p[i-1].T()).S.todense()))
        alpha = Core.norm2(r[i-1])**2/(p[i-1]*A*p[i-1].T())
        x.append(x[i-1]+p[i-1].T()*alpha)
        #print("------")
        #print((A*p[i-1].T())*alpha)
        #print("*******")
        #print(r[i-1].S.shape)
        #print(((A*p[i-1].T())*alpha).T().S.shape)
        r.append(r[i-1]-((A*p[i-1].T())*alpha).T())
        #print(r)
        beta = Core.norm2(r[i])**2/Core.norm2(r[i-1])**2
        p.append(r[i]+p[i-1]*beta)
        i = i+1
        if A.Sparse == False:
            if (i > maxit):
                if(sequence):
                    return [Matrix(ToList(x[i-1][0,:])[0],onlyNP=A.onlyNP,onlyMP=A.onlyMP),i-1]
                else:
                    return [x[i-1][0,:],i-1]
            if (beta <= tol):
                if(sequence):
                    return [x,i-1]
                else:
                    return [Matrix(ToList(x[i-1][0,:])[0],onlyNP=A.onlyNP,onlyMP=A.onlyMP),i-1]
        else:
            if (i > maxit):
                if(sequence):
                    return [x,i-1]
                else:
                    return [x[i-1],i-1]
            if (beta <= tol):
                if(sequence):
                    return [x,i-1]
                else:
                    return [x[i-1],i-1]
def PreConditionCG(A,b,maxit,tol,option,sequence=False,BetaBreak=False):
    """
    ! Not implemented for high precision matrix
    """
    #print(A)
    if option == "Jacobi":
        if A.Sparse == True:
            B = A.S.todense()
        else:
            B = A.NP
        D = np.diag(A.S.todense())
        D = np.diag(1/D)
        P = Matrix(D.tolist(),Sparse=A.Sparse)
    x = []
    r = []
    z = []
    p = []
    s = []
    Beta = []
    B = 0
    x.append(Matrix([[0] for _ in range(A.rows)],onlyNP=A.onlyNP,onlyMP=A.onlyMP,Sparse=A.Sparse))
    r.append(b)
    #print(b.S.shape)
    #print(P.S.shape)
    p.append(P*r[0])
    z.append(p[0])
    Beta.append(-1)
    i = 1
    while(True):
        if (p[i-1].T()*A*p[i-1]==0):
            #print(p[i-1].T()*A*p[i-1])
            if A.Sparse == False:
                if(sequence):
                    return [Matrix(ToList(x[i-1][0,:])[0],onlyNP=A.onlyNP,onlyMP=A.onlyMP),i-1]
                else:
                    return [x[i-1][0,:],i-1]
            else:
                if(sequence):
                    return [x,i-1]
                else:
                    return [x[i-1],i-1]
        #print((r[i-1].T()*z[i-1]))
        #print((p[i-1].T()*A*p[i-1]))
        alpha = (r[i-1].T()*z[i-1])/(p[i-1].T()*A*p[i-1])
        x.append(x[i-1]+p[i-1]*alpha)
        r.append(r[i-1]-(A*p[i-1]*alpha))
        z.append(P*r[i])
        beta = (r[i].T()*z[i])/(r[i-1].T()*z[i-1])
        Beta.append(beta)
        if Beta[i] == Beta[i-1]:
            B = B+1
        p.append(z[i]+p[i-1]*beta)
        i = i+1

        if A.Sparse == False:
            if (i > maxit):
                if(sequence):
                    return [Matrix(ToList(x[i-1][0,:])[0],onlyNP=A.onlyNP,onlyMP=A.onlyMP),i-1]
                else:
                    return [x[i-1][0,:],i-1]
            if (beta <= tol):
                if(sequence):
                    return [x,i-1]
                else:
                    return [Matrix(ToList(x[i-1][0,:])[0],onlyNP=A.onlyNP,onlyMP=A.onlyMP),i-1]
            if B > maxit/10:
                if(sequence):
                    return [x,maxit]
                else:
                    return [Matrix(ToList(x[i-1][0,:])[0],onlyNP=A.onlyNP,onlyMP=A.onlyMP),maxit]
        else:
            if (i > maxit):
                if(sequence):
                    return [x,i-1]
                else:
                    return [x[i-1],i-1]
            if (beta <= tol):
                if(sequence):
                    return [x,i-1]
                else:
                    return [x[i-1],i-1]
            if B > maxit/10:
                if(sequence):
                    return [x,maxit]
                else:
                    return [x[i-1],maxit]

def Arnoldi (A,N):
    """
    ! ONLY WORK WITH NP MATRIX !
    """
    A = A.NP
    H = np.zeros((N,N))
    Q = np.zeros(A.shape)
    n = A.shape[0]
    b = np.random.rand(n)
    Q[:,0] = b/np.linalg.norm(b)
    for i in range(0,N):
        for j in range(0,i+1):
            H[j,i] = np.dot(Q[:,j],np.dot(A,Q[:,i]))
        if i!=N-1:
            q = np.dot(A,Q[:,i]) - np.dot(Q[:,0:i+1],H[0:i+1,i])
            H[i+1,i] = np.linalg.norm(q)
            Q[:,i+1] = q/H[i+1,i]
    Q = Matrix(Q)
    H = Matrix(H)
    return [Q,H]
def RitzSpectrum(A,n):
    """
    ! ONLY WORK WITH NP MATRIX !
    """
    Q, H = Arnoldi(A,n)
    return np.linalg.eig(H.NP)[0]
def Redi(A,x):
    A = x.T()*A*x
    B = x.T()*x
    return (A/B)
def PowerIteration(A,Shift=False,mu=0,ShiftUpdate=False,maxit=1000000,Succession=False):
        if A.rows != A.cols:
            print("Power Iteration only with square matrix (I'm lazy not yet implemented not square)!")
        if A.Sparse == True:
            print("Power Itaration only works with dense matrix please convert it first")
        if A.onlyMP == True:
            print("Power Iteration with shift not yet implemented for high precision")

        if ShiftUpdate==True:
            #print("A")
            v=[]
            l=[]
            vec = Core.MatrixGallery("rvector",A.rows)
            vec = InRow(vec)
            v.append(vec/Core.norm2(vec))
            l.append(v[0].T()*A*v[0])
            it = 1
            for k in tqdm(range(1,maxit)):
                B = (A-Core.MatrixGallery("I",A.rows)*l[k-1])
                Bi = Matrix(np.linalg.inv(B.NP).tolist(),onlyNP=True)
                w = Bi*v[k-1]
                v.append(w/Core.norm2(w))
                l.append((1/Redi(Bi,v[-1]))+l[k-1])
                it=it+1
            if (Succession):
                return [v,l,it]
            else:
                return [v[-1],l[-1],it]
        elif Shift==True:
            #print("B")
            v=[]
            vec = Core.MatrixGallery("rvector",A.rows)
            vec = InRow(vec)
            v.append(vec/Core.norm2(vec))
            l = []
            B = (A-Core.MatrixGallery("I",A.rows)*mu)
            #print(B)
            B = np.linalg.inv(B.NP)
            #print(B)
            B = Matrix(B.tolist(),onlyNP=True)
            it = 1
            for k in range(1,maxit):
                w = B*v[k-1]
                v.append(w*(1/Core.norm2(w)))
                l.append((1/Redi(B,v[-1]))+mu)
                it=it+1
            if (Succession):
                return [v,l,it]
            else:
                return [v[-1],l[-1],it]
        elif Shift == False:
            #print("C")
            v=[]
            vec = Core.MatrixGallery("rvector",A.rows)
            vec = InRow(vec)
            v.append(vec/Core.norm2(vec))
            #print(Core.norm2(v[0]))
            l = []
            it = 1
            for k in range(1,maxit):
                w = A*v[k-1]
                v.append(w/Core.norm2(w))
                l.append(v[k].T()*A*v[k])
                it=it+1
            if (Succession):
                return [v,l,it]
            else:
                return [v[-1],l[-1],it]
