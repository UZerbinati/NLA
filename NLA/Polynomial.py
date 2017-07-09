from __init__ import *
from Core import *

def Poly(K):
    return np.polynomial.polynomial.Polynomial(K)
def PolyVal(x,K):
    return np.polyval(K,x)
def PolyFit(x,y,n):
    return np.polyfit(x, y, n)
def Domain(a,b,n=100):
    return np.linspace(a,b,n)
def Roots(k):
    return np.roots(k)
def ChPoly(H):
    return np.poly(H.NP)
def ChPolyVal(H,x):
    I = Core.MatrixGallery("I",H.rows,H.cols)
    A = H - I*x
    return np.linalg.det(A.NP)
def ArnoldiLamniscates(H,A):
    print(H.NP.shape)
    B=A[:H.rows,:H.cols]
    print(B.NP.shape)
    q=Core.MatrixGallery("rvector",H.NP.shape[0])
    print(q)
    C = Core.norm2(q*ChPolyVal(H,B))/Core.norm2(q)
    print("C",C)
    print(ChPoly(H))
    return np.polysub(ChPoly(H),[-C])
@np.vectorize
def wilkinson(x,order):
    p = np.prod(np.array([x - i for i in range(1, order)]))
    return p
