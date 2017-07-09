import numpy as np #Import Numeric Library
import mpmath #Import Floating Point Arrithmetic
import sys
import scipy.sparse as sparse
from copy import *

class InfixProduct(np.ndarray):
    def __new__(cls, function):
        obj = np.ndarray.__new__(cls, 0)
        obj.function = function
        return obj
    def __array_finalize__(self, obj):
        if obj is None: return
        self.function = getattr(obj, 'function', None)
    def __rmul__(self, other):
        return InfixProduct(lambda x,self=self,
            other=other:self.function(other,x))
    def __mul__(self, other):
        return self.function(other)
    def __call__(self, value1, value2):
        return self.function(value1, value2)

class Infix(object):
    def __init__(self, func):
        self.func = func
    def __or__(self, other):
        return self.func(other)
    def __ror__(self, other):
        return Infix(partial(self.func, other))
    def __call__(self, v1, v2):
        return self.func(v1, v2)
class Matrix:
    def __init__(self,Array,onlyMP=False,onlyNP=True,Sparse=False):
        self.Array=Array
        if type(Array[0])== int or type(Array[0])== float:
            self.rows=len(Array)
            self.cols=1
        elif type(Array[0]) == np.float64:
            self.rows=len(Array)
            self.cols=1
        elif Array.__class__.__name__ == 'csr_matrix':
            #print("INN")
            self.rows = Array.shape[0]
            self.cols = Array.shape[1]
        else:
            self.rows=len(Array)
            self.cols=len(Array[0])
        if Sparse == True:
            self.onlyNP=False
            self.onlyMP=False
            self.Sparse=True
            self.S=sparse.csr_matrix(Array)
            self.NP=None
            self.MP=None
        elif onlyMP == True:
            self.onlyMP=True
            self.onlyNP=False
            self.Sparse=False
            self.MP=mpmath.mp.matrix(Array)
            self.NP = None
        elif onlyNP == True:
            self.onlyNP=True
            self.onlyMP=False
            self.Sparse=False
            self.NP=np.array(Array)
            self.MP=None
        else:
            self.onlyMP=False
            self.onlyNP=False
            self.NP=np.array(Array)
            self.MP=mpmath.mp.matrix(Array)
    def __str__(self):
        if self.onlyMP == True:
            return str(self.MP)
        elif self.onlyNP == True:
            return str(self.NP)
        elif self.Sparse==True:
            return str(self.S)
        else:
            return str(self.MP)+'\n'+str(self.NP)
    def __len__(self):
        return self.cols
    def __getitem__(self,key):
        if type(key[0])==int and type(key[1])==int:
            return self.Array[key[0]][key[1]]
        elif type(key[0])==int and type(key[1])!=int:
            return Matrix([self.Array[key[0]][key[1]]],self.onlyMP,self.onlyNP)
        elif type(key[0])!=int and type(key[1])==int:
            return Matrix([row[key[1]] for row in self.Array[key[0]]],self.onlyMP,self.onlyNP)
        else:
            return Matrix([row[key[1]] for row in self.Array[key[0]]],self.onlyMP,self.onlyNP)
    def __setitem__(self,key,value):
        if value.__class__.__name__ == 'Matrix':
            if type(key[0])==int and type(key[1])==int:
                Array = self.Array
                Array[key[0]][key[1]] = value.Array
                self.Array = Array
                if self.onlyMP == True:
                    self.onlyMP=True
                    self.onlyNP=False
                    self.Sparse=False
                    self.MP=mpmath.mp.matrix(Array)
                    self.NP = None
                elif self.onlyNP == True:
                    self.onlyNP=True
                    self.onlyMP=False
                    self.Sparse=False
                    self.NP=np.array(Array)
                    self.MP=None
                else:
                    self.onlyMP=False
                    self.onlyNP=False
                    self.Sparse=False
                    self.NP=np.array(Array)
                    self.MP=mpmath.mp.matrix(Array)
            elif type(key[0])==int and type(key[1])!=int:
                Array = self.Array
                Array[key[0]][key[1]] = value.Array
                self.Array = Array
                if self.onlyMP == True:
                    self.onlyMP=True
                    self.onlyNP=False
                    self.Sparse=False
                    self.MP=mpmath.mp.matrix(Array)
                    self.NP = None
                elif self.onlyNP == True:
                    self.onlyNP=True
                    self.onlyMP=False
                    self.Sparse=False
                    self.NP=np.array(Array)
                    self.MP=None
                else:
                    self.onlyMP=False
                    self.onlyNP=False
                    self.Sparse=False
                    self.NP=np.array(Array)
                    self.MP=mpmath.mp.matrix(Array)
            elif type(key[0])!=int and type(key[1])==int:
                Array = self.Array
                i = 0
                data = []
                for row in Array[key[0]]:
                    row[key[1]] = value.Array[i]
                    i = i+1
                    data.append(copy(row))
                self.Array = data
                if self.onlyMP == True:
                    self.onlyMP=True
                    self.onlyNP=False
                    self.Sparse=False
                    self.MP=mpmath.mp.matrix(Array)
                    self.NP = None
                elif self.onlyNP == True:
                    self.onlyNP=True
                    self.onlyMP=False
                    self.Sparse=False
                    self.NP=np.array(Array)
                    self.MP=None
                else:
                    self.onlyMP=False
                    self.onlyNP=False
                    self.Sparse=False
                    self.NP=np.array(Array)
                    self.MP=mpmath.mp.matrix(Array)
                #return [row[key[1]] for row in self.Array[key[0]]]
            else:
                raise(NLANotImplemented('Set block value of a matrix has\'t been implemented.'))
        else:
            if type(key[0])==int and type(key[1])==int:
                Array = self.Array
                Array[key[0]][key[1]] = value
                self.Array = Array
                if self.onlyMP == True:
                    self.onlyMP=True
                    self.onlyNP=False
                    self.MP=mpmath.mp.matrix(Array)
                    self.NP = None
                elif self.onlyNP == True:
                    self.onlyNP=True
                    self.onlyMP=False
                    self.NP=np.array(Array)
                    self.MP=None
                else:
                    self.onlyMP=False
                    self.onlyNP=False
                    self.NP=np.array(Array)
                    self.MP=mpmath.mp.matrix(Array)
            elif type(key[0])==int and type(key[1])!=int:
                Array = self.Array
                Array[key[0]][key[1]] = value
                self.Array = Array
                if self.onlyMP == True:
                    self.onlyMP=True
                    self.onlyNP=False
                    self.MP=mpmath.mp.matrix(Array)
                    self.NP = None
                elif self.onlyNP == True:
                    self.onlyNP=True
                    self.onlyMP=False
                    self.NP=np.array(Array)
                    self.MP=None
                else:
                    self.onlyMP=False
                    self.onlyNP=False
                    self.NP=np.array(Array)
                    self.MP=mpmath.mp.matrix(Array)
            elif type(key[0])!=int and type(key[1])==int:
                Array = self.Array
                data = [[0] * self.cols for i in range(self.rows)]
                i = 0
                for row in Array[key[0]]:
                    row[key[1]] = value[i]
                    i = i+1
                self.Array = Array
                if self.onlyMP == True:
                    self.onlyMP=True
                    self.onlyNP=False
                    self.MP=mpmath.mp.matrix(Array)
                    self.NP = None
                elif self.onlyNP == True:
                    self.onlyNP=True
                    self.onlyMP=False
                    self.NP=np.array(Array)
                    self.MP=None
                else:
                    self.onlyMP=False
                    self.onlyNP=False
                    self.NP=np.array(Array)
                    self.MP=mpmath.mp.matrix(Array)
                #return [row[key[1]] for row in self.Array[key[0]]]
            else:
                raise(NLANotImplemented('Set block value of a matrix has\'t been implemented.'))
    def __mul__(self,B):
        #print(type(B))
        if type(B)==float or type(B)==int or type(B)==np.float64 or type(B) == np.complex128:
            if self.onlyNP == True:
                return Matrix(np.dot(self.NP,B).tolist(),onlyNP=True)
            elif self.Sparse == True:
                return Matrix(self.S*B,Sparse=True)
        elif self.onlyNP == True and B.onlyNP == True:
            if(type(np.dot(self.NP,B.NP).tolist())==float or type(np.dot(self.NP,B.NP).tolist())==int):
                return(np.dot(self.NP,B.NP).tolist())
            else:
                return Matrix(np.dot(self.NP,B.NP).tolist(),onlyNP=True)
        elif self.onlyMP == True and B.onlyMP == True:
            return Matrix(mpmath.fdot(self.MP,B.MP).tolist(),onlyMP=True)
        elif self.Sparse == True and B.Sparse == True:
            #print(self.S.shape)
            #print(B.S.shape)
            M = self.S*B.S
            if M.shape == (1,1):
                return float(M.todense())
            else:
                return Matrix(self.S*B.S,Sparse=True)
    def __div__(self,B):
        if self.onlyNP == True and (type(B)== int or type(B) == float):
            return Matrix(np.divide(self.NP,B).tolist(),onlyNP=True)
        elif self.onlyNP == True and (type(B)== int or isinstance(B, np.float64)):
            return Matrix(np.divide(self.NP,B).tolist(),onlyNP=True)
        elif self.onlyNP == True and (type(B)== int or type(B) == float):
            f = mpmath.fraction(1,B)
            return Matrix(mpmath.fdot(self.MP,f).tolist(),onlyMP=True)
    def __sub__(self,B):
        if B.__class__.__name__ == "Matrix":
            if self.onlyNP == True and B.onlyNP == True:
                return Matrix(self.NP-B.NP,onlyNP=True)
            elif self.Sparse == True and B.Sparse == True:
                return Matrix(self.S-B.S,Sparse=True)
        elif type(B) == np.float64 or type(B) == int or type(B)== float or type(B) == np.complex128:
            #print("test")
            #print(B)
            if self.onlyNP == True:
                return Matrix(np.subtract(self.NP,B).tolist(),onlyNP=True)
            elif self.Sparse == True:
                print("Not Supported operetion between scalar and sparse matrix.")
    def __add__(self,B):
        #print(type(B))
        #print(self.S.shape)
        #print(B.S.shape)
        if B.__class__.__name__ == "Matrix":
            if self.onlyNP == True and B.onlyNP == True:
                return Matrix(self.NP+B.NP,onlyNP=True)
            elif self.Sparse == True and B.Sparse == True:
                return Matrix(self.S+B.S,Sparse=True)
        elif type(B) == np.float64 or type(B) == int or type(B)== float or type(B) == np.complex128:
            #print("test")
            #print(B)
            if self.onlyNP == True:
                return Matrix(np.add(self.NP,B).tolist(),onlyNP=True)
            elif self.Sparse == True:
                print("Not Supported operetion between scalar and sparse matrix.")
    def T(self):
        if self.Sparse == True:
            return Matrix(sparse.csr_matrix(self.S.transpose()),Sparse=True)
        elif (self.cols==1):
            return Matrix(self.Array,self.onlyMP,self.onlyNP)
        else:
            Array = list(map(list, zip(*self.Array)))
            #print(self.Array)
            return Matrix(Array,self.onlyMP,self.onlyNP)
    def tolist(self):
            Array = []
            for element in self.Array:
                if type(element) == float or type(element)==int:
                    Array.append(element)
                else:
                    Array.append(element[0])
            return Array
    def todense(self,onlyMP=False):
        if self.Sparse == False:
            print ("To be converted into dense must be sparse")
        else:
            A=(self.S.todense()).tolist()
            return Matrix(A,onlyMP)
def floatPrecision(prec):
    mpmath.dps=prec
def mphveceq(v,u):
    m=v.cols
    n=u.cols
    if (m!=n):
        return False
    else:
        for i in range(0,m):
            if v[0,i] != u[0,i]:
                return False
        return True
""" |OLD PRDOUCT AND DIVISON OPERTOR|
Those were the old product and division operetor when those operation
weren't yet included in the Matrix class, using __mul__ and __div__
def Product(A,B):
    if type(A)==float or type(A)==int:
         return Matrix(np.dot(A,B.NP).tolist(),onlyNP=True)
    elif type(B)==float or type(B)==int:
         return Matrix(np.dot(A.NP,B).tolist(),onlyNP=True)
    elif A.onlyNP == True and B.onlyNP == True:
        if(type(np.dot(A.NP,B.NP).tolist())==float or type(np.dot(A.NP,B.NP).tolist())==int):
            return(np.dot(A.NP,B.NP).tolist())
        else:
            return Matrix(np.dot(A.NP,B.NP).tolist(),onlyNP=True)
    elif A.onlyMP == True and B.onlyMP == True:
        return Matrix(mpmath.fdot(A.MP,B.MP).tolist(),onlyMP=True)
def Divide(A,B):
    if A.onlyNP == True and (type(B)== int or type(B) == float):
        return Matrix(np.divide(A.NP,B).tolist(),onlyNP=True)
    elif A.onlyNP == True and (type(B)== int or isinstance(B, np.float64)):
        return Matrix(np.divide(A.NP,B).tolist(),onlyNP=True)
    elif A.onlyNP == True and (type(B)== int or type(B) == float):
        f = mpmath.fraction(1,B)
        return Matrix(mpmath.fdot(A.MP,f).tolist(),onlyMP=True)
dot = InfixProduct(Product)
div = InfixProduct(Divide)
"""
def Slice(start,end,ExtIncluded=False):
    if ExtIncluded == True:
        return range(start, end+1)
    else:
        return range(start, end)
def InColum(v,ONP=True,OMP=False):
    u = Matrix([[0] for _ in range(len(v))],onlyNP=ONP,onlyMP=OMP)
    i=0
    for element in v:
        u[i,0] = element
        i = i+1
    return u
def InRow(v,onlyNP=True,onlyMP=False,Sparse=False):
    #print(v.rows);
    #u = Matrix([0 for _ in range(v.rows)],onlyNP=ONP,onlyMP=OMP)
    #print(v.Array)
    if Sparse ==True:
        #print("-----")
        #print(v)
        #print(v.S.transpose().__class__.__name__)
        tmp = sparse.csr_matrix(v.S.transpose())
        #print("-----")
        #print(tmp)
        return Matrix(tmp,Sparse=True)
    u = [];
    i=0
    for element in v.Array:
        #print(element)
        u.append(element[0])
    u = Matrix(u,onlyNP=onlyNP,onlyMP=onlyMP,)
    return u
from .Core import *
from .Errors import *
from .Factorization import *
from .NOpt import *
from .Polynomial import *
from .Iteractive import *
from .Graph import *
from wathen import *
