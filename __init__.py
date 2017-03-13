import numpy as np #Import Numeric Library
import mpmath #Import Floating Point Arrithmetic
import sys

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
    def __init__(self,Array,onlyMP=False,onlyNP=True):
        if onlyMP == False and onlyNP == False:
            raise NLALogicError('Logic Error Only one between onlyMP and onlyNP could be false.')
        elif onlyMP == True:
            self.onlyMP=True
            self.onlyNP=False
            self.MP=mpmath.mp.matrix(Array)
            self.NP = None
        elif onlyNP == True:
            self.onlyNP=True
            self.onlyMP=False
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
        elif onlyNP == True:
            return str(self.NP)
        else:
            return str(self.MP)+'\n'+str(self.NP)

dot = InfixProduct(np.dot)

from .Core import *
from .Errors import *
from .Factorization import *
