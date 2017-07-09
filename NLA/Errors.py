class NLAMathError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
class NLALogicError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
class NLATypeError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
class NLANotImplemented(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
