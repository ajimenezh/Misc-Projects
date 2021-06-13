import numpy as np
from numpy import linalg as LA

def norm(v):
    return LA.norm(v)

def toUnit(v):
    if norm(v) == 0:
        return v
    return v / norm(v)

