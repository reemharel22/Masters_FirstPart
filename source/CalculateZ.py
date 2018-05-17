import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from sympy.solvers import solve
from sympm import Symbol

def f(x):
    z = Symbol('z');
    return solve((cosh(z) - x - 1/z));
#fname = ""
#ff = open(fname,"w");
def main(argv=None):
    if argv is None:
        argv = sys.argv
def main(argv=sys.argv):
    # etc.
