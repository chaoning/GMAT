import numpy as np
from scipy import linalg


def is_positive(x):
    return np.all(linalg.eigvals(x) > 0)
