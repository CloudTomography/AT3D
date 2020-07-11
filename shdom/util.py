"""
Utility functions
"""
import numpy as np

def float_round(x):
    """Round a float or np.float32 to a 3 digits float"""
    if type(x) == np.float32:
        x = x.item()
    return round(x,3)

def int_round(x):
    """Round a float or np.float32 to a 3 digits integer by 1000x scaling"""
    return int(np.round(x*1000))

def find_nearest(array, value):
    """Find the nearest element index in an array"""
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def set_pyshdom_path():
    """set path to pyshdom parent directory"""
    import os, shdom
    from pathlib import Path
    os.chdir(str(Path(shdom.__path__[0]).parent))