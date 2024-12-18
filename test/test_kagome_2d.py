import matplotlib.pyplot as plt
import numpy as np
from ansys.mapdl.core import launch_mapdl
from ansys.mapdl import reader as pymapdl_reader
import pyvista as pv

from functions import lattice_tensile_test


# Start an MAPDL instance
try:
    mapdl = launch_mapdl(start_instance=False, port=50052)
except:
    mapdl = launch_mapdl(start_instance=True, port=50052)

lattice_tensile_test(mapdl, mesh_size=0.3, n_cell=5, relative_density=0.5, three_dim=False, filename='2d')
mapdl.exit()