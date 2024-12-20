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

for n_cell in range(2, 8):
    for relative_density in np.arange(0.15, 0.55, 0.1):
        lattice_tensile_test(mapdl, lattice_type='kagome_sine_h', mesh_size=0.3, n_cell=n_cell, relative_density=relative_density,
                             amplitude_ratio=0.5, three_dim=False, filename=f'n{n_cell}_r{int(relative_density*100)}', external_wall=False)


mapdl.exit()