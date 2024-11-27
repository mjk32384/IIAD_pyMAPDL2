import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from ansys.mapdl.core import launch_mapdl
from ansys.mapdl import reader as pymapdl_reader
import pyvista as pv
import time
import traceback
import os

from create_geometry import kagome_horizontal

def kagome_tensile_test(mapdl, mesh_size):
    mapdl.clear()
    mapdl.prep7()

    print(f"Simulation Start. Mesh size: {mesh_size}")
    design_width = 60/1000
    design_height = 25/1000
    total_width = 144/1000
    total_z = 5/1000
    wall_thickness = kagome_horizontal(mapdl, n_cell=5, relative_density=0.5)

    # 물성치 입력
    stress_strain = pd.read_excel('PLA stress-strain.xlsx')
    stress_strain = stress_strain.iloc[:8]  # due to limitations of MISO model
    stress_strain['y'] = stress_strain['y (MPa)']*1e6
    stress_strain['x'] = stress_strain['x']/100
    UTS = stress_strain['y'].max()
    Youngs_Modulus = stress_strain['y'][1]/stress_strain['x'][1]
    max_plastic_strain = stress_strain['x'].iloc[-1]

    mapdl.mp('EX', 1, Youngs_Modulus)   # Young's Modulus in Pascals
    mapdl.mp('PRXY', 1, 0.33)   # Poisson's Ratio

    # Non linear 부분 정의
    stress_strain['x'] = stress_strain['x'] - stress_strain['y']/Youngs_Modulus
    mapdl.tb('PLASTIC', 1, '', 'MISO')
    for i in range(len(stress_strain)):
        mapdl.tbpt(i, stress_strain['x'][i], stress_strain['y'][i])

    # Meshing
    line_design_part = mapdl.lsel("S", "LOC", "X", 0, design_width)
    mapdl.lsel("NONE")
    for line in line_design_part:
        mapdl.lsel("A", "LINE", vmin=line, vmax=line)
    mapdl.lesize("ALL", wall_thickness*mesh_size, kforc=1)
    mapdl.lsel("ALL")

    # Define mesh size
    mapdl.esize(design_height/3)  # Default size

    mapdl.mopt("EXPND", 0.7)  # Decrease the area mesh expansion.
    mapdl.et(1, 'SOLID187')

    # Mesh the volume
    mapdl.vmesh("all")
    print("Preprocessing Finished")


    # Solution
    mapdl.solution()
    mapdl.antype('STATIC')
    # mapdl.nlgeom('ON')  # Enable large deformation effects

    total_displacement = 0.02  # 20 mm total displacement
    number_of_steps = 200
    displacement_increment = total_displacement / number_of_steps

    # mapdl.nropt('UNSYM')  # Use unsymmetric matrix if convergence issues occur
    mapdl.autots('ON')    # Enable automatic time stepping
    # mapdl.nliter(50, 0, 0, 0, 0)  # Set maximum number of equilibrium iterations
    # mapdl.cnvtol('F', 0.1)  # Force convergence tolerance
    # mapdl.nsubst(20, 100, 1)
    # mapdl.deltim(dtime = "AUTOTS", dtmin = 0, dtmax = 1)

    # Initialize arrays to store results
    reaction_forces = []
    applied_displacements = []

    for step in range(1, number_of_steps + 1):
        time_start = time.time()
        current_disp = displacement_increment * step
        mapdl.ddele('ALL', 'UX')  # Remove previous displacement BC
        mapdl.nsel('S', 'LOC', 'X', (total_width + design_width)/2)
        mapdl.d('ALL', 'UX', current_disp)
        mapdl.d("ALL", "UY", 0)
        mapdl.d("ALL", "UZ", 0)
        mapdl.nsel('S', 'LOC', 'X', (-total_width + design_width)/2)
        mapdl.d("ALL", "ALL", 0)
        mapdl.nsel('ALL')
        
        time_solutionstart = time.time()
        # Solve current load step
        try:
            mapdl.solve()
        except:
            with open("result/kagome/error_log.txt", "a") as file:
                file.write(f"An error occurred at mesh size:{mesh_size}, iteration:{step}\n")
                file.write(traceback.format_exc())  # Writes the full stack trace
                file.write("\n")
            print("Solution Error!\n\n")
            break
        time_solutionend = time.time()
        # print(f"Solution time = {int(time_solutionend-time_solutionstart):4}s")

        # Post-process to check for failure criterion
        mapdl.post1()
        mapdl.set('LAST')
        
        # Get nodal principal stress
        mapdl.nsel('ALL')
        principal_stress = np.maximum.reduce([mapdl.post_processing.nodal_principal_stress("1"),
                    mapdl.post_processing.nodal_principal_stress("2"),
                    mapdl.post_processing.nodal_principal_stress("3")])


        mapdl.nsel('ALL')
        mapdl.nsel('S', 'LOC', 'X', (-total_width + design_width)/2)

        reaction_forces_selected = mapdl.prrsol('FX').to_dataframe(columns=["NODE", "FX"])
        reaction_force = abs(reaction_forces_selected['FX'].sum())

        # Get the reaction forces from the result object
        # result = mapdl.result
        # # 이거를 result 안쓰고 하는 법 없나? result 문법을 잘 모르겠음
        # nnum, reaction_force = result.nodal_static_forces(0)
        # mapdl.nsel('ALL')
        # mapdl.nsel('S', 'LOC', 'X', 0)
        # selected_nodes = mapdl.mesh.nnum  # Node numbers in the current selection
        # reaction_force = reaction_force[selected_nodes-1].sum(axis=0)[0]
        
        reaction_forces.append(reaction_force)
        applied_displacements.append(current_disp)

        print(f"Step{step:4} completed. Current displacement : {current_disp*1000:6.2f} Elapsed Time : {int(time.time()-time_start):4}s")

        # Principal stress / UTS

        if (principal_stress > UTS).sum() > 0:
            print("Specimen failed!\n\n")
            break

        mapdl.solution()

    current_dir = os.getcwd()
    df = pd.DataFrame({'displacement': applied_displacements, 'force': reaction_forces})
    df.to_csv(f'{current_dir}/result/kagome/stress_strain_{mesh_size:1.2f}'.replace(".", "_")+'.csv')
    mapdl.save(fname=f"{current_dir}/result/kagome/db_{mesh_size:1.2f}".replace(".", "_"), ext="db")
    mapdl.reswrite(fname = f"{current_dir}/result/kagome/res_{mesh_size:1.2f}".replace(".", "_"))

    mapdl.finish()


# Start an MAPDL instance
try:
    mapdl = launch_mapdl(start_instance=False, port=50053)
except:
    mapdl = launch_mapdl(start_instance=True, port=50053)

for mesh_size in [2, 1.5, 1, 0.7, 0.5, 0.3]:
    kagome_tensile_test(mapdl, mesh_size)

mapdl.exit()