import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import traceback
import os

from ansys.mapdl.core import launch_mapdl


def lattice_tensile_test(mapdl, lattice_type, mesh_size, n_cell, relative_density, three_dim=False, filename='', external_wall=True):
    mapdl.clear()
    mapdl.prep7()

    lattice_type_list = ['kagome_h', 'kagome_v', 'kagome_sine_h', 'kagome_sine_v', 'triangle_h', 'triangle_v']
    if not lattice_type in lattice_type_list:
        raise Exception("Invalid lattice type")

    print(f"{'3D' if three_dim else '2D'} {lattice_type} Simulation Start. Mesh Size: {mesh_size}")
    print(f"Name: {filename}")
    design_width = 60/1000
    design_height = 25/1000
    total_width = 144/1000
    total_z = 5/1000
    if lattice_type == 'kagome_h':
        wall_thickness = kagome_horizontal(mapdl, n_cell, relative_density, three_dim=three_dim, external_wall=external_wall)
    elif lattice_type == 'kagome_v':
        pass
    elif lattice_type == 'kagome_sine_h':
        pass
    elif lattice_type == 'kagome_sine_v':
        pass
    elif lattice_type == 'triagnle_h':
        pass
    elif lattice_type == 'triangle_v':
        pass

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
    mapdl.lesize("ALL", wall_thickness*mesh_size, kforc=1)  # 격자 부분의 mesh를 더 작게 함
    mapdl.lsel("ALL")

    # Define default mesh size
    if three_dim:
        mapdl.esize(design_height/3)
    else:
        mapdl.esize(design_height/6)

    mapdl.mopt("EXPND", 0.7)  # Decrease the area mesh expansion.
    mapdl.et(1, 'PLANE')

    if three_dim:
        mapdl.et(1, 'SOLID187')
        mapdl.vmesh("all")
    else:
        # define a PLANE183 element type with thickness
        mapdl.et(1, "PLANE183", kop1=1, kop3=3)
        mapdl.r(1, total_z)  # thickness of 5mm
        # Mesh the volume
        mapdl.amesh("all")
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
        if three_dim:
            mapdl.d("ALL", "UZ", 0)
        mapdl.nsel('S', 'LOC', 'X', (-total_width + design_width)/2)
        mapdl.d("ALL", "ALL", 0)
        mapdl.nsel('ALL')
        
        time_solutionstart = time.time()
        # Solve current load step
        try:
            mapdl.solve()
        except:
            with open(f"result/{lattice_type}_{'wall' if external_wall else 'nowall'}/error_log.txt", "a") as file:
                file.write(f"An error occurred at mesh size:{mesh_size}, iteration:{step}, filename: {filename}\n")
                file.write(traceback.format_exc())  # Writes the full stack trace
                file.write("\n")
            print("Solution Error!\n\n")
            try:
                mapdl = launch_mapdl(start_instance=False, port=50052)
            except:
                mapdl = launch_mapdl(start_instance=True, port=50052)
                return        
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

        print(f"Step{step:3} completed.  Current displacement : {current_disp*1000:5.2f}  Elapsed Time : {int(time.time()-time_start):4}s")

        # Principal stress / UTS

        if (principal_stress > UTS).sum() > 0:
            print("Specimen failed!\n\n")
            break

        mapdl.solution()

    current_dir = os.getcwd()
    df = pd.DataFrame({'displacement': applied_displacements, 'force': reaction_forces})
    df.to_csv(f'{current_dir}/result/{lattice_type}_{"wall" if external_wall else "nowall"}/stress_strain_{filename}_{mesh_size:1.2f}'.replace(".", "_")+'.csv')
    mapdl.save(fname=f'{current_dir}/result/{lattice_type}_{"wall" if external_wall else "nowall"}/db_{filename}_{mesh_size:1.2f}'.replace(".", "_"), ext="db")
    mapdl.reswrite(fname = f'{current_dir}/result/{lattice_type}_{"wall" if external_wall else "nowall"}/res_{filename}_{mesh_size:1.2f}'.replace(".", "_"))

    mapdl.finish()

def triangle_horizontal(mapdl, n_cell, wall_thickness, external_wall):
    design_width = 60/1000
    design_height = 25/1000
    total_width = 144/1000
    total_z = 5/1000
    external_wall_thickness = wall_thickness

    cell_size = (design_height - 2*5/6*np.sqrt(3)*wall_thickness
              - (n_cell - 1)*np.sqrt(3)*wall_thickness)/n_cell + np.sqrt(3)*wall_thickness
    
    i = 0
    design_part = mapdl.rectng(0, design_width, 0, design_height)
    while True:
        if i%2 == 0:
            for j in range(n_cell):
                k1 = mapdl.k(x = i*cell_size*np.sqrt(3)/2, y = j*cell_size + 5/6*np.sqrt(3)*wall_thickness)
                k2 = mapdl.k(x = i*cell_size*np.sqrt(3)/2, y = (j+1)*cell_size - 1/6*np.sqrt(3)*wall_thickness)
                k3 = mapdl.k(x = (i+1)*cell_size*np.sqrt(3)/2 - 3/2*wall_thickness, y = (j+1/2)*cell_size + 1/3*np.sqrt(3)*wall_thickness)
                hole = mapdl.a(k1, k2, k3)
                try:
                    design_part = mapdl.asba(na1 = design_part, na2 = hole)
                except:
                    pass
                
            for j in range(n_cell-1):
                k1 = mapdl.k(x = (i+1)*cell_size*np.sqrt(3)/2 - wall_thickness, y = (j+1/2)*cell_size + 5/6*np.sqrt(3)*wall_thickness)
                k2 = mapdl.k(x = (i+1)*cell_size*np.sqrt(3)/2 - wall_thickness, y = (j+3/2)*cell_size - 1/6*np.sqrt(3)*wall_thickness)
                k3 = mapdl.k(x = i*cell_size*np.sqrt(3)/2 + 1/2*wall_thickness, y = (j+1)*cell_size + 1/3*np.sqrt(3)*wall_thickness)
                hole = mapdl.a(k1, k2, k3)
                try:
                    design_part = mapdl.asba(na1 = design_part, na2 = hole)
                except:
                    pass
            k1 = mapdl.k(x = i*cell_size*np.sqrt(3)/2 - wall_thickness/2, y = 0)
            k2 = mapdl.k(x = (i+2)*cell_size*np.sqrt(3)/2 - wall_thickness/2, y = 0)
            k3 = mapdl.k(x = (i+1)*cell_size*np.sqrt(3)/2 - wall_thickness/2, y = 1/2*cell_size)
            hole = mapdl.a(k1, k2, k3)
            try:
                design_part = mapdl.asba(na1 = design_part, na2 = hole)
            except:
                pass
            k1 = mapdl.k(x = i*cell_size*np.sqrt(3)/2 - wall_thickness/2, y = design_height)
            k2 = mapdl.k(x = (i+2)*cell_size*np.sqrt(3)/2 - wall_thickness/2, y = design_height)
            k3 = mapdl.k(x = (i+1)*cell_size*np.sqrt(3)/2 - wall_thickness/2, y = design_height - 1/2*cell_size)
            hole = mapdl.a(k1, k2, k3)
            try:
                design_part = mapdl.asba(na1 = design_part, na2 = hole)
            except:
                pass
        else:
            for j in range(n_cell):
                k1 = mapdl.k(x = (i+1)*cell_size*np.sqrt(3)/2 - wall_thickness, y = j*cell_size + 5/6*np.sqrt(3)*wall_thickness)
                k2 = mapdl.k(x = (i+1)*cell_size*np.sqrt(3)/2 - wall_thickness, y = (j+1)*cell_size - 1/6*np.sqrt(3)*wall_thickness)
                k3 = mapdl.k(x = i*cell_size*np.sqrt(3)/2 + 1/2*wall_thickness, y = (j+1/2)*cell_size + 1/3*np.sqrt(3)*wall_thickness)
                hole = mapdl.a(k1, k2, k3)
                try:
                    design_part = mapdl.asba(na1 = design_part, na2 = hole)
                except:
                    pass
            for j in range(n_cell-1):
                k1 = mapdl.k(x = i*cell_size*np.sqrt(3)/2, y = (j+1/2)*cell_size + 5/6*np.sqrt(3)*wall_thickness)
                k2 = mapdl.k(x = i*cell_size*np.sqrt(3)/2, y = (j+3/2)*cell_size - 1/6*np.sqrt(3)*wall_thickness)
                k3 = mapdl.k(x = (i+1)*cell_size*np.sqrt(3)/2 - 3/2*wall_thickness, y = (j+1)*cell_size + 1/3*np.sqrt(3)*wall_thickness)
                hole = mapdl.a(k1, k2, k3)
                try:
                    design_part = mapdl.asba(na1 = design_part, na2 = hole)
                except:
                    pass
        if (i+1)*cell_size*np.sqrt(3)/2 > design_width:
            break
        i+=1

        mapdl.rectng(x1=-total_width/2+design_width/2, x2=0, y1=0, y2=design_height)
        mapdl.rectng(design_width, total_width/2+design_width/2, 0, design_height)
        if external_wall:
            mapdl.rectng(x1=0, x2=design_width, y1=0, y2=external_wall_thickness)
            mapdl.rectng(x1=0, x2=design_width, y1=design_height, y2=design_height-external_wall_thickness)
        mapdl.aadd("all")
        mapdl.vext("ALL", dz = total_z)

def kagome_horizontal(mapdl, n_cell, relative_density, external_wall=False, three_dim=True):
    design_width = 60/1000
    design_height = 25/1000
    total_width = 144/1000
    total_z = 5/1000

    t_over_l = np.sqrt(3)/2 - np.sqrt(3 - 4*relative_density)/2
    cell_size = design_height/(np.sqrt(3)*n_cell + t_over_l)
    wall_thickness = cell_size*t_over_l
    external_wall_thickness = wall_thickness

    side_length = cell_size - 1/np.sqrt(3)*wall_thickness
    design_part = mapdl.rectng(-1/1000, design_width, 0, design_height)

    for i in range(n_cell):
        j = 0
        while True:
            center = [cell_size/2 - np.sqrt(3)/2*wall_thickness + 2*j*cell_size - (i%2)*cell_size,
                    (cell_size + 1/np.sqrt(3)*wall_thickness)*np.sqrt(3)/2 + i*np.sqrt(3)*cell_size]
            # 조건문 넣을 필요 없이, design part 늘리면 해결되긴 함
            if center[0] - side_length < design_width:
                k_list = []
                for k in range(6):
                    k_list.append(mapdl.k(x = center[0] + side_length*np.cos(np.pi/3*k), y = center[1] + side_length*np.sin(np.pi/3*k)))
                hole = mapdl.a(*k_list)
                try:
                    design_part = mapdl.asba(na1 = design_part, na2 = hole)
                except:
                    print(f"error2 at i={i}, j={j}")
            # 조건문 넣을 필요 없이, design part 늘리면 해결되긴 함
            if cell_size + 2*j*cell_size - (i%2)*cell_size < design_width:
                k1 = mapdl.k(x = cell_size + 2*j*cell_size - (i%2)*cell_size, y = wall_thickness + i*np.sqrt(3)*cell_size)
                k2 = mapdl.k(x = 2*cell_size - np.sqrt(3)*wall_thickness + 2*j*cell_size - (i%2)*cell_size, y = wall_thickness + i*np.sqrt(3)*cell_size)
                k3 = mapdl.k(x = (3/2)*cell_size - 1/2*np.sqrt(3)*wall_thickness + 2*j*cell_size - (i%2)*cell_size, y = np.sqrt(3)/2*cell_size - 1/2*wall_thickness + i*np.sqrt(3)*cell_size)
                hole = mapdl.a(k1, k2, k3)
                try:
                    design_part = mapdl.asba(na1 = design_part, na2 = hole)
                except Exception as e:
                    print(f"error2 at i={i}, j={j}")
                    print(e)

                k1 = mapdl.k(x = cell_size + 2*j*cell_size - (i%2)*cell_size, y = np.sqrt(3)*cell_size + i*np.sqrt(3)*cell_size)
                k2 = mapdl.k(x = 2*cell_size - np.sqrt(3)*wall_thickness + 2*j*cell_size - (i%2)*cell_size, y = np.sqrt(3)*cell_size + i*np.sqrt(3)*cell_size)
                k3 = mapdl.k(x = (3/2)*cell_size - 1/2*np.sqrt(3)*wall_thickness + 2*j*cell_size - (i%2)*cell_size, y = np.sqrt(3)/2*cell_size + 3/2*wall_thickness + i*np.sqrt(3)*cell_size)
                hole = mapdl.a(k1, k2, k3)
                try:
                    design_part = mapdl.asba(na1 = design_part, na2 = hole)
                except Exception as e:
                    print(f"error3 at i={i}, j={j}")
                    print(e)

            if cell_size/2 - 1/np.sqrt(3)*wall_thickness + 2*j*cell_size > design_width:
                break

            j += 1

    mapdl.rectng(x1=-total_width/2+design_width/2, x2=0, y1=0, y2=design_height)
    mapdl.rectng(design_width, total_width/2+design_width/2, 0, design_height)
    mapdl.aadd("all")

    if three_dim:
        mapdl.vext("ALL", dz = total_z)
    
    return wall_thickness