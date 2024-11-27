import numpy as np

def triangle_parallel(mapdl, n_cell, wall_thickness, external_wall):
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