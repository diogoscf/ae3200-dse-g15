# type: ignore # TODO: remove this when this file works


def MOI(y, wbox=WINGBOX):
    stringers_top, stringers_bottom = stringers(y, wbox)
    chord = ((WING["taper_ratio"] - 1) * abs(y) + 1) * WING["root_chord"]
    spars = sorted([wbox["front_spar"], wbox["rear_spar"], *[s[0] for s in wbox["other_spars"] if s[1] >= abs(y)]])
    centroid_y = centroid(y, stringers_top, stringers_bottom, wbox)
    Ixx = 0
    Izz = 0

    # for top stringers
    for stringer in stringers_top:
        l_spar_idx = bisect.bisect_left(spars, stringer) - 1
        l_spar, r_spar = spars[l_spar_idx], spars[l_spar_idx + 1]
        z_coord = np.interp(stringer, (l_spar, r_spar), (airfoil_info(l_spar)[2], airfoil_info(r_spar)[2]))
        position = (chord * stringer, chord * z_coord)  # coordinates converted to meters
        rel_position = (position[0] - centroid_y[0], position[1] - centroid_y[1])  # relative position to centroid
        Ixx += wbox["stringer_area"] * (rel_position[1] ** 2)
        Izz += wbox["stringer_area"] * (rel_position[0] ** 2)

    # for bottom stringers
    for stringer in stringers_bottom:
        l_spar_idx = bisect.bisect_left(spars, stringer) - 1
        l_spar, r_spar = spars[l_spar_idx], spars[l_spar_idx + 1]
        z_coord = np.interp(stringer, (l_spar, r_spar), (airfoil_info(l_spar)[3], airfoil_info(r_spar)[3]))
        position = (chord * stringer, chord * z_coord)  # coordinates converted to meters
        rel_position = (position[0] - centroid_y[0], position[1] - centroid_y[1])  # relative position to centroid
        Ixx += wbox["stringer_area"] * (rel_position[1] ** 2)
        Izz += wbox["stringer_area"] * (rel_position[0] ** 2)
        # print(stringer, l_spar, r_spar, rel_position)

    # for spars
    for x in spars:
        position = (chord * x, chord * airfoil_info(x)[1])  # coordinates converted to meters
        rel_position = (position[0] - centroid_y[0], position[1] - centroid_y[1])  # relative position to centroid
        t = thickness_y(y, *wbox["spar_thickness"])
        h = airfoil_info(x)[0] * chord
        Ixx += t * (h**3) / 12 + t * h * (rel_position[1] ** 2)
        Izz += (t**3) * h / 12 + t * h * (rel_position[0] ** 2)

    # for skin
    for i in range(len(spars) - 1):
        left_spar = spars[i]
        right_spar = spars[i + 1]
        t = wbox["skin_thickness"]

        # top element
        center = (
            (left_spar + right_spar) / 2,
            (airfoil_info(left_spar)[2] + airfoil_info(right_spar)[2]) / 2,
        )  # center of skin element
        theta = np.arctan(
            (airfoil_info(right_spar)[2] - airfoil_info(left_spar)[2]) / (right_spar - left_spar)
        )  # angle to x-axis
        position = (chord * center[0], chord * center[1])  # coordinates converted to meters
        rel_position = (position[0] - centroid_y[0], position[1] - centroid_y[1])  # relative position to centroid
        w = (
            np.sqrt((left_spar - right_spar) ** 2 + (airfoil_info(left_spar)[2] - airfoil_info(right_spar)[2]) ** 2)
        ) * chord
        Ixx += (w**3) * t * (np.sin(theta) ** 2) / 12 + w * t * (rel_position[1] ** 2)
        Izz += (w**3) * t * (np.cos(theta) ** 2) / 12 + w * t * (rel_position[0] ** 2)

        # bottom element
        center = ((left_spar + right_spar) / 2, (airfoil_info(left_spar)[3] + airfoil_info(right_spar)[3]) / 2)
        theta = np.arctan((airfoil_info(right_spar)[3] - airfoil_info(left_spar)[3]) / (right_spar - left_spar))
        position = (chord * center[0], chord * center[1])
        rel_position = (position[0] - centroid_y[0], position[1] - centroid_y[1])
        w = (
            np.sqrt((left_spar - right_spar) ** 2 + (airfoil_info(left_spar)[3] - airfoil_info(right_spar)[3]) ** 2)
        ) * chord
        Ixx += (w**3) * t * (np.sin(theta) ** 2) / 12 + w * t * (rel_position[1] ** 2)
        Izz += (w**3) * t * (np.cos(theta) ** 2) / 12 + w * t * (rel_position[0] ** 2)

    return Ixx, Izz
