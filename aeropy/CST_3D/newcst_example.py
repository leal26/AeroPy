import aeropy.CST_3D as cst
from aeropy.filehandling.vtk import generate_surface
import panairwrapper.mesh_tools as meshtools

import numpy as np
import matplotlib.pyplot as plt


# wing parameters
span = 2.
eta_cp = [0., 1.]
chord = [1.5, 0.4]
sweep = [0., -1.5]
dihedral = [0., .1]

# fuselage parameters
length = 4.

f_sx = cst.piecewise_linear(eta_cp, chord)
f_sy_upper = cst.BernstienPolynomial(5, [0.172802, 0.167353, 0.130747,
                                         0.172053, 0.112797, 0.168891])
f_sy_lower = cst.BernstienPolynomial(5, [0.163339, 0.175407, 0.134176,
                                         0.152834, 0.133240, 0.161677])
f_etashear = cst.piecewise_linear(eta_cp, sweep)
f_zetashear_u = cst.piecewise_linear(eta_cp, np.array(dihedral))
f_zetashear_l = cst.piecewise_linear(eta_cp, -np.array(dihedral))

wing_upper = cst.CST3D(rotation=(0., 0., 90.),
                       location=(1.5, 0., 0.),
                       XYZ=(span/2., 1., 1.),
                       ref=(0., 1., 0.),
                       sx=f_sx,
                       nx=(0., 0.),
                       sy=f_sy_upper,
                       ny=(1., 1.),
                       etashear=f_etashear,
                       zetashear=f_zetashear_u)

wing_lower = cst.CST3D(rotation=(0., 0., 90.),
                       location=(1.5, 0., 0.),
                       XYZ=(span/2., 1., -1.),
                       ref=(0., 1., 0.),
                       sx=f_sx,
                       nx=(0., 0.),
                       sy=f_sy_lower,
                       ny=(1., 1.),
                       etashear=f_etashear,
                       zetashear=f_zetashear_l)

# psi_t = 0.5
# eta_t = 0.5
# x_t, y_t, z_t = wing_upper(psi_t, eta_t)
# psi_t, eta_t = wing_upper.inverse(x_t, y_t, z_t)
# print(psi_t, eta_t)

fuselage = cst.CST3D(rotation=(-90., 0., 0.),
                     XYZ=(length, 1.5, 1.5),
                     nx=(1., 1.),
                     ref=(0., 0.5, 0.),
                     ny=(.5, .5))


# generate mesh
N_chord = 5
N_span = 10

eta_spacing_w = meshtools.cosine_spacing(0., 1., N_chord)

# upper wing intersection and mesh
p_intersect_u = cst.intersection(wing_upper, fuselage, eta_spacing_w, 0.3)
np.savetxt("intersection_upper.csv", p_intersect_u)

psi_wu_intrsct, eta_wu_intrsct = wing_upper.inverse(p_intersect_u[:, 0],
                                                    p_intersect_u[:, 1],
                                                    p_intersect_u[:, 2])
psi_limit_wu = np.array([psi_wu_intrsct, eta_spacing_w]).T
psi_wu, eta_wu = meshtools.meshparameterspace((N_span, N_chord), flip=False, cos_spacing=True,
                                              psi_limits=(psi_limit_wu, None))

mesh_wu = wing_upper(psi_wu, eta_wu)

# lower wing intersection and mesh
p_intersect_l = cst.intersection(wing_lower, fuselage, eta_spacing_w, 0.3)
np.savetxt("intersection_lower.csv", p_intersect_l)

psi_wl_intrsct, eta_wl_intrsct = wing_lower.inverse(p_intersect_l[:, 0],
                                                    p_intersect_l[:, 1],
                                                    p_intersect_l[:, 2])
psi_limit_wl = np.array([psi_wl_intrsct, eta_spacing_w]).T
psi_wl, eta_wl = meshtools.meshparameterspace((N_span, N_chord), flip=False, cos_spacing=True,
                                              psi_limits=(psi_limit_wl, None))

mesh_wl = wing_lower(psi_wl, eta_wl)

# fuselage intersections and mesh
N_nose = 3
N_tail = 3
N_circ = 5
psi_fu_intersect, eta_fu_intersect = fuselage.inverse(p_intersect_u[:, 0],
                                                      p_intersect_u[:, 1],
                                                      p_intersect_u[:, 2])
psi_fl_intersect, eta_fl_intersect = fuselage.inverse(p_intersect_l[:, 0],
                                                      p_intersect_l[:, 1],
                                                      p_intersect_l[:, 2])
psi_frontpoint = psi_fu_intersect[-1]
eta_frontpoint = eta_fu_intersect[-1]
psi_rearpoint = psi_fu_intersect[0]
eta_rearpoint = eta_fu_intersect[0]
# print("front point", psi_frontpoint, eta_frontpoint)
# print("rear point", psi_rearpoint, eta_rearpoint)

front_section_fuse = np.full((N_nose, 2,), 0.)
rear_section_fuse = np.full((N_tail, 2,), 0.)

front_section_fuse[:, 0] = np.flipud(np.linspace(0., psi_frontpoint, N_nose))
front_section_fuse[:, 1] = eta_frontpoint
rear_section_fuse[:, 1] = eta_rearpoint
# rear_section_fuse[:, 1] = np.linspace(rear_point_fuse[1], 1., N_tail)
rear_section_fuse[:, 0] = np.flipud(meshtools.cosine_spacing(psi_rearpoint, 1., N_tail))

int_upperfuse_p = np.concatenate((rear_section_fuse[:-1], np.array([psi_fu_intersect, eta_fu_intersect]).T, front_section_fuse[1:]))
int_lowerfuse_p = np.concatenate((rear_section_fuse[:-1], np.array([psi_fl_intersect, eta_fl_intersect]).T, front_section_fuse[1:]))

psi_fu, eta_fu = meshtools.meshparameterspace((N_chord+N_nose+N_tail-2, N_circ), flip=True,
                                              eta_limits=(None, np.flipud(int_upperfuse_p)),
                                              cos_spacing=False)
psi_fl, eta_fl = meshtools.meshparameterspace((N_chord+N_nose+N_tail-2, N_circ), flip=True,
                                              eta_limits=(np.flipud(int_lowerfuse_p), None),
                                              cos_spacing=False)

# print(psi_fu)
# print(eta_fu)
# 
# print(int_upperfuse_p[:, 0])
for j in range(N_circ):
    psi_fu[:, j] = int_upperfuse_p[:, 0]
    psi_fl[:, j] = int_lowerfuse_p[:, 0]
# print(psi_fu)
# print(eta_fu)
    # pmesh_upperfuse[i, :, 0] = int_upperfuse_p[:, 1]
    # pmesh_lowerfuse[i, :, 0] = int_lowerfuse_p[:, 1]
mesh_fu = fuselage(psi_fu, eta_fu)
mesh_fl = fuselage(psi_fl, eta_fl)

network_wu = np.dstack(mesh_wu)
network_wl = np.dstack(mesh_wl)
network_fu = np.dstack(mesh_fu)
network_fl = np.dstack(mesh_fl)

# generate cap
network_cap = np.zeros((N_chord, 2, 3))
network_cap[:, 0, :] = network_wl[-1, :, :]
network_cap[:, 1, :] = network_wu[-1, :, :]

# calculate wake
wing_trailing_edge = network_wu[:, 0, :]
fuselage_wake_boundary = fuselage(rear_section_fuse[:, 0], rear_section_fuse[:, 1])

fuselage_wake_boundary = np.flipud(np.array([fuselage_wake_boundary[0],
                                   fuselage_wake_boundary[1],
                                   fuselage_wake_boundary[2]]).T)

inner_endpoint = np.copy(fuselage_wake_boundary[-1])
aoa = 0.
n_wake_streamwise = len(fuselage_wake_boundary)
wake = meshtools.generate_wake(wing_trailing_edge, inner_endpoint[0],
                               n_wake_streamwise, aoa, cos_spacing=True)

wingbody_wake = np.zeros((n_wake_streamwise, 2, 3))
wingbody_wake[:, 0] = fuselage_wake_boundary
wingbody_wake[:, 1] = wake[:, 0]

# output_w = {'upper': wing_upper.mesh_surface,
#             'lower': wing_lower.mesh_surface}
# output_f = fuselage.mesh_surface

# Generating vtk files
generate_surface(network_wu, "wingupper")
generate_surface(network_wl, "winglower")
generate_surface(network_fu, "fuselageupper")
generate_surface(network_fl, "fuselagelower")
generate_surface(network_cap, "wingcap")
generate_surface(wake, "wake")
generate_surface(wingbody_wake, "wingbody_wake")
