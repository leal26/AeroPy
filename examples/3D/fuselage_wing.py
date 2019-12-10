import aeropy.CST_3D as cst
from aeropy.filehandling.vtk import generate_surface
import aeropy.CST_3D.mesh_tools as meshtools

import numpy as np
import matplotlib.pyplot as plt


# wing parameters
span = 2.
eta_cp = [0., 1.]
taper = [1.0, 0.3]
sweep = [0., 1.5]
dihedral = [0., .0]
twist = [0., 1.]
TE_thickness = [0.01, 0.01]
# fuselage parameters
length = 4.

f_sy = cst.piecewise_linear(eta_cp, taper)
f_sx_upper = cst.BernsteinPolynomial(5, [0.172802, 0.167353, 0.130747,
                                         0.172053, 0.112797, 0.168891])
f_sx_lower = cst.BernsteinPolynomial(5, [0.163339, 0.175407, 0.134176,
                                         0.152834, 0.133240, 0.161677])
f_xshear = cst.piecewise_linear(eta_cp, sweep)
f_zshear = cst.piecewise_linear(eta_cp, dihedral)
f_twist = cst.piecewise_linear(eta_cp, twist)
f_TE = cst.piecewise_linear(eta_cp, TE_thickness)

wing_upper = cst.CST3D(rotation=(0., 0., 0.),
                       location=(1.5, 0., 0.),
                       XYZ=(1.5, span/2., .2),
                       ref=(0., 0., 0.),
                       sx=f_sx_upper,
                       nx=(1., 1.),
                       sy=f_sy,
                       ny=(0., 0.),
                       xshear=f_xshear,
                       zshear=f_zshear,
                       twist=f_twist,
                       TE_thickness=f_TE)

wing_lower = cst.CST3D(rotation=(0., 0., 0.),
                       location=(1.5, 0., 0.),
                       XYZ=(1.5, span/2., -.2),
                       ref=(0., 0., 0.),
                       sx=f_sx_lower,
                       nx=(1., 1.),
                       sy=f_sy,
                       ny=(0., 0.),
                       xshear=f_xshear,
                       zshear=f_zshear,
                       twist=f_twist,
                       TE_thickness=cst.piecewise_linear(eta_cp, [0, 0]))

psi_w, eta_w = meshtools.meshparameterspace((20, 20),
                                            psi_spacing='cosine',
                                            eta_spacing='cosine')

fuselage = cst.CST3D(rotation=(0., -90., -90.),
                     XYZ=(.4, length, .2),
                     nx=(.5, .5),
                     ref=(0.5, 0., 0.),
                     ny=(1., 1.),
                     TE_thickness=f_TE)

# generate mesh
N_chord = 5
N_span = 10
N_nose = 3
N_tail = 3
N_circ = 5

# wing and fuselage intersection
psi_spacing_w = meshtools.cosine_spacing()(0., 1., N_chord)
intersection_f_wu = cst.intersection(wing_upper, fuselage, psi_spacing_w, 0.3)
intersection_f_wl = cst.intersection(wing_lower, fuselage, psi_spacing_w, 0.3)

# upper wing network
np.savetxt("intersection_upper.csv", np.array(intersection_f_wu).T)
edge_wu = meshtools.gen_network_edge(wing_upper, intersection_f_wu)
psi_wu, eta_wu = meshtools.meshparameterspace((N_chord, N_span),
                                              psi_spacing='cosine',
                                              eta_spacing='cosine',
                                              eta_limits=(edge_wu, None))

mesh_wu = wing_upper(psi_wu, eta_wu)
# lower wing network
# np.savetxt("intersection_lower.csv", p_intersect_l)
edge_wl = meshtools.gen_network_edge(wing_lower, intersection_f_wl)

psi_wl, eta_wl = meshtools.meshparameterspace((N_chord, N_span),
                                              psi_spacing='cosine',
                                              eta_spacing='cosine',
                                              eta_limits=(edge_wl, None))

mesh_wl = wing_lower(psi_wl, eta_wl)

# fuselage mesh
edge_fu = meshtools.gen_network_edge(fuselage, intersection_f_wu, N_b=N_nose,
                                     N_t=N_tail, vertical=True)
edge_fl = meshtools.gen_network_edge(fuselage, intersection_f_wl, N_b=N_nose,
                                     N_t=N_tail, vertical=True)

psi_fu, eta_fu = meshtools.meshparameterspace((N_circ, N_chord+N_nose+N_tail-2),
                                              psi_spacing='cosine',
                                              eta_spacing='uniform',
                                              psi_limits=(edge_fu, None))


psi_fl, eta_fl = meshtools.meshparameterspace((N_circ, N_chord+N_nose+N_tail-2),
                                              psi_spacing='cosine',
                                              eta_spacing='uniform',
                                              psi_limits=(None, edge_fl))

mesh_fu = fuselage(psi_fu, eta_fu)
mesh_fl = fuselage(psi_fl, eta_fl)

# format as networks
network_wu = np.dstack(mesh_wu)
network_wl = np.dstack(mesh_wl)
network_fu = np.dstack(mesh_fu)
network_fl = np.dstack(mesh_fl)


# generate cap
network_cap = np.zeros((N_chord, 2, 3))
network_cap[:, 0, :] = network_wu[:, -1, :]
network_cap[:, 1, :] = network_wl[:, -1, :]

# calculate wake
wing_trailing_edge = network_wu[-1, :, :]

fuselage_wake_boundary = network_fu[0, -N_tail:]
inner_endpoint = np.copy(fuselage_wake_boundary[-1])
aoa = 0.
n_wake_streamwise = len(fuselage_wake_boundary)
spacing = meshtools.cosine_spacing()
wake = meshtools.generate_wake(wing_trailing_edge, inner_endpoint[0],
                               n_wake_streamwise, aoa, user_spacing=spacing)

wingbody_wake = np.zeros((n_wake_streamwise, 2, 3))
wingbody_wake[:, 0] = fuselage_wake_boundary
wingbody_wake[:, 1] = wake[:, 0]

# Generating vtk files
generate_surface(network_wu, "wingupper")
generate_surface(network_wl, "winglower")
generate_surface(network_fu, "fuselageupper")
generate_surface(network_fl, "fuselagelower")
generate_surface(network_cap, "wingcap")
generate_surface(wake, "wake")
generate_surface(wingbody_wake, "wingbody_wake")
