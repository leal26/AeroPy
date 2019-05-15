# import aeropy.CST_lib as cst
import aeropy.CST_3D as cst
import aeropy.CST_3D.mesh_tools as meshtools
import panairwrapper
from aeropy.filehandling.vtk import generate_surface

import numpy as np
import matplotlib.pyplot as plt


# wing_points = stl_processing.points_from_stl('./wing_right.stl')
# wing_points = np.genfromtxt('./wing_upper_points.csv', skip_header=1, delimiter=',')
aoa = 2.3067
# wing parameters
span = 4.578*2.  # meters
eta_cp = [0., 0.276004, .552007, 1.]
taper = [1.0, 0.375317, .204817, .080035]
chord_root = 21.43
sweep = [0., 14.4476, 19.1612, 24.79405]
# sweep = [0., 0.]
dihedral = [0., 0.095311, 0.172374]  # , 0., 0.]
eta_dihedral = [0., .552007, 1.]
twist = [0.1, -0.05, -0.1, -0.1]


f_sy = cst.piecewise_linear(eta_cp, taper)

# shape coefficients upper wing
a_mat_upper = np.array([[0.01263563, 0.02351035, 0.01410948, 0.02870071, 0.00672103, 0.01423143],
                        [0.01493258, 0.02318954, 0.0256106, 0.02404792, 0.0188881, 0.01735003],
                        [0.01451941, 0.03214948, 0.03073803, 0.03342343, 0.03707097, 0.03925866],
                        [0.01687441, 0.01973094, 0.05100545, 0.01154137, 0.05669724, 0.02297795]]).T

A0_u = cst.piecewise_linear(eta_cp, a_mat_upper[0])
A1_u = cst.piecewise_linear(eta_cp, a_mat_upper[1])
A2_u = cst.piecewise_linear(eta_cp, a_mat_upper[2])
A3_u = cst.piecewise_linear(eta_cp, a_mat_upper[3])
A4_u = cst.piecewise_linear(eta_cp, a_mat_upper[4])
A5_u = cst.piecewise_linear(eta_cp, a_mat_upper[5])

# shape coefficients lower wing
a_mat_lower = np.array([[-0.01115817, -0.00544143, -0.0136072, -0.01172747, -0.01432724, -0.01362163],
                        [-0.01205803, 0.00105121, -0.01297246, -0.00643726, -0.0084237, -0.00878171],
                        [-0.0115331, 0.00720371, 0.00013933, -0.0095523, 0.02167392, 0.00636361],
                        [-0.00905945, 0.0107793, -0.01631087, 0.0092086, 0.01001222, 0.01334978]]).T

A0_l = cst.piecewise_linear(eta_cp, a_mat_lower[0])
A1_l = cst.piecewise_linear(eta_cp, a_mat_lower[1])
A2_l = cst.piecewise_linear(eta_cp, a_mat_lower[2])
A3_l = cst.piecewise_linear(eta_cp, a_mat_lower[3])
A4_l = cst.piecewise_linear(eta_cp, a_mat_lower[4])
A5_l = cst.piecewise_linear(eta_cp, a_mat_lower[5])

A_upper = [A0_u, A1_u, A2_u, A3_u, A4_u, A5_u]
A_lower = [A0_l, A1_l, A2_l, A3_l, A4_l, A5_l]

f_sx_upper = cst.BernsteinPolynomial(5, A_upper)
f_sx_lower = cst.BernsteinPolynomial(5, A_lower)

f_xshear = cst.piecewise_linear(eta_cp, sweep)
f_zshear = cst.piecewise_linear(eta_dihedral, dihedral)
f_twist = cst.piecewise_linear(eta_cp, twist)

wing_upper = cst.CST3D(rotation=(0., aoa, 0.),
                       location=(12., 0., -0.535),
                       XYZ=(chord_root, span/2., chord_root),
                       ref=(0., 0., 0.),
                       sx=f_sx_upper,
                       nx=(0.5, 1.),
                       sy=f_sy,
                       ny=(0., 0.),
                       xshear=f_xshear,
                       zshear=f_zshear,
                       twist=f_twist)

wing_lower = cst.CST3D(rotation=(0., aoa, 0.),
                       location=(12., 0., -0.535),
                       XYZ=(chord_root, span/2., chord_root),
                       ref=(0., 0., 0.),
                       sx=f_sx_lower,
                       nx=(0.5, 1.),
                       sy=f_sy,
                       ny=(0., 0.),
                       xshear=f_xshear,
                       zshear=f_zshear,
                       twist=f_twist)

# psi_w, eta_w = meshtools.meshparameterspace((20, 50),
#                                             psi_spacing='cosine',
#                                             eta_spacing='cosine')
# psi_stl, eta_stl = wing_upper.inverse(np.array([wing_points[:, 1]]),
#                                       np.array([wing_points[:, 2]]),
#                                       np.array([wing_points[:, 3]]))
# plt.scatter(psi_stl, eta_stl)
# plt.show()

# wing_mesh_u = wing_upper(psi_w, eta_w)
# wing_mesh_l = wing_lower(psi_w, eta_w)
# network_wu = np.dstack(wing_mesh_u)
# network_wl = np.dstack(wing_mesh_l)
# generate_surface(network_wu, "wing_upper")
# generate_surface(network_wl, "wing_lower")


# fuselage parameters
fuse_data = np.genfromtxt('./fuselage_raw.txt', skip_header=1)
i_sort = np.argsort(fuse_data[:, 0])
fuse_data = fuse_data[i_sort]

y_section, x_LE, c, N1, N2, A0, A1, A2, A3, A4, A5, error = fuse_data.T
length = y_section[-1]
width = np.max(c)
eta_f_cp = y_section/length
fuse_xshear = cst.piecewise_linear(eta_f_cp, x_LE)
fuse_xtaper = cst.piecewise_linear(eta_f_cp, 1.001*c/width)
fuse_nx1 = cst.piecewise_linear(eta_f_cp, N1)
fuse_nx2 = cst.piecewise_linear(eta_f_cp, N2)
A0_f = cst.piecewise_linear(eta_f_cp, A0)
A1_f = cst.piecewise_linear(eta_f_cp, A1)
A2_f = cst.piecewise_linear(eta_f_cp, A2)
A3_f = cst.piecewise_linear(eta_f_cp, A3)
A4_f = cst.piecewise_linear(eta_f_cp, A4)
A5_f = cst.piecewise_linear(eta_f_cp, A5)
f_sx_fuse = cst.BernsteinPolynomial(5, [A0_f, A1_f, A2_f, A3_f, A4_f, A5_f])


fuselage = cst.CST3D(rotation=(-90., -90.+aoa, 0.),
                     XYZ=(width, length, width),
                     ref=(0., 0., 0.),
                     sx=f_sx_fuse,
                     nx=(fuse_nx1, fuse_nx2),
                     sy=fuse_xtaper,
                     ny=(0., 0.),
                     xshear=fuse_xshear)

# fuse_mesh = fuselage(psi_w, eta_w)
# network_f = np.dstack(fuse_mesh)
# generate_surface(network_f, "fuselage")

# generate mesh
N_chord = 20
N_span = 10
N_nose = 20
N_tail = 50
N_circ = 20

# wing and fuselage intersection
cos_space = meshtools.cosine_spacing()
psi_spacing_w = cos_space(0., 1., N_chord)
intersection_f_wu = cst.intersection(wing_upper, fuselage, psi_spacing_w, 0.3)
intersection_f_wl = cst.intersection(wing_lower, fuselage, psi_spacing_w, 0.3)

# upper wing network
np.savetxt("intersection_upper.csv", np.array(intersection_f_wu).T)
edge_wf = meshtools.gen_network_edge(wing_upper, intersection_f_wu)
# edge_w1 = meshtools.gen_network_edge(wing_upper, intersection_f_wu)
edge_1 = meshtools.constant_eta_edge(eta_cp[1], N_chord)
edge_2 = meshtools.constant_eta_edge(eta_cp[2], N_chord)
edge_3 = meshtools.constant_eta_edge(eta_cp[3], N_chord)
psi_wu1, eta_wu1 = meshtools.meshparameterspace((N_chord, N_span),
                                                psi_spacing='cosine',
                                                eta_spacing='cosine',
                                                eta_limits=(edge_wf, edge_1))
psi_wu2, eta_wu2 = meshtools.meshparameterspace((N_chord, N_span),
                                                psi_spacing='cosine',
                                                eta_spacing='cosine',
                                                eta_limits=(edge_1, edge_2))
psi_wu3, eta_wu3 = meshtools.meshparameterspace((N_chord, N_span),
                                                psi_spacing='cosine',
                                                eta_spacing='cosine',
                                                eta_limits=(edge_2, edge_3))

mesh_wu1 = wing_upper(psi_wu1, eta_wu1)
mesh_wu2 = wing_upper(psi_wu2, eta_wu2)
mesh_wu3 = wing_upper(psi_wu3, eta_wu3)
# lower wing network
# np.savetxt("intersection_lower.csv", p_intersect_l)
edge_wl = meshtools.gen_network_edge(wing_lower, intersection_f_wl)

psi_wl1, eta_wl1 = meshtools.meshparameterspace((N_chord, N_span),
                                                psi_spacing='cosine',
                                                eta_spacing='cosine',
                                                eta_limits=(edge_wl, edge_1))
psi_wl2, eta_wl2 = meshtools.meshparameterspace((N_chord, N_span),
                                                psi_spacing='cosine',
                                                eta_spacing='cosine',
                                                eta_limits=(edge_1, edge_2))
psi_wl3, eta_wl3 = meshtools.meshparameterspace((N_chord, N_span),
                                                psi_spacing='cosine',
                                                eta_spacing='cosine',
                                                eta_limits=(edge_2, edge_3))

mesh_wl1 = wing_lower(psi_wl1, eta_wl1)
mesh_wl2 = wing_lower(psi_wl2, eta_wl2)
mesh_wl3 = wing_lower(psi_wl3, eta_wl3)

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
network_wu1 = np.dstack(mesh_wu1)
network_wu2 = np.dstack(mesh_wu2)
network_wu3 = np.dstack(mesh_wu3)
network_wl1 = np.dstack(mesh_wl1)
network_wl2 = np.dstack(mesh_wl2)
network_wl3 = np.dstack(mesh_wl3)
network_fu = np.dstack(mesh_fu)
network_fl = np.dstack(mesh_fl)


# generate  wing cap
network_wingcap = np.zeros((N_chord, 2, 3))
network_wingcap[:, 0, :] = network_wu3[:, -1, :]
network_wingcap[:, 1, :] = network_wl3[:, -1, :]

# generate tail cap
# print(network_fl[:, -1, 2])
# print(network_fu[:, -1, 2])
network_fusecap = np.zeros((2*N_circ-1, 2, 3))
network_fusecap[:, 1, :] = np.concatenate((network_fl[:-1, -1],
                                           network_fu[:, -1]))
network_fusecap[:, 0, 0] = network_fusecap[:, 1, 0]
network_fusecap[:, 0, 2] = network_fusecap[:, 1, 2]
# np.savetxt('tail_cap0.csv', network_fuse_cap[:, 0, :])
# np.savetxt('tail_cap1.csv', network_fuse_cap[:, 1, :])
# exit()


# calculate wake
wing_trailing_edge = np.concatenate((network_wu1[-1, :-1, :],
                                     network_wu2[-1, :-1, :],
                                     network_wu3[-1, :, :]))

# print(wing_trailing_edge)
# exit()

fuselage_wake_boundary = network_fu[0, -N_tail:]
inner_endpoint = np.copy(fuselage_wake_boundary[-1])
n_wake_streamwise = len(fuselage_wake_boundary+1)
cos_space = meshtools.cosine_spacing()

body_wake_l = meshtools.generate_wake(network_fl[:, -1], inner_endpoint[0]+.05,
                                      2, aoa, user_spacing=cos_space)
body_wake_u = meshtools.generate_wake(network_fu[:, -1], inner_endpoint[0]+.05,
                                      2, aoa, user_spacing=cos_space)

wing_wake = meshtools.generate_wake(wing_trailing_edge, inner_endpoint[0]+.05,
                                    n_wake_streamwise+1, aoa, user_spacing=cos_space)


wingbody_wake = np.zeros((n_wake_streamwise+1, 2, 3))
wingbody_wake[:-1, 0] = fuselage_wake_boundary
wingbody_wake[-1, 0] = body_wake_u[-1, 0]
wingbody_wake[:, 1] = wing_wake[:, 0]

# Generating vtk files
generate_surface(network_wu1, "wingupper1")
generate_surface(network_wu2, "wingupper2")
generate_surface(network_wu3, "wingupper3")
generate_surface(network_wl1, "winglower1")
generate_surface(network_wl2, "winglower2")
generate_surface(network_wl3, "winglower3")
generate_surface(network_fu, "fuselageupper")
generate_surface(network_fl, "fuselagelower")
generate_surface(network_wingcap, "wingcap")
generate_surface(network_fusecap, "fusecap")
generate_surface(wing_wake, "wake")
generate_surface(body_wake_l, "body_wake_l")
generate_surface(body_wake_u, "body_wake_u")
generate_surface(wingbody_wake, "wingbody_wake")

# run in Panair
gamma = 1.4
MACH = 1.6
panair = panairwrapper.PanairWrapper('wingbody', exe='panair.exe')
panair.set_aero_state(MACH, aoa)
panair.set_symmetry(1, 0)
panair.add_network("wing_u1", np.flipud(network_wu1))
panair.add_network("wing_u2", np.flipud(network_wu2))
panair.add_network("wing_u3", np.flipud(network_wu3))
panair.add_network("wing_l1", network_wl1)
panair.add_network("wing_l2", network_wl2)
panair.add_network("wing_l3", network_wl3)
panair.add_network("fuselage_u", np.flipud(network_fu))
panair.add_network("fuselage_l", np.flipud(network_fl))
panair.add_network("wing_cap", np.flipud(network_wingcap))
panair.add_network("fuse_cap", network_fusecap, 5.)
panair.add_network("wake", wing_wake, 18.)
panair.add_network("body_wake_l", body_wake_l, 20.)
panair.add_network("body_wake_u", body_wake_u, 20.)
panair.add_network("wingbody_wake", wingbody_wake, 20.)

panair.set_sensor(MACH, aoa, 3, 32.92, 2.)

results = panair.run()
results.write_vtk()


offbody_data = results.get_offbody_data()
distance_along_sensor = offbody_data[:, 2]
dp_over_p = 0.5*gamma*MACH**2*offbody_data[:, -2]
nearfield_sig = np.array([distance_along_sensor, dp_over_p]).T

plt.plot(nearfield_sig[:, 0], nearfield_sig[:, 1])
plt.title("nearfield signature")
plt.show()
np.savetxt('nearfield_sig', nearfield_sig)
