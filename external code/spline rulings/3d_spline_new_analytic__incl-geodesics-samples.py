import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D

#noise = 0.1
noise = 0.005
# noise = 0

# smoothing parameter (0 means fit the points exactly, and otherwise might wanna understand the thing written here : https://stackoverflow.com/questions/8719754/scipy-interpolate-univariatespline-not-smoothing-regardless-of-parameters?rq=1)
if noise == 0:
	smooth_p = 1 # smooth_p = 0 means pass through the points exactly
else:
	smooth_p = 2 # might need to play with this, we just want smooth derivatives


num_sample_pts = 80
total_rad = 10
s_true = np.linspace(0, total_rad, num_sample_pts)

"""
# arc
s_true = np.linspace(-1, 2, num_sample_pts)
z_factor = 3
x_true = np.cos(s_true)*100
y_true = np.sin(s_true)*100
z_true = s_true/z_factor
"""

"""
# helix
z_factor = 1
x_true = np.cos(s_true)
y_true = np.sin(s_true)
z_true = s_true/z_factor
"""


"""
# vivians curve
a = 1
x_true = a*(1+np.cos(s_true))
y_true = a*np.sin(s_true)
z_true = 2*a*np.sin(s_true/2)
"""


"""
# conical spiral
a = 2
h = 0.75
r = 0.5
x_true = s_true*r*np.cos(a*s_true)
y_true = s_true*r*np.sin(a*s_true)
z_true = s_true
"""

"""
# cycloid
R = 1
r = 1
x_true = R*np.cos((r/R)*(s_true-np.sin(s_true)))
y_true = R*np.sin((r/R)*(s_true-np.sin(s_true)))
z_true = r-r*np.cos(s_true)
"""

# # bunny_v1097 -- jump 3:
# x_true = [7.9088, 7.47324, 6.87891, 6.26775, 5.61109, 4.89023, 4.42728, 3.88281, 3.08641, 2.54517, 1.77788, 1.25005, 0.674907, 0.305249, -0.195438, -0.558647, -0.860882, -1.22209, ]
# y_true = [-2.34168, -1.12291, 0.437691, 1.76631, 2.66145, 3.39475, 3.74999, 3.97927, 4.08408, 4.03875, 3.84606, 3.55298, 2.9185, 2.40365, 1.40497, 0.514436, -0.289832, -1.3698, ]
# z_true = [-3.62766, -3.18741, -2.63774, -2.02581, -1.26969, -0.393336, 0.205758, 0.974225, 2.13637, 2.94928, 4.12349, 4.95085, 5.88726, 6.49153, 7.26929, 7.8084, 8.23049, 8.68012, ]
# num_sample_pts = 18

# # bunny_v1097 -- jump 2:
# x_true = [7.9088, 7.60553, 7.20775, 6.87891, 6.45422, 6.01838, 5.61109, 5.19026, 4.81619, 4.42728, 4.00048, 3.54515, 3.08641, 2.63805, 2.20208, 1.77788, 1.37173, 1.00747, 0.674907, 0.355147, 0.0456219, -0.195438, -0.472208, -0.70796, -0.860882, -1.21354, -1.22723, ]
# y_true = [-2.34168, -1.43953, -0.509434, 0.437691, 1.44723, 2.1498, 2.66145, 3.09606, 3.45767, 3.74999, 3.94455, 4.05568, 4.08408, 4.05512, 3.96694, 3.84606, 3.65771, 3.3303, 2.9185, 2.48113, 2.00386, 1.40497, 0.743547, 0.121844, -0.289832, -1.34393, -1.38841, ]
# z_true = [-3.62766, -3.3239, -2.90482, -2.63774, -2.23076, -1.74559, -1.26969, -0.759332, -0.299596, 0.205758, 0.805901, 1.46108, 2.13637, 2.80796, 3.47301, 4.12349, 4.75422, 5.34423, 5.88726, 6.41051, 6.91349, 7.26929, 7.68735, 8.01907, 8.23049, 8.66954, 8.68455, ]
# num_sample_pts = 27

# # bunny_v1097:
# x_true = [7.9088, 7.75218, 7.60553, 7.47324, 7.20775, 6.96528, 6.87891, 6.65759, 6.45422, 6.26775, 6.01838, 5.84274, 5.61109, 5.37457, 5.19026, 4.89023, 4.81619, 4.43753, 4.42728, 4.37987, 4.00048, 3.88281, 3.54515, 3.28175, 3.08641, 2.67488, 2.63805, 2.54517, 2.20208, 2.06273, 1.77788, 1.44714, 1.37173, 1.25005, 1.00747, 0.865375, 0.674907, 0.419963, 0.355147, 0.305249, 0.0456219, -0.185802, -0.195438, -0.205226, -0.472208, -0.558647, -0.70796, -0.791774, -0.860882, -1.00895, -1.21354, -1.22209, -1.22723, -1.3989]
# y_true = [-2.34168, -1.73037, -1.43953, -1.12291, -0.509434, 0.172459, 0.437691, 0.98269, 1.44723, 1.76631, 2.1498, 2.39239, 2.66145, 2.91853, 3.09606, 3.39475, 3.45767, 3.74325, 3.74999, 3.77329, 3.94455, 3.97927, 4.05568, 4.07774, 4.08408, 4.05834, 4.05512, 4.03875, 3.96694, 3.93041, 3.84606, 3.69839, 3.65771, 3.55298, 3.3303, 3.1568, 2.9185, 2.5756, 2.48113, 2.40365, 2.00386, 1.43, 1.40497, 1.38084, 0.743547, 0.514436, 0.121844, -0.103161, -0.289832, -0.752355, -1.34393, -1.3698, -1.38841, -2.07147]
# z_true = [-3.62766, -3.49224, -3.3239, -3.18741, -2.90482, -2.70031, -2.63774, -2.43075, -2.23076, -2.02581, -1.74559, -1.54414, -1.26969, -0.98579, -0.759332, -0.393336, -0.299596, 0.192095, 0.205758, 0.271983, 0.805901, 0.974225, 1.46108, 1.8477, 2.13637, 2.75262, 2.80796, 2.94928, 3.47301, 3.68637, 4.12349, 4.63632, 4.75422, 4.95085, 5.34423, 5.57618, 5.88726, 6.30462, 6.41051, 6.49153, 6.91349, 7.25533, 7.26929, 7.28384, 7.68735, 7.8084, 8.01907, 8.13521, 8.23049, 8.39176, 8.66954, 8.68012, 8.68455, 8.77841]
# num_sample_pts = 54

# # bunny_v1382:
# x_true = [-6.64047, -6.35423, -6.34084, -6.32511, -6.03546, -5.68388, -5.67895, -5.67173, -5.2996, -4.78956, -4.64375, -4.43602, -3.91939, -3.54665, -3.18667, -2.60741, -2.45997, -2.36123, -1.72924, -1.13008, -0.997596, -0.828209, -0.270246, 0.145889, 0.461473, 0.897415, 1.20194, 1.40737, 1.9445, 2.62785, 2.69456, 2.8077, 3.44619, 3.80227, 4.2, 4.94858, 5.6719, 6.07315, 6.37674, 6.82813, 7.05728, 7.25455, 7.71763]
# y_true = [6.17752, 5.59635, 5.56313, 5.52791, 4.94257, 4.35125, 4.34292, 4.33178, 3.71869, 3.53769, 3.46895, 3.44041, 3.32644, 3.26538, 3.1991, 3.05336, 3.00236, 3.01012, 3.04201, 3.19644, 3.22964, 3.28935, 3.45811, 3.58205, 3.67021, 3.76653, 3.82517, 3.86457, 3.95019, 4.05234, 4.05957, 4.07332, 4.13047, 4.1417, 4.13943, 4.06902, 3.8795, 3.72896, 3.60738, 3.39468, 3.28108, 3.17588, 2.91252]
# z_true = [3.6831, 3.56444, 3.55725, 3.54975, 3.42653, 3.30312, 3.30138, 3.29907, 3.17201, 3.12096, 3.10172, 3.09919, 3.08164, 3.07264, 3.06255, 3.03855, 3.03127, 3.02658, 2.99566, 2.97674, 2.97246, 2.9689, 2.95529, 2.94512, 2.93723, 2.92538, 2.91705, 2.91144, 2.89686, 2.87822, 2.87629, 2.87312, 2.85468, 2.84471, 2.83355, 2.81136, 2.78796, 2.77436, 2.76389, 2.74733, 2.73886, 2.73165, 2.71456]
# num_sample_pts = 43

# bunny_v1308:
x_true = [14.8203, 14.4652, 14.0698, 13.8179, 13.4, 13.0331, 12.649, 12.3205, 12.299, 12.2161, 12.1715, 12.1589, 12.1383, 12.0885, 12.0596, 12.0117, 11.9068, 11.9034, 11.895, 11.6984, 11.5364, 11.4166, 11.2227, 11.1203, 11.0463, 10.758, 10.4203, 10.3449, 10.2017, 9.90885, 9.74073, 9.42757, 8.99342, 8.89979, 8.69739, 8.33388, 8.12858, 7.75604, 7.22832, 7.12713, 7.01182, 6.47033, 6.14126, 5.79395, 5.07607, 4.33594, 3.96553, 3.58166, 2.82789, 2.82384, 2.82158, 2.06752, 1.62788, 1.31492, 0.816607, 0.573275, 0.409678, -0.159005, -0.812278, -0.886539, -1.00559, -1.62013, -2.00277, -2.35766, -2.95502, -3.10476, -3.19524, -3.84757, -4.64567, -4.66205, -4.68108, -5.31847, -5.57541, -5.71115, -5.89032, -6.10097, -6.22751, -6.55773, -6.82403, -6.92094, -7.0454, -7.14107, -7.16502, -7.39643]
y_true = [-8.1747, -7.64119, -7.33365, -7.21889, -7.17387, -7.14692, -7.04164, -6.91437, -6.69792, -6.20006, -5.57308, -5.45378, -5.23209, -4.68701, -4.39449, -3.92051, -3.18132, -3.16146, -3.13072, -2.42673, -2.00848, -1.71738, -1.24918, -1.01727, -0.876916, -0.349629, 0.162067, 0.278546, 0.484881, 0.88803, 1.09432, 1.45747, 1.88726, 1.97624, 2.14725, 2.44783, 2.5954, 2.87207, 3.16411, 3.22057, 3.27674, 3.50883, 3.65845, 3.79003, 4.00349, 4.07674, 4.08378, 4.07393, 4.0163, 4.01594, 4.01578, 3.92786, 3.86464, 3.80986, 3.70917, 3.6481, 3.59902, 3.41632, 3.20178, 3.17619, 3.14718, 2.95893, 2.94502, 2.91499, 3.03073, 3.04399, 3.06472, 3.15952, 3.22532, 3.22485, 3.23549, 3.53741, 3.93103, 4.10768, 4.38084, 4.67443, 4.83214, 5.19654, 5.63657, 5.7782, 6.14287, 6.4022, 6.50949, 6.98713]
z_true = [2.26258, 2.27374, 2.29263, 2.3147, 2.38596, 2.45073, 2.50799, 2.55618, 2.58084, 2.63914, 2.71256, 2.72643, 2.75234, 2.81608, 2.85017, 2.90537, 2.9899, 2.99214, 2.99534, 3.06856, 3.10959, 3.13771, 3.18286, 3.20492, 3.21803, 3.26705, 3.31291, 3.32338, 3.34179, 3.3779, 3.39704, 3.43114, 3.47234, 3.48089, 3.49734, 3.52634, 3.54119, 3.56865, 3.6027, 3.60925, 3.61627, 3.64586, 3.66518, 3.68244, 3.71304, 3.72718, 3.72973, 3.73007, 3.72641, 3.72638, 3.72637, 3.71932, 3.71406, 3.70954, 3.7015, 3.69615, 3.69114, 3.67171, 3.64868, 3.64589, 3.64293, 3.62014, 3.63238, 3.64027, 3.67847, 3.68368, 3.69222, 3.73493, 3.77299, 3.77301, 3.77919, 3.9615, 4.14755, 4.23078, 4.3616, 4.50177, 4.57746, 4.75248, 4.96913, 5.03787, 5.23265, 5.37025, 5.42876, 5.65848]
num_sample_pts = 84

# # bunny_v1306:
# x_true = [14.4346, 14.0956, 14.0944, 14.0916, 13.4208, 12.9693, 12.7437, 12.2936, 12.2661, 12.2641, 12.203, 12.2024, 12.167, 12.1396, 12.1147, 12.045, 12.0355, 12.0301, 11.9224, 11.8331, 11.7192, 11.4669, 11.4359, 11.3654, 11.1266, 11.0042, 10.7603, 10.4372, 10.3417, 10.0789, 9.90686, 9.83882, 9.41861, 9.17621, 8.88153, 8.45391, 8.31711, 7.98621, 7.76632, 7.64895, 7.16492, 6.61993, 6.51832, 6.37958, 5.86816, 5.55924, 5.16262, 4.42791, 3.68377, 3.20968, 2.92677, 2.49212, 2.17226, 1.94776, 1.41966, 0.711312, 0.682614, 0.660152, -0.0390715, -0.622947, -0.755002, -0.934326, -1.47985, -1.93032, -2.21938, -2.55921, -2.97243, -3.34824, -3.71251, -4.05428, -4.55256, -4.9956, -5.31766, -5.77569, -5.85846, -6.24011, -6.24066, -6.24076, -6.80173, -6.86384, -7.43628, -7.58427, -8.15913, -8.38257, -8.90879, -9.09353, -9.46626, -9.64281, -9.9658, -10.1602, -10.6326, -10.8788, -11.4519, -11.7506, -12.2548, -12.66, -12.9627]
# y_true = [-8.56874, -7.88781, -7.88556, -7.88366, -7.50797, -7.36736, -7.3284, -6.89997, -6.86674, -6.85056, -6.17304, -5.9361, -5.44349, -4.94089, -4.68333, -4.01468, -3.92258, -3.88743, -3.1636, -2.82094, -2.42496, -1.78836, -1.71513, -1.54974, -1.01686, -0.787734, -0.34928, 0.132998, 0.278456, 0.652517, 0.88885, 0.972198, 1.4554, 1.69346, 1.96994, 2.32677, 2.43104, 2.69387, 2.85666, 2.9213, 3.13492, 3.36963, 3.41175, 3.47289, 3.6847, 3.77879, 3.879, 3.97774, 3.97455, 3.94823, 3.92093, 3.88768, 3.85251, 3.82181, 3.72677, 3.54822, 3.53912, 3.53108, 3.27544, 3.06287, 3.01175, 2.95465, 2.74949, 2.63373, 2.58402, 2.60745, 2.62686, 2.71618, 2.77826, 2.73317, 2.67453, 2.56452, 2.62751, 2.80514, 2.8534, 3.13431, 3.13469, 3.13475, 3.40651, 3.43807, 3.65305, 3.70133, 3.85322, 3.90627, 4.04756, 4.13899, 4.35065, 4.44417, 4.66411, 4.75327, 4.93403, 4.99891, 5.09521, 5.11624, 5.10393, 5.03594, 4.96542]
# z_true = [1.09075, 1.10759, 1.10769, 1.10873, 1.38015, 1.62208, 1.75609, 1.93364, 1.94435, 1.94775, 2.08835, 2.13571, 2.23557, 2.33742, 2.38938, 2.52405, 2.5426, 2.5496, 2.69394, 2.76146, 2.83923, 2.96239, 2.9765, 3.00835, 3.11079, 3.15522, 3.24038, 3.33494, 3.3634, 3.4372, 3.48419, 3.50143, 3.60343, 3.65628, 3.71858, 3.80209, 3.82712, 3.88971, 3.92956, 3.9487, 4.01976, 4.09854, 4.1129, 4.1332, 4.20511, 4.24026, 4.28084, 4.33879, 4.37675, 4.39607, 4.40556, 4.42142, 4.43111, 4.43655, 4.44488, 4.44742, 4.44708, 4.44652, 4.42772, 4.41229, 4.40795, 4.40528, 4.38699, 4.39248, 4.40743, 4.47037, 4.54162, 4.65023, 4.74264, 4.78365, 4.84645, 4.87824, 5.05516, 5.49603, 5.61137, 6.28357, 6.28447, 6.28462, 6.84847, 6.91468, 7.3508, 7.44807, 7.75282, 7.85956, 8.14511, 8.32314, 8.72885, 8.90899, 9.33259, 9.49238, 9.79884, 9.88932, 9.97764, 9.96303, 9.84664, 9.65039, 9.46908]
# num_sample_pts = 97

# # synthetic example: 2D, arc line arc, regular sampling
# x_true = [-1, -0.996917, -0.987688, -0.97237, -0.951057, -0.92388, -0.891007, -0.85264, -0.809017, -0.760406, -0.707107, -0.649448, -0.587785, -0.522499, -0.45399, -0.382683, -0.309017, -0.233445, -0.156434, -0.0784591, 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.07846, 1.15643, 1.23345, 1.30902, 1.38268, 1.45399, 1.5225, 1.58779, 1.64945, 1.70711, 1.76041, 1.80902, 1.85264, 1.89101, 1.92388, 1.95106, 1.97237, 1.98769, 1.99692]
# y_true = [-1, -0.921541, -0.843566, -0.766555, -0.690983, -0.617317, -0.54601, -0.477501, -0.412215, -0.350552, -0.292893, -0.239594, -0.190983, -0.14736, -0.108993, -0.0761205, -0.0489435, -0.0276301, -0.0123117, -0.00308267, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00308267, 0.0123117, 0.0276301, 0.0489435, 0.0761205, 0.108993, 0.14736, 0.190983, 0.239594, 0.292893, 0.350552, 0.412215, 0.477501, 0.54601, 0.617317, 0.690983, 0.766555, 0.843566, 0.921541]
# z_true = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# num_sample_pts = 60


# #noise = 0.1
# #smooth_p = 10
noise = 0.0
smooth_p = 5

num_spline_pts = 300

# add noise possibly (set noise to 0.1)
x_sample = x_true + noise * np.random.randn(num_sample_pts)
y_sample = y_true + noise * np.random.randn(num_sample_pts)
z_sample = z_true + noise * np.random.randn(num_sample_pts)



tck, u = interpolate.splprep([x_sample,y_sample,z_sample], k=5, s=smooth_p) # change s to something bigger for instance 2e3 with noise = 0.1
x_knots, y_knots, z_knots = interpolate.splev(tck[0], tck)
u_fine = np.linspace(0,1,num_spline_pts)
x_fine, y_fine, z_fine = interpolate.splev(u_fine, tck)

# Compute ruling directions as the cross product of the second derivative with the third derivative, normalized
x_der1, y_der1, z_der1 = interpolate.splev(u_fine, tck,1)
der1 = np.vstack((x_der1,y_der1,z_der1)).transpose()
der1_norm = np.sqrt((der1 * der1).sum(axis=1))
x_der2, y_der2, z_der2 = interpolate.splev(u_fine, tck,2)
x_der3, y_der3, z_der3 = interpolate.splev(u_fine, tck,3)
der2 = np.vstack((x_der2,y_der2,z_der2)).transpose()

der2_norm = np.sqrt((der2 * der2).sum(axis=1)) # cross products norm
der2_normalized = der2/der2_norm.reshape(len(der2_norm),1)
der3 = np.vstack((x_der3,y_der3,z_der3)).transpose()

der1_normalized = der1/der1_norm.reshape(len(der1_norm),1)
tan = der1_normalized
print("tan ", tan);

der2_tan_proj_l = (der2*tan).sum(axis = 1)
prin_n = der2-tan*der2_tan_proj_l.reshape(len(der2_tan_proj_l),1)
prin_n_norm = np.sqrt((prin_n*prin_n).sum(axis=1))
prin_n = prin_n/prin_n_norm.reshape(len(prin_n_norm),1)
print("prin_n", prin_n);


# Not normalized than: prin_n' = (g''-<g'/|g'|,g''>g'/|g'|)' = g'''- (<g'/|g'|,g''>g'/|g'|)'        					      (1)
# using the multiplication rule we get: (<g'/|g'|,g''>g'/|g'|)' = (<g'/|g'|,g''>)'g'/|g'| + (<g'/|g'|,g''>)(g'/|g'|)'         (2)
# now (<g'/|g'|,g''>)' = <(g'/|g'|)',g''> + <(g'/|g'|,g'''>				         											  (3) 
# and  (g'/|g'|)' = (g''*|g'|-g'*|g'|') / |g'|^2 																			  (4)
# lastly |g'|' = <g'(t),g''(t)> / |g'(t)|																					  (5)

# start compute things
eq_5 = (der1 * der2).sum(axis=1) / der1_norm
eq_5 = eq_5.reshape(len(eq_5),1)

eq_4 = (der2*der1_norm.reshape(len(der1_norm),1)-der1*eq_5)/((der1 * der1).sum(axis=1)).reshape(len(der1_norm),1)
eq_3 = (eq_4*der2).sum(axis=1) + (tan*der3).sum(axis=1)
#print 'eq_3 = ', eq_3
eq_2 = eq_3.reshape(len(eq_3),1)*der1_normalized + (der1_normalized*der2).sum(axis=1).reshape(len(der1_norm),1)*eq_4
#print 'eq_2 = ', eq_2.shape
eq_1 = der3-eq_2


old_r = np.cross(der2,der3)
l2norm = np.sqrt((old_r * old_r).sum(axis=1)) # cross products norm
old_r = old_r/l2norm.reshape(len(l2norm),1) # normalize the ruling direction

r = np.cross(prin_n,eq_1)
l2norm = np.sqrt((r * r).sum(axis=1)) # cross products norm
r = r/l2norm.reshape(len(l2norm),1) # normalize the ruling direction


# Third way with curvature and torsion
binormal = np.cross(tan,prin_n)
print("binormal ", binormal);

# compute curvature in every parameterization (local expression at https://en.wikipedia.org/wiki/Curvature)
der1_der2_cross = np.cross(der1,der2)
print("der1_der2_cross ", der1_der2_cross);
der1_der2_cross_norm = np.sqrt((der1_der2_cross * der1_der2_cross).sum(axis=1))
print("der1_der2_cross_norm ", der1_der2_cross_norm);

k = der1_der2_cross_norm/np.power(der1_norm,3)
print("k ", k);

torsion_up = (der1_der2_cross*der3).sum(axis = 1)
print("torsion_up ", torsion_up);

torsion_down = (der1_der2_cross * der1_der2_cross).sum(axis=1)
print("torsion_down ", torsion_down);

torsion = torsion_up/torsion_down
print("torsion ", torsion);


k_div_torsion = k/torsion
print("k_div_torsion ", k_div_torsion);

#tan_angle = np.arctan(k_div_torsion)
#r_curv_torsion = tan*np.cos(tan_angle).reshape(len(tan_angle),1)+binormal*np.sin(tan_angle).reshape(len(tan_angle),1)
tan_angle = np.arctan(k_div_torsion)
print("tan_angle ", tan_angle);

r_curv_torsion = tan*np.cos(tan_angle).reshape(len(tan_angle),1)+binormal*np.sin(tan_angle).reshape(len(tan_angle),1)
print("r_curv_torsion ", r_curv_torsion);

scale = 1 # plot rulings at length of "scale"
r = scale * r
old_r = scale * old_r
r_curv_torsion = scale * r_curv_torsion 


# principal normal up to +- sign (change the order to np.cross(r,der1_normalized) to flip the sign)
# print(der1_norm)

fig2 = plt.figure(2)
ax3d = fig2.add_subplot(111, projection='3d')
# ax3d.plot(x_true, y_true, z_true, 'b') # curve points, without noise
ax3d.plot(x_sample, y_sample, z_sample, 'r*') # sample points (possiblyy with noise, depending on noise parameter)
ax3d.plot(x_knots, y_knots, z_knots, 'go') # knots
ax3d.plot(x_fine, y_fine, z_fine, 'g') # spline points
ax3d.quiver(x_fine-old_r[:,0], y_fine-old_r[:,1], z_fine-old_r[:,2], 2*old_r[:,0], 2*old_r[:,1], 2*old_r[:,2], arrow_length_ratio=0,color='y') # rulings
ax3d.quiver(x_fine-r[:,0], y_fine-r[:,1], z_fine-r[:,2], 2*r[:,0], 2*r[:,1], 2*r[:,2], arrow_length_ratio=0,color='c') # rulings
ax3d.quiver(x_fine-r_curv_torsion[:,0], y_fine - r_curv_torsion[:,1], z_fine - r_curv_torsion[:,2], 2*r_curv_torsion[:,0], 2*r_curv_torsion[:,1], 2*r_curv_torsion[:,2], arrow_length_ratio=0, color='m') # normals

fig2.show()



"""
# derivative plotting (for debugging, to make sure the derivative themselves are smooth and have a reasonable variation)
x_der1, y_der1, z_der1 = interpolate.splev(u_fine, tck,1)
fig3 = plt.figure(3)
ax3d = fig3.add_subplot(111, projection='3d')
ax3d.plot(x_der1, y_der1, z_der1, 'r')
fig3.show()

x_der2, y_der2, z_der2 = interpolate.splev(u_fine, tck,2)
fig4 = plt.figure(4)
ax3d = fig4.add_subplot(111, projection='3d')
ax3d.plot(x_der2, y_der2, z_der2, 'r')
fig4.show()

x_der3, y_der3, z_der3 = interpolate.splev(u_fine, tck,3)
fig5 = plt.figure(5)
ax3d = fig5.add_subplot(111, projection='3d')
ax3d.plot(x_der3, y_der3, z_der3, 'r')
fig5.show()
"""

plt.show()


"""
lim t-> 0 n(x) cross n(x+eps)

lim t-> 0 bla(x)*n(x) cross bla(x+t)n(x+t)

bla(x)n(x) cross (bla(x)n(x))' = bla(x)n(x) cross (bla(x)n(x)'+bla(x)'n(x)') which is parallel to n(x) cross n(x)'

The latter is since cross(n(x),n(x)) = 0

so bla(x)n(x) = g''(t)-<g'(t)/norm(g'(t)),g''(t)>g'(t)   (1)

deriving this is:
g'''(t) - derivative of other expression: (2)

(<g'(t)/norm(g'(t)),g''(t)>g'(t))' = (<g'(t)/norm(g'(t)),g''(t)>)'g'(t)+(<g'(t)|g'(t)|,g''(t)>)g''(t) (3)

with (<g'(t)/|g'(t)|,g''(t)>)' = <(g'(t)/|g'(t)|)',g''(t)> + <g'(t)/|g'(t)},g'''(t)> (4)
 

 By the quotient rule: (g'(t)/|g'(t)|)' = (g''(t)*|g'(t)|-g'(t)*|g'(t)|') / |g'(t)|^2 (5) 
 with:
 |g'(t)|' = <g'(t),g''(t)> / |g'(t)| (6)

 Wrapping it up:
 Plug in (6) in (5) we get:
 (g'(t)/|g'(t)|)' = (g''(t)*|g'(t)|-g'(t)*|g'(t)|') / |g'(t)|^2 = (g''(t)*|g'(t)|-g'(t)*<g'(t),g''(t)> / |g'(t)|) / |g'(t)|^2 (7)

 Plug in (7) in (4) we get:
(<g'(t)/|g'(t)|,g''(t)>)' = <(g''(t)*|g'(t)|-g'(t)*<g'(t),g''(t)> / |g'(t)|) / |g'(t)|^2,g''(t)> + <g'(t)/|g'(t)|,g'''(t)> (8)

Plug in (8) in (3) we get:
(<g'(t)/norm(g'(t)),g''(t)>g'(t))' = (<g'(t)/norm(g'(t)),g''(t)>)'g'(t)+(<g'(t)|g'(t)|,g''(t)>)g''(t) = 
		(<(g''(t)*|g'(t)|-g'(t)*<g'(t),g''(t)> / |g'(t)|) / |g'(t)|^2,g''(t)> + <g'(t)/|g'(t)|,g'''(t)>)g'(t) + (<g'(t)|g'(t)|,g''(t)>)g''(t) (9)

Plug in (9) in (2) we get:

g'''(t) - ((<(g''(t)*|g'(t)|-g'(t)*<g'(t),g''(t)> / |g'(t)|) / |g'(t)|^2,g''(t)> + <g'(t)/|g'(t)|,g'''(t)>)g'(t) + (<g'(t)|g'(t)|,g''(t)>)g''(t))
"""