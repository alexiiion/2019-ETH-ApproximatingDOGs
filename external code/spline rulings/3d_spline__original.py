import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D

#noise = 0.1
noise = 0.05
#noise = 0

# smoothing parameter (0 means fit the points exactly, and otherwise might wanna understand the thing written here : https://stackoverflow.com/questions/8719754/scipy-interpolate-univariatespline-not-smoothing-regardless-of-parameters?rq=1)
if noise == 0:
	smooth_p = 1 # smooth_p = 0 means pass through the points exactly
else:
	smooth_p = 5 # might need to play with this, we just want smooth derivatives


num_sample_pts = 80
total_rad = 10
s_true = np.linspace(0, total_rad, num_sample_pts)


"""
# helix
z_factor = 3
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

# conical spiral
a = 2
h = 0.75
r = 0.5
x_true = s_true*r*np.cos(a*s_true)
y_true = s_true*r*np.sin(a*s_true)
z_true = s_true


"""
# cycloid
R = 1
r = 1
x_true = R*np.cos((r/R)*(s_true-np.sin(s_true)))
y_true = R*np.sin((r/R)*(s_true-np.sin(s_true)))
z_true = r-r*np.cos(s_true)
"""



# add noise possibly (set noise to 0.1)
x_sample = x_true + noise * np.random.randn(num_sample_pts)
y_sample = y_true + noise * np.random.randn(num_sample_pts)
z_sample = z_true + noise * np.random.randn(num_sample_pts)

tck, u = interpolate.splprep([x_sample,y_sample,z_sample], k=5, s=smooth_p) # change s to something bigger for instance 2e3 with noise = 0.1
x_knots, y_knots, z_knots = interpolate.splev(tck[0], tck)
u_fine = np.linspace(0,1,num_sample_pts)
x_fine, y_fine, z_fine = interpolate.splev(u_fine, tck)

# Compute ruling directions as the cross product of the second derivative with the third derivative, normalized
x_der2, y_der2, z_der2 = interpolate.splev(u_fine, tck,2)
x_der3, y_der3, z_der3 = interpolate.splev(u_fine, tck,3)
der2 = np.vstack((x_der2,y_der2,z_der2)).transpose()
der3 = np.vstack((x_der3,y_der3,z_der3)).transpose()
r = np.cross(der2,der3)
l2norm = np.sqrt((r * r).sum(axis=1)) # cross products norm
r = r/l2norm.reshape(len(l2norm),1) # normalize the ruling direction
scale = 2 # plot rulings at length of "scale"
r = scale * r

fig2 = plt.figure(2)
ax3d = fig2.add_subplot(111, projection='3d')
ax3d.plot(x_true, y_true, z_true, 'b') # curve points, without noise
ax3d.plot(x_sample, y_sample, z_sample, 'r*') # sample points (possiblyy with noise, depending on noise parameter)
ax3d.plot(x_knots, y_knots, z_knots, 'go') # knots
ax3d.plot(x_fine, y_fine, z_fine, 'g') # spline points
ax3d.quiver(x_fine, y_fine, z_fine, r[:,0], r[:,1], r[:,2], arrow_length_ratio=0) # rulings
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
