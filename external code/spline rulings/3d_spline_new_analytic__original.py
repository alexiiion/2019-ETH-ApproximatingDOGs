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
der2_tan_proj_l = (der2*tan).sum(axis = 1)
prin_n = der2-tan*der2_tan_proj_l.reshape(len(der2_tan_proj_l),1)
prin_n_norm = np.sqrt((prin_n*prin_n).sum(axis=1))
prin_n = prin_n/prin_n_norm.reshape(len(prin_n_norm),1)


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

# compute curvature in every parameterization (local expression at https://en.wikipedia.org/wiki/Curvature)
der1_der2_cross = np.cross(der1,der2)
der1_der2_cross_norm = np.sqrt((der1_der2_cross * der1_der2_cross).sum(axis=1))
k = der1_der2_cross_norm/np.power(der1_norm,3)

torsion_up = (der1_der2_cross*der3).sum(axis = 1)
torsion_down = (der1_der2_cross * der1_der2_cross).sum(axis=1)
torsion = torsion_up/torsion_down

k_div_torsion = k/torsion
#tan_angle = np.arctan(k_div_torsion)
#r_curv_torsion = tan*np.cos(tan_angle).reshape(len(tan_angle),1)+binormal*np.sin(tan_angle).reshape(len(tan_angle),1)
tan_angle = np.arctan(k_div_torsion)
r_curv_torsion = tan*np.cos(tan_angle).reshape(len(tan_angle),1)+binormal*np.sin(tan_angle).reshape(len(tan_angle),1)

scale = 2 # plot rulings at length of "scale"
r = scale * r
old_r = scale * old_r
r_curv_torsion = scale * r_curv_torsion 


# principal normal up to +- sign (change the order to np.cross(r,der1_normalized) to flip the sign)
print der1_norm

fig2 = plt.figure(2)
ax3d = fig2.add_subplot(111, projection='3d')
ax3d.plot(x_true, y_true, z_true, 'b') # curve points, without noise
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