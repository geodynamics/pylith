from numpy import *

p0 = array([-1.3, -1.4, -0.8], dtype=float64)
p1 = array([1.2, -1.5, -0.9], dtype=float64)
p2 = array([1.4, 0.7, -1.2], dtype=float64)
p3 = array([-1.6, 0.4, -0.5], dtype=float64)
p4 = array([-1.7, -0.8, 1.8], dtype=float64)
p5 = array([2.1, -1.7, 0.6], dtype=float64)
p6 = array([2.3, 0.2, 1.9], dtype=float64)
p7 = array([-1.8, 0.3, 2.2], dtype=float64)

j000 = zeros( (3,3), dtype=float64)
j000[:,0] = p1-p0
j000[:,1] = p3-p0
j000[:,2] = p4-p0
print 'j000\n',j000

j100 = zeros( (3,3), dtype=float64)
j100[:,0] = p1-p0
j100[:,1] = p2-p1
j100[:,2] = p5-p1
print 'j100\n',j100

j010 = zeros( (3,3), dtype=float64)
j010[:,0] = p2-p3
j010[:,1] = p3-p0
j010[:,2] = p7-p3
print 'j010\n',j010

j110 = zeros( (3,3), dtype=float64)
j110[:,0] = p2-p3
j110[:,1] = p2-p1
j110[:,2] = p6-p2
print 'j110\n',j110

j001 = zeros( (3,3), dtype=float64)
j001[:,0] = p5-p4
j001[:,1] = p7-p4
j001[:,2] = p4-p0
print 'j001\n',j001

j101 = zeros( (3,3), dtype=float64)
j101[:,0] = p5-p4
j101[:,1] = p6-p5
j101[:,2] = p5-p1
print 'j101\n',j101

j011 = zeros( (3,3), dtype=float64)
j011[:,0] = p6-p7
j011[:,1] = p7-p4
j011[:,2] = p7-p3
print 'j011\n',j011

j111 = zeros( (3,3), dtype=float64)
j111[:,0] = p6-p7
j111[:,1] = p6-p5
j111[:,2] = p6-p2
print 'j111\n',j111

(x0,y0,z0) = p0
(x1,y1,z1) = p1
(x2,y2,z2) = p2
(x3,y3,z3) = p3
(x4,y4,z4) = p4
(x5,y5,z5) = p5
(x6,y6,z6) = p6
(x7,y7,z7) = p7
x = 0.2
y = 0.8
z = 0.7

f_xy = x2-x1-x3+x0
g_xy = y2-y1-y3+y0
h_xy = z2-z1-z3+z0

f_yz = x7-x3-x4+x0
g_yz = y7-y3-y4+y0
h_yz = z7-z3-z4+z0

f_xz = x5-x1-x4+x0
g_xz = y5-y1-y4+y0
h_xz = z5-z1-z4+z0

f_xyz = x6-x0+x1-x2+x3+x4-x5-x7
g_xyz = y6-y0+y1-y2+y3+y4-y5-y7
h_xyz = z6-z0+z1-z2+z3+z4-z5-z7

j00 = x1-x0 + f_xy*y + f_xz*z + f_xyz*y*z
j10 = y1-y0 + g_xy*y + g_xz*z + g_xyz*y*z
j20 = z1-z0 + h_xy*y + h_xz*z + h_xyz*y*z

j01 = x3-x0 + f_xy*x + f_yz*z + f_xyz*x*z
j11 = y3-y0 + g_xy*x + g_yz*z + g_xyz*x*z
j21 = z3-z0 + h_xy*x + h_yz*z + h_xyz*x*z

j02 = x4-x0 + f_xz*x + f_yz*y + f_xyz*x*y
j12 = y4-y0 + g_xz*x + g_yz*y + g_xyz*x*y
j22 = z4-z0 + h_xz*x + h_yz*y + h_xyz*x*y

j = array([ [j00, j01, j02],
            [j10, j11, j12],
            [j20, j21, j22] ], dtype=float64)
print 'j\n',j

