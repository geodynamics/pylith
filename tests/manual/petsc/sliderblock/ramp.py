import matplotlib.pyplot as pyplot
import numpy

df = 3.0
ta = 0.5
tr = 4.0
dt = 0.01


amax = 0.5 * df * ta / (1/6*tr**3 -1/3*(tr-0.5*ta)**3 +0.5*tr*(tr-0.5*ta)**2 -0.5*tr**2*(tr-0.5*ta) +0.5*(tr-ta)*(tr-0.5*ta)**2 -0.5*(tr-ta)**2*(tr-0.5*ta) +1/4*ta**2*(tr-0.5*ta) +1/6*(tr-ta)**3 -1/8*ta**3)
c = 2.0 * amax / ta;

t = numpy.arange(0.0, 1.5*tr, dt)


def acc(t):
    mask1 = t <= 0.5*ta
    mask2 = numpy.logical_and(t > 0.5*ta, t <= ta)
    mask3 = numpy.logical_and(t > ta, t <= tr-ta)
    mask4 = numpy.logical_and(t > tr-ta, t <= tr-0.5*ta)
    mask5 = numpy.logical_and(t > tr-0.5*ta, t <= tr)
    mask6 = t > tr
    a = mask1 * c*t
    a += mask2 * c*(ta-t)
    a += mask3 * 0.0
    a += mask4 * -c*(t-(tr-ta))
    a += mask5 * c*(t-tr)
    a += mask6 * 0.0
    return a

def vel(t):
    mask1 = t <= 0.5*ta
    mask2 = numpy.logical_and(t > 0.5*ta, t <= ta)
    mask3 = numpy.logical_and(t > ta, t <= tr-ta)
    mask4 = numpy.logical_and(t > tr-ta, t <= tr-0.5*ta)
    mask5 = numpy.logical_and(t > tr-0.5*ta, t <= tr)
    mask6 = t > tr
    v = mask1 * c * 1/2*t**2
    v += mask2 * c*(ta*t-1/2*t**2-1/4*ta**2)
    v += mask3 * c*1/4*ta**2
    v += mask4 * -c*(1/2*t**2-(tr-ta)*t+1/2*(tr-ta)**2-1/4*ta**2)
    v += mask5 * c*(1/2*t**2-tr*t+1/2*tr**2)
    return v

def disp(t):
    mask1 = t <= 0.5*ta
    mask2 = numpy.logical_and(t > 0.5*ta, t <= ta)
    mask3 = numpy.logical_and(t > ta, t <= tr-ta)
    mask4 = numpy.logical_and(t > tr-ta, t <= tr-0.5*ta)
    mask5 = numpy.logical_and(t > tr-0.5*ta, t <= tr)
    mask6 = t > tr
    d = mask1 * c * 1/6*t**3
    d += mask2 * c*(0.5*ta*t**2 -1/6.0*t**3 -0.25*ta**2*t +1/24.0*ta**3)
    d += mask3 * c*(1/4*ta**2*t-1/8*ta**3)
    d += mask4 * -c*(1/6*t**3 -0.5*(tr-ta)*t**2 +0.5*(tr-ta)**2*t -1/4*ta**2*t -1/6*(tr-ta)**3 +1/8*ta**3)
    d += mask5 * c*(1/6*t**3 -0.5*tr*t**2 +0.5*tr**2*t -1/3*(tr-0.5*ta)**3 +0.5*tr*(tr-0.5*ta)**2 -0.5*tr**2*(tr-0.5*ta) +0.5*(tr-ta)*(tr-0.5*ta)**2 -0.5*(tr-ta)**2*(tr-0.5*ta) +1/4*ta**2*(tr-0.5*ta) +1/6*(tr-ta)**3 -1/8*ta**3)
    d += mask6 * c*(1/6*tr**3 -1/3*(tr-0.5*ta)**3 +0.5*tr*(tr-0.5*ta)**2 -0.5*tr**2*(tr-0.5*ta) +0.5*(tr-ta)*(tr-0.5*ta)**2 -0.5*(tr-ta)**2*(tr-0.5*ta) +1/4*ta**2*(tr-0.5*ta) +1/6*(tr-ta)**3 -1/8*ta**3)
    return d

a = acc(t)
pyplot.plot(t, a)
pyplot.show()

v = vel(t)
pyplot.plot(t, v, t, dt*numpy.cumsum(a), '--')
pyplot.show()

d = disp(t)
pyplot.plot(t, d, t, dt*numpy.cumsum(v), '--')
pyplot.show()
