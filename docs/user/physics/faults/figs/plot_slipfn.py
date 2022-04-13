#!/usr/bin/env python3

import numpy
import matplotlib.pyplot as pyplot

FINAL_SLIP = 1.0
RISE_TIME = 2.0
IMPULSE_DURATION = RISE_TIME / 4.0

class Step():
    NAME = "step"
    
    def ticks(self):
        return ((0, "$t_r$"),)
    
    def slipfn(self, t):
        return FINAL_SLIP * (t > 0)

class ConstRate():
    NAME = "constrate"
    
    def ticks(self):
        return ((0, "$t_r$"),)
    
    def slipfn(self, t):
        slip_rate = FINAL_SLIP / numpy.max(t)
        return slip_rate * t * (t > 0)

class Ramp():
    NAME = "ramp"
    
    def _acc(self):
        cacc = 1.0/6.0 * RISE_TIME**3 \
          -1.0/3.0 * (RISE_TIME-0.5*IMPULSE_DURATION)**3 \
          +0.5*RISE_TIME * (RISE_TIME-0.5*IMPULSE_DURATION)**2 \
          -0.5*RISE_TIME**2 * (RISE_TIME-0.5*IMPULSE_DURATION) \
          +0.5*(RISE_TIME-IMPULSE_DURATION) * (RISE_TIME-0.5*IMPULSE_DURATION)**2 \
          -0.5*(RISE_TIME-IMPULSE_DURATION)**2 * (RISE_TIME-0.5*IMPULSE_DURATION) \
          +0.25*IMPULSE_DURATION**2 * (RISE_TIME-0.5*IMPULSE_DURATION) \
          +1.0/6.0*(RISE_TIME-IMPULSE_DURATION)**3 \
          -0.125*IMPULSE_DURATION**3
        return 0.5 * FINAL_SLIP * IMPULSE_DURATION / cacc;

    def ticks(self):
        return (
            (0, "$t_r$"),
            (IMPULSE_DURATION, "$t_r+t_\mathit{acc}$"),
            (RISE_TIME-IMPULSE_DURATION, "$t_\mathit{rise}-t_\mathit{acc}$"),
            (RISE_TIME, "$t_\mathit{rise}$"),
            )
    
    def slipfn(self, t):
        mask2 = numpy.logical_and(t > 0.0, t <= 0.5*IMPULSE_DURATION)
        fn2 = 1.0/6.0*t**3

        mask3 = numpy.logical_and(t > 0.5*IMPULSE_DURATION, t <= IMPULSE_DURATION)
        fn3 = 0.5*IMPULSE_DURATION*t**2 \
            -1/6.0*t**3 \
            -0.25*IMPULSE_DURATION**2*t \
            +1.0/24.0*IMPULSE_DURATION**3

        mask4 = numpy.logical_and(t > IMPULSE_DURATION, t <= RISE_TIME-IMPULSE_DURATION)
        fn4 = 0.25*IMPULSE_DURATION**2*t \
            -0.125*IMPULSE_DURATION**3

        mask5 = numpy.logical_and(t > RISE_TIME-IMPULSE_DURATION, t <= RISE_TIME-0.5*IMPULSE_DURATION)
        fn5 = -1.0/6.0*t**3 \
            +0.5*(RISE_TIME-IMPULSE_DURATION)*t**2 \
            -0.5*(RISE_TIME-IMPULSE_DURATION)**2*t \
            +0.25*IMPULSE_DURATION**2*t \
            +1.0/6.0*(RISE_TIME-IMPULSE_DURATION)**3 \
            -0.125*IMPULSE_DURATION**3

        mask6 = numpy.logical_and(t > RISE_TIME-0.5*IMPULSE_DURATION, t <= RISE_TIME)
        ti = RISE_TIME - 0.5*IMPULSE_DURATION
        fn6 = 1.0/6.0*t**3 \
            -0.5*RISE_TIME*t**2 \
            +0.5*RISE_TIME**2*t \
            -1.0/3.0*ti**3 \
            +0.5*RISE_TIME*ti**2 \
            -0.5*RISE_TIME**2*ti \
            +0.5*(RISE_TIME-IMPULSE_DURATION)*ti**2 \
            -0.5*(RISE_TIME-IMPULSE_DURATION)**2*ti \
            +0.25*IMPULSE_DURATION**2*ti \
            +1.0/6.0*(RISE_TIME-IMPULSE_DURATION)**3 \
            -0.125*IMPULSE_DURATION**3

        mask7 = t > RISE_TIME
        fn7 = 1.0/6.0*RISE_TIME**3 \
            -1.0/3.0*ti**3 \
            +0.5*RISE_TIME*ti**2 \
            -0.5*RISE_TIME**2*ti \
            +0.5*(RISE_TIME-IMPULSE_DURATION)*ti**2 \
            -0.5*(RISE_TIME-IMPULSE_DURATION)**2*ti \
            +0.25*IMPULSE_DURATION**2*ti \
            +1.0/6.0*(RISE_TIME-IMPULSE_DURATION)**3 \
            -0.125*IMPULSE_DURATION**3

        slip = mask2*fn2
        slip += mask3*fn3
        slip += mask4*fn4
        slip += mask5*fn5
        slip += mask6*fn6
        slip += mask7*fn7
        slip *= 2.0 * self._acc() / IMPULSE_DURATION
        return slip
    
class Brune():
    NAME = "brune"
    
    def ticks(self):
        return ((0, "$t_r$"), (RISE_TIME, "$t_\mathit{rise}$"))
    
    def slipfn(self, t):
        tau = 0.21081916 * RISE_TIME
        slip = FINAL_SLIP * (1.0 - numpy.exp(-t/tau)*(1.0+t/tau)) * (t > 0)
        return slip

class LiuCos():
    NAME = "liucos"
    
    def ticks(self):
        return ((0, "$t_r$"), (RISE_TIME, "$t_\mathit{rise}$"))
    
    def slipfn(self, t):
        tau = RISE_TIME
        tau1 = 0.13 * tau
        tau2 = tau - tau1
        Cn = numpy.pi /  (1.4 * numpy.pi * tau1 + 1.2 * tau1 + 0.3 * numpy.pi * tau2)

        mask1 = numpy.logical_and(t > 0.0, t <= tau1)
        fn1 = 0.7*t \
          - 0.7*tau1/numpy.pi*numpy.sin(numpy.pi*t/tau1) \
          - 0.6*tau1/(0.5*numpy.pi)*(numpy.cos(0.5*numpy.pi*t/tau1) - 1.0)

        mask2 = numpy.logical_and(t > tau1, t <= 2.0*tau1)
        fn2 = 1.0*t \
          - 0.7*tau1/numpy.pi*numpy.sin(numpy.pi*t/tau1) \
          + 0.3*tau2/numpy.pi*numpy.sin(numpy.pi*(t-tau1)/tau2) \
          + 1.2*tau1/numpy.pi - 0.3*tau1

        mask3 = numpy.logical_and(t > 2.0*tau1, t <= tau)
        fn3 = 0.3*t \
          + 0.3*tau2/numpy.pi*numpy.sin(numpy.pi*(t-tau1)/tau2) \
          + 1.1*tau1 + 1.2*tau1/numpy.pi;

        mask4 = t > tau
        fn4 = 1.0

        slip = Cn*mask1*fn1 + Cn*mask2*fn2 + Cn*mask3*fn3 + mask4*fn4
        return slip

class App():

    def main(self):
        FIG_WIDTH = 7.0
        FIG_HEIGHT = 2.5
        SUBPLOTS = {
            "left": 0.08,
            "bottom": 0.18,
            "right": 0.99,
            "top": 0.96,
            "wspace": 0.3,
            "hspace": 0,
            }
        NROWS = 1
        NCOLS = 2
    
        t = numpy.arange(-0.2*RISE_TIME, 1.2*RISE_TIME, 0.01)
        for fn_class in [Step, ConstRate, Ramp, Brune, LiuCos]:
            fn = fn_class()

            figure = pyplot.figure(figsize=(FIG_WIDTH, FIG_HEIGHT), dpi=300)
            figure.subplots_adjust(**SUBPLOTS)
            ticks = fn.ticks()

            ax = figure.add_subplot(NROWS, NCOLS, 1)
            slip = fn.slipfn(t)
            ax.plot(t, slip, 'b-', lw=2)
            ax.set_xlabel("Time")
            ax.set_ylabel("Slip")
            ax.set_xticks([t[0] for t in ticks])
            ax.set_xticklabels([t[1] for t in ticks])
            ax.set_xlim(numpy.min(t), numpy.max(t))
            ax.set_ylim(0, 1.0)
            
            ax = figure.add_subplot(NROWS, NCOLS, 2)
            slip_rate = numpy.gradient(slip, t)
            ax.plot(t, slip_rate, 'b-', lw=2)
            ax.set_xlabel("Time")
            ax.set_ylabel("Slip Rate")
            ax.set_xticks([t[0] for t in ticks])
            ax.set_xticklabels([t[1] for t in ticks])
            ax.set_xlim(numpy.min(t), numpy.max(t))
            ax.set_ylim(0, numpy.max(slip_rate))

            figure.savefig(f"slipfn-{fn.NAME}.pdf", pad_inches=0.02)
            pyplot.close(figure)
    

if __name__ == "__main__":
    App().main()


# End of file
