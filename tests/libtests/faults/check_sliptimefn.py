#!/usr/bin/env python3

import numpy
import matplotlib.pyplot as pyplot

FINAL_SLIP = 2.5
RISE_TIME = 3.0
T0 = 0.5
DT = 0.02

t = numpy.arange(0.0, 2*RISE_TIME, DT)


def brune(field):
    TAU = 0.21081916 * RISE_TIME
    mask = t >= T0
    if field == "slip":
        time_history = 1.0 - (numpy.exp(-(t-T0)/TAU) * (1.0 + (t-T0)/TAU))
    elif field == "slip_rate":
        time_history = (t-T0)/(TAU*TAU) * numpy.exp(-(t-T0)/TAU)
    elif field == "slip_acc":
        time_history = 1.0 / (TAU*TAU) * \
            (1.0 - (t-T0)/TAU) * numpy.exp(-(t-T0) / TAU)
    return FINAL_SLIP * mask * time_history


def liucos(field):
    TAU = RISE_TIME * 1.525
    TAU1 = 0.13 * TAU
    TAU2 = TAU - TAU1
    PI = numpy.pi
    CN = PI / (1.4 * PI * TAU1 + 1.2 * TAU1 + 0.3 * PI * TAU2)
    tR = t - T0

    if field == "slip":
        mask1 = numpy.logical_and(tR >= 0, tR <= TAU1)
        time_history1 = 0.7*tR - 0.7*TAU1/PI * \
            numpy.sin(PI*tR/TAU1) - 0.6*TAU1/(0.5*PI) * \
            (numpy.cos(0.5*PI*tR/TAU1) - 1.0)
        time_history1 *= CN

        mask2 = numpy.logical_and(tR > TAU1, tR <= 2*TAU1)
        time_history2 = 1.0*tR - 0.7*TAU1/PI * \
            numpy.sin(PI*tR/TAU1) + 0.3*TAU2/PI*numpy.sin(PI *
                                                          (tR-TAU1)/TAU2) + 1.2*TAU1/PI - 0.3*TAU1
        time_history2 *= CN

        mask3 = numpy.logical_and(tR > 2.0*TAU1, tR <= TAU)
        time_history3 = 0.3*tR + 0.3*TAU2/PI * \
            numpy.sin(PI*(tR-TAU1)/TAU2) + 1.1*TAU1 + 1.2*TAU1/PI
        time_history3 *= CN

        mask4 = tR > TAU
        time_history4 = 1.0

        time_history = mask1 * time_history1 + mask2 * time_history2 + \
            mask3 * time_history3 + mask4 * time_history4

    elif field == "slip_rate":
        mask1 = numpy.logical_and(tR >= 0, tR <= TAU1)
        time_history1 = 0.7 - 0.7 * \
            numpy.cos(PI*tR/TAU1) + 0.6*numpy.sin(0.5*PI*tR/TAU1)
        time_history1 *= CN

        mask2 = numpy.logical_and(tR > TAU1, tR <= 2*TAU1)
        time_history2 = 1.0 - 0.7 * \
            numpy.cos(PI*tR/TAU1) + 0.3*numpy.cos(PI*(tR-TAU1)/TAU2)
        time_history2 *= CN

        mask3 = numpy.logical_and(tR > 2.0*TAU1, tR <= TAU)
        time_history3 = 0.3 + 0.3*numpy.cos(PI*(tR-TAU1)/TAU2)
        time_history3 *= CN

        mask4 = tR > TAU
        time_history4 = 0.0

        time_history = mask1 * time_history1 + mask2 * time_history2 + \
            mask3 * time_history3 + mask4 * time_history4
    elif field == "slip_acc":
        mask1 = numpy.logical_and(tR >= 0, tR <= TAU1)
        time_history1 = PI/TAU1 * \
            (0.7*numpy.sin(PI*tR/TAU1) + 0.3*numpy.cos(0.5*PI*tR/TAU1))
        time_history1 *= CN

        mask2 = numpy.logical_and(tR > TAU1, tR <= 2*TAU1)
        time_history2 = PI/TAU1*0.7 * \
            numpy.sin(PI*tR/TAU1) - PI/TAU2*0.3*numpy.sin(PI*(tR-TAU1)/TAU2)
        time_history2 *= CN

        mask3 = numpy.logical_and(tR > 2.0*TAU1, tR <= TAU)
        time_history3 = -PI/TAU2*0.3*numpy.sin(PI*(tR-TAU1)/TAU2)
        time_history3 *= CN

        mask4 = tR > TAU
        time_history4 = 0.0

        time_history = mask1 * time_history1 + mask2 * time_history2 + \
            mask3 * time_history3 + mask4 * time_history4
    return FINAL_SLIP * time_history


def finite_diff(time_history):
    return numpy.gradient(time_history, DT)


pyplot.plot(t, brune("slip"), 'r-',
            t, liucos("slip_acc"), 'g-',
            t, finite_diff(liucos("slip_rate")), 'k--')
pyplot.show()
