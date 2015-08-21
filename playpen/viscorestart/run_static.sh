#!/bin/bash
# Run all static simulations (no initial state variables).
pylith quad.cfg grav_static_genmaxps.cfg
pylith quad.cfg grav_static_maxps.cfg
pylith quad.cfg grav_static_powerlawps
pylith hex grav_static_genmax3d.cfg
pylith hex grav_static_max3d.cfg
pylith hex grav_static_powerlaw3d
