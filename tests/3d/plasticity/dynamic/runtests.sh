#!/bin/bash

pylith hex8.cfg swave.cfg hex8_swave.cfg
plot_data.py hex8 swave &

pylith hex8.cfg pwave.cfg hex8_pwave.cfg
plot_data.py hex8 pwave &

pylith tet4.cfg swave.cfg tet4_swave.cfg
plot_data.py tet4 swave &

pylith tet4.cfg pwave.cfg tet4_pwave.cfg
plot_data.py tet4 pwave &
