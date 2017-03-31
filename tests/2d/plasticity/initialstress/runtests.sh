#!/bin/bash
pylith elastic.cfg elastic_dt01.cfg
pylith elastic.cfg elastic_dt02.cfg
pylith elastic.cfg elastic_dt05.cfg
pylith elastic.cfg elastic_dt10.cfg

pylith plastic.cfg plastic_dt01.cfg
pylith plastic.cfg plastic_dt02.cfg
pylith plastic.cfg plastic_dt05.cfg
pylith plastic.cfg plastic_dt10.cfg

./plot_invars.py
