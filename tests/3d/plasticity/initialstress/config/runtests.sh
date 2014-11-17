#!/bin/bash
pylith elastic_sshear_dt01.cfg
pylith elastic_sshear_dt02.cfg
pylith elastic_sshear_dt05.cfg
pylith elastic_sshear_dt10.cfg
pylith plastic_sshear_dt01.cfg
pylith plastic_sshear_dt02.cfg
pylith plastic_sshear_dt05.cfg
pylith plastic_sshear_dt10.cfg
python ./plot_invars.py
