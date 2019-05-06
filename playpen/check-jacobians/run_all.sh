#!/bin/bash
./run_all_2d.sh
./run_all_3d.sh
./check_jacobians.py 2>&1 | tee check_jacobians.log
