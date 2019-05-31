#!/bin/bash
# Run all models.
pylith quad.cfg power_law_2d.cfg 2>&1 | tee power_law_2d.log
pylith quad.cfg power_law_2d_initstress.cfg 2>&1 | tee power_law_2d_initstress.log
pylith hex.cfg power_law_3d.cfg 2>&1 | tee power_law_3d.log
pylith hex.cfg power_law_3d_initstress.cfg 2>&1 | tee power_law_3d_initstress.log
