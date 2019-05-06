#!/bin/bash
# Run all 3D models.
pylith hex.cfg drucker_prager_3d_solve.cfg 2>&1 | tee drucker_prager_3d_solve.log
pylith hex.cfg drucker_prager_3d_jacobian.cfg 2>&1 | tee drucker_prager_3d_jacobian.log
# pylith hex.cfg drucker_prager_3d_tensile_solve.cfg 2>&1 | tee drucker_prager_3d_tensile_solve.log
# pylith hex.cfg drucker_prager_3d_tensile_jacobian.cfg 2>&1 | tee drucker_prager_3d_tensile_jacobian.log
pylith hex.cfg drucker_prager_3d_nonassoc_solve.cfg 2>&1 | tee drucker_prager_3d_nonassoc_solve.log
pylith hex.cfg drucker_prager_3d_nonassoc_jacobian.cfg 2>&1 | tee drucker_prager_3d_nonassoc_jacobian.log
pylith hex.cfg elastic_3d_solve.cfg 2>&1 | tee elastic_3d_solve.log
pylith hex.cfg elastic_3d_jacobian.cfg 2>&1 | tee elastic_3d_jacobian.log
pylith hex.cfg gen_max_3d_solve.cfg 2>&1 | tee gen_max_3d_solve.log
pylith hex.cfg gen_max_3d_jacobian.cfg 2>&1 | tee gen_max_3d_jacobian.log
pylith hex.cfg maxwell_3d_solve.cfg 2>&1 | tee maxwell_3d_solve.log
pylith hex.cfg maxwell_3d_jacobian.cfg 2>&1 | tee maxwell_3d_jacobian.log
pylith hex.cfg power_law_3d_solve.cfg 2>&1 | tee power_law_3d_solve.log
pylith hex.cfg power_law_3d_jacobian.cfg 2>&1 | tee power_law_3d_jacobian.log

# Extract Jacobian info and output to separate files.
grep -F -A 66 'Explicit preconditioning Jacobian' drucker_prager_3d_jacobian.log > pylith_dp_3d_jacobian.txt
# grep -F -A 66 'Explicit preconditioning Jacobian' drucker_prager_3d_tensile_jacobian.log > pylith_dp_3d_tensile_jacobian.txt
grep -F -A 66 'Explicit preconditioning Jacobian' drucker_prager_3d_nonassoc_jacobian.log > pylith_dp_3d_nonassoc_jacobian.txt
grep -F -A 66 'Explicit preconditioning Jacobian' elastic_3d_jacobian.log > pylith_el_3d_jacobian.txt
grep -F -A 66 'Explicit preconditioning Jacobian' gen_max_3d_jacobian.log > pylith_gm_3d_jacobian.txt
grep -F -A 66 'Explicit preconditioning Jacobian' maxwell_3d_jacobian.log > pylith_mx_3d_jacobian.txt
grep -F -A 66 'Explicit preconditioning Jacobian' power_law_3d_jacobian.log > pylith_pl_3d_jacobian.txt
grep -F -A 66 'Finite difference Jacobian' drucker_prager_3d_jacobian.log > finite_diff_dp_3d_jacobian.txt
# grep -F -A 66 'Finite difference Jacobian' drucker_prager_3d_tensile_jacobian.log > finite_diff_dp_3d_tensile_jacobian.txt
grep -F -A 66 'Finite difference Jacobian' drucker_prager_3d_nonassoc_jacobian.log > finite_diff_dp_3d_nonassoc_jacobian.txt
grep -F -A 66 'Finite difference Jacobian' elastic_3d_jacobian.log > finite_diff_el_3d_jacobian.txt
grep -F -A 66 'Finite difference Jacobian' gen_max_3d_jacobian.log > finite_diff_gm_3d_jacobian.txt
grep -F -A 66 'Finite difference Jacobian' maxwell_3d_jacobian.log > finite_diff_mx_3d_jacobian.txt
grep -F -A 66 'Finite difference Jacobian' power_law_3d_jacobian.log > finite_diff_pl_3d_jacobian.txt
grep -F -A 66 'minus finite difference Jacobian' drucker_prager_3d_jacobian.log > diff_dp_3d_jacobian.txt
# grep -F -A 66 'minus finite difference Jacobian' drucker_prager_3d_tensile_jacobian.log > diff_dp_3d_tensile_jacobian.txt
grep -F -A 66 'minus finite difference Jacobian' drucker_prager_3d_nonassoc_jacobian.log > diff_dp_3d_nonassoc_jacobian.txt
grep -F -A 66 'minus finite difference Jacobian' elastic_3d_jacobian.log > diff_el_3d_jacobian.txt
grep -F -A 66 'minus finite difference Jacobian' gen_max_3d_jacobian.log > diff_gm_3d_jacobian.txt
grep -F -A 66 'minus finite difference Jacobian' maxwell_3d_jacobian.log > diff_mx_3d_jacobian.txt
grep -F -A 66 'minus finite difference Jacobian' power_law_3d_jacobian.log > diff_pl_3d_jacobian.txt
