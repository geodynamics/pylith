#!/bin/bash
# Run all 3D models.
pylith quad.cfg drucker_prager_2d_solve.cfg 2>&1 | tee drucker_prager_2d_solve.log
pylith quad.cfg drucker_prager_2d_jacobian.cfg 2>&1 | tee drucker_prager_2d_jacobian.log
pylith quad.cfg drucker_prager_2d_tensile_solve.cfg 2>&1 | tee drucker_prager_2d_tensile_solve.log
pylith quad.cfg drucker_prager_2d_tensile_jacobian.cfg 2>&1 | tee drucker_prager_2d_tensile_jacobian.log
pylith quad.cfg drucker_prager_2d_nonassoc_solve.cfg 2>&1 | tee drucker_prager_2d_nonassoc_solve.log
pylith quad.cfg drucker_prager_2d_nonassoc_jacobian.cfg 2>&1 | tee drucker_prager_2d_nonassoc_jacobian.log
pylith quad.cfg elastic_2d_solve.cfg 2>&1 | tee elastic_2d_solve.log
pylith quad.cfg elastic_2d_jacobian.cfg 2>&1 | tee elastic_2d_jacobian.log
pylith quad.cfg gen_max_2d_solve.cfg 2>&1 | tee gen_max_2d_solve.log
pylith quad.cfg gen_max_2d_jacobian.cfg 2>&1 | tee gen_max_2d_jacobian.log
pylith quad.cfg maxwell_2d_solve.cfg 2>&1 | tee maxwell_2d_solve.log
pylith quad.cfg maxwell_2d_jacobian.cfg 2>&1 | tee maxwell_2d_jacobian.log
pylith quad.cfg power_law_2d_solve.cfg 2>&1 | tee power_law_2d_solve.log
pylith quad.cfg power_law_2d_jacobian.cfg 2>&1 | tee power_law_2d_jacobian.log

# Extract Jacobian info and output to separate files.
grep -F -A 18 'Explicit preconditioning Jacobian' drucker_prager_2d_jacobian.log > pylith_dp_2d_jacobian.txt
grep -F -A 18 'Explicit preconditioning Jacobian' drucker_prager_2d_tensile_jacobian.log > pylith_dp_2d_tensile_jacobian.txt
grep -F -A 18 'Explicit preconditioning Jacobian' drucker_prager_2d_nonassoc_jacobian.log > pylith_dp_2d_nonassoc_jacobian.txt
grep -F -A 18 'Explicit preconditioning Jacobian' elastic_2d_jacobian.log > pylith_el_2d_jacobian.txt
grep -F -A 18 'Explicit preconditioning Jacobian' gen_max_2d_jacobian.log > pylith_gm_2d_jacobian.txt
grep -F -A 18 'Explicit preconditioning Jacobian' maxwell_2d_jacobian.log > pylith_mx_2d_jacobian.txt
grep -F -A 18 'Explicit preconditioning Jacobian' power_law_2d_jacobian.log > pylith_pl_2d_jacobian.txt
grep -F -A 18 'Finite difference Jacobian' drucker_prager_2d_jacobian.log > finite_diff_dp_2d_jacobian.txt
grep -F -A 18 'Finite difference Jacobian' drucker_prager_2d_tensile_jacobian.log > finite_diff_dp_2d_tensile_jacobian.txt
grep -F -A 18 'Finite difference Jacobian' drucker_prager_2d_nonassoc_jacobian.log > finite_diff_dp_2d_nonassoc_jacobian.txt
grep -F -A 18 'Finite difference Jacobian' elastic_2d_jacobian.log > finite_diff_el_2d_jacobian.txt
grep -F -A 18 'Finite difference Jacobian' gen_max_2d_jacobian.log > finite_diff_gm_2d_jacobian.txt
grep -F -A 18 'Finite difference Jacobian' maxwell_2d_jacobian.log > finite_diff_mx_2d_jacobian.txt
grep -F -A 18 'Finite difference Jacobian' power_law_2d_jacobian.log > finite_diff_pl_2d_jacobian.txt
grep -F -A 18 'minus finite difference Jacobian' drucker_prager_2d_jacobian.log > diff_dp_2d_jacobian.txt
grep -F -A 18 'minus finite difference Jacobian' drucker_prager_2d_tensile_jacobian.log > diff_dp_2d_tensile_jacobian.txt
grep -F -A 18 'minus finite difference Jacobian' drucker_prager_2d_nonassoc_jacobian.log > diff_dp_2d_nonassoc_jacobian.txt
grep -F -A 18 'minus finite difference Jacobian' elastic_2d_jacobian.log > diff_el_2d_jacobian.txt
grep -F -A 18 'minus finite difference Jacobian' gen_max_2d_jacobian.log > diff_gm_2d_jacobian.txt
grep -F -A 18 'minus finite difference Jacobian' maxwell_2d_jacobian.log > diff_mx_2d_jacobian.txt
grep -F -A 18 'minus finite difference Jacobian' power_law_2d_jacobian.log > diff_pl_2d_jacobian.txt
