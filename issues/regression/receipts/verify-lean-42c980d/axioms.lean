import Enthymemata.ActiveSubspace

-- ProjectionMean
#print axioms ActiveSubspace.frobSq_eq_trace
#print axioms ActiveSubspace.IsOrthProj.isHermitian
#print axioms ActiveSubspace.IsOrthProj.one_sub
#print axioms ActiveSubspace.IsOrthProj.posSemidef
#print axioms ActiveSubspace.IsOrthProj.diag_mem
#print axioms ActiveSubspace.IsOrthProj.frobSq_sub
#print axioms ActiveSubspace.projMean_posSemidef
#print axioms ActiveSubspace.one_sub_projMean_posSemidef
#print axioms ActiveSubspace.trace_projMean
#print axioms ActiveSubspace.projLoss_eq
#print axioms ActiveSubspace.eigU_transpose_mul
#print axioms ActiveSubspace.eigU_mul_transpose
#print axioms ActiveSubspace.eq_eigU
#print axioms ActiveSubspace.trace_sub_two_trace_mul_eq
#print axioms ActiveSubspace.trace_sub_two_trace_mul_ge
#print axioms ActiveSubspace.summaryProj_isOrthProj
#print axioms ActiveSubspace.summaryProj_trace_sub
#print axioms ActiveSubspace.summaryLoss_le
#print axioms ActiveSubspace.IsOrthProj.isHermitian_compress
#print axioms ActiveSubspace.compress_eigenvalues_mem
#print axioms ActiveSubspace.half_frobSq_cross_rank

-- MatrixPoincare
#print axioms ActiveSubspace.loewner_iff
#print axioms ActiveSubspace.meanOf_add
#print axioms ActiveSubspace.meanOf_sub
#print axioms ActiveSubspace.meanOf_const_mul
#print axioms ActiveSubspace.meanOf_sum
#print axioms ActiveSubspace.meanOf_mono
#print axioms ActiveSubspace.covOf_isHermitian
#print axioms ActiveSubspace.energyOf_isHermitian
#print axioms ActiveSubspace.quadForm_eq
#print axioms ActiveSubspace.meanOf_sq_lin
#print axioms ActiveSubspace.covForm_eq_var
#print axioms ActiveSubspace.energyForm_eq
#print axioms ActiveSubspace.meanOf_mul_sub_coordMean
#print axioms ActiveSubspace.var_le_of_coord_poincare
#print axioms ActiveSubspace.coordwise_le_scalar
#print axioms ActiveSubspace.matrix_poincare
#print axioms ActiveSubspace.energy_box

-- DirectionMoves
#print axioms ActiveSubspace.mhRatio_of_symm
#print axioms ActiveSubspace.powerSpherical_symm
#print axioms ActiveSubspace.powerSpherical_antipode
#print axioms ActiveSubspace.change_ratio
#print axioms ActiveSubspace.change_detailed_balance
#print axioms ActiveSubspace.change_ratio_powerSpherical
#print axioms ActiveSubspace.prod_insertAt
#print axioms ActiveSubspace.birth_green_ratio
#print axioms ActiveSubspace.poisson_count_ratio
#print axioms ActiveSubspace.double_count_ratio
#print axioms ActiveSubspace.prior_draw_cancel
#print axioms ActiveSubspace.birth_death_balance

-- Refit
#print axioms ActiveSubspace.ridgeGrid_det
#print axioms ActiveSubspace.ridgeHinge_not_separable
#print axioms ActiveSubspace.hinge_isHingeSum
#print axioms ActiveSubspace.hinge_comp_not_hingeSum
#print axioms ActiveSubspace.hingeSum_not_closed
#print axioms ActiveSubspace.hingeFactor_neg
#print axioms ActiveSubspace.hingeFactor_scale
#print axioms ActiveSubspace.outer_pullback
#print axioms ActiveSubspace.activeMatrix_pullback
#print axioms ActiveSubspace.pullback_eigen_of_orthogonal
#print axioms ActiveSubspace.conformal_pullback
#print axioms ActiveSubspace.pullback_ker
#print axioms ActiveSubspace.pullback_range
#print axioms ActiveSubspace.pullback_ex
#print axioms ActiveSubspace.leading_C_coord
#print axioms ActiveSubspace.leading_pullback_coord
#print axioms ActiveSubspace.leading_moves
#print axioms ActiveSubspace.hat_invariant
#print axioms ActiveSubspace.rawRidge_moves
#print axioms ActiveSubspace.trainingRange_factor
#print axioms ActiveSubspace.trainingRange_hinge
#print axioms ActiveSubspace.kappa_moves
#print axioms ActiveSubspace.kappa_ratio
