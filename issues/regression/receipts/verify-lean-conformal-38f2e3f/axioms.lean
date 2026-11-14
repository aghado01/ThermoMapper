import Enthymemata.ConformalValidity
open Conformal

-- Slice 14's follow-up 1: the weighted predictive distribution's PIT
#print axioms pit_weighted
#print axioms sum_mul_wPitCdf
#print axioms wPit_estimated

-- Slice 14's follow-up 2: the cross-fitted forms under weights
#print axioms sum_weight_strange_le
#print axioms strangeWeight_swap_ge
#print axioms crossFit_weighted_swap_tv
#print axioms crossFit_weighted_exchangeable
#print axioms strangeWeight_of_jackRefusesW
#print axioms jackknifePlus_weighted
#print axioms tournament_weighted
#print axioms strangeWeight_comp
#print axioms crossFit_likelihoodRatio
#print axioms jackknifePlus_likelihoodRatio
#print axioms tournament_likelihoodRatio
#print axioms fold_refused_of_refused
#print axioms cvPlus_oneFold_likelihoodRatio
