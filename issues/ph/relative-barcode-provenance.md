# Is the relative barcode table new?

**Date:** 2026-11-02. **Asked by the user:** whether the fixed-K₀ table is new or folklore.

**Answer: folklore in substance. Don't present it as new.** I didn't find it written down as a standalone statement, but it
is a short specialization of published results, and software has computed the same bars for years. A machine-checked proof
in Lean may be the only new part, and that search wasn't exhaustive either.

## The statement

With K₀ strictly first in the filtration order, one absolute reduction gives the barcode of H(K_δ, K₀):

| Absolute bar | Relative bar |
|---|---|
| Born outside K₀ | The same bar |
| Born in K₀, killed outside at δ | [δ, ∞), one degree up |
| Born in K₀, killed within it, or essential | None |

## Where it comes from

| Source | What it has | What it lacks |
|---|---|---|
| Cohen-Steiner, Edelsbrunner, Harer, Morozov, *Persistent homology for kernels, images, and cokernels* (SODA 2009) | Rules for the fixed-subcomplex setting L_t = K_t ∩ L, from reductions with L's rows first. | States kernels and cokernels, not the relative module. |
| Lentz, Henselman-Petrusek, Ziegelmeier, [*A U-match algorithm for persistent relative homology*](https://arxiv.org/abs/2602.03163) (2026), §6 | For "lag filtrations" G_t = F_{t-l}, one R = DV gives bases for the relative cycle and boundary spaces, by the same mechanism as here. | Doesn't state a constant first subcomplex, or a table against the absolute barcode. |
| de Silva, Morozov, Vejdemo-Johansson, *Dualities in persistent (co)homology* (2011), Prop 2.4 | The analogous table for the moving module H(Xₙ, Xᵢ). | The fixed case. |
| Edelsbrunner, Morozov, *Persistent homology* (handbook chapter) | The same moving-case table, with generators. | The fixed case. |
| [Dionysus 2](https://mrzv.org/software/dionysus2/tutorial/basics.html) `relative=`, [TamerOp](https://www.tamerop.com/implementation/relative_persistence.html) "graded quotient" | Compute relative persistence against a fixed subcomplex by deleting its cells and reducing the rest. So `finrank_quotient_persistence` is standard practice. | Neither relates the output to the absolute barcode. |
| Blaser, Brun, [*Relative persistent homology*](https://arxiv.org/abs/1911.07484) (2019) | Relative persistence of Čech-type pairs. | Not this statement, as far as its abstract shows. |

**The specialization.** I checked this against the 2009 paper's rules. Put L first, so the reordered matrix they call R_im is just R_f:

- *Kernel birth:* σ ∈ K−L, negative in R_f, with its low in L. That is the second row.
- *Kernel death:* needs τ ∈ L that is negative in R_g and positive in R_f. With L first that can't happen, so these bars are essential.
- *Cokernel birth:* positive in R_f and in K−L. The rule's other case, σ ∈ L negative in R_g, can't also be positive in R_f once L is first. That is the first row.

**Why these two pieces make up the relative module.** The long exact sequence of the pair gives

0 → coker(H_p(K₀) → H_p(K_δ)) → H_p(K_δ, K₀) → ker(H_{p-1}(K₀) → H_{p-1}(K_δ)) → 0.

**The one point that takes thought.** The sequence alone doesn't say it splits as persistence modules. A cokernel bar [b, δ)
could in principle join a kernel bar [δ, ∞) into [b, ∞). The table says they never do, because no cell both kills and is born
(`column_zero_of_isLow`). I didn't find this remark written down, but it follows at once from standard facts, so I wouldn't
treat it as a result.

## Formalization

| Work | Prover | What it proves |
|---|---|---|
| Heras, Coquand, Mörtberg, Siles, [*Computing persistent homology within Coq/SSReflect*](https://dl.acm.org/doi/10.1145/2528929) (2013) | Coq | Persistent Betti numbers and the Fundamental Lemma, with certified programs. From what I read, not the pairing of a column reduction, and not relative persistence. |
| `lean/Lemmas/RelativeHomology.lean` and `RelativePersistence.lean` | Lean 4 | This repo: reduction computes persistent homology, absolute and relative, over any field and for ℤ/2 simplicial chains. |

I found no other formalization, and nothing in Lean.

## What the search didn't reach

So "not found" is weak evidence. These could state it outright:

- Edelsbrunner and Harer's *Computational Topology*, whose text I didn't read.
- Oineus's [relative homology page](https://oineus.readthedocs.io/en/stable/topics/relative_homology.html). The fetch timed out.
- Morozov's thesis, PHAT's documentation, and Pokorny et al. (2016) on trajectory classification, which use fixed-subcomplex relative persistence.
- Filtered mapping-cone results in the Floer literature. One search turned up Usher and Zhang's Floer-type complexes, with no such formula on the pages I saw.
