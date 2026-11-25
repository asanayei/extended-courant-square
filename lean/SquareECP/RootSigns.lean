import Mathlib

/-!
# Roots in `(-1, 1)` and corner signs

The elementary ingredient of the corner cap in the paper (Lemma "Roots and
corner signs"). Let `p` be a real polynomial that does not vanish at `±1`. If the sign of
`p 1` differs from the sign of the leading coefficient, then `p` has a root in `(1, ∞)`.
Likewise at `-1`, comparing with the sign of `p` at `-∞`. So each sign mismatch costs one
root inside `(-1, 1)`.
-/

open Polynomial Filter

namespace SquareECP.RootSigns

/-- If `p 1` and the leading coefficient have opposite signs, `p` has a root `> 1`. -/
theorem exists_root_gt_one (p : ℝ[X]) (h : p.eval 1 * p.leadingCoeff < 0) :
    ∃ x, 1 < x ∧ p.IsRoot x := by
  set c := p.leadingCoeff
  have hc : c ≠ 0 := by
    rintro hc; rw [hc, mul_zero] at h; exact lt_irrefl _ h
  have hdeg : 0 < p.natDegree := by
    by_contra h0
    push Not at h0
    have : p = C c := by
      rw [eq_C_of_natDegree_eq_zero (Nat.le_zero.mp h0)]
      simp [c, leadingCoeff, Nat.le_zero.mp h0]
    rw [this, eval_C] at h
    nlinarith [mul_self_nonneg c]
  set q := C c * p
  have hq_eval : ∀ x, q.eval x = c * p.eval x := fun x => by simp [q]
  have hq_deg : 0 < q.degree := by
    rw [degree_C_mul hc]
    exact natDegree_pos_iff_degree_pos.mp hdeg
  have hq_lc : 0 ≤ q.leadingCoeff := by
    rw [leadingCoeff_mul, leadingCoeff_C]; exact mul_self_nonneg c
  have hlim := tendsto_atTop_of_leadingCoeff_nonneg (P := q) hq_deg hq_lc
  obtain ⟨x, hxpos, hx1⟩ :
      ∃ x, 0 < q.eval x ∧ 2 ≤ x :=
    ((hlim.eventually (eventually_gt_atTop 0)).and (eventually_ge_atTop 2)).exists
  have hq1 : q.eval 1 < 0 := by rw [hq_eval]; linarith [mul_comm c (p.eval 1)]
  have hcont : ContinuousOn (fun t => q.eval t) (Set.Icc 1 x) :=
    q.continuous.continuousOn
  obtain ⟨y, hy, hy0⟩ :=
    intermediate_value_Ioo (by linarith) hcont ⟨hq1, hxpos⟩
  refine ⟨y, hy.1, ?_⟩
  have : c * p.eval y = 0 := by rw [← hq_eval]; exact hy0
  exact (mul_eq_zero.mp this).resolve_left hc

/-- If `p (-1)` and the sign of `p` at `-∞` differ, `p` has a root `< -1`. -/
theorem exists_root_lt_neg_one (p : ℝ[X])
    (h : p.eval (-1) * ((-1) ^ p.natDegree * p.leadingCoeff) < 0) :
    ∃ x, x < -1 ∧ p.IsRoot x := by
  set q := p.comp (-X)
  have hX : (-X : ℝ[X]).natDegree = 1 := by simp
  have hlc : q.leadingCoeff = (-1) ^ p.natDegree * p.leadingCoeff := by
    rw [leadingCoeff_comp (by rw [hX]; exact one_ne_zero)]
    simp [mul_comm]
  have hq : q.eval 1 * q.leadingCoeff < 0 := by
    rw [hlc]; simpa [q] using h
  obtain ⟨x, hx, hroot⟩ := exists_root_gt_one q hq
  refine ⟨-x, by linarith, ?_⟩
  simpa [q, IsRoot] using hroot

/-- **Root lemma.** Roots in `(-1, 1)` (with multiplicity), plus one for each corner-sign
mismatch, is at most the degree. -/
theorem card_roots_Ioo_add_le (p : ℝ[X]) (hp : p ≠ 0) :
    (p.roots.filter (fun x => -1 < x ∧ x < 1)).card
      + (if p.eval 1 * p.leadingCoeff < 0 then 1 else 0)
      + (if p.eval (-1) * ((-1) ^ p.natDegree * p.leadingCoeff) < 0 then 1 else 0)
      ≤ p.natDegree := by
  classical
  set s := p.roots
  set A := s.filter (fun x => -1 < x ∧ x < 1)
  set R := s.filter (fun x => 1 < x)
  set L := s.filter (fun x => x < -1)
  -- the three filters are disjoint, so their cards add up to at most `card s`
  have hsum : A.card + R.card + L.card ≤ s.card := by
    have h1 := Multiset.filter_add_filter (fun x : ℝ => -1 < x ∧ x < 1) (fun x => 1 < x) s
    have h2 := Multiset.filter_add_filter
      (fun x : ℝ => (-1 < x ∧ x < 1) ∨ 1 < x) (fun x => x < -1) s
    have e1 : s.filter (fun x : ℝ => (-1 < x ∧ x < 1) ∧ 1 < x) = 0 :=
      Multiset.filter_eq_nil.mpr fun x _ hx => by linarith [hx.1.2, hx.2]
    have e2 : s.filter (fun x : ℝ => ((-1 < x ∧ x < 1) ∨ 1 < x) ∧ x < -1) = 0 :=
      Multiset.filter_eq_nil.mpr fun x _ hx => by
        rcases hx.1 with h | h <;> linarith [hx.2]
    rw [e1, add_zero] at h1
    rw [e2, add_zero] at h2
    have hle := Multiset.card_le_card
      (Multiset.filter_le (fun x : ℝ => ((-1 < x ∧ x < 1) ∨ 1 < x) ∨ x < -1) s)
    have c1 := congrArg Multiset.card h1
    have c2 := congrArg Multiset.card h2
    simp only [Multiset.card_add] at c1 c2
    simp only [A, R, L]
    omega
  have hR : (if p.eval 1 * p.leadingCoeff < 0 then 1 else 0) ≤ R.card := by
    split_ifs with h
    · obtain ⟨x, hx, hroot⟩ := exists_root_gt_one p h
      exact Multiset.card_pos_iff_exists_mem.mpr
        ⟨x, Multiset.mem_filter.mpr ⟨(mem_roots hp).mpr hroot, hx⟩⟩
    · exact Nat.zero_le _
  have hL : (if p.eval (-1) * ((-1) ^ p.natDegree * p.leadingCoeff) < 0 then 1 else 0)
      ≤ L.card := by
    split_ifs with h
    · obtain ⟨x, hx, hroot⟩ := exists_root_lt_neg_one p h
      exact Multiset.card_pos_iff_exists_mem.mpr
        ⟨x, Multiset.mem_filter.mpr ⟨(mem_roots hp).mpr hroot, hx⟩⟩
    · exact Nat.zero_le _
  have : s.card ≤ p.natDegree := card_roots' p
  omega

end SquareECP.RootSigns

/-! ## The ends cap for cubic sides (levels λ = 17, 18)

On `S(17)` and `S(18)` the only monomial with `Y`-degree 3 is `Y^3`, and the only one with
`X`-degree 3 is `X^3`. So the side restrictions `R = P(1, ·)`, `L = P(-1, ·)` share their
`t^3` coefficient, and so do `T = P(·, 1)` and `B = P(·, -1)`. Corner values agree. We
prove that the four sides together have at most 10 roots in the open sides, not 12.
-/

namespace SquareECP.RootSigns

/-- Number of roots in `(-1, 1)`, with multiplicity. -/
noncomputable def nIoo (p : ℝ[X]) : ℕ := (p.roots.filter (fun x => -1 < x ∧ x < 1)).card

theorem nIoo_le_natDegree (p : ℝ[X]) : nIoo p ≤ p.natDegree :=
  (Multiset.card_le_card (Multiset.filter_le _ _)).trans (card_roots' p)

/-- One side, cubic with `t^3` coefficient `c ≠ 0`. -/
theorem side_bound (p : ℝ[X]) (c : ℝ) (hdeg : p.natDegree ≤ 3) (hc : p.coeff 3 = c)
    (hc0 : c ≠ 0) :
    nIoo p + (if p.eval 1 * c < 0 then 1 else 0)
      + (if p.eval (-1) * (-c) < 0 then 1 else 0) ≤ 3 := by
  have hp : p ≠ 0 := by rintro rfl; simp at hc; exact hc0 hc.symm
  have hnd : p.natDegree = 3 := le_antisymm hdeg (le_natDegree_of_ne_zero (hc ▸ hc0))
  have hlc : p.leadingCoeff = c := by rw [leadingCoeff, hnd, hc]
  have := card_roots_Ioo_add_le p hp
  rw [hnd, hlc, show ((-1 : ℝ) ^ 3 * c) = -c by ring] at this
  simpa only [nIoo] using this

/-- One side whose `t^3` coefficient vanishes: at most 2 roots. -/
theorem side_bound_zero (p : ℝ[X]) (hdeg : p.natDegree ≤ 3) (hc : p.coeff 3 = 0) :
    nIoo p ≤ 2 := by
  have : p.natDegree ≤ 2 := by
    by_contra h
    push Not at h
    have h3 : p.natDegree = 3 := by omega
    have : p.leadingCoeff ≠ 0 := by
      intro h0; rw [leadingCoeff_eq_zero] at h0; simp [h0] at h3
    rw [leadingCoeff, h3, hc] at this; exact this rfl
  exact (nIoo_le_natDegree p).trans this

/-- A nonzero corner value contributes at least one mismatch for two leading signs that agree. -/
theorem corner_same (e c c' : ℝ) (he : e ≠ 0) (h : 0 < c * c') :
    1 ≤ (if e * c < 0 then 1 else 0) + (if e * (-c') < 0 then 1 else 0) := by
  split_ifs with h1 h2 h2 <;> try omega
  exfalso
  push Not at h1 h2
  have := mul_nonneg h1 h2
  nlinarith [mul_pos (mul_self_pos.mpr he) h]

/-- Same, for leading signs of opposite sign. -/
theorem corner_opp (e c c' : ℝ) (he : e ≠ 0) (h : c * c' < 0) :
    1 ≤ (if e * c < 0 then 1 else 0) + (if e * c' < 0 then 1 else 0) := by
  split_ifs with h1 h2 h2 <;> try omega
  exfalso
  push Not at h1 h2
  have := mul_nonneg h1 h2
  nlinarith [mul_self_pos.mpr he]

theorem corner_opp' (e c c' : ℝ) (he : e ≠ 0) (h : c * c' < 0) :
    1 ≤ (if e * (-c) < 0 then 1 else 0) + (if e * (-c') < 0 then 1 else 0) := by
  have := corner_opp e (-c) (-c') he (by nlinarith)
  exact this

/-- **Ends cap for cubic sides.** Four side cubics with shared leading coefficients on
opposite sides and consistent nonzero corner values have at most 10 roots in total on
the open sides. (Corners: `R 1 = T 1`, `R (-1) = B 1`, `L 1 = T (-1)`, `L (-1) = B (-1)`.) -/
theorem ends_cap_cubic (R L T B : ℝ[X]) (c c' : ℝ)
    (hR : R.natDegree ≤ 3) (hL : L.natDegree ≤ 3) (hT : T.natDegree ≤ 3) (hB : B.natDegree ≤ 3)
    (hRc : R.coeff 3 = c) (hLc : L.coeff 3 = c) (hTc : T.coeff 3 = c') (hBc : B.coeff 3 = c')
    (k1 : R.eval 1 = T.eval 1) (k2 : R.eval (-1) = B.eval 1)
    (k3 : L.eval 1 = T.eval (-1)) (k4 : L.eval (-1) = B.eval (-1))
    (n1 : R.eval 1 ≠ 0) (n2 : R.eval (-1) ≠ 0) (n3 : L.eval 1 ≠ 0) (n4 : L.eval (-1) ≠ 0) :
    nIoo R + nIoo L + nIoo T + nIoo B ≤ 10 := by
  by_cases hc0 : c = 0
  · have := side_bound_zero R hR (hRc.trans hc0)
    have := side_bound_zero L hL (hLc.trans hc0)
    have := (nIoo_le_natDegree T).trans hT
    have := (nIoo_le_natDegree B).trans hB
    omega
  by_cases hc0' : c' = 0
  · have := side_bound_zero T hT (hTc.trans hc0')
    have := side_bound_zero B hB (hBc.trans hc0')
    have := (nIoo_le_natDegree R).trans hR
    have := (nIoo_le_natDegree L).trans hL
    omega
  have sR := side_bound R c hR hRc hc0
  have sL := side_bound L c hL hLc hc0
  have sT := side_bound T c' hT hTc hc0'
  have sB := side_bound B c' hB hBc hc0'
  rw [← k1, ← k3] at sT
  rw [← k2, ← k4] at sB
  rcases lt_or_gt_of_ne (mul_ne_zero hc0 hc0') with hneg | hpos
  · -- opposite signs: corners (1,1) and (-1,-1) each force a mismatch
    have a1 := corner_opp (R.eval 1) c c' n1 hneg
    have a2 := corner_opp' (L.eval (-1)) c c' n4 hneg
    omega
  · -- same sign: corners (1,-1) and (-1,1) each force a mismatch
    have a1 := corner_same (R.eval (-1)) c' c n2 (by linarith [mul_comm c c'])
    have a2 := corner_same (L.eval 1) c c' n3 hpos
    omega

end SquareECP.RootSigns
