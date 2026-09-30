# Graphene s–p (and s–p–d) Wannier models

Companion to `../Graphene_wannier/`, which holds the 2-orbital p_z model.
Same cell, same SCF, same 15×15×1 k-mesh — only the projection stage differs,
so the three models are directly comparable.

| model | dir / seed | num_wann | status |
|---|---|---|---|
| p_z | `../Graphene_wannier/graphene` | 2 | excellent — exact at K, 2.8 meV within 4 eV of E_F |
| s,p | `graphene_sp` | 8 | usable; see caveats below |
| s,p,d | `graphene_spd` | 18 | best σ description |

Orbital ordering is **sublattice A first, then B**, in every basis:

```
sp :  WF 1-4  = A (s, pz, px, py)            WF 5-8   = B (same)
spd:  WF 1-9  = A (s, pz, px, py, dz2,       WF 10-18 = B (same)
                   dxz, dyz, dx2-y2, dxy)
```

so a site-resolved on-site potential is a block-diagonal perturbation and
`TBModel.with_staggered_onsite(delta)` works unchanged on all three.

Because the same Δ is applied to *every* orbital of a site, the perturbation is
Δ·I on that site's block, which is invariant under any unitary mixing within the
site. It therefore means the same thing whether the on-site orbitals came out as
clean s/p/d or as hybrids of them — only the **atom-centring** matters, not the
orbital purity.

---

## Why the sp model is hard (and the p_z model is not)

The failure is specific and diagnosable. In the projection gauge the eight sp
Wannier functions come out as:

```
A_pz    spread  0.985 Ang^2   <- perfect, sits exactly on the atom
A_s     spread 17.1   Ang^2
A_px    spread 27.5   Ang^2
A_py    spread 22.2   Ang^2   <- smeared over several unit cells
```

p_z works because its antibonding partner π\* is a real, reachable state at
+9.218 eV. s/p_x/p_y fail because *their* antibonding partners — the σ\* bands —
lie **above the vacuum level** and are shredded across the vacuum/interlayer
continuum of a 19.7 Å cell. You cannot build a localised atomic orbital when its
antibonding partner is not in the subspace.

This is why `exclude_bands` and energy windows cannot rescue it: σ\* and the
vacuum states *overlap in energy*, so no choice of `dis_froz_*` admits one
without the other.

## Things that do not work (checked, so they don't get retried)

| approach | result |
|---|---|
| maximal localisation (`num_iter > 0`) | rotates s,p into bond-pointing hybrids: 2 WFs land on C–C bonds, site assignment destroyed, **and** occupied bands get *worse* (26 → 34 meV) |
| `site_symmetry` (SAWF) | aborts: `symmetrize_ukirr: not converged` at every `symmetrize_eps`; also refuses a frozen window in this build |
| `slwf_constrain` | silently inert — wannier90 enables selective localisation only when `slwf_num < num_wann` (`wannier90_readwrite.F90:757`), so constraining all centres does nothing; every λ gave byte-identical output |
| `dis_mix_ratio = 1.0` | disentanglement oscillates (Ω_I bounced 22.09 ↔ 22.44 for 1400+ iterations). Use the 0.5 default for these larger manifolds |

## THE ANSWER: projectability to discard + energy window to pin

The working recipe is the **combination**. Neither piece works alone:

| configuration | spreads (Å²) | occupied-band error |
|---|---|---|
| energy window only, `dis_win_max = 11` | 0.99 – 27.5 | 26.2 meV |
| energy window only, `dis_win_max = 22` | 0.92 – 6.03 | 16.7 meV |
| projectability only (`dis_proj_min/max`) | **0.62 – 0.75** | 545 meV ✗ |
| **both together** | **0.62 – 0.94** | **0.52 meV** ✓ |

Projectability alone localises beautifully but *un-pins* the bands; the energy
window alone pins the bands but cannot exclude the vacuum continuum. Together:

```
dis_win_max       =  40.0     ! generous outer bound
dis_proj_min      =  0.01     ! DISCARD vacuum states (no atomic character)
dis_proj_max      =  1.0      ! freeze nothing on projectability alone
dis_froz_min      = -27.0     ! PIN the occupied bands by energy
dis_froz_max      =  0.73
dis_num_iter      =  600
dis_conv_tol      =  1.0e-9   ! needed: it converges to ~5e-10 and would
dis_mix_ratio     =  0.5      ! otherwise spin to dis_num_iter and get killed
num_iter          =  0        ! projection gauge -- keeps atomic character
```

### Final 8-orbital model

Ω_I = 3.554, Ω_D = 0.0159, Ω_OD = 2.292, **Ω_total = 5.862 Å²**

```
WF          x         y         z     spread   off-site
A_s   -0.00002   1.42116  -0.00000   0.74105    0.00088
A_pz  -0.00000   1.42028  -0.00000   0.94382    0.00000
A_px   0.00001   1.49123   0.00000   0.62693    0.07094
A_py  -0.00000   1.34806   0.00000   0.61907    0.07222
B_s    1.23002   0.70926   0.00000   0.74105    0.00088
B_pz   1.23000   0.71014   0.00000   0.94382    0.00000
B_px   1.22999   0.63920  -0.00000   0.62693    0.07094
B_py   1.23000   0.78236  -0.00000   0.61907    0.07222
```

- A/B spreads identical to **0.00e+00** for all four orbital types
- on-site energies degenerate between sublattices to **0.00e+00**
  (ε_s = −2.3211, ε_pz = −1.9145, ε_p∥ = +1.0814 / +0.8359 eV)
- occupied bands **0.52 meV mean, 11.3 meV max**; 5.16 meV within 6 eV of E_F
- eigenvalues at K reproduce DFT exactly on all five occupied states
- Dirac gap opens as exactly 2Δ (0.200000, 0.500000, 1.000000 eV for
  Δ = 0.1, 0.25, 0.5)

**d orbitals turned out to be unnecessary.** With projectability doing the
state selection, the plain 8-orbital s+p basis is both better localised and far
more accurate than the 18-orbital spd model (which gave 64.9 meV on occupied
bands and spreads out to 23.9 Å²). The d shell was compensating for a
disentanglement problem that is better solved directly.

### Known caveat: the in-plane p doublet

ε_px and ε_py come out split by 0.246 eV, and the on-site block carries a small
s–p∥ element (0.053 eV). Under D3h site symmetry that doublet should be exactly
degenerate, so this is a genuine (if small, ~2.5 % of the p bandwidth) symmetry
violation left by the unconstrained disentanglement — the two in-plane p
functions are not exactly p_x and p_y but two slightly inequivalent in-plane
combinations (spreads 0.62693 vs 0.61907).

**This does not affect a sublattice potential.** Applying the same Δ to all four
orbitals of a site is Δ·I on that block, invariant under any mixing within the
site, and the A↔B relation is exact to machine precision. It matters only if you
want orbital-resolved on-site energies that distinguish p_x from p_y.

## What does work

**1. Projection gauge, `num_iter = 0`.** Keeps exact atomic character on each
site and, for this manifold, gives *better* bands than maximal localisation.

**2. A wide outer window.** Pushing `dis_win_max` 11 → 22 eV let σ\* states into
the subspace: σ spreads fell 17–27 → 4.3–6.0 Å², occupied-band error 26 → 16.7
meV, A/B pairs matching to 0.01 Å². This is capped by `nbnd`, hence the 60-band
NSCF (`graphene_nscf60.in`).

**3. Projectability disentanglement** (`dis_proj_min` / `dis_proj_max`),
Qiao et al., [arXiv:2303.07877](https://arxiv.org/abs/2303.07877), supported in
wannier90 3.1. Instead of selecting states by *energy*, it selects them by how
much atomic-orbital *character* they carry — exactly the right criterion when
σ\* and vacuum states are energetically interleaved. This is the single biggest
improvement found: Ω_I dropped from 23.9 (oscillating) to 15.9 (converging
monotonically) for the spd model.

**4. d orbitals — helpful diagnostically, but not needed in the end.** With
s+p+d the projections reach atomic-character weight out to 26 eV (bands 1–28)
versus 7.8 eV for s+p alone, which confirmed that the antibonding manifold was
the missing ingredient. But once projectability disentanglement was doing the
state selection, the plain s+p basis outperformed spd on every measure (see
above), so the shipped model is 8 orbitals. Keep `graphene_spd` only if you
specifically want d character in the basis.

## Calibration against the literature

Worth knowing before chasing further accuracy: the NIST **JARVIS-WTB** database
([Sci. Data 8, 92 (2021)](https://www.nature.com/articles/s41597-021-00885-z),
1406 3D + 365 2D materials) treats a Wannier Hamiltonian as "useful for most
applications" at a **maximum band error below 0.1 eV**, and reports that only
**64 %** of their materials meet that on high-symmetry paths. Their protocol is
the same SMV disentanglement with a default ±2 eV frozen window, widened to
cover the valence bands when that fails.

By that yardstick the sp model here (≈17–29 meV mean on occupied bands) is
comfortably inside the accepted range, despite the alarming-looking spreads.

If you want an independent reference rather than another local run:
- **JARVIS-WTB**, <https://jarvis.nist.gov/jarviswtb> (login required); the
  `_hr.dat` files are distributed via Figshare and load with `jarvis-tools`'
  `WanHam` class.
- **Materials Cloud** high-throughput Wannierisation, Vitale et al.,
  [npj Comput. Mater. 6, 66 (2020)](https://arxiv.org/abs/1909.00433),
  DOI `10.24435/materialscloud:2019.0044/v2`, and the 21,737-material PDWF
  dataset from arXiv:2303.07877.

## Pipeline

The SCF and the 15×15×1 k-mesh are reused from `../Graphene_wannier`; `tmp` here
is a symlink to that scratch directory. Only the 60-band NSCF is new.

```bash
L=/home/haney411/wrk/qe/fix_snprintf.so
W90=/home/haney411/wrk/wannier90_source/wannier90/wannier90.x

LD_PRELOAD=$L mpirun -np 4 -x LD_PRELOAD=$L pw.x -in graphene_nscf60.in > graphene_nscf60.out

python3 make_basis.py spd 60 ./tmp60      # or: sp 60 ./tmp60
$W90 -pp graphene_spd
LD_PRELOAD=$L mpirun -np 4 -x LD_PRELOAD=$L pw2wannier90.x -in graphene_spd.pw2wan > graphene_spd_pw2wan.out
python3 analyze_windows_sp.py graphene_spd 18 graphene_nscf60.out
$W90 graphene_spd

python3 eval_sp.py graphene_spd           # centring + spreads
python3 compare_sp.py graphene_spd        # bands vs DFT
```

Note `mpirun -np 4` — OpenMPI offers only 4 slots on this machine, and `-np 6`
fails with a slot error that looks like a physics crash but isn't.

The two bases share their `.mmn` and `.eig`: the b-vectors are byte-identical
(verified by md5 of the `nnkpts` block), so only the `.amn` needs recomputing
when the projections change. That turns a ~4 min pw2wannier90 run into seconds.
