# Action difference in swaps: decomposition into $P_1$, $P_2$

This note proves that the action difference of a swap between two replicas can be written as a combination of **only two
sums per configuration**, $P_1$ and $P_2$, which depend neither on the replica pair nor on the replica index.
It is the basis of `swap_plaquette_sums` (host path), `swap_plaquette_sums_replica` (GPU path; its implementation is
described in [swap_translation_gpu.md](swap_translation_gpu.md)) and `swap_group` in `lib/gauge_conf_paral_temp.c`.

## 1. Definitions (as in the code)

- $\mu$ = `param->d_defect_dir`, the defect direction.
- Set of defect links (depends only on the geometry, `geo->d_defect_mask`, `is_on_defect`):
  $$\mathcal D = \{(x,\mu) \;:\; x_\mu = L_\mu - 1,\ \ x_\nu < L^{\rm def}_\nu \ \ \forall \nu \neq \mu\}.$$
  Only links in direction $\mu$ can lie on the defect.
- Link factor for the coefficient $c$ (`link_defect_coeff`):
  $$k_c(x,\nu) = \begin{cases} c & (x,\nu) \in \mathcal D \\ 1 & \text{otherwise.}\end{cases}$$
- Plaquette factor (`plaquette_defect_coeff`): for the plaquette $p = (x;\, \alpha, \beta)$ with links
  $(x,\alpha)$, $(x+\hat\alpha,\beta)$, $(x+\hat\beta,\alpha)$, $(x,\beta)$,
  $$K_p(c) = k_c(x,\alpha)\, k_c(x+\hat\alpha,\beta)\, k_c(x+\hat\beta,\alpha)\, k_c(x,\beta).$$
- Plaquette weight of the configuration $X$ (links $U_X$ and twist factors $Z_X$, which travel with the configuration),
  computed by `plaquettep`:
  $$w_p(X) = \frac{1}{N_c}\,\mathrm{Re\,Tr}\big(Z_{X,p}\, \Pi_p[U_X]\big).$$
- Action of the configuration $X$ when it has the coefficient $c$:
  $$S_c[X] = \beta \sum_p K_p(c)\,\big(1 - w_p(X)\big).$$

The replicas $a = 0,\dots,N_r-1$ have coefficients $c_a$ (`d_pt_bound_cond_coeff[a]`). The swap of the pair $(a, b)$
exchanges the replica indices of the configurations $A$ (index $a$) and $B$ (index $b$); the configurations do not
move.

## 2. Action difference of a swap

$$\Delta S = \big(S_{c_a}[B] + S_{c_b}[A]\big) - \big(S_{c_a}[A] + S_{c_b}[B]\big).$$

The constant terms $\beta \sum_p K_p(c)$ depend only on $c$: $c_a$ and $c_b$ each appear once with the $+$ sign and
once with the $-$ sign, so they cancel. What remains is
$$\Delta S = -\beta \sum_p \Big[K_p(c_a)\,w_p(B) + K_p(c_b)\,w_p(A) - K_p(c_a)\,w_p(A) - K_p(c_b)\,w_p(B)\Big]
= \beta\,\big(D_A - D_B\big),$$
$$D_X(c_a, c_b) = \sum_p \big[K_p(c_a) - K_p(c_b)\big]\, w_p(X).$$

This is the formula in the comment to `swap_group`. Only plaquettes with $K_p(c_a) \neq K_p(c_b)$ contribute, i.e.
those with at least one link in $\mathcal D$. They all lie in the $(\mu, j)$ planes, $j \neq \mu$, with base site $r$
in `swap_rect` (the rectangle with a border of 1 around the defect), and the code counts each of them once as
`plaquettep(r, mu, j)`. The sum over `swap_rect` is therefore the same before and after the change.

## 3. Lemma: $K_p(c) = c^{\,n_p}$ with $n_p \in \{0, 1, 2\}$ independent of $c$

Let $n_p$ be the number of links of the plaquette $p$ that lie in $\mathcal D$ (`plaquette_defect_links`). Each of the
four factors of $K_p(c)$ equals $c$ if its link lies in $\mathcal D$ and $1$ otherwise, hence
$$K_p(c) = c^{\,n_p}.$$
$n_p$ depends only on $p$ and $\mathcal D$, i.e. on the geometry, not on $c$ or on the configuration.

Bound $n_p \le 2$: $\mathcal D$ contains only links in direction $\mu$. A plaquette in the $(\alpha,\beta)$ plane,
with $\alpha \neq \beta$, has two links in direction $\alpha$ and two in direction $\beta$, and at most one of the two
directions is $\mu$. Therefore:
- planes without $\mu$: $n_p = 0$;
- $(\mu, j)$ planes: only $(r,\mu)$ and $(r+\hat j,\mu)$ can lie on the defect, hence $n_p \le 2$. $\square$

## 4. Theorem: decomposition

Grouping the plaquettes by $n_p$ and using the lemma:
$$D_X(c_a, c_b) = \sum_p \big(c_a^{\,n_p} - c_b^{\,n_p}\big)\, w_p(X)
= \sum_{n=0}^{2} \big(c_a^{\,n} - c_b^{\,n}\big) \sum_{p:\, n_p = n} w_p(X).$$
The $n = 0$ term vanishes ($c^0 = 1$). Defining
$$P_n(X) = \sum_{p:\, n_p = n} w_p(X), \qquad n = 1, 2,$$
one obtains, **for every pair of coefficients**:
$$\boxed{\,D_X(c_a, c_b) = (c_a - c_b)\,P_1(X) + (c_a^2 - c_b^2)\,P_2(X)\,}$$
and therefore
$$\Delta S = \beta\Big[(c_a - c_b)\big(P_1(A) - P_1(B)\big) + (c_a^2 - c_b^2)\big(P_2(A) - P_2(B)\big)\Big]. \qquad\square$$

## 5. Consequences for the code

1. **Independence of the pair.** $P_1(X)$ and $P_2(X)$ depend only on the configuration $X$ (links and $Z$) and on the
   geometry. They do not depend on $c_a$, $c_b$, the pair, or the replica index of $X$. Two numbers per configuration
   (per slot `g`) are therefore enough, computed with a reduction on the device of the configuration.
2. **A single computation for both groups.** During `swap()` the configurations do not change: only
   `replica_index` and `replica_to_gc` change. The second group, which uses the indices permuted by the first, uses the
   same $P_n$: `P1[replica_to_gc[a]]` is read with the updated permutation. The device does not need the
   `replica_index` values during the swap.
3. **Same Markov chain.** $\Delta S$ is the same real number as before, so the Metropolis probability
   $\min(1, e^{-\Delta S})$, detailed balance and the `rng_host` sequence are unchanged. The multicanonical term
   `delta_topo_potential_swap` is added as before.
4. **Rounding.** The $K_p$ are the same floating-point numbers as before: $c\cdot 1\cdot c\cdot 1$ in
   `plaquette_defect_coeff` is exactly $c\cdot c$, because multiplying by 1 is exact. Only the summation order changes,
   together with the distributivity $(c_a - c_b)\sum w$ instead of $\sum (K_a - K_b)\, w$. $\Delta S$ therefore differs
   by a relative $O(\varepsilon)$ ($\varepsilon \approx 10^{-16}$). A Metropolis decision changes only if the random
   number falls within $O(\varepsilon)$ of the threshold, which is negligible in practice. The previous OpenMP
   reduction did not fix a "canonical" summation order either. The summation order of the GPU block sums is described
   in [swap_translation_gpu.md](swap_translation_gpu.md), section 2.5.
5. **Cost.** Previously each configuration was summed over `swap_rect` once for every group it took part in, about
   twice per swap. Now it is summed only once, and only plaquettes with $n_p \ge 1$ are accumulated.

## 6. Remark: defect extended to the whole plane

If $L^{\rm def}_\nu = L_\nu$ for every $\nu \neq \mu$, every $(\mu, j)$ plaquette with $r_\mu = L_\mu - 1$ has both
links in direction $\mu$ on the defect, even with periodic boundary conditions in $j$. Hence $n_p = 2$ always,
$P_1 = 0$ and $\Delta S = \beta (c_a^2 - c_b^2)\big(P_2(A) - P_2(B)\big)$. With a partial defect, $n_p = 1$ appears
only on the border of the defect ($r_j = L^{\rm def}_j - 1$ or $r_j = L_j - 1$).

## 7. Numerical check

`src/debug_device_data.c`, `test_swap_sums`: for each slot, with random links copied to the device, it compares
$(c_a - c_b)P_1 + (c_a^2 - c_b^2)P_2$ (sums computed on the device) with the direct host sum of `swap_action_site`
over `swap_rect`. The pairs used are $(c_0, c_1)$ from the input and an arbitrary pair $(0.3, 1.7)$; the relative
tolerance is $10^{-10}$.
