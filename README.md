# Five Honest Questions About Sums of Three Cubes

Hi — I'm Nihar. Before you scroll any further, I want to tell you what this
repository actually is, because I think that's more useful than a fancy
tagline: it's five small, self-contained C++ programs, each one asking a
single honest question about James Maynard's 2026 JLMS survey *"Sums of
three positive cubes,"* and each one answering that question with real
numbers instead of hand-waving. Some of the answers are exciting. One of
them is a rigorous "no," and I've kept it exactly as it came out, because a
well-measured "no" is worth more than a hopeful "maybe."

Nothing here reimplements Davenport's, Vaughan's, or Wooley's proofs. I
didn't want to walk the same road they already paved so beautifully — I
wanted to walk five *different* roads that happen to arrive next to theirs,
using tools number theory doesn't usually reach for: extreme value theory,
equidistribution, elliptic curves, ecology's species-counting statistics,
and formal model-selection regression. If that sounds like an unusual guest
list, stick around — I think you'll enjoy meeting all five.

## The one paragraph you need before any of this makes sense

Maynard's paper studies `R(N)` — how many whole numbers up to `N` can be
written as `a³+b³+c³` for positive whole numbers `a,b,c`. Every proof in the
paper (Davenport 1950 → Vaughan 1985/1989 → Wooley 2015) is an *existence*
argument: it shows `R(N)` is *at least* `N` raised to some exponent, using a
Cauchy-Schwarz trick applied to the "recipe count" `ρ₁(n)` (how many ways a
given `n` can be written as such a sum). What none of those proofs ever do
is hand you an actual, checkable number, or explain *why* some numbers get
extra recipes while most get none. That gap is where all five extensions
live.

## How the dataset was chosen, and why

I wanted real data, not a synthetic stand-in, and I wanted it to be the
*exact same object* the paper studies — not something merely "related."
So every program downloads, automatically, for free, with no login and no
API key:

- **OEIS A003072** — "Numbers that are the sum of 3 positive cubes." This
  sequence's own reference list cites Davenport's original 1950 paper and
  Wooley's 2015 paper — literally the same two sources Maynard's survey
  leans on. I didn't have to reach for it; it was already the community's
  agreed-upon ledger of exactly this object.
- **OEIS A003325** — "Numbers that are the sum of 2 positive cubes," pulled
  in specifically for Extension 3, since that program's elliptic curves live
  on top of two-cube sums, not three.

Each program parses the downloaded file and cross-checks it against its own
from-scratch sieve. If your network can't reach `oeis.org` (sandboxed CI,
airplane wifi, corporate firewall — I've hit all three), the program
notices, says so honestly in the output, and keeps going on its own
computation. I'd rather a program tell you it's flying without a co-pilot
than pretend the co-pilot is there.

---

## Extension 1 — Does the paper's biggest hidden assumption actually hold?

**Why this extension.** Every one of the paper's proofs quietly assumes
that recipe-rich numbers are rare — that `ρ₁(n)` doesn't typically explode
whenever it's nonzero. The paper calls this the "enemy scenario" and admits,
in plain text, that nobody has ever ruled it out unconditionally. An
assumption that important deserves to be tested, not just trusted.

**The math.** Extreme value theory (Fisher & Tippett, 1928; Gnedenko, 1943)
is the branch of statistics built to answer "how large can the largest few
outcomes of a process plausibly get" — the same machinery insurers use for
100-year floods. I fit a Gumbel distribution to block maxima of `ρ₁(n)`, and
separately regressed the *running record* recipe-count against `N` on a
log-log scale to get a growth exponent.

**What we found.** A perfect match against the paper's own Table 1 on every
row, and a clean, independent match against 119,980 OEIS entries with zero
disagreements. The measured growth exponent came out to **0.178** (R²=0.96)
— comfortably above Mahler's proven floor of 1/12≈0.083, and clearly
inconsistent with Hardy-Littlewood's original (and long-disproved)
guess of no growth at all. I also computed, for the first time at a real,
finite N, exactly how tight the paper's core Cauchy-Schwarz inequality is:
**86.5%** of the true answer, a number that literally cannot be extracted
from an asymptotic theorem.

**What this means.** The enemy scenario didn't show up at this scale, and
we now have an honest, numeric feel for how loose the paper's central
inequality really is — something the paper itself never states, because its
language (`N^{β-ε}`) structurally can't.

## Extension 2 — What if we sieve with a completely different comb?

**Why this extension.** Every restriction set in the paper — Davenport's,
Vaughan's — is built from *divisibility* (prime factors, in one form or
another). That's one legitimate way to filter triples, but it's not the
only one, and it's also why the paper can never state an explicit constant:
divisibility-based arguments live entirely in asymptotic land.

**The math.** I built a sieve based on *equidistribution* instead — the
three-distance (Steinhaus) theorem and Beatty sequences (Van Ravenstein,
1988; O'Bryant, 2002), using the golden ratio conjugate `φ=(√5−1)/2`, the
irrational number famous for spreading itself as evenly as possible.

**What we found.** A fully explicit, hand-computed lower bound: **at least
249,812** makeable numbers below 4,000,000 (true answer: 411,818), captured
with every constant literally written down — 60.7% of the truth. Honest
footnote: the best setting we tested was the widest one on our grid, so we
proved you *can* get a real, checkable number this way, more than we proved
we'd found a clever sweet spot.

**What this means.** The paper can tell you an exponent for `N→∞`; it can
never tell you "at N=4,000,000, you're guaranteed at least 249,812." Now
there's a number sitting right there for anyone who wants one.

## Extension 3 — Can we crank out new recipes algebraically?

**Why this extension.** The paper never explains why a handful of numbers
get bonus recipes (the way 1729 famously does) while most get none or one.
Elliptic curves give a real algebraic "crank" — the same family of tools
Elkies (2000) and Elsenhans & Jahnel (2009) used to computationally crack
this exact equation for numbers like 33.

**The math.** If `u³+v³=m`, the point `(X,Y) = (12m/(u+v), 36m(u−v)/(u+v))`
sits on the elliptic curve `Y²=X³−432m²`. I derived this by hand, verified
it against the famous taxicab identity `1³+12³=1729=9³+10³`, and then
turned the group-law "doubling" crank on hundreds of real recipes, checking
every result with exact 128-bit integer arithmetic — never a floating-point
guess.

**What we found.** The self-test passed perfectly on both known recipes for
1729. Then, across 963 real seed recipes, the crank produced **zero** new
whole-number recipes — every single doubled point was a mathematically
valid point on the curve, just not one that lands on the integer lattice
(mean denominator size: under 2 digits, so genuinely *close*, just not
whole).

**What this means.** This is a rigorous negative result, and it's a good
one: it's exactly what Siegel's 1929 finiteness theorem predicts (integer
points on an elliptic curve are always scarce), and it tells us something
the paper never says out loud — Mahler's bonus-recipe construction *had* to
be a clever explicit polynomial identity, because simple algebraic
cranking, on its own, structurally can't manufacture the extra recipes.

## Extension 4 — Can ecology's "how many species are out there" trick count our numbers?

**Why this extension.** Ecologists estimate a rainforest's total species
count from one limited sample (Chao, 1984; Good, 1953 — the same statistics
behind "how many words did Shakespeare know but never use," Efron &
Thisted, 1976). `R(N)` is structurally the same question — "how many
distinct outcomes exist, given only a partial look" — asked about an
equation instead of a jungle.

**The math.** I searched only part of the triple-space (10%, 20%, …, 100%
of the way out) and applied the classical richness estimators to each
partial slice, comparing every guess against the true, fully-computed
answer.

**What we found.** Two honest results, not one. First, a clean proof (then
confirmed by direct computation, in all 10 rows): the "doubletons" count
these formulas divide by is **provably, always, exactly zero** for this
problem — recipes only ever arrive in bundles of 1, 3, or 6, and no
combination of those ever totals exactly 2. Second, an accuracy curve: at
10% sampled, the classical guess was off by **−99.8%**; by 100% sampled, it
was within about 1%.

**What this means.** The formulas assume something like random sampling;
our search is a deterministic march outward, which leaves an entire
unexplored region invisible to the statistics. That's not the formulas
failing — it's a genuine, quantified answer to "how far does ecology's
trick travel outside its home turf," which nobody had asked of this problem
before.

## Extension 5 — Where is that percentage actually heading?

**Why this extension.** The paper's own Table 1 shows `R(N)/N` drifting
down toward a conjectured 9.99425%, but the paper never fits a curve to
that drift or predicts responsibly past its own table, which stops at
`10^11`.

**The math.** A parallel, bit-packed sieve (in the engineering spirit of
Lagarias-Miller-Odlyzko's and Deléglise-Rivat's large-scale counting-
function work, applied here to a different function) computed real values
out to a billion, and then three competing curve shapes were fit and
ranked using Akaike's (1974) and Schwarz's (1978) formal model-selection
criteria — the same statistical toolkit used to choose between competing
weather-forecast models.

**What we found.** Eight-for-eight exact matches against the paper's table,
including the billion row. The winning model's predicted final destination:
**0.099709**. The paper's independently-conjectured value: **0.099943** — a
gap of about **0.23%**, reached with pure curve-fitting and zero cube-sum
theory whatsoever.

**What this means.** That's a genuinely pleasing independent confirmation —
but it's an honest 0.23% gap, not a bullseye, and I'd rather tell you that
plainly than round it up to sound more impressive than it is.

---

## Running it yourself

Each program is one self-contained `.cpp` file. Build and run any of them
with:

```bash
g++ -O3 -std=c++17 -pthread -o ext1 ext1_extreme_value_cube_multiplicity.cpp && ./ext1
```

Or use `run.sh` to build, run, and clean up all five in one shot, leaving
only the `results_*.txt` files behind:

```bash
chmod +x run.sh
./run.sh
```

## Citations

Every claim of "this technique comes from field X" above is backed by a
real, peer-reviewed or foundational source, cited in full inside each
program's header comment:

- R. A. Fisher & L. H. C. Tippett, *Limiting forms of the frequency
  distribution of the largest or smallest member of a sample*, Proc.
  Cambridge Phil. Soc. 24 (1928), 180–190.
- B. V. Gnedenko, *Sur la distribution limite du terme maximum d'une série
  aléatoire*, Annals of Mathematics 44(3) (1943), 423–453.
- L. de Haan & A. Ferreira, *Extreme Value Theory: An Introduction*,
  Springer, 2006.
- T. van Ravenstein, *The three gap theorem (Steinhaus conjecture)*, J.
  Austral. Math. Soc. Ser. A 45(3) (1988), 360–370.
- K. O'Bryant, *A generating function technique for Beatty sequences and
  other step sequences*, Journal of Number Theory 94 (2002), 299–319.
- V. Beresnevich & N. Leong, *Sums of reciprocals and the three distance
  theorem*, arXiv:1712.03758, 2017.
- N. D. Elkies, *Rational points near curves and small nonzero |x³−y²| via
  lattice reduction*, ANTS-IV, LNCS 1838, Springer, 2000, 33–63.
- A.-S. Elsenhans & J. Jahnel, *New sums of three cubes*, Mathematics of
  Computation 78(266) (2009), 1227–1230.
- D. R. Heath-Brown, W. M. Lioen, H. J. J. te Riele, *On solving the
  Diophantine equation x³+y³+z³=k on a vector computer*, Mathematics of
  Computation 61(203) (1993), 235–244.
- C. L. Siegel, *Über einige Anwendungen diophantischer Approximationen*,
  Abh. Preuss. Akad. Wiss. Phys.-Math. Kl. (1929), 41–69.
- A. Chao, *Nonparametric estimation of the number of classes in a
  population*, Scandinavian Journal of Statistics 11(4) (1984), 265–270.
- I. J. Good, *The population frequencies of species and the estimation of
  population parameters*, Biometrika 40(3-4) (1953), 237–264.
- B. Efron & R. Thisted, *Estimating the number of unseen species: How many
  words did Shakespeare know?*, Biometrika 63(3) (1976), 435–447.
- H. Akaike, *A new look at the statistical model identification*, IEEE
  Transactions on Automatic Control 19(6) (1974), 716–723.
- G. Schwarz, *Estimating the dimension of a model*, Annals of Statistics
  6(2) (1978), 461–464.
- M. Deléglise & J. Rivat, *Computing π(x): the Meissel, Lehmer, Lagarias,
  Miller, Odlyzko method*, Mathematics of Computation 65(213) (1996),
  235–245.

And, of course, the survey that started this whole conversation:
J. Maynard, *Sums of three positive cubes*, J. London Math. Soc. (2)
113 (2026), e70554.

## Disclaimer

This repository is an independent, self-contained research exploration. It
does not reproduce, depend on, or require any component of Maynard's
survey or the proofs of Davenport, Vaughan, or Wooley to function — every
program stands on its own, built from scratch. Nothing here claims to
improve on, or compete with, those theorems; the five results above are
computational and statistical complements sitting *alongside* deep proofs,
not replacements for them. Where a finding is a rigorous "no" (Extension 3)
or an honest, non-zero gap (Extensions 2 and 5), I've reported it exactly
as measured, because that's the only kind of result worth publishing.

Created by **Nihar Mahesh Jani** — reach me at
**niharmaheshjani@gmail.com**.

If you found a bug, a sharper citation, or a better idea for a sixth
question worth asking — I'd genuinely love to hear it.
