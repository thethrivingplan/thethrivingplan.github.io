# Three Stars, Four Stars — Real Newtonian Runs

*Run November 24, 2026, in answer to the @satisfactman TikTok ("Three stars.
Real Newtonian gravity. Which star gets ejected first?") — Tristan's order:
run the same test with three, and with four. These are real computed runs:
the integrator was written for this test and executed here. Script:
`nbody.py` in this folder; raw numbers in `results_raw.json`.*

## The setup, on the record

- **Physics:** pure Newtonian gravity, no softening, no tricks. G = 1,
  every star has mass 1, distances of order 1.
- **Integrator:** leapfrog (kick–drift–kick). Runs A and B use a fixed
  timestep; Run C uses the same leapfrog with an adaptive timestep (see
  the honesty notes for why).
- **Ejection means:** a star is more than 20 distance units from the
  barycentre of the others, its energy relative to them is positive
  (unbound), and it is still receding.
- **The honesty check:** total energy must not drift. Every run prints its
  energy drift; a run whose energy blows up is a broken run, not a result.

## Run A — three stars, generic (the video's test)

Initial conditions (centre-of-mass frame; a rough triangle, small
tangential velocities):

| Star | Position | Velocity |
|---|---|---|
| 0 | (+1.000, +0.010) | (+0.053, +0.307) |
| 1 | (−0.550, +0.860) | (−0.277, −0.113) |
| 2 | (−0.450, −0.870) | (+0.223, −0.193) |

Timestep dt = 0.0002.

**Outcome: Star 0 — the one that started at the right-hand point of the
triangle — was ejected at t = 7.36.** It left fast: speed 3.46 relative to
the remaining pair's centre as it crossed the 20-unit line, settling to an
asymptotic departure speed of **3.43**. The triple did not dance for long;
it had one hard encounter and threw a star.

Stars 1 and 2 (the two left-hand starters) stayed together as a hard,
tight binary: semi-major axis **0.091**, orbital period **0.122** time
units. The ejected star's speed and the pair's tightness are the same
event seen from two sides — the pair fell deep into its own well and the
kick that bound it is the kick that threw Star 0 out.

Energy drift: final **0.025%**, worst moment **0.057%** — inside the 0.1%
bar.

## Run B — three stars, the figure-eight (the control)

The Chenciner–Montgomery figure-eight, standard published initial
conditions: positions (−0.97000436, 0.24308753), (0.97000436, −0.24308753),
(0, 0); velocities (0.4662036850, 0.4323657300) for stars 1 and 2 and
(−0.9324073700, −0.8647314600) for star 3. Timestep dt = 0.0005, three
full periods (period ≈ 6.3259, total t ≈ 18.98).

**Outcome: it held.** Over all three periods, no star ever strayed more
than **0.00068** distance units from the shared figure-eight track — the
three chase each other around the same 8 forever, and nobody is ejected.
Energy drift: **0.00001%** worst moment — the integrator is not inventing
anything here.

This run is the control, and it matters: *three* does not automatically
mean chaos and ejection. A three can be a perfect, repeating weave when
the arrangement is exactly balanced. The generic three falls apart; the
balanced three holds. The difference is the arrangement, not the number.

## Run C — four stars, generic

Initial conditions: two loose pairs on a slow mutual approach (pair A =
stars 0,1 around (−1.8, 0); pair B = stars 2,3 around (1.8, 0.3); each
pair internally near-circular, pair B turning the opposite way), shifted
to the centre-of-mass frame:

| Star | Position | Velocity |
|---|---|---|
| 0 | (−2.200, −0.150) | (0.000, −1.241) |
| 1 | (−1.400, −0.150) | (0.000, +0.341) |
| 2 | (+1.350, +0.150) | (0.000, +1.195) |
| 3 | (+2.250, +0.150) | (0.000, −0.295) |

Adaptive leapfrog: dt = 0.002 × the closest pair's dynamical time
(η = 0.002, capped at dt = 0.01); 728,596 steps.

**Outcome: the four resolved as two ejections plus one exchanged pair.**

- **Star 3** (pair B's outer star) was ejected first, at **t = 72.57** —
  gently, departure speed settling to **0.52**.
- **Star 1** (pair A's second star) was ejected at **t = 110.14**, even
  more gently — asymptotic speed **0.16**, barely unbound.
- The two that stayed, **Star 0 and Star 2 — one from each original
  pair** — ended bound to each other: an exchange, not a survival of
  either starting pair. Their binary: semi-major axis **0.256**, period
  **0.577**.

So the four did not end as two pairs, and not as a triple-plus-one. It
ended as **one pair and two singles**, and the pair that survived was
made of strangers.

Energy drift: final **0.052%**; worst moment **0.121%**, a brief excursion
during the single closest encounter — stated plainly rather than rounded
away.

## Honesty notes

- **Run C is chaotic, and its details are sensitive.** The same initial
  conditions integrated at η = 0.005 instead of 0.002 gave the same *kind*
  of ending — two ejections and one exchanged bound pair — but with a
  different individual star leaving first (star 2 at t = 78, then star 1
  at t = 129; final pair 0–3). Which star goes first in a chaotic system
  is not a fact about the stars; it is a fact about the exact numbers.
  The robust result is the shape of the resolution, not the name of the
  first leaver. Run A is chaotic too, though its outcome (star 0 out at
  t ≈ 7.4, hard pair left behind) was identical at two different fixed
  timesteps.
- **Fixed-step attempts at Run C failed honestly before this run:** a
  collapsing four-body system throws encounters too close for a fixed
  step, and the integrator answered with energy errors of thousands of
  percent — broken runs, discarded, which is why Run C uses the adaptive
  step. The drift figures above are the ones from the runs reported.
- In the video's terms — *which star gets ejected first?* — our three
  answered: the right-hand one, almost immediately. Our four answered:
  nobody for a long while, and then the system traded partners before it
  let anyone go.
