#!/usr/bin/env python3
"""
Newtonian N-body runs for Tristan's three-star / four-star test.
Units: G = 1, every star mass = 1, distances ~ 1.
Integrator: leapfrog (kick-drift-kick), fixed timestep, no softening.
Energy drift is printed for every run as the honesty check.
"""
import json
import numpy as np

G = 1.0


def accel(pos, mass):
    """Newtonian acceleration, no softening. pos: (N,2)."""
    n = len(mass)
    a = np.zeros_like(pos)
    for i in range(n):
        d = pos - pos[i]                 # vectors from i to j
        r2 = np.sum(d * d, axis=1)
        r2[i] = np.inf                   # skip self
        inv_r3 = r2 ** -1.5
        inv_r3[i] = 0.0
        a[i] = G * np.sum((d * (mass * inv_r3)[:, None]), axis=0)
    return a


def energy(pos, vel, mass):
    ke = 0.5 * np.sum(mass[:, None] * vel * vel)
    pe = 0.0
    n = len(mass)
    for i in range(n):
        for j in range(i + 1, n):
            r = np.linalg.norm(pos[i] - pos[j])
            pe -= G * mass[i] * mass[j] / r
    return ke + pe


def integrate(pos, vel, mass, dt, t_end, check=None, check_every=100):
    """Leapfrog KDK. check(t, pos, vel) may return a dict to stop early."""
    pos = pos.copy(); vel = vel.copy()
    a = accel(pos, mass)
    E0 = energy(pos, vel, mass)
    max_drift = 0.0
    t = 0.0
    step = 0
    event = None
    nsteps = int(round(t_end / dt))
    while step < nsteps:
        vel += 0.5 * dt * a
        pos += dt * vel
        a = accel(pos, mass)
        vel += 0.5 * dt * a
        t += dt; step += 1
        if step % check_every == 0:
            E = energy(pos, vel, mass)
            drift = abs(E - E0) / abs(E0)
            max_drift = max(max_drift, drift)
            if check is not None:
                event = check(t, pos, vel)
                if event is not None:
                    break
    E1 = energy(pos, vel, mass)
    return dict(t=t, pos=pos, vel=vel, E0=E0, E1=E1,
                drift_final=abs(E1 - E0) / abs(E0), drift_max=max_drift,
                event=event)


def to_com(pos, vel, mass):
    """Shift to centre-of-mass frame (zero total momentum)."""
    M = mass.sum()
    com = (mass[:, None] * pos).sum(axis=0) / M
    vcom = (mass[:, None] * vel).sum(axis=0) / M
    return pos - com, vel - vcom


# ---------------------------------------------------------------- RUN A
def run_A():
    print("=== RUN A: three stars, generic chaotic triple ===")
    mass = np.ones(3)
    pos = np.array([[1.00, 0.00],
                    [-0.55, 0.85],
                    [-0.45, -0.88]])
    vel = np.array([[0.05, 0.30],
                    [-0.28, -0.12],
                    [0.22, -0.20]])
    pos, vel = to_com(pos, vel, mass)
    print("ICs (COM frame):")
    for i in range(3):
        print(f"  star {i}: pos=({pos[i,0]:+.6f},{pos[i,1]:+.6f}) "
              f"vel=({vel[i,0]:+.6f},{vel[i,1]:+.6f})")
    dt = 0.0002

    def eject_check(t, p, v):
        for i in range(3):
            others = [j for j in range(3) if j != i]
            bcen = p[others].mean(axis=0)
            bvel = v[others].mean(axis=0)
            rel_p = p[i] - bcen
            rel_v = v[i] - bvel
            dist = np.linalg.norm(rel_p)
            if dist > 20.0:
                espec = 0.5 * np.dot(rel_v, rel_v) - G * 2.0 / dist
                radial = np.dot(rel_p, rel_v) / dist
                if espec > 0 and radial > 0:
                    return dict(star=i, t=t, dist=dist,
                                speed=float(np.linalg.norm(rel_v)),
                                v_inf=float(np.sqrt(2 * espec)),
                                others=others)
        return None

    res = integrate(pos, vel, mass, dt, 2000.0, check=eject_check)
    out = dict(run="A", dt=dt, E0=res["E0"], E1=res["E1"],
               drift_final=res["drift_final"], drift_max=res["drift_max"],
               t_end=res["t"])
    ev = res["event"]
    if ev:
        i = ev["star"]; j, k = ev["others"]
        # remaining pair orbital elements from relative state at ejection
        rp = res["pos"][j] - res["pos"][k]
        rv = res["vel"][j] - res["vel"][k]
        r = np.linalg.norm(rp)
        mu = G * 2.0
        e_rel = 0.5 * np.dot(rv, rv) - mu / r      # per unit reduced mass
        a_sem = -mu / (2 * e_rel)
        period = 2 * np.pi * np.sqrt(a_sem ** 3 / mu)
        out.update(ejected=True, ejected_star=i, t_eject=ev["t"],
                   pair=[j, k], pair_sep_now=float(r),
                   pair_semi_major=float(a_sem), pair_period=float(period),
                   eject_speed=float(ev["speed"]), eject_v_inf=ev["v_inf"])
        print(f"EJECTION: star {i} ejected at t = {ev['t']:.2f}")
        print(f"  remaining pair = stars {j},{k}; separation now {r:.3f}, "
              f"semi-major axis {a_sem:.3f}, period {period:.3f}")
        print(f"  ejected star speed vs pair barycentre {ev['speed']:.4f}, "
              f"asymptotic v_inf {ev['v_inf']:.4f}")
    else:
        print(f"No ejection by t = {res['t']:.1f}")
        out.update(ejected=False)
    print(f"Energy: E0={res['E0']:.8f} E1={res['E1']:.8f} "
          f"drift_final={res['drift_final']*100:.4f}% "
          f"drift_max={res['drift_max']*100:.4f}%")
    return out


# ---------------------------------------------------------------- RUN B
def run_B():
    print("=== RUN B: figure-eight (Chenciner-Montgomery), control ===")
    mass = np.ones(3)
    pos = np.array([[-0.97000436, 0.24308753],
                    [0.97000436, -0.24308753],
                    [0.0, 0.0]])
    v1 = np.array([0.4662036850, 0.4323657300])
    vel = np.array([v1, v1, -2 * v1])
    dt = 0.0005
    T = 6.32591398                      # published period
    # First period: sample the shared track.
    r1 = integrate(pos, vel, mass, dt, T)
    track_run = integrate(pos, vel, mass, dt, T)  # placeholder, replaced below
    # Re-run sampling the track densely:
    p = pos.copy(); v = vel.copy(); a = accel(p, mass)
    E0 = energy(p, v, mass)
    track = []
    steps_per = int(round(T / dt))
    for s in range(steps_per):
        v += 0.5 * dt * a; p += dt * v; a = accel(p, mass); v += 0.5 * dt * a
        if s % 5 == 0:
            track.append(p.copy())
    track = np.array(track)             # (samples, 3, 2)
    track_pts = track.reshape(-1, 2)    # union of all stars' positions = the 8
    # Continue for two more periods, measuring distance of each star to track.
    max_dev = 0.0
    max_drift = 0.0
    for s in range(2 * steps_per):
        v += 0.5 * dt * a; p += dt * v; a = accel(p, mass); v += 0.5 * dt * a
        if s % 20 == 0:
            d = np.linalg.norm(p[:, None, :] - track_pts[None, :, :], axis=2)
            max_dev = max(max_dev, float(d.min(axis=1).max()))
            E = energy(p, v, mass)
            max_drift = max(max_drift, abs(E - E0) / abs(E0))
    E1 = energy(p, v, mass)
    print(f"Integrated 3 full periods (t = {3*T:.4f})")
    print(f"Max deviation of any star from the shared figure-eight track: "
          f"{max_dev:.6f}")
    print(f"Energy: E0={E0:.8f} E1={E1:.8f} "
          f"drift_final={abs(E1-E0)/abs(E0)*100:.5f}% "
          f"drift_max={max_drift*100:.5f}%")
    return dict(run="B", dt=dt, period=T, periods=3, max_dev=max_dev,
                E0=E0, E1=E1, drift_final=abs(E1 - E0) / abs(E0),
                drift_max=max_drift, held=bool(max_dev < 0.05))


def integrate_adaptive(pos, vel, mass, eta, t_end, check=None,
                       check_dt=0.05, dt_max=0.01, max_steps=60_000_000):
    """Leapfrog KDK with a timestep set by the closest pair's dynamical
    time: dt = eta * min_ij sqrt(r_ij^3 / (G (m_i + m_j))), capped at
    dt_max. Used for Run C, where a collapsing four-body system produces
    encounters too close for any fixed step to resolve honestly."""
    pos = pos.copy(); vel = vel.copy()
    a = accel(pos, mass)
    E0 = energy(pos, vel, mass)
    max_drift = 0.0
    t = 0.0; step = 0; next_check = 0.0
    event = None
    n = len(mass)
    while t < t_end and step < max_steps:
        dt = dt_max
        for i in range(n):
            for j in range(i + 1, n):
                r = np.linalg.norm(pos[i] - pos[j])
                dt = min(dt, eta * np.sqrt(r ** 3 / (G * (mass[i] + mass[j]))))
        vel += 0.5 * dt * a
        pos += dt * vel
        a = accel(pos, mass)
        vel += 0.5 * dt * a
        t += dt; step += 1
        if t >= next_check:
            next_check = t + check_dt
            E = energy(pos, vel, mass)
            max_drift = max(max_drift, abs(E - E0) / abs(E0))
            if check is not None:
                event = check(t, pos, vel)
                if event is not None:
                    break
    E1 = energy(pos, vel, mass)
    return dict(t=t, pos=pos, vel=vel, E0=E0, E1=E1,
                drift_final=abs(E1 - E0) / abs(E0), drift_max=max_drift,
                event=event, steps=step)


# ---------------------------------------------------------------- RUN C
def run_C():
    print("=== RUN C: four stars, generic (two loose pairs) ===")
    mass = np.ones(4)
    # Two loose pairs on a slow mutual approach, built by hand:
    # pair A = stars 0,1 (centre (-1.8, 0), separation 0.8, circular
    # relative speed sqrt(2/0.8)=1.581, split evenly); pair B = stars 2,3
    # (centre (1.8, 0.3), separation 0.9, opposite internal sense).
    # Pair centres get a slightly sub-circular mutual tangential drift.
    pos = np.array([[-2.20, 0.00],
                    [-1.40, 0.00],
                    [1.35, 0.30],
                    [2.25, 0.30]])
    vel = np.array([[0.00, -0.7906 - 0.45],
                    [0.00, 0.7906 - 0.45],
                    [0.00, 0.7454 + 0.45],
                    [0.00, -0.7454 + 0.45]])
    pos, vel = to_com(pos, vel, mass)
    print("ICs (COM frame):")
    for i in range(4):
        print(f"  star {i}: pos=({pos[i,0]:+.6f},{pos[i,1]:+.6f}) "
              f"vel=({vel[i,0]:+.6f},{vel[i,1]:+.6f})")
    eta = 0.002
    ejected = {}

    def check(t, p, v):
        # A star counts as ejected when it is >20 units from the barycentre
        # of the not-yet-ejected rest, unbound from them, and receding.
        for i in range(4):
            if i in ejected:
                continue
            others = [j for j in range(4) if j != i and j not in ejected]
            bcen = p[others].mean(axis=0); bvel = v[others].mean(axis=0)
            rel_p = p[i] - bcen; rel_v = v[i] - bvel
            dist = np.linalg.norm(rel_p)
            if dist > 20.0:
                espec = 0.5 * np.dot(rel_v, rel_v) - G * len(others) / dist
                if espec > 0 and np.dot(rel_p, rel_v) > 0:
                    ejected[i] = dict(t=t, speed=float(np.linalg.norm(rel_v)),
                                      v_inf=float(np.sqrt(2 * espec)))
        if len(ejected) >= 2:
            return dict(done=True)
        return None

    res = integrate_adaptive(pos, vel, mass, eta, 2000.0, check=check)
    print(f"  (adaptive leapfrog: eta={eta}, steps taken={res['steps']})")
    p, v = res["pos"], res["vel"]
    print(f"Run ended at t = {res['t']:.2f}")
    for i, ev in ejected.items():
        print(f"  star {i} ejected at t = {ev['t']:.2f}, "
              f"speed {ev['speed']:.4f}, v_inf {ev['v_inf']:.4f}")
    remaining = [i for i in range(4) if i not in ejected]
    print(f"  remaining stars: {remaining}")
    # pairwise bound analysis among remaining
    pairs = []
    for a_ in range(len(remaining)):
        for b_ in range(a_ + 1, len(remaining)):
            i, j = remaining[a_], remaining[b_]
            rp = p[i] - p[j]; rv = v[i] - v[j]
            r = np.linalg.norm(rp)
            espec = 0.5 * np.dot(rv, rv) - G * 2.0 / r
            pairs.append((i, j, float(r), float(espec)))
            print(f"  pair {i}-{j}: separation {r:.3f}, "
                  f"specific rel. energy {espec:+.4f} "
                  f"({'BOUND' if espec < 0 else 'unbound'})")
    print(f"Energy: E0={res['E0']:.8f} E1={res['E1']:.8f} "
          f"drift_final={res['drift_final']*100:.4f}% "
          f"drift_max={res['drift_max']*100:.4f}%")
    return dict(run="C", eta=eta, steps=res["steps"],
                E0=res["E0"], E1=res["E1"],
                drift_final=res["drift_final"], drift_max=res["drift_max"],
                t_end=res["t"],
                ejected={int(k): vv for k, vv in ejected.items()},
                remaining=remaining,
                pairs=[list(x) for x in pairs],
                final_pos=p.tolist(), final_vel=v.tolist())


if __name__ == "__main__":
    results = {}
    results["A"] = run_A()
    print()
    results["B"] = run_B()
    print()
    results["C"] = run_C()
    with open("results_raw.json", "w") as f:
        json.dump(results, f, indent=2, default=float)
    print("\nRaw results written to results_raw.json")
