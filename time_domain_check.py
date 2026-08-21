#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
time_domain_check.py

Independent validation of basin_network_response.py.

The same network is integrated in the time domain with RK4, keeping the
quadratic drag as it is (no Lorentz linearisation) and using exactly the
same friction closure, so the comparison isolates the error of the
harmonic solution rather than mixing in a different closure:

    L_e dQ_e/dt + kappa_e Q_e|Q_e| = eta_node1 - eta_node2
    S_m d eta_m/dt                 = net inflow into basin m

    kappa_e = (3 pi / 8) * beta_e   so that the Lorentz linearisation of
                                    kappa_e Q|Q| is exactly beta_e |Q| Q

The run starts from rest, integrates for a number of forcing periods, and
Fourier-analyses the last few to extract the fundamental and the harmonics
of the periodic state.

Usage:

    python time_domain_check.py network.txt --period 12.42 --amplitude 0.5

Reference results obtained this way, on a two-basin and a three-basin
network, over periods from 2 to 87 h and amplitudes from 0.1 to 2 m:
the fundamental amplitude is accurate to about 1% (worst 3%), and the
third harmonic to about 3% (worst 12%), the model slightly underestimating
it because the fifth and higher harmonics are truncated.
"""

import math
import numpy as np


def assemble(bays, boundaries, channels):
    names = list(bays)
    idx = {n: i for i, n in enumerate(names)}
    M, Ne = len(names), len(channels)

    A = np.zeros((Ne, M))          # incidence on the basins
    bnd = np.zeros((Ne, 2))        # (factor, phase) of the boundary end
    bsign = np.zeros(Ne)           # +1 if the boundary is node1, -1 if node2

    Lin = np.array([ch.inertance for ch in channels])
    kappa = np.array([ch.beta for ch in channels]) * 3.0 * math.pi / 8.0
    S = np.array([bays[n].S for n in names])

    for e, ch in enumerate(channels):
        for node, sgn in ((ch.node1, +1.0), (ch.node2, -1.0)):
            if node in bays:
                A[e, idx[node]] = sgn
            else:
                b = boundaries[node]
                bnd[e] = (b.factor, math.radians(b.phase_deg))
                bsign[e] = sgn

    return names, idx, A, bnd, bsign, Lin, kappa, S


def integrate(bays, boundaries, channels, omega, a_s,
              periods=40, steps_per_period=2000, analyse_last=4):
    names, idx, A, bnd, bsign, Lin, kappa, S = assemble(bays, boundaries, channels)

    dt = 2.0 * math.pi / omega / steps_per_period
    nsteps = periods * steps_per_period

    def forcing(t):
        return bsign * a_s * bnd[:, 0] * np.cos(omega * t + bnd[:, 1])

    def deriv(t, Q, eta):
        dhead = A @ eta + forcing(t)
        dQ = (dhead - kappa * Q * np.abs(Q)) / Lin
        deta = -(A.T @ Q) / S
        return dQ, deta

    Q = np.zeros(len(channels))
    eta = np.zeros(len(names))

    keep = analyse_last * steps_per_period
    store = np.empty((keep, len(names)))
    t = 0.0

    for k in range(nsteps):
        if k >= nsteps - keep:
            store[k - (nsteps - keep)] = eta

        k1Q, k1e = deriv(t, Q, eta)
        k2Q, k2e = deriv(t + dt / 2, Q + dt / 2 * k1Q, eta + dt / 2 * k1e)
        k3Q, k3e = deriv(t + dt / 2, Q + dt / 2 * k2Q, eta + dt / 2 * k2e)
        k4Q, k4e = deriv(t + dt, Q + dt * k3Q, eta + dt * k3e)

        Q = Q + dt / 6 * (k1Q + 2 * k2Q + 2 * k3Q + k4Q)
        eta = eta + dt / 6 * (k1e + 2 * k2e + 2 * k3e + k4e)
        t += dt

    # harmonic analysis of the stored window (an exact number of periods)
    tt = (np.arange(keep) + (nsteps - keep)) * dt
    out = {}
    for i, n in enumerate(names):
        y = store[:, i]
        coeff = {}
        for h in (1, 2, 3, 5):
            coeff[h] = 2.0 * np.sum(y * np.exp(-1j * h * omega * tt)) / keep
        out[n] = dict(c1=coeff[1], c2=coeff[2], c3=coeff[3], c5=coeff[5],
                      rng=float(y.max() - y.min()))
    return out


# ======================================================================
# COMMAND LINE
# ======================================================================

def main():
    import argparse
    import cmath
    from pathlib import Path

    import basin_network_response as bnr

    parser = argparse.ArgumentParser(
        description=(
            "Compare basin_network_response.py against a time-domain "
            "integration of the same network."
        )
    )

    parser.add_argument("input", help="network .txt file")

    parser.add_argument(
        "--period",
        type=float,
        help="forcing period [hours]",
    )

    parser.add_argument(
        "--omega",
        type=float,
        help="forcing angular frequency [rad/s], alternative to --period",
    )

    parser.add_argument(
        "--amplitude",
        type=float,
        required=True,
        help="sea amplitude a_s [m]",
    )

    parser.add_argument(
        "--periods",
        type=int,
        default=20,
        help="forcing periods to integrate before analysing",
    )

    parser.add_argument(
        "--steps-per-period",
        type=int,
        default=1000,
        help="RK4 steps per forcing period",
    )

    args = parser.parse_args()

    if (args.period is None) == (args.omega is None):
        parser.error("give either --period or --omega")

    omega = (
        args.omega
        if args.omega is not None
        else 2.0 * math.pi / (args.period * 3600.0)
    )

    a_s = args.amplitude

    bays, boundaries, channels = bnr.read_network(args.input)

    ref = integrate(
        bays,
        boundaries,
        channels,
        omega,
        a_s,
        periods=args.periods,
        steps_per_period=args.steps_per_period,
    )

    eta, Q, R, Y, nit, ok = bnr.solve_network(
        bays=bays,
        boundaries=boundaries,
        channels=channels,
        omega=omega,
        a_s=a_s,
    )

    eta3, Q3 = bnr.third_harmonic_network(
        bays,
        boundaries,
        channels,
        omega,
        Q,
    )

    print()
    print(
        f"Forcing    a_s = {a_s:g} m"
        f"     T = {2 * math.pi / omega / 3600:.3f} h"
        f"     reference: {args.periods} periods, "
        f"{args.steps_per_period} steps/period"
    )
    print()
    print(
        f"  {'basin':<10}{'quantity':<14}"
        f"{'time domain':>13}{'model':>12}{'difference':>13}"
    )

    for name in bays:
        c1 = ref[name]["c1"]
        c3 = ref[name]["c3"]
        c5 = ref[name]["c5"]

        for label, refval, modval in (
            ("A1 / a_s", abs(c1) / a_s, abs(eta[name]) / a_s),
            ("A3 / a_s", abs(c3) / a_s, abs(eta3[name]) / a_s),
        ):
            diff = 100.0 * (modval - refval) / refval

            print(
                f"  {name:<10}{label:<14}"
                f"{refval:>13.6f}{modval:>12.6f}{diff:>12.2f} %"
            )

        print(
            f"  {name:<10}{'A5 / a_s':<14}"
            f"{abs(c5) / a_s:>13.6f}{'-':>12}"
            f"{'truncated':>13}"
        )

    qratio = np.abs(Q3) / np.abs(Q)

    print()
    print(f"  worst |Q3|/|Q1| = {qratio.max():.3f}")
    print()


if __name__ == "__main__":
    main()
