#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
multi_inlet_response.py

Compute the gain G of a single semi-enclosed basin connected to the sea by
N parallel inlets, using the fully self-consistent Lorentz linearisation,
sweeping:

    omega : from 1e-5 to 1e-3 rad/s
    a_s   : from 0.05 to 5 m

The input file must follow this format:

    N = 2
    L = 3200 4500
    D = 1.2 3.4
    B = 150 255.4
    kse = 30 30
    S = 10000000

Assumptions:
    A_j  = B_j D_j
    Rh_j = B_j D_j / (B_j + 2 D_j)
    n_j  = 1 / kse_j
    fc   = 0

Any a_s value found in the input file is ignored, because here the
forcing amplitude is swept.
"""

from __future__ import annotations

import argparse
import cmath
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


G_GRAV = 9.81

# Quality thresholds for the third harmonic, on the per-inlet |Q3|/|Q1|.
Q3_CAUTION = 0.20
Q3_WARN = 0.50

# ------------------------------------------------------------
# DEFAULT RANGES
# ------------------------------------------------------------
OMEGA_MIN = 1.0e-5       # rad/s
OMEGA_MAX = 1.0e-3       # rad/s
N_OMEGA = 700

# Forcing amplitudes to be compared [m]
AMPLITUDES = np.array([0.05, 0.10, 0.20, 0.50, 1.0, 2.0, 5.0])

TOL = 1.0e-12
ITMAX = 20000
RELAX = 0.5


# ============================================================
# INPUT
# ============================================================

def read_input(path):
    data = {}

    with open(path, "r", encoding="utf-8") as f:
        for lineno, raw in enumerate(f, 1):
            line = raw.split("#", 1)[0].strip()

            if not line:
                continue

            if "=" not in line:
                raise ValueError(f"Line {lineno}: missing '='")

            key, value = line.split("=", 1)
            key = key.strip().lower()
            vals = value.split()

            if key == "n":
                data[key] = int(vals[0])

            elif key == "s":
                data[key] = float(vals[0])

            elif key in ("l", "d", "b", "kse"):
                data[key] = np.array(
                    [float(x) for x in vals],
                    dtype=float
                )

            elif key == "a_s":
                # Optional: used as the default amplitude of a
                # single-point run, ignored by the sweep.
                data[key] = float(vals[0])

            else:
                raise ValueError(
                    f"Line {lineno}: unknown parameter '{key}'"
                )

    required = ("n", "l", "d", "b", "kse", "s")

    missing = [key for key in required if key not in data]
    if missing:
        raise ValueError(
            "Missing parameters: " + ", ".join(missing)
        )

    N = data["n"]

    for key in ("l", "d", "b", "kse"):
        if len(data[key]) != N:
            raise ValueError(
                f"{key.upper()}: {N} values expected, "
                f"{len(data[key])} found"
            )

        if np.any(data[key] <= 0):
            raise ValueError(
                f"{key.upper()}: all values must be positive"
            )

    if data["s"] <= 0:
        raise ValueError("S must be positive")

    return data


# ============================================================
# GEOMETRIA
# ============================================================

def prepare_geometry(data, a_s):
    """
    Rectangular cross-section:

        A  = B D
        Rh = A / (B + 2D)

    Manning-Strickler:

        n = 1/kse
    """

    Lch = data["l"]
    D = data["d"]
    B = data["b"]
    kse = data["kse"]
    S = data["s"]

    A = B * D

    wetted_perimeter = B + 2.0 * D
    Rh = A / wetted_perimeter

    n = 1.0 / kse

    # Inertance
    inertance = Lch / (G_GRAV * A)

    # Reduced geometric resistance as defined in the paper
    # fc = 0 -> nr = n
    rgeo = (
        8.0
        * Lch
        * n**2
        * S
        * a_s
        /
        (
            3.0
            * math.pi
            * A**2
            * Rh**(4.0 / 3.0)
        )
    )

    return A, Rh, inertance, rgeo


# ============================================================
# CLOSED-FORM INITIAL GUESS
# ============================================================

def closed_form_initial_guess(S, omega, inertance, rgeo):

    weights = rgeo**(-0.5)

    s = np.sum(weights)

    phi = weights / s

    req = 1.0 / s**2

    leq = (
        np.sum(inertance / rgeo)
        / s**2
    )

    Ddet = 1.0 - omega**2 * S * leq

    theta = (
        omega**4
        * S**2
        * req**2
    )

    G = math.sqrt(
        2.0
        /
        (
            Ddet**2
            + math.sqrt(
                Ddet**4 + 4.0 * theta
            )
        )
    )

    return G, phi


# ============================================================
# PARALLEL REDUCTION
# ============================================================

def equivalent_impedance(
    G,
    phi,
    inertance,
    rgeo,
    omega
):

    r = phi * rgeo

    Z = omega * (
        G * r
        + 1j * inertance
    )

    Y = 1.0 / Z

    Ysum = np.sum(Y)

    Zeq = 1.0 / Ysum

    Req = Zeq.real
    Leq = Zeq.imag / omega

    return Req, Leq, Y, Ysum


def gain_from_equivalent(
    Req,
    Leq,
    omega,
    S
):

    denominator = complex(
        1.0 - omega**2 * S * Leq,
        omega * S * Req
    )

    return 1.0 / abs(denominator)


# ============================================================
# FIRST SECONDARY HARMONIC
# ============================================================

def complex_response(
    S,
    omega,
    a_s,
    inertance,
    rgeo,
    G,
    phi
):
    """
    Complex basin elevation and channel discharges of the fundamental,
    with the sea level taken as the real reference a_s + 0j:

        eta1 = Ysum a_s / (Ysum + i omega S)
        Q1_j = Y_j a_s i omega S / (Ysum + i omega S)

    so |eta1|/a_s reproduces the gain G and arg(eta1) is the phase lag.
    """

    _, _, Y, Ysum = equivalent_impedance(
        G,
        phi,
        inertance,
        rgeo,
        omega
    )

    denom = Ysum + 1j * omega * S

    eta1 = Ysum * a_s / denom
    Q1 = Y * a_s * 1j * omega * S / denom

    return eta1, Q1


def third_harmonic(
    S,
    omega,
    a_s,
    inertance,
    rgeo,
    G,
    phi
):
    """
    Perturbative 3*omega response, given the converged fundamental.

    Expanding the quadratic drag for Q = |Q| cos(theta) gives odd
    harmonics only. The n = 1 term is the Lorentz resistance already
    carried by the fundamental; the n = 3 term is a known source once
    the fundamental has converged, so the 3*omega problem is linear:

        R_j^(3) = (3/2) beta_j |Q_j|
        S_j     = (1/5) beta_j Q_j^3 / |Q_j|
        beta_j  = rgeo_j / (S a_s)

    with eta3 = 0 imposed at the sea. Continuity then gives directly

        eta3 = - sum_j Y3_j S_j / (3 i omega S + sum_j Y3_j)

    Returns the complex eta3 and the complex fundamental eta1.
    """

    eta1, Q1 = complex_response(
        S,
        omega,
        a_s,
        inertance,
        rgeo,
        G,
        phi
    )

    beta = rgeo / (S * a_s)

    absQ = np.abs(Q1)
    safe = np.where(absQ > 0.0, absQ, 1.0)

    R3 = 1.5 * beta * absQ
    source = 0.2 * beta * Q1**3 / safe

    Y3 = 1.0 / (
        R3
        + 3j * omega * inertance
    )

    eta3 = (
        -np.sum(Y3 * source)
        / (3j * omega * S + np.sum(Y3))
    )

    Q3 = Y3 * (-eta3 - source)

    return eta3, eta1, Q3


def third_harmonic_quality(ratio_max):
    """Return (level, lines) describing how far the perturbation is trusted.

    The test is the per-inlet ratio |Q3|/|Q1|. For a single basin it stays
    around 0.1 - 0.15 (exactly 2/15 in the friction-dominated limit of a
    single inlet), and validation against a time-domain integration of the
    same equations puts the error on the 3*omega amplitude within about 10%
    there. Larger values mean 3*omega is close to a resonance.
    The fundamental is not affected in either case.
    """

    if ratio_max > Q3_WARN:
        return "warn", [
            "3*omega is at or near a resonance, so the perturbation has "
            "broken down:",
            "treat the 3w amplitude as indicative only. "
            "The fundamental is not affected.",
        ]

    if ratio_max > Q3_CAUTION:
        return "caution", [
            "3*omega is approaching a resonance: the 3w amplitude may be "
            "in error by",
            "more than 10%. The fundamental is not affected.",
        ]

    return "ok", []


# ============================================================
# SELF-CONSISTENT SOLUTION
# ============================================================

def self_consistent_gain(
    S,
    omega,
    inertance,
    rgeo,
    tol=TOL,
    itmax=ITMAX,
    relax=RELAX
):

    G, phi = closed_form_initial_guess(
        S,
        omega,
        inertance,
        rgeo
    )

    for iteration in range(1, itmax + 1):

        Req, Leq, Y, Ysum = (
            equivalent_impedance(
                G,
                phi,
                inertance,
                rgeo,
                omega
            )
        )

        G_new = gain_from_equivalent(
            Req,
            Leq,
            omega,
            S
        )

        phi_new = np.abs(
            Y / Ysum
        )

        err_G = (
            abs(G_new - G)
            / max(abs(G), 1.0e-30)
        )

        err_phi = np.max(
            np.abs(phi_new - phi)
        )

        if (
            err_G <= tol
            and err_phi <= tol
        ):

            return (
                G_new,
                phi_new,
                iteration,
                True
            )

        G = (
            (1.0 - relax) * G
            + relax * G_new
        )

        phi = (
            (1.0 - relax) * phi
            + relax * phi_new
        )

    return (
        G,
        phi,
        itmax,
        False
    )


# ============================================================
# MAIN
# ============================================================

def build_parser():

    parser = argparse.ArgumentParser(
        description=(
            "Tidal gain of a single basin connected to the sea "
            "by N inlets in parallel."
        ),
        epilog=(
            "Examples:\n"
            "  one point   : %(prog)s inlets.txt --period 12.42 --amplitude 0.8\n"
            "  full sweep  : %(prog)s inlets.txt\n"
            "  custom sweep: %(prog)s inlets.txt --amplitudes 0.2 1 --nfreq 800"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "input",
        help="input .txt file"
    )

    point = parser.add_argument_group(
        "single point (give --omega or --period to use it)"
    )

    point.add_argument(
        "--omega",
        type=float,
        help="forcing angular frequency [rad/s]"
    )

    point.add_argument(
        "--period",
        type=float,
        help="forcing period [hours], alternative to --omega"
    )

    point.add_argument(
        "--amplitude",
        type=float,
        help=(
            "sea amplitude a_s [m]; "
            "defaults to the a_s of the input file"
        )
    )

    point.add_argument(
        "--inlets",
        action="store_true",
        help="also list the discharge carried by each inlet"
    )

    sweep = parser.add_argument_group("sweep")

    sweep.add_argument(
        "--amplitudes",
        type=float,
        nargs="+",
        default=list(AMPLITUDES),
        help="sea amplitudes to compare [m]"
    )

    sweep.add_argument(
        "--omega-min",
        type=float,
        default=OMEGA_MIN,
        help="minimum angular frequency [rad/s]"
    )

    sweep.add_argument(
        "--omega-max",
        type=float,
        default=OMEGA_MAX,
        help="maximum angular frequency [rad/s]"
    )

    sweep.add_argument(
        "--nfreq",
        type=int,
        default=N_OMEGA,
        help="number of log-spaced frequencies"
    )

    sweep.add_argument(
        "--show",
        action="store_true",
        help="display the plot"
    )

    return parser


def solve_point(data, omega, a_s):
    """Fundamental response at one frequency and one amplitude."""

    S = data["s"]

    A, Rh, inertance, rgeo = prepare_geometry(data, a_s)

    G, phi, nit, ok = self_consistent_gain(
        S,
        omega,
        inertance,
        rgeo
    )

    eta1, Q1 = complex_response(
        S,
        omega,
        a_s,
        inertance,
        rgeo,
        G,
        phi
    )

    return G, phi, eta1, Q1, inertance, rgeo, ok


def report_point(data, omega, a_s, show_inlets):

    G, phi, eta1, Q1, inertance, rgeo, ok = solve_point(
        data,
        omega,
        a_s
    )

    S = data["s"]

    eta3, _, Q3 = third_harmonic(
        S,
        omega,
        a_s,
        inertance,
        rgeo,
        G,
        phi
    )

    qratio = np.abs(Q3) / np.abs(Q1)

    period_h = 2.0 * math.pi / omega / 3600.0
    lag_deg = -math.degrees(cmath.phase(eta1))
    lag_h = lag_deg / 360.0 * period_h

    print()
    print(
        f"Forcing    a_s = {a_s:g} m"
        f"     omega = {omega:.4e} rad/s"
        f"     T = {period_h:.3f} h"
    )
    print()
    print(f"  gain             G = {G:.4f}")
    print(f"  basin amplitude    = {G * a_s:.4f} m")
    print(f"  basin tidal range  = {2.0 * G * a_s:.4f} m")
    print(f"  phase lag          = {lag_deg:.2f} deg  ({lag_h:.3f} h)")
    print(f"  3w amplitude       = {abs(eta3):.4f} m")

    if not ok:
        print("  [WARNING] the fixed point did not converge")

    if show_inlets:
        print()
        print("  inlet    |Q| [m3/s]     share   |Q3|/|Q1|")

        absQ = np.abs(Q1)
        total = absQ.sum()

        for j in range(data["n"]):
            print(
                f"  {j + 1:>5}  {absQ[j]:12.1f}"
                f"  {100.0 * absQ[j] / total:8.1f} %"
                f"  {qratio[j]:10.3f}"
            )

    level, lines = third_harmonic_quality(qratio.max())

    if level != "ok":
        tag = "WARNING" if level == "warn" else "NOTE"

        print()
        print(
            f"  [{tag}] worst |Q3|/|Q1| = {qratio.max():.2f} "
            f"(inlet {int(np.argmax(qratio)) + 1})"
        )

        for line in lines:
            print(f"           {line}")

    print()


def main():

    parser = build_parser()
    args = parser.parse_args()

    data = read_input(args.input)

    N = data["n"]
    S = data["s"]

    # ------------------------------------------------------------
    # Single point
    # ------------------------------------------------------------

    if args.omega is not None or args.period is not None:

        if args.omega is not None and args.period is not None:
            parser.error("give either --omega or --period, not both")

        omega = (
            args.omega
            if args.omega is not None
            else 2.0 * math.pi / (args.period * 3600.0)
        )

        if omega <= 0:
            parser.error("the forcing frequency must be positive")

        a_s = args.amplitude

        if a_s is None:
            a_s = data.get("a_s")

        if a_s is None:
            parser.error(
                "no amplitude: pass --amplitude, "
                "or set a_s in the input file"
            )

        if a_s <= 0:
            parser.error("the amplitude must be positive")

        report_point(data, omega, a_s, args.inlets)
        return

    # ------------------------------------------------------------
    # Sweep
    # ------------------------------------------------------------

    amplitudes = np.array(args.amplitudes, dtype=float)

    if np.any(amplitudes <= 0):
        parser.error("all amplitudes must be positive")

    if not (0 < args.omega_min < args.omega_max):
        parser.error("0 < omega-min < omega-max is required")

    if args.nfreq < 2:
        parser.error("nfreq must be >= 2")

    omega_values = np.logspace(
        np.log10(args.omega_min),
        np.log10(args.omega_max),
        args.nfreq
    )

    results = {}
    results3 = {}
    worst_q = 0.0
    worst_where = None

    print(
        f"{N} inlets   S = {S:.4g} m2   "
        f"omega {args.omega_min:.1e} - {args.omega_max:.1e} rad/s   "
        f"({args.nfreq} points)"
    )
    print()

    for a_s in amplitudes:

        A, Rh, inertance, rgeo = prepare_geometry(data, a_s)

        gain_values = np.empty_like(omega_values)
        third_values = np.zeros_like(omega_values)

        all_converged = True

        for i, omega in enumerate(omega_values):

            G, phi, nit, ok = self_consistent_gain(
                S,
                omega,
                inertance,
                rgeo
            )

            gain_values[i] = G

            eta3, eta1, Q3 = third_harmonic(
                S,
                omega,
                a_s,
                inertance,
                rgeo,
                G,
                phi
            )

            third_values[i] = abs(eta3) / abs(eta1)

            _, Q1 = complex_response(
                S,
                omega,
                a_s,
                inertance,
                rgeo,
                G,
                phi
            )

            qratio = np.abs(Q3) / np.abs(Q1)

            if qratio.max() > worst_q:
                worst_q = float(qratio.max())
                worst_where = (
                    int(np.argmax(qratio)) + 1,
                    omega,
                    a_s,
                )

            if not ok:
                all_converged = False

        results[a_s] = gain_values
        results3[a_s] = third_values

        imax = int(np.argmax(gain_values))

        line = (
            f"  a_s = {a_s:6g} m :  Gmax = {gain_values[imax]:.4f}"
            f"  at T = "
            f"{2 * math.pi / omega_values[imax] / 3600:8.3f} h"
        )

        line += f"   3w/1w up to {third_values.max():.4f}"

        if not all_converged:
            line += "   [WARNING: not converged]"

        print(line)

    level, lines = third_harmonic_quality(worst_q)

    if level != "ok":
        inlet, omega_bad, a_bad = worst_where
        tag = "WARNING" if level == "warn" else "NOTE"

        print()
        print(
            f"  [{tag}] worst |Q3|/|Q1| = {worst_q:.2f} in inlet {inlet} "
            f"at T = {2 * math.pi / omega_bad / 3600:.2f} h, "
            f"a_s = {a_bad:g} m"
        )

        for text in lines:
            print(f"           {text}")

    # ------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------

    omega_diurnal = 2.0 * math.pi / (24.0 * 3600.0)
    omega_semidiurnal = 2.0 * math.pi / (12.42 * 3600.0)

    fig, ax = plt.subplots(figsize=(10.5, 6.0))

    for a_s in amplitudes:
        ax.plot(
            omega_values,
            results[a_s],
            linewidth=1.8,
            label=rf"$a_s={a_s:g}$ m"
        )

    ax.set_xscale("log")
    ax.set_xlim(args.omega_min, args.omega_max)

    ax.set_xlabel(
        r"Angular frequency $\omega$ [rad/s]",
        fontsize=13
    )

    ax.set_ylabel(
        r"Gain $G=a_b/a_s$",
        fontsize=13
    )

    ax.grid(True, which="both", alpha=0.25)

    for omega_ref, label in (
        (omega_diurnal, "24 h"),
        (omega_semidiurnal, "M2"),
    ):
        if args.omega_min <= omega_ref <= args.omega_max:
            ax.axvline(
                omega_ref,
                linestyle="--",
                linewidth=1.0,
                alpha=0.65
            )

            ax.text(
                omega_ref,
                0.97 * ax.get_ylim()[1],
                label,
                rotation=90,
                va="top",
                ha="right"
            )

    ax.legend(title="Sea amplitude", fontsize=9, ncol=2)

    fig.tight_layout()

    input_path = Path(args.input)

    png_path = input_path.with_name(
        input_path.stem + "_gain.png"
    )

    fig.savefig(png_path, dpi=220)

    if args.show:
        plt.show()
    else:
        plt.close(fig)

    # ------------------------------------------------------------
    # CSV
    # ------------------------------------------------------------

    columns = [omega_values]
    headers = ["omega_rad_s"]

    for a_s in amplitudes:
        columns.append(results[a_s])
        headers.append(f"G_as_{a_s:g}m")

        columns.append(results3[a_s])
        headers.append(f"eta3_over_eta1_as_{a_s:g}m")

    csv_path = input_path.with_name(
        input_path.stem + "_gain.csv"
    )

    np.savetxt(
        csv_path,
        np.column_stack(columns),
        delimiter=",",
        header=",".join(headers),
        comments=""
    )

    print()
    print(f"  {png_path}")
    print(f"  {csv_path}")



if __name__ == "__main__":
    main()
