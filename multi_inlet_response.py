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
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


G_GRAV = 9.81

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
                # Ignored: a_s is swept in this program.
                pass

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

    _, _, Y, Ysum = equivalent_impedance(
        G,
        phi,
        inertance,
        rgeo,
        omega
    )

    # Complex fundamental, with the sea level taken as a_s + 0j
    denom = Ysum + 1j * omega * S

    eta1 = Ysum * a_s / denom
    Q1 = Y * a_s * 1j * omega * S / denom

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

    return eta3, eta1


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

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Basin gain as a function of omega "
            "for several forcing amplitudes."
        )
    )

    parser.add_argument(
        "input",
        help="input .txt file"
    )

    parser.add_argument(
        "--no-third-harmonic",
        action="store_true",
        help=(
            "skip the perturbative 3*omega correction "
            "(it costs a few percent of the run time)"
        )
    )

    parser.add_argument(
        "--show",
        action="store_true",
        help="display the plot"
    )

    args = parser.parse_args()

    data = read_input(args.input)

    N = data["n"]
    S = data["s"]

    # Angular frequencies [rad/s]
    omega_values = np.logspace(
        np.log10(OMEGA_MIN),
        np.log10(OMEGA_MAX),
        N_OMEGA
    )

    third = not args.no_third_harmonic

    results = {}
    results3 = {}

    print("=" * 72)
    print("SINGLE BASIN, MULTIPLE INLETS — FREQUENCY / AMPLITUDE RESPONSE")
    print("=" * 72)

    print(f"N = {N}")
    print(f"S = {S:.6g} m^2")

    print()
    print(
        f"omega range = "
        f"{OMEGA_MIN:.1e} - "
        f"{OMEGA_MAX:.1e} rad/s"
    )

    print(
        "a_s = "
        + ", ".join(
            f"{x:g}" for x in AMPLITUDES
        )
        + " m"
    )

    print()

    # --------------------------------------------------------
    # Compute every curve
    # --------------------------------------------------------

    for a_s in AMPLITUDES:

        A, Rh, inertance, rgeo = (
            prepare_geometry(
                data,
                a_s
            )
        )

        gain_values = np.empty_like(
            omega_values
        )

        third_values = np.zeros_like(
            omega_values
        )

        max_iterations = 0
        all_converged = True

        for i, omega in enumerate(
            omega_values
        ):

            G, phi, nit, ok = self_consistent_gain(
                S,
                omega,
                inertance,
                rgeo
            )

            gain_values[i] = G

            if third:
                eta3, eta1 = third_harmonic(
                    S,
                    omega,
                    a_s,
                    inertance,
                    rgeo,
                    G,
                    phi
                )

                third_values[i] = (
                    abs(eta3) / abs(eta1)
                )

            max_iterations = max(
                max_iterations,
                nit
            )

            if not ok:
                all_converged = False

        results[a_s] = gain_values
        results3[a_s] = third_values

        imax = int(
            np.argmax(gain_values)
        )

        line = (
            f"a_s = {a_s:5.2f} m : "
            f"Gmax = {gain_values[imax]:.5f} "
            f"at omega = "
            f"{omega_values[imax]:.4e} rad/s"
            f"   max iter = {max_iterations}"
        )

        if third:
            line += (
                f"   max|eta3/eta1| = "
                f"{third_values.max():.4f}"
            )

            if third_values.max() > 0.10:
                line += "  [PERTURBATION SUSPECT]"

        if not all_converged:
            line += "  [WARNING]"

        print(line)

    # --------------------------------------------------------
    # Reference tidal frequencies
    # --------------------------------------------------------

    # 24 h: diurnal reference
    omega_diurnal = (
        2.0 * math.pi
        / (24.0 * 3600.0)
    )

    # M2: 12.42 h
    omega_semidiurnal = (
        2.0 * math.pi
        / (12.42 * 3600.0)
    )

    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(10.5, 6.0)
    )

    for a_s in AMPLITUDES:

        ax.plot(
            omega_values,
            results[a_s],
            linewidth=1.8,
            label=rf"$a_s={a_s:g}$ m"
        )

    ax.set_xscale("log")

    ax.set_xlim(
        OMEGA_MIN,
        OMEGA_MAX
    )

    ax.set_xlabel(
        r"Angular frequency $\omega$ [rad/s]",
        fontsize=13
    )

    ax.set_ylabel(
        r"Gain $G=a_b/a_s$",
        fontsize=13
    )

    ax.grid(
        True,
        which="both",
        alpha=0.25
    )

    # Tidal reference lines
    ax.axvline(
        omega_diurnal,
        linestyle="--",
        linewidth=1.0,
        alpha=0.65
    )

    ax.axvline(
        omega_semidiurnal,
        linestyle="--",
        linewidth=1.0,
        alpha=0.65
    )

    ymax = ax.get_ylim()[1]

    ax.text(
        omega_diurnal,
        0.97 * ymax,
        "24 h",
        rotation=90,
        va="top",
        ha="right"
    )

    ax.text(
        omega_semidiurnal,
        0.97 * ymax,
        "M2",
        rotation=90,
        va="top",
        ha="right"
    )

    ax.legend(
        title="Sea amplitude",
        fontsize=9,
        ncol=2
    )

    fig.tight_layout()

    input_path = Path(args.input)

    png_path = input_path.with_name(
        input_path.stem
        + "_response.png"
    )

    pdf_path = input_path.with_name(
        input_path.stem
        + "_response.pdf"
    )

    fig.savefig(
        png_path,
        dpi=220
    )

    fig.savefig(
        pdf_path
    )

    # --------------------------------------------------------
    # Third-harmonic figure
    # --------------------------------------------------------

    png3_path = None

    if third:

        fig3, ax3 = plt.subplots(
            figsize=(10.5, 5.2)
        )

        for a_s in AMPLITUDES:

            ax3.plot(
                omega_values,
                results3[a_s],
                linewidth=1.6,
                label=rf"$a_s={a_s:g}$ m"
            )

        ax3.axhline(
            2.0 / 45.0,
            color="k",
            linestyle=":",
            linewidth=1.0,
            label=r"bound $2/45$"
        )

        ax3.set_xscale("log")

        ax3.set_xlim(
            OMEGA_MIN,
            OMEGA_MAX
        )

        ax3.set_xlabel(
            r"Angular frequency $\omega$ [rad/s]",
            fontsize=13
        )

        ax3.set_ylabel(
            r"$|\hat{\eta}^{(3)}|/|\hat{\eta}^{(1)}|$",
            fontsize=13
        )

        ax3.grid(
            True,
            which="both",
            alpha=0.25
        )

        ax3.legend(
            fontsize=9,
            ncol=2
        )

        fig3.tight_layout()

        png3_path = input_path.with_name(
            input_path.stem
            + "_third_harmonic.png"
        )

        fig3.savefig(
            png3_path,
            dpi=220
        )

        if not args.show:
            plt.close(fig3)

    # --------------------------------------------------------
    # CSV
    # --------------------------------------------------------

    csv_path = input_path.with_name(
        input_path.stem
        + "_response.csv"
    )

    columns = [omega_values]

    headers = ["omega_rad_s"]

    for a_s in AMPLITUDES:

        columns.append(
            results[a_s]
        )

        headers.append(
            f"G_as_{a_s:g}m"
        )

        if third:
            columns.append(
                results3[a_s]
            )

            headers.append(
                f"eta3_over_eta1_as_{a_s:g}m"
            )

    table = np.column_stack(
        columns
    )

    np.savetxt(
        csv_path,
        table,
        delimiter=",",
        header=",".join(headers),
        comments=""
    )

    print()
    print(
        f"PNG: {png_path}"
    )

    print(
        f"PDF: {pdf_path}"
    )

    print(
        f"CSV: {csv_path}"
    )

    if png3_path is not None:
        print(
            f"PNG: {png3_path}"
        )

    if args.show:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
