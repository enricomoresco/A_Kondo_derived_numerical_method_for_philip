#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
basin_network_response.py

Frequency-domain tidal response of a network of M bays and N channels.

Every bay is a dynamic node with compliance C_m = S_m.
Every external boundary has a prescribed harmonic water level.
Every channel connects two nodes and has impedance

    Z_e = R_e + i*omega*L_e

with

    L_e = l_e/(g*A_e)

and a Lorentz resistance updated self-consistently:

    R_e = beta_e * |Q_e|

    beta_e = 8*l_e*n_e^2 /
             (3*pi*A_e^2*Rh_e^(4/3))

For a rectangular cross-section:

    A_e  = B_e*D_e
    Rh_e = A_e/(B_e + 2*D_e)
    n_e  = 1/kse_e

At fixed R_e the nodal system is linear and complex:

    A(omega,R) * eta = b

The fixed point then updates R_e from the resulting harmonic discharge.

-----------------------------------------------------------------------
INPUT FORMAT
-----------------------------------------------------------------------

N_bays = 3
N_channels = 7

bay , B1 , 10000000
bay , B2 , 25000000
bay , B3 , 8000000

# boundary,name[,relative_amplitude,phase_deg]
boundary , sea
# equivalent to: boundary , sea , 1.0 , 0.0

# channel,name,node1,node2,L[m],D[m],B[m],kse
channel , C1 , sea , B1 , 3200 , 1.2 , 150.0 , 30
channel , C2 , sea , B1 , 4500 , 3.4 , 255.4 , 30
channel , C3 , sea , B2 , 1800 , 2.5 , 220.0 , 35
channel , C4 , sea , B3 , 2700 , 1.8 , 130.0 , 32
channel , C5 , B1  , B2 , 5000 , 2.0 , 180.0 , 30
channel , C6 , B2  , B3 , 3500 , 1.5 , 120.0 , 28
channel , C7 , B1  , B3 , 6200 , 1.1 , 90.0  , 25

-----------------------------------------------------------------------
EXAMPLE RUN
-----------------------------------------------------------------------

python basin_network_response.py network.txt --show

or

python basin_network_response.py network.txt \
    --omega-min 1e-5 --omega-max 1e-3 --nfreq 500 \
    --amplitudes 0.05 0.1 0.2 0.5 1 2 5 --show

OUTPUT
-----------------------------------------------------------------------

1) Long-format CSV:
       amplitude, omega, period_h,
       G_B1, phase_B1_deg, G3_B1, phase3_B1_deg, ...
   where G3 / phase3 are the perturbative first secondary
   harmonic (3*omega), see third_harmonic_network().

2) One PNG plot per bay:
       G_m(omega) for every requested amplitude.

3) An NPZ file that also stores the full complex eta and Q.
"""

from __future__ import annotations

import argparse
import cmath
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


G_GRAV = 9.81

DEFAULT_OMEGA_MIN = 1.0e-5
DEFAULT_OMEGA_MAX = 1.0e-3
DEFAULT_NFREQ = 400
DEFAULT_AMPLITUDES = (0.05, 0.10, 0.20, 0.50, 1.0, 2.0, 5.0)

DEFAULT_TOL = 1.0e-10
DEFAULT_ITMAX = 5000
DEFAULT_RELAX = 0.5


@dataclass
class Bay:
    name: str
    S: float


@dataclass
class Boundary:
    name: str
    factor: float = 1.0
    phase_deg: float = 0.0

    def eta(self, a_s: float) -> complex:
        return a_s * self.factor * cmath.exp(
            1j * math.radians(self.phase_deg)
        )


@dataclass
class Channel:
    name: str
    node1: str
    node2: str
    length: float
    depth: float
    width: float
    kse: float

    A: float = 0.0
    Rh: float = 0.0
    n: float = 0.0
    inertance: float = 0.0
    beta: float = 0.0

    def prepare(self):
        self.A = self.width * self.depth
        self.Rh = self.A / (self.width + 2.0 * self.depth)
        self.n = 1.0 / self.kse

        self.inertance = self.length / (G_GRAV * self.A)

        # R = beta * |Q|
        self.beta = (
            8.0
            * self.length
            * self.n**2
            / (
                3.0
                * math.pi
                * self.A**2
                * self.Rh**(4.0 / 3.0)
            )
        )


# ======================================================================
# INPUT
# ======================================================================

def read_network(path):
    bays = {}
    boundaries = {}
    channels = []

    declared_n_bays = None
    declared_n_channels = None

    with open(path, "r", encoding="utf-8") as f:
        for lineno, raw in enumerate(f, 1):
            line = raw.split("#", 1)[0].strip()

            if not line:
                continue

            if "=" in line and "," not in line:
                key, value = line.split("=", 1)
                key = key.strip().lower()
                value = value.strip()

                if key == "n_bays":
                    declared_n_bays = int(value)
                elif key == "n_channels":
                    declared_n_channels = int(value)
                else:
                    raise ValueError(
                        f"Line {lineno}: unknown parameter {key!r}"
                    )
                continue

            p = [x.strip() for x in line.split(",")]
            kind = p[0].lower()

            if kind == "bay":
                if len(p) != 3:
                    raise ValueError(
                        f"Line {lineno}: required format: bay,name,S"
                    )

                name = p[1]
                S = float(p[2])

                if S <= 0:
                    raise ValueError(
                        f"Line {lineno}: S must be positive"
                    )

                if name in bays or name in boundaries:
                    raise ValueError(
                        f"Line {lineno}: duplicate node {name!r}"
                    )

                bays[name] = Bay(name=name, S=S)

            elif kind == "boundary":
                if len(p) not in (2, 3, 4):
                    raise ValueError(
                        f"Line {lineno}: required format: "
                        "boundary,name[,factor,phase_deg]"
                    )

                name = p[1]
                factor = float(p[2]) if len(p) >= 3 else 1.0
                phase_deg = float(p[3]) if len(p) == 4 else 0.0

                if factor < 0:
                    raise ValueError(
                        f"Line {lineno}: factor must be >= 0"
                    )

                if name in bays or name in boundaries:
                    raise ValueError(
                        f"Line {lineno}: duplicate node {name!r}"
                    )

                boundaries[name] = Boundary(
                    name=name,
                    factor=factor,
                    phase_deg=phase_deg,
                )

            elif kind == "channel":
                # channel,name,node1,node2,L,D,B,kse
                if len(p) != 8:
                    raise ValueError(
                        f"Line {lineno}: required format: "
                        "channel,name,node1,node2,L,D,B,kse"
                    )

                ch = Channel(
                    name=p[1],
                    node1=p[2],
                    node2=p[3],
                    length=float(p[4]),
                    depth=float(p[5]),
                    width=float(p[6]),
                    kse=float(p[7]),
                )

                for field_name, value in (
                    ("L", ch.length),
                    ("D", ch.depth),
                    ("B", ch.width),
                    ("kse", ch.kse),
                ):
                    if value <= 0:
                        raise ValueError(
                            f"Line {lineno}: {field_name} must be positive"
                        )

                ch.prepare()
                channels.append(ch)

            else:
                raise ValueError(
                    f"Line {lineno}: unknown entry type {kind!r}"
                )

    if not bays:
        raise ValueError("No bay defined")

    if not boundaries:
        raise ValueError("No external boundary defined")

    if not channels:
        raise ValueError("No channel defined")

    nodes = set(bays) | set(boundaries)

    channel_names = set()

    for ch in channels:
        if ch.name in channel_names:
            raise ValueError(
                f"Duplicate channel name: {ch.name!r}"
            )
        channel_names.add(ch.name)

        if ch.node1 not in nodes:
            raise ValueError(
                f"Channel {ch.name}: node {ch.node1!r} is not defined"
            )
        if ch.node2 not in nodes:
            raise ValueError(
                f"Channel {ch.name}: node {ch.node2!r} is not defined"
            )
        if ch.node1 == ch.node2:
            raise ValueError(
                f"Channel {ch.name}: the two nodes coincide"
            )

        # A boundary-to-boundary channel is useless in the bay model.
        if ch.node1 in boundaries and ch.node2 in boundaries:
            raise ValueError(
                f"Channel {ch.name}: it connects two external boundaries; "
                "it does not contribute to the bay dynamics"
            )

    if declared_n_bays is not None and declared_n_bays != len(bays):
        raise ValueError(
            f"N_bays={declared_n_bays}, but {len(bays)} bays are defined"
        )

    if (
        declared_n_channels is not None
        and declared_n_channels != len(channels)
    ):
        raise ValueError(
            f"N_channels={declared_n_channels}, "
            f"but {len(channels)} channels are defined"
        )

    return bays, boundaries, channels


# ======================================================================
# LINEAR NETWORK SOLVE FOR FIXED RESISTANCES
# ======================================================================

def solve_linear_network(
    bays,
    boundaries,
    channels,
    omega,
    a_s,
    resistances,
):
    """
    Solve for the complex bay levels at fixed R_e.

    KCL at bay m:

        i*omega*S_m*eta_m
        = sum_e Y_e (eta_neighbor - eta_m)

    hence:

        (i*omega*S_m + sum Y_e) eta_m
        - sum_internal Y_e eta_n
        = sum_boundary Y_e eta_boundary
    """

    bay_names = list(bays.keys())
    idx = {name: i for i, name in enumerate(bay_names)}
    M = len(bay_names)

    A = np.zeros((M, M), dtype=complex)
    rhs = np.zeros(M, dtype=complex)

    # Storage/compliance
    for name, bay in bays.items():
        i = idx[name]
        A[i, i] += 1j * omega * bay.S

    # Boundary elevations
    eta_boundary = {
        name: bc.eta(a_s)
        for name, bc in boundaries.items()
    }

    admittances = np.empty(len(channels), dtype=complex)

    for e, ch in enumerate(channels):
        Z = complex(
            resistances[e],
            omega * ch.inertance,
        )

        Y = 1.0 / Z
        admittances[e] = Y

        n1_is_bay = ch.node1 in bays
        n2_is_bay = ch.node2 in bays

        if n1_is_bay and n2_is_bay:
            i = idx[ch.node1]
            j = idx[ch.node2]

            A[i, i] += Y
            A[j, j] += Y
            A[i, j] -= Y
            A[j, i] -= Y

        elif n1_is_bay:
            i = idx[ch.node1]
            eta_ext = eta_boundary[ch.node2]

            A[i, i] += Y
            rhs[i] += Y * eta_ext

        elif n2_is_bay:
            j = idx[ch.node2]
            eta_ext = eta_boundary[ch.node1]

            A[j, j] += Y
            rhs[j] += Y * eta_ext

        else:
            raise RuntimeError(
                f"Channel {ch.name}: neither end is a bay"
            )

    try:
        eta_vec = np.linalg.solve(A, rhs)
    except np.linalg.LinAlgError as exc:
        raise RuntimeError(
            "Singular or ill-conditioned nodal system. "
            "Check the topology and the parameters."
        ) from exc

    eta = {
        name: eta_vec[idx[name]]
        for name in bay_names
    }

    # All nodal levels, boundary conditions included
    eta_all = dict(eta_boundary)
    eta_all.update(eta)

    # Q is positive from node1 towards node2
    Q = np.empty(len(channels), dtype=complex)

    for e, ch in enumerate(channels):
        Q[e] = (
            admittances[e]
            * (
                eta_all[ch.node1]
                - eta_all[ch.node2]
            )
        )

    return eta, Q, admittances


# ======================================================================
# SELF-CONSISTENT FIXED POINT
# ======================================================================

def initial_resistances(
    bays,
    boundaries,
    channels,
    omega,
    a_s,
):
    """
    Simple but physically scaled guess.

    The discharge scale is omega*S_tot*a_s divided by the number of
    external connections. It is only used to start the fixed point.
    """

    S_tot = sum(bay.S for bay in bays.values())

    n_external = sum(
        (ch.node1 in boundaries) or (ch.node2 in boundaries)
        for ch in channels
    )

    q_scale = (
        omega * S_tot * a_s
        / max(n_external, 1)
    )

    # Avoid an exactly vanishing R.
    q_scale = max(q_scale, 1.0e-12)

    return np.array(
        [ch.beta * q_scale for ch in channels],
        dtype=float,
    )


def solve_network(
    bays,
    boundaries,
    channels,
    omega,
    a_s,
    R0=None,
    tol=DEFAULT_TOL,
    itmax=DEFAULT_ITMAX,
    relax=DEFAULT_RELAX,
):
    """
    Self-consistent fixed point:

        R_e -> Z_e -> eta_m -> Q_e -> beta_e |Q_e| -> R_e

    with damping/relaxation.
    """

    if R0 is None:
        R = initial_resistances(
            bays,
            boundaries,
            channels,
            omega,
            a_s,
        )
    else:
        R = np.array(R0, dtype=float, copy=True)

    beta = np.array(
        [ch.beta for ch in channels],
        dtype=float,
    )

    eta_old = None

    for iteration in range(1, itmax + 1):
        eta, Q, _ = solve_linear_network(
            bays,
            boundaries,
            channels,
            omega,
            a_s,
            R,
        )

        R_target = beta * np.abs(Q)

        scale_R = np.maximum(
            np.maximum(np.abs(R), np.abs(R_target)),
            1.0e-14,
        )

        err_R = np.max(
            np.abs(R_target - R) / scale_R
        )

        if eta_old is None:
            err_eta = np.inf
        else:
            e_new = np.array(
                [eta[name] for name in bays],
                dtype=complex,
            )
            e_old = np.array(
                [eta_old[name] for name in bays],
                dtype=complex,
            )

            scale_eta = np.maximum(
                np.maximum(np.abs(e_new), np.abs(e_old)),
                max(a_s, 1.0e-14) * 1.0e-12,
            )

            err_eta = np.max(
                np.abs(e_new - e_old) / scale_eta
            )

        if err_R <= tol and err_eta <= tol:
            # One final evaluation with the self-consistent R.
            R = R_target

            eta, Q, Y = solve_linear_network(
                bays,
                boundaries,
                channels,
                omega,
                a_s,
                R,
            )

            return eta, Q, R, Y, iteration, True

        eta_old = eta

        R = (
            (1.0 - relax) * R
            + relax * R_target
        )

    eta, Q, Y = solve_linear_network(
        bays,
        boundaries,
        channels,
        omega,
        a_s,
        R,
    )

    return eta, Q, R, Y, itmax, False


# ======================================================================
# THIRD HARMONIC (PERTURBATIVE)
# ======================================================================

# Quality of the third harmonic, judged on the per-channel ratio |Q3|/|Q1|.
# For a single basin that ratio stays around 0.1 - 0.15 (it is exactly 2/15 in
# the friction-dominated limit of a single inlet), and validation against a
# time-domain integration of the same equations puts the error on the 3*omega
# amplitude within about 10% there. Larger values mean 3*omega is close to a
# resonance of the network, where the drag linearisation - calibrated on the
# fundamental discharge - underestimates its own damping.
# The fundamental is not affected in either case.
Q3_CAUTION = 0.20
Q3_WARN = 0.50


def third_harmonic_quality(ratio_max):
    """Return (level, lines) describing how far the perturbation is trusted."""

    if ratio_max > Q3_WARN:
        return "warn", [
            "3*omega is at or near a resonance of the network, so the "
            "perturbation has",
            "broken down: treat the 3w amplitude as indicative only. "
            "The fundamental is",
            "not affected.",
        ]

    if ratio_max > Q3_CAUTION:
        return "caution", [
            "3*omega is approaching a resonance of the network: the 3w "
            "amplitude may be",
            "in error by more than 10%. The fundamental is not affected.",
        ]

    return "ok", []


def third_harmonic_network(
    bays,
    boundaries,
    channels,
    omega,
    Q1,
):
    """
    First secondary harmonic, given the converged fundamental Q1.

    The quadratic drag expands, for Q_e = |Q_e| cos(theta_e), as

        Q|Q| = |Q|^2 * sum_{n odd} 8 (-1)^((n-1)/2) cos(n theta)
                                   / (pi n (4 - n^2))

    The n = 1 term is the Lorentz resistance already carried by the
    fundamental. The n = 3 term is a known source once the fundamental is
    known, and the perturbation Q3 sees the linearised drag

        d/dQ (Q|Q|) = 2|Q|   ->   <2|Q| |cos|> = (4/pi)|Q|

    so that, relative to the Lorentz coefficient (8/3pi)|Q|,

        R_e^(3) = (3/2) * beta_e * |Q_e|
        S_e     = (1/5) * beta_e * Q_e^3 / |Q_e|      (branch e.m.f. at 3w)

    The resulting 3w problem is LINEAR: one nodal solve, no fixed point.
    The external boundaries are monochromatic, hence eta3 = 0 there.
    """

    beta = np.array(
        [ch.beta for ch in channels],
        dtype=float,
    )

    absQ = np.abs(Q1)
    safe = np.where(absQ > 0.0, absQ, 1.0)

    R3 = 1.5 * beta * absQ
    S = 0.2 * beta * Q1**3 / safe

    bay_names = list(bays.keys())
    idx = {name: i for i, name in enumerate(bay_names)}
    M = len(bay_names)

    A = np.zeros((M, M), dtype=complex)
    rhs = np.zeros(M, dtype=complex)

    for name, bay in bays.items():
        A[idx[name], idx[name]] += 3j * omega * bay.S

    Y3 = np.empty(len(channels), dtype=complex)

    for e, ch in enumerate(channels):
        Y3[e] = 1.0 / complex(
            R3[e],
            3.0 * omega * ch.inertance,
        )

        n1_is_bay = ch.node1 in bays
        n2_is_bay = ch.node2 in bays

        if n1_is_bay and n2_is_bay:
            i = idx[ch.node1]
            j = idx[ch.node2]

            A[i, i] += Y3[e]
            A[j, j] += Y3[e]
            A[i, j] -= Y3[e]
            A[j, i] -= Y3[e]

            rhs[i] += Y3[e] * S[e]
            rhs[j] -= Y3[e] * S[e]

        elif n1_is_bay:
            i = idx[ch.node1]
            A[i, i] += Y3[e]
            rhs[i] += Y3[e] * S[e]

        else:
            j = idx[ch.node2]
            A[j, j] += Y3[e]
            rhs[j] -= Y3[e] * S[e]

    try:
        eta3_vec = np.linalg.solve(A, rhs)
    except np.linalg.LinAlgError as exc:
        raise RuntimeError(
            "Singular nodal system at 3*omega."
        ) from exc

    eta3 = {
        name: eta3_vec[idx[name]]
        for name in bay_names
    }

    eta3_all = {name: 0.0 + 0.0j for name in boundaries}
    eta3_all.update(eta3)

    Q3 = np.array(
        [
            Y3[e]
            * (
                eta3_all[ch.node1]
                - eta3_all[ch.node2]
                - S[e]
            )
            for e, ch in enumerate(channels)
        ],
        dtype=complex,
    )

    return eta3, Q3


# ======================================================================
# SWEEP
# ======================================================================

def run_sweep(
    bays,
    boundaries,
    channels,
    omega_values,
    amplitudes,
    tol,
    itmax,
    relax,
):
    bay_names = list(bays.keys())
    channel_names = [ch.name for ch in channels]

    Na = len(amplitudes)
    Nw = len(omega_values)
    M = len(bay_names)
    Ne = len(channels)

    eta_arr = np.empty((Na, Nw, M), dtype=complex)
    Q_arr = np.empty((Na, Nw, Ne), dtype=complex)
    R_arr = np.empty((Na, Nw, Ne), dtype=float)
    iterations = np.empty((Na, Nw), dtype=int)
    converged = np.empty((Na, Nw), dtype=bool)

    eta3_arr = np.zeros((Na, Nw, M), dtype=complex)
    Q3_arr = np.zeros((Na, Nw, Ne), dtype=complex)

    for ia, a_s in enumerate(amplitudes):
        # Continuation in frequency:
        # the converged R at one omega is a good initial guess for the next.
        R_guess = None

        for iw, omega in enumerate(omega_values):
            eta, Q, R, _, nit, ok = solve_network(
                bays=bays,
                boundaries=boundaries,
                channels=channels,
                omega=omega,
                a_s=a_s,
                R0=R_guess,
                tol=tol,
                itmax=itmax,
                relax=relax,
            )

            eta_arr[ia, iw, :] = [
                eta[name] for name in bay_names
            ]
            Q_arr[ia, iw, :] = Q
            R_arr[ia, iw, :] = R
            iterations[ia, iw] = nit
            converged[ia, iw] = ok

            eta3, Q3 = third_harmonic_network(
                bays,
                boundaries,
                channels,
                omega,
                Q,
            )

            eta3_arr[ia, iw, :] = [
                eta3[name] for name in bay_names
            ]
            Q3_arr[ia, iw, :] = Q3

            R_guess = R

    return {
        "bay_names": bay_names,
        "channel_names": channel_names,
        "eta": eta_arr,
        "Q": Q_arr,
        "R": R_arr,
        "iterations": iterations,
        "converged": converged,
        "eta3": eta3_arr,
        "Q3": Q3_arr,
    }


# ======================================================================
# OUTPUT
# ======================================================================

def save_csv(
    path,
    result,
    omega_values,
    amplitudes,
):
    bay_names = result["bay_names"]
    eta = result["eta"]
    iterations = result["iterations"]
    converged = result["converged"]

    headers = [
        "a_s_m",
        "omega_rad_s",
        "period_h",
    ]

    eta3 = result.get("eta3")

    for name in bay_names:
        headers.extend([
            f"G_{name}",
            f"phase_{name}_deg",
        ])

        if eta3 is not None:
            headers.extend([
                f"G3_{name}",
                f"phase3_{name}_deg",
            ])

    headers.extend([
        "iterations",
        "converged",
    ])

    rows = []

    for ia, a_s in enumerate(amplitudes):
        for iw, omega in enumerate(omega_values):
            period_h = (
                2.0 * math.pi / omega / 3600.0
            )

            row = [
                a_s,
                omega,
                period_h,
            ]

            for im, _ in enumerate(bay_names):
                z = eta[ia, iw, im]
                G = abs(z) / a_s
                phase = math.degrees(cmath.phase(z))

                row.extend([
                    G,
                    phase,
                ])

                if eta3 is not None:
                    z3 = eta3[ia, iw, im]

                    row.extend([
                        abs(z3) / a_s,
                        math.degrees(cmath.phase(z3)),
                    ])

            row.extend([
                iterations[ia, iw],
                int(converged[ia, iw]),
            ])

            rows.append(row)

    np.savetxt(
        path,
        np.asarray(rows, dtype=float),
        delimiter=",",
        header=",".join(headers),
        comments="",
    )


def save_plots(
    out_dir,
    stem,
    result,
    omega_values,
    amplitudes,
    show=False,
):
    eta = result["eta"]
    bay_names = result["bay_names"]

    generated = []

    # Reference tidal frequencies
    omega_24h = 2.0 * math.pi / (24.0 * 3600.0)
    omega_M2 = 2.0 * math.pi / (12.42 * 3600.0)

    for im, bay_name in enumerate(bay_names):
        fig, ax = plt.subplots(
            figsize=(9.5, 5.8)
        )

        for ia, a_s in enumerate(amplitudes):
            G = (
                np.abs(eta[ia, :, im])
                / a_s
            )

            ax.plot(
                omega_values,
                G,
                linewidth=1.8,
                label=rf"$a_s={a_s:g}$ m",
            )

        ax.set_xscale("log")
        ax.set_xlim(
            omega_values.min(),
            omega_values.max(),
        )

        ax.set_xlabel(
            r"Angular frequency $\omega$ [rad/s]"
        )
        ax.set_ylabel(
            rf"Gain $G_{{{bay_name}}}=|\hat{{\eta}}_{{{bay_name}}}|/a_s$"
        )

        ax.grid(
            True,
            which="both",
            alpha=0.25,
        )

        if (
            omega_values.min()
            <= omega_24h
            <= omega_values.max()
        ):
            ax.axvline(
                omega_24h,
                linestyle="--",
                linewidth=1.0,
                alpha=0.55,
            )

        if (
            omega_values.min()
            <= omega_M2
            <= omega_values.max()
        ):
            ax.axvline(
                omega_M2,
                linestyle="--",
                linewidth=1.0,
                alpha=0.55,
            )

        ax.legend(
            title="Incoming amplitude",
            fontsize=8.5,
            ncol=2,
        )

        fig.tight_layout()

        safe_name = "".join(
            c if c.isalnum() or c in "-_" else "_"
            for c in bay_name
        )

        png = out_dir / (
            f"{stem}_gain_{safe_name}.png"
        )

        fig.savefig(
            png,
            dpi=220,
        )

        generated.append(png)

        if show:
            plt.show()
        else:
            plt.close(fig)

    return generated


# ======================================================================
# MAIN
# ======================================================================

def build_parser():

    parser = argparse.ArgumentParser(
        description=(
            "Tidal gain of a network of semi-enclosed basins "
            "connected by frictional channels."
        ),
        epilog=(
            "Examples:\n"
            "  one point   : %(prog)s network.txt --period 12.42 "
            "--amplitude 0.5\n"
            "  full sweep  : %(prog)s network.txt\n"
            "  custom sweep: %(prog)s network.txt --amplitudes 0.2 1 "
            "--nfreq 800"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "input",
        help="input .txt file describing the network",
    )

    point = parser.add_argument_group(
        "single point (give --omega or --period to use it)"
    )

    point.add_argument(
        "--omega",
        type=float,
        help="forcing angular frequency [rad/s]",
    )

    point.add_argument(
        "--period",
        type=float,
        help="forcing period [hours], alternative to --omega",
    )

    point.add_argument(
        "--amplitude",
        type=float,
        help="sea amplitude a_s [m]",
    )

    point.add_argument(
        "--channels",
        action="store_true",
        help="also list the discharge carried by each channel",
    )

    sweep = parser.add_argument_group("sweep")

    sweep.add_argument(
        "--amplitudes",
        type=float,
        nargs="+",
        default=list(DEFAULT_AMPLITUDES),
        help="sea amplitudes to compare [m]",
    )

    sweep.add_argument(
        "--omega-min",
        type=float,
        default=DEFAULT_OMEGA_MIN,
        help="minimum angular frequency [rad/s]",
    )

    sweep.add_argument(
        "--omega-max",
        type=float,
        default=DEFAULT_OMEGA_MAX,
        help="maximum angular frequency [rad/s]",
    )

    sweep.add_argument(
        "--nfreq",
        type=int,
        default=DEFAULT_NFREQ,
        help="number of log-spaced frequencies",
    )

    sweep.add_argument(
        "--npz",
        action="store_true",
        help="also save the full complex fields to a .npz archive",
    )

    sweep.add_argument(
        "--show",
        action="store_true",
        help="display the plots in addition to saving them",
    )

    numerics = parser.add_argument_group("numerics")

    numerics.add_argument(
        "--tol",
        type=float,
        default=DEFAULT_TOL,
        help="fixed-point tolerance",
    )

    numerics.add_argument(
        "--itmax",
        type=int,
        default=DEFAULT_ITMAX,
        help="maximum number of iterations",
    )

    numerics.add_argument(
        "--relax",
        type=float,
        default=DEFAULT_RELAX,
        help="relaxation factor, 0 < relax <= 1",
    )

    return parser


def report_point(
    bays,
    boundaries,
    channels,
    omega,
    a_s,
    show_channels,
    tol,
    itmax,
    relax,
):
    """Print the response of every basin at one frequency, one amplitude."""

    eta, Q, R, Y, nit, ok = solve_network(
        bays=bays,
        boundaries=boundaries,
        channels=channels,
        omega=omega,
        a_s=a_s,
        tol=tol,
        itmax=itmax,
        relax=relax,
    )

    eta3, Q3 = third_harmonic_network(
        bays,
        boundaries,
        channels,
        omega,
        Q,
    )

    qratio = np.abs(Q3) / np.abs(Q)

    period_h = 2.0 * math.pi / omega / 3600.0

    print()
    print(
        f"Forcing    a_s = {a_s:g} m"
        f"     omega = {omega:.4e} rad/s"
        f"     T = {period_h:.3f} h"
    )
    print()

    width = max(
        [len("basin")]
        + [len(name) for name in bays]
    )

    print(
        f"  {'basin':<{width}}"
        f"    gain   amplitude [m]   range [m]"
        f"   lag [deg]   lag [h]   3w amp [m]"
    )

    for name in bays:
        z = eta[name]
        G = abs(z) / a_s
        lag_deg = -math.degrees(cmath.phase(z))
        lag_h = lag_deg / 360.0 * period_h

        print(
            f"  {name:<{width}}"
            f"  {G:7.4f}"
            f"  {abs(z):13.4f}"
            f"  {2.0 * abs(z):10.4f}"
            f"  {lag_deg:10.2f}"
            f"  {lag_h:8.3f}"
            f"  {abs(eta3[name]):11.4f}"
        )

    if not ok:
        print()
        print("  [WARNING] the fixed point did not converge")

    if show_channels:
        cwidth = max(
            [len("channel")]
            + [len(ch.name) for ch in channels]
        )

        nwidth = max(
            [len(ch.node1) for ch in channels]
            + [len(ch.node2) for ch in channels]
        )

        link_width = max(
            2 * nwidth + 4,
            len("from -> to"),
        )

        print()
        print(
            f"  {'channel':<{cwidth}}"
            f"  {'from -> to':<{link_width}}"
            f"  {'|Q| [m3/s]':>12}"
            f"  {'lag [deg]':>10}"
            f"  {'|Q3|/|Q1|':>10}"
        )

        for e, ch in enumerate(channels):
            link = f"{ch.node1:>{nwidth}} -> {ch.node2:<{nwidth}}"

            print(
                f"  {ch.name:<{cwidth}}"
                f"  {link:<{link_width}}"
                f"  {abs(Q[e]):12.1f}"
                f"  {-math.degrees(cmath.phase(Q[e])):10.2f}"
                f"  {qratio[e]:10.3f}"
            )

    level, lines = third_harmonic_quality(qratio.max())

    if level != "ok":
        worst = channels[int(np.argmax(qratio))].name
        tag = "WARNING" if level == "warn" else "NOTE"

        print()
        print(
            f"  [{tag}] worst |Q3|/|Q1| = {qratio.max():.2f} "
            f"in channel {worst}"
        )

        for line in lines:
            print(f"           {line}")

    print()


def main():

    parser = build_parser()
    args = parser.parse_args()

    if not (0 < args.relax <= 1):
        parser.error("relax must lie between 0 and 1")

    bays, boundaries, channels = read_network(args.input)

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

        if args.amplitude is None:
            parser.error(
                "a single-point run needs --amplitude"
            )

        if args.amplitude <= 0:
            parser.error("the amplitude must be positive")

        report_point(
            bays,
            boundaries,
            channels,
            omega,
            args.amplitude,
            args.channels,
            args.tol,
            args.itmax,
            args.relax,
        )
        return

    # ------------------------------------------------------------
    # Sweep
    # ------------------------------------------------------------

    if not (0 < args.omega_min < args.omega_max):
        parser.error("0 < omega-min < omega-max is required")

    if args.nfreq < 2:
        parser.error("nfreq must be >= 2")

    amplitudes = np.array(args.amplitudes, dtype=float)

    if np.any(amplitudes <= 0):
        parser.error("all amplitudes must be positive")

    omega_values = np.logspace(
        math.log10(args.omega_min),
        math.log10(args.omega_max),
        args.nfreq,
    )

    print(
        f"{len(bays)} basins, {len(channels)} channels, "
        f"{len(boundaries)} boundaries   "
        f"omega {args.omega_min:.1e} - {args.omega_max:.1e} rad/s   "
        f"({args.nfreq} points)"
    )
    print()

    result = run_sweep(
        bays=bays,
        boundaries=boundaries,
        channels=channels,
        omega_values=omega_values,
        amplitudes=amplitudes,
        tol=args.tol,
        itmax=args.itmax,
        relax=args.relax,
    )

    n_bad = int(
        result["converged"].size
        - np.count_nonzero(result["converged"])
    )

    for im, bay_name in enumerate(result["bay_names"]):
        print(f"  [{bay_name}]")

        for ia, a_s in enumerate(amplitudes):
            G = np.abs(result["eta"][ia, :, im]) / a_s
            j = int(np.argmax(G))

            line = (
                f"    a_s = {a_s:6g} m :  Gmax = {G[j]:.4f}"
                f"  at T = "
                f"{2 * math.pi / omega_values[j] / 3600:8.3f} h"
            )

            ratio = (
                np.abs(result["eta3"][ia, :, im])
                / np.abs(result["eta"][ia, :, im])
            )

            line += f"   3w/1w up to {ratio.max():.4f}"

            print(line)

    if n_bad:
        print()
        print(
            f"  [WARNING] {n_bad} of "
            f"{result['converged'].size} points did not converge"
        )

    # Is the third harmonic still a perturbation everywhere?
    qratio = (
        np.abs(result["Q3"])
        / np.abs(result["Q"])
    )

    level, lines = third_harmonic_quality(qratio.max())

    if level != "ok":
        ia, iw, ie = np.unravel_index(
            int(np.argmax(qratio)),
            qratio.shape,
        )

        tag = "WARNING" if level == "warn" else "NOTE"

        print()
        print(
            f"  [{tag}] worst |Q3|/|Q1| = {qratio.max():.2f} "
            f"in channel {result['channel_names'][ie]} "
            f"at T = {2 * math.pi / omega_values[iw] / 3600:.2f} h, "
            f"a_s = {amplitudes[ia]:g} m"
        )

        for line in lines:
            print(f"           {line}")

    input_path = Path(args.input)
    out_dir = input_path.parent
    stem = input_path.stem

    csv_path = out_dir / (stem + "_gain.csv")

    save_csv(
        csv_path,
        result,
        omega_values,
        amplitudes,
    )

    plots = save_plots(
        out_dir=out_dir,
        stem=stem,
        result=result,
        omega_values=omega_values,
        amplitudes=amplitudes,
        show=args.show,
    )

    print()
    print(f"  {csv_path}")

    for path in plots:
        print(f"  {path}")

    if args.npz:
        npz_path = out_dir / (stem + "_fields.npz")

        np.savez_compressed(
            npz_path,
            omega=omega_values,
            amplitudes=amplitudes,
            bay_names=np.array(result["bay_names"], dtype=str),
            channel_names=np.array(
                result["channel_names"], dtype=str
            ),
            eta=result["eta"],
            Q=result["Q"],
            R=result["R"],
            eta3=(
                result["eta3"]
                if result["eta3"] is not None
                else np.zeros(0, dtype=complex)
            ),
            Q3=(
                result["Q3"]
                if result["Q3"] is not None
                else np.zeros(0, dtype=complex)
            ),
            iterations=result["iterations"],
            converged=result["converged"],
        )

        print(f"  {npz_path}")



if __name__ == "__main__":
    main()
