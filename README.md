# Harmonic Response of Semi-Enclosed Basins

Two dependency-light Python programs that compute the **tidal gain**

$$G = \frac{a_b}{a_s}$$

of one or more semi-enclosed basins connected to the open sea through frictional
inlets, using a lumped **resistance–inertance (R–L) network** with the quadratic
drag treated by Lorentz equivalent linearisation, closed **self-consistently** on
the forcing amplitude. A perturbative first secondary harmonic (`3 omega`) is
computed on top of the fundamental.

The formulation generalises the classical single-inlet treatment to (i) an
arbitrary number of parallel inlets and (ii) an arbitrary network of
interconnected basins and channels.

---

## Physical model

Each channel `e` is a series impedance

```
Z_e = R_e + i * omega * L_e
```

with inertance

```
L_e = l_e / (g * A_e)
```

and a Lorentz-equivalent linear resistance proportional to the discharge amplitude

```
R_e    = beta_e * |Q_e|
beta_e = 8 * l_e * n_e^2 / (3 * pi * A_e^2 * Rh_e^(4/3))
```

Each basin `m` is a storage element of surface area `S_m`, so continuity reads

```
i * omega * S_m * eta_m = sum_e Y_e * (eta_neighbour - eta_m),      Y_e = 1 / Z_e
```

For a rectangular cross-section the geometry follows from the input as

```
A_e  = B_e * D_e
Rh_e = A_e / (B_e + 2 * D_e)
n_e  = 1 / kse_e            (Manning from the Strickler coefficient)
```

Because `R_e` depends on `|Q_e|`, which in turn depends on `R_e`, the system is
nonlinear. Both programs close it with a **damped fixed-point iteration**

```
R_e -> Z_e -> eta_m -> Q_e -> beta_e |Q_e| -> R_e
```

until resistances and basin elevations stop changing. The response is therefore
amplitude-dependent: the gain curves are **not** the same for different `a_s`.

### First secondary harmonic

Equivalent linearisation retains only the fundamental of the quadratic drag.
Expanding `Q|Q|` for `Q_e = |Q_e| cos(theta_e)` gives odd harmonics only,

```
Q|Q| = |Q|^2 * sum_{n odd} 8 (-1)^((n-1)/2) cos(n theta) / (pi n (4 - n^2))
```

The `n = 1` term is the Lorentz resistance. Once the fundamental has converged the
`n = 3` term is a *known* source, so the `3 omega` problem is **linear**: the same
nodal system, reassembled at `3 omega`, with

```
R_e^(3) = (3/2) * beta_e * |Q_e|          (from <2|Q| |cos|> = (4/pi)|Q|)
S_e     = (1/5) * beta_e * Q_e^3 / |Q_e|  (branch e.m.f. at 3 omega)
```

and `eta3 = 0` at the external boundaries, the forcing being monochromatic. This
costs one extra linear solve per point — a few percent of the run time — and is
always computed.

For a single basin with a single inlet the result reduces to the closed-form
third-harmonic amplitude of the accompanying manuscript, to machine precision.

### When the third harmonic can be trusted

The drag linearisation is calibrated on the *fundamental* discharge, so the
perturbation is self-consistent only while the `3 omega` discharge stays small
compared with it. Both programs therefore report the per-channel ratio
`|Q3|/|Q1|` and warn on it:

| worst `|Q3|/|Q1|` | meaning |
|---|---|
| below 0.20 | normal; the `3 omega` amplitude is accurate to roughly 10% |
| 0.20 – 0.50 | `[NOTE]` — `3 omega` is approaching a network resonance, error may exceed 10% |
| above 0.50 | `[WARNING]` — `3 omega` is at a resonance, the perturbation has broken down |

In the last case the third harmonic can be amplified beyond the fundamental,
which happens when a natural frequency of the network sits near `3 omega`. **The
fundamental is unaffected in all three cases** — the warning concerns the
harmonic correction only.

For a single basin the elevation ratio `|eta3|/|eta1|` never exceeds
`2/45 ~ 0.0444`, a bound that is saturated but not crossed; for a network of two
or more basins it can be exceeded, so it is not used as the diagnostic.

---

## Contents

| File | Purpose |
|---|---|
| `multi_inlet_response.py` | One basin, `N` inlets **in parallel**. Sweeps `omega` and `a_s`, produces `G(omega)` curves. |
| `basin_network_response.py` | General network of `M` basins and `N` channels, arbitrary topology, multiple open boundaries. |
| `time_domain_check.py` | Independent validation: integrates the same network in the time domain and compares. |
| `example_multi_inlet.txt` | Example input for the parallel-inlet program (17 inlets). |
| `example_network.txt` | Example input for the network program (3 basins, 7 channels). |

---

## Requirements

- Python >= 3.8
- `numpy`
- `matplotlib`

```bash
pip install numpy matplotlib
```

---

## 1. Parallel inlets — `multi_inlet_response.py`

A single basin of surface area `S` connected to the sea by `N` inlets in parallel.

### Input format

Whitespace-separated values, one keyword per line; `#` starts a comment.

```
N   = 2
L   = 3200 4500          # inlet lengths [m]
D   = 1.2 3.4            # inlet depths [m]
B   = 150 255.4          # inlet widths [m]
kse = 30 30              # Strickler coefficients [m^(1/3)/s]
S   = 10000000           # basin surface area [m^2]
```

`L`, `D`, `B` and `kse` must each contain exactly `N` positive values.

### One operating point

Give a frequency (`--period` in hours, or `--omega` in rad/s) and an amplitude,
and the program prints the answer instead of running a sweep:

```bash
python multi_inlet_response.py example_multi_inlet.txt --period 12.42 --amplitude 0.6
```

```
Forcing    a_s = 0.6 m     omega = 1.4053e-04 rad/s     T = 12.420 h

  gain             G = 0.4070
  basin amplitude    = 0.2442 m
  basin tidal range  = 0.4884 m
  phase lag          = 73.01 deg  (2.519 h)
```

Add `--inlets` for the discharge carried by each inlet and its share of the total.

### Sweep

With no frequency given, the program sweeps `omega` for a set of amplitudes:

```bash
python multi_inlet_response.py example_multi_inlet.txt
python multi_inlet_response.py example_multi_inlet.txt --amplitudes 0.2 1 --nfreq 800
```

Output, next to the input file: `<stem>_gain.png` and `<stem>_gain.csv`
(`omega_rad_s` plus one gain column per amplitude).

---

## 2. Networks — `basin_network_response.py`

An arbitrary graph of basins (dynamic nodes with storage) and external boundaries
(nodes with a prescribed harmonic level), connected by channels. Channels may link
a boundary to a basin or two basins to each other; boundary-to-boundary channels
are rejected, since they do not affect the basin dynamics.

### Input format

Comma-separated records; `#` starts a comment. `N_bays` and `N_channels`, if
given, are checked against the number of records actually declared.

```
N_bays     = 3
N_channels = 7

# bay , name , S [m2]
bay , B1 , 10000000
bay , B2 , 25000000
bay , B3 ,  8000000

# boundary , name [, relative_amplitude , phase_deg]
boundary , sea                 # equivalent to: boundary , sea , 1.0 , 0.0

# channel , name , node1 , node2 , L [m] , D [m] , B [m] , kse
channel , C1 , sea , B1 , 3200 , 1.2 , 150.0 , 30
channel , C5 , B1  , B2 , 5000 , 2.0 , 180.0 , 30
```

Several boundaries can be declared, each with its own amplitude factor and phase:
the imposed elevation is `a_s * factor * exp(i * phase)`. This allows, for
instance, a strait open at both ends with a phase lag between the two seas.

### One operating point

```bash
python basin_network_response.py example_network.txt --period 12.42 --amplitude 0.5 --channels
```

```
Forcing    a_s = 0.5 m     omega = 1.4053e-04 rad/s     T = 12.420 h

  basin    gain   amplitude [m]   range [m]   lag [deg]   lag [h]
  B1      0.7237         0.3619      0.7237       37.85     1.306
  B2      0.4065         0.2033      0.4065       73.34     2.530
  B3      0.3429         0.1715      0.3429       82.72     2.854

  channel  from -> to    |Q| [m3/s]   lag [deg]
  C1       sea -> B1           64.3      -42.98
  C2       sea -> B1          518.2      -38.86
  ...
```

`--channels` is optional; without it only the basin table is printed.

### Sweep

```bash
python basin_network_response.py example_network.txt
python basin_network_response.py example_network.txt --amplitudes 0.2 1 --nfreq 800
```

Output: `<stem>_gain.csv` (`a_s_m`, `omega_rad_s`, `period_h`, then `G_<basin>`
and `phase_<basin>_deg` for every basin) and one `<stem>_gain_<basin>.png` per basin.

### Options

| Option | Meaning |
|---|---|
| `--omega` / `--period` | single operating point: frequency [rad/s] or period [h] |
| `--amplitude` | single operating point: sea amplitude `a_s` [m] |
| `--channels` | also list the discharge in each channel |
| `--amplitudes` | sweep: sea amplitudes to compare [m] |
| `--omega-min` / `--omega-max` / `--nfreq` | sweep: frequency range and resolution |
| `--npz` | also save the full complex fields to a `.npz` archive |
| `--show` | display the plots in addition to saving them |
| `--tol` / `--itmax` / `--relax` | fixed-point controls |

---

## 3. Validation — `time_domain_check.py`

The harmonic solution can be checked against a direct time-domain integration of
the same network. `time_domain_check.py` integrates

```
L_e dQ_e/dt + kappa_e Q_e|Q_e| = eta_node1 - eta_node2
S_m d eta_m/dt                 = net inflow into basin m
```

with RK4, keeping the quadratic drag as it is and using **the same friction
closure** as the harmonic solver, so the comparison isolates the error of the
Lorentz linearisation instead of mixing in a different closure. It starts from
rest, runs for a number of forcing periods, and Fourier-analyses the last few.

```bash
python time_domain_check.py example_network.txt --period 12.42 --amplitude 0.5
```

```
  basin     quantity        time domain       model   difference
  B1        A1 / a_s           0.731718    0.723705       -1.10 %
  B1        A3 / a_s           0.045205    0.042776       -5.37 %
  B1        A5 / a_s           0.010296           -    truncated
  ...
  worst |Q3|/|Q1| = 0.139
```

Over a two-basin and a three-basin network, periods from 2 to 87 h and amplitudes
from 0.1 to 2 m, the fundamental amplitude comes out accurate to about 1% (worst
3%) and the third harmonic to about 3% (worst 12%). The bias is systematic and
negative: the model slightly underestimates the third harmonic, because the fifth
and higher harmonics — of order 1% of the forcing amplitude — are truncated and
their feedback is lost. The error grows with distance from the open boundary:
basins facing the sea come out within a few percent, strongly choked inner basins
at the bottom of the range.

Use `--periods` and `--steps-per-period` to check that the reference itself has
converged; the defaults (20 and 1000) are enough for the supplied examples.

---

## Numerical notes

- **Continuation in frequency.** In the network solver the converged resistances at
  one frequency seed the next one, which makes the sweep cheaper and more robust.
- **Relaxation.** The fixed point is damped (`relax = 0.5` by default). Strongly
  frictional or strongly resonant configurations may need a smaller value.
- **Convergence.** Non-converged points are flagged in the CSV (`converged = 0`)
  and summarised on screen; results there should not be trusted blindly.
- **Cost.** A single operating point is instantaneous; a sweep of 400 frequencies by
  7 amplitudes on a three-basin network runs in a few seconds, and the
  third-harmonic correction adds roughly 4% to that.
- The initial guess for the parallel-inlet program comes from a closed-form estimate
  of the gain and of the discharge shares, which is why it converges in a few tens
  of iterations.

---

## References

Kondo, H. (1975). *Depth of Maximum Velocity and Minimum Flow Area of Tidal Entrances*.
Coastal Engineering in Japan, 18(1), 167–183.

Lorentz, H. A. (1922). *Het in rekening brengen van den weerstand bij schommelende
vloeistofbewegingen*. De Ingenieur, 37(36), 695–696.
