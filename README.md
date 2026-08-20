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
costs **one extra linear solve per point** — a few percent of the run time — and
is enabled by default (`--no-third-harmonic` disables it).

For a single basin with a single inlet the result reduces to the closed-form
third-harmonic amplitude of the accompanying manuscript, to machine precision.

The perturbation is meaningful only while `|eta3| << |eta1|`. Points where the
ratio exceeds 0.10 are flagged `[PERTURBATION SUSPECT]` in the console report.
For reference, the single-basin analysis bounds the ratio by `2/45 ~ 0.0444`, and
that bound is drawn on the third-harmonic figures.

---

## Contents

| File | Purpose |
|---|---|
| `multi_inlet_response.py` | One basin, `N` inlets **in parallel**. Sweeps `omega` and `a_s`, produces `G(omega)` curves. |
| `basin_network_response.py` | General network of `M` basins and `N` channels, arbitrary topology, multiple open boundaries. |
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
The program computes `G(omega)` for a set of forcing amplitudes and marks the 24 h
and M2 tidal frequencies on the frequency-response plot.

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

`L`, `D`, `B` and `kse` must each contain exactly `N` positive values. An `a_s`
entry may be present but is **ignored**: the forcing amplitude is swept.

### Usage

```bash
python multi_inlet_response.py example_multi_inlet.txt --show
```

Swept ranges are set at the top of the file (`OMEGA_MIN`, `OMEGA_MAX`, `N_OMEGA`,
`AMPLITUDES`).

### Output

Written next to the input file, using its stem as prefix:

- `<stem>_response.png` / `.pdf` — `G(omega)`, one curve per amplitude
- `<stem>_response.csv` — `omega_rad_s`, then `G_as_<a>m` and
  `eta3_over_eta1_as_<a>m` for each amplitude
- `<stem>_third_harmonic.png` — relative amplitude of the `3 omega` harmonic

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

### Usage

```bash
python basin_network_response.py example_network.txt --show

python basin_network_response.py example_network.txt \
    --omega-min 1e-5 --omega-max 1e-3 --nfreq 500 \
    --amplitudes 0.05 0.1 0.2 0.5 1 2 5 --show
```

| Option | Default | Meaning |
|---|---|---|
| `--omega-min` | `1e-5` | minimum angular frequency [rad/s] |
| `--omega-max` | `1e-3` | maximum angular frequency [rad/s] |
| `--nfreq` | `400` | number of log-spaced frequencies |
| `--amplitudes` | `0.05 ... 5` | incoming amplitudes `a_s` [m] |
| `--tol` | `1e-10` | fixed-point tolerance |
| `--itmax` | `5000` | maximum number of iterations |
| `--relax` | `0.5` | relaxation factor, `0 < relax <= 1` |
| `--no-third-harmonic` | off | skip the `3 omega` correction |
| `--show` | off | display the plots in addition to saving them |

### Output

- `<stem>_response.csv` — long format: `a_s_m`, `omega_rad_s`, `period_h`, then
  `G_<basin>`, `phase_<basin>_deg`, `G3_<basin>`, `phase3_<basin>_deg` for every
  basin, plus `iterations` and `converged`
- `<stem>_response.npz` — full complex fields: `eta`, `eta3` (basins), `Q`, `Q3`
  and `R` (channels), with shapes `(n_amplitudes, n_frequencies, n_nodes/n_channels)`
- `<stem>_gain_<basin>.png` — one frequency-response plot per basin
- `<stem>_third_<basin>.png` — relative amplitude of the `3 omega` harmonic per basin

The console report lists, for every basin and amplitude, the peak gain, the
corresponding angular frequency and period, and the largest third-harmonic ratio
over the sweep.

---

## Numerical notes

- **Continuation in frequency.** In the network solver the converged resistances at
  one frequency seed the next one, which makes the sweep cheaper and more robust.
- **Relaxation.** The fixed point is damped (`relax = 0.5` by default). Strongly
  frictional or strongly resonant configurations may need a smaller value.
- **Convergence.** Non-converged points are flagged in the CSV (`converged = 0`)
  and summarised on screen; results there should not be trusted blindly.
- **Cost.** A sweep of 400 frequencies by 7 amplitudes on a three-basin network runs
  in a few seconds; the third-harmonic correction adds roughly 4% to that.
- The initial guess for the parallel-inlet program comes from a closed-form estimate
  of the gain and of the discharge shares, which is why it converges in a few tens
  of iterations.

---

## References

Kondo, H. (1975). *Depth of Maximum Velocity and Minimum Flow Area of Tidal Entrances*.
Coastal Engineering in Japan, 18(1), 167–183.

Lorentz, H. A. (1922). *Het in rekening brengen van den weerstand bij schommelende
vloeistofbewegingen*. De Ingenieur, 37(36), 695–696.
