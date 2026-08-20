# Tidal Gain of Bay–Channel Systems 

Two small, dependency-light Python programs that compute the **tidal gain**

$$G = \frac{a_b}{a_s}$$

of one or more bays connected to the open sea through frictional channels, using
a lumped **resistance–inertance (R–L) network** with a Lorentz-linearised,
amplitude-dependent friction term that is solved **self-consistently** .

The model generalises the classical single-inlet treatment of
Kondo (1975) to (i) an arbitrary number of parallel inlets and
(ii) an arbitrary network of interconnected bays and channels.

---

## Physical model

Each channel `e` is represented as a series impedance

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

Each bay `m` acts as a storage element (compliance) of surface area `S_m`, so that
continuity at the bay reads

```
i * omega * S_m * eta_m = sum_e Y_e * (eta_neighbour - eta_m),      Y_e = 1 / Z_e
```

For a rectangular cross-section the geometry is derived from the input as

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

until the resistances and the bay elevations stop changing. This is the
"fully self-consistent"; the amplitude of the forcing therefore matters,
and the gain curves are **not** the same for different `a_s`.

---

## Contents

| File | Purpose |
|---|---|
| `gain_method4_frequency_amplitude.py` | Single bay, `N` channels **in parallel**. Sweeps `omega` and `a_s`, produces `G(omega)` curves. |
| `network_method4.py` | General network of `M` bays and `N` channels with arbitrary topology and multiple external boundaries. |
| `esempio_canali_method4.txt` | Example input for the parallel-channel program (17 channels). |
| `network_example.txt` | Example input for the network program (3 bays, 7 channels). |

---

## Requirements

- Python ≥ 3.8
- `numpy`
- `matplotlib`

```bash
pip install numpy matplotlib
```

---

## 1. Parallel channels — `gain_method4_frequency_amplitude.py`

A single bay of surface area `S` is connected to the sea by `N` channels in parallel.
The program computes `G(omega)` for a set of forcing amplitudes and produces a
frequency-response plot with the 24 h and M2 tidal frequencies marked.

### Input format

Whitespace-separated values, one keyword per line; `#` starts a comment.

```
N   = 2
L   = 3200 4500          # channel lengths [m]
D   = 1.2 3.4            # channel depths [m]
B   = 150 255.4          # channel widths [m]
kse = 30 30              # Strickler coefficients [m^(1/3)/s]
S   = 10000000           # bay surface area [m^2]
```

`L`, `D`, `B` and `kse` must each contain exactly `N` positive values.
An `a_s` entry may be present but is **ignored**: the forcing amplitude is swept
by the program.

### Usage

```bash
python gain_method4_frequency_amplitude.py esempio_canali_method4.txt --show
```

The swept ranges are set at the top of the file
(`OMEGA_MIN`, `OMEGA_MAX`, `N_OMEGA`, `AMPLITUDES`).

### Output

Written next to the input file, using its stem as prefix:

- `<stem>_gain_method4_frequency_amplitude.png` / `.pdf` — `G(omega)`, one curve per amplitude
- `<stem>_gain_method4_frequency_amplitude.csv` — columns `omega_rad_s`, `G_as_<a>m`, ...

---

## 2. Networks — `network_method4.py`

An arbitrary graph of bays (dynamic nodes with storage) and external boundaries
(nodes with a prescribed harmonic level), connected by channels. Channels may link
a boundary to a bay or two bays to each other; boundary-to-boundary channels are
rejected because they do not affect the bay dynamics.

### Input format

Comma-separated records; `#` starts a comment. `N_bays` and `N_channels`, if given,
are checked against the number of records actually declared.

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
the imposed elevation is `a_s * factor * exp(i * phase)`. This allows, for instance,
a strait open at both ends with a phase lag between the two seas.

### Usage

```bash
python network_method4.py network_example.txt --show

python network_method4.py network_example.txt \
    --omega-min 1e-5 --omega-max 1e-3 --nfreq 500 \
    --amplitudes 0.05 0.1 0.2 0.5 1 2 5 --show
```

| Option | Default | Meaning |
|---|---|---|
| `--omega-min` | `1e-5` | minimum angular frequency [rad/s] |
| `--omega-max` | `1e-3` | maximum angular frequency [rad/s] |
| `--nfreq` | `400` | number of log-spaced frequencies |
| `--amplitudes` | `0.05 … 5` | incoming amplitudes `a_s` [m] |
| `--tol` | `1e-10` | fixed-point tolerance |
| `--itmax` | `5000` | maximum number of iterations |
| `--relax` | `0.5` | relaxation factor, `0 < relax ≤ 1` |
| `--show` | off | display the plots in addition to saving them |

### Output

- `<stem>_network_method4.csv` — long format: `a_s_m`, `omega_rad_s`, `period_h`,
  then `G_<bay>` and `phase_<bay>_deg` for every bay, plus `iterations` and `converged`
- `<stem>_network_method4.npz` — full complex fields: `eta` (bays), `Q` and `R` (channels),
  with shapes `(n_amplitudes, n_frequencies, n_nodes/n_channels)`
- `<stem>_gain_<bay>.png` — one frequency-response plot per bay

The console report lists, for every bay and amplitude, the peak gain and the
corresponding angular frequency and period.

---

## Numerical notes

- **Continuation in frequency.** In the network solver the converged resistances at
  one frequency are reused as the initial guess for the next one, which makes the
  sweep considerably cheaper and more robust.
- **Relaxation.** The fixed point is damped (`relax = 0.5` by default). Strongly
  frictional or strongly resonant configurations may need a smaller value.
- **Convergence.** Non-converged points are flagged in the CSV (`converged = 0`) and
  summarised on screen; results at those points should not be trusted blindly.
- The initial guess for the parallel-channel program comes from a closed-form
  ("Method 2") estimate of the gain and of the discharge shares, which is why it
  typically converges in a few tens of iterations.

---

## Reference

Kondo, H. (1975). *Depth of Maximum Velocity and Minimum Flow Area of Tidal Entrances*.
Coastal Engineering in Japan, 18(1), 167–183.
