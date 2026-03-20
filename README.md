# spring_grid

A real-time 2D spring-mass lattice simulation written in C17 with SDL2.

An N×M grid of identical point masses is connected by Hookean springs.
Perturb the grid and watch transverse waves, standing modes, and shock fronts
propagate through the lattice.  Toggle shear (diagonal) springs on and off
mid-run to feel the difference in rigidity.

![spring grid screenshot placeholder](screenshot.png)

---

## Features

- Configurable grid size, mass, spring constant, and damping via `#define`
- Cardinal springs (horizontal + vertical) with optional diagonal shear springs
- Velocity Verlet integration — symplectic, stable for oscillatory systems
- 16 physics sub-steps per rendered frame for accuracy without visual lag
- Springs colour-coded by strain (tension → red, compression → blue)
- Masses colour-coded by speed (cold → warm)
- Interactive: grab and drag any mass with the mouse

---

## Dependencies

| Platform | Command |
|---|---|
| macOS (Homebrew) | `brew install sdl2` |
| Debian / Ubuntu | `apt install libsdl2-dev` |

A C17-capable compiler (`gcc` or `clang`) and `make` are also required.

---

## Build

```sh
make          # optimised release build  →  ./spring_grid
make debug    # debug build with ASan + UBSan  →  ./spring_grid_debug
make clean    # remove built binaries
```

---

## Controls

| Input | Action |
|---|---|
| Left-click + drag | Grab and pull a mass |
| Right-click | Random impulse on all masses |
| Space | Pluck the centre mass |
| D | Toggle diagonal springs on / off |
| R | Reset grid to equilibrium |
| + / − | Increase / decrease damping |
| G | Chaotic global impulse |
| Q / Escape | Quit |

---

## Configuration

All physical and display parameters are `#define` constants at the top of
`spring_grid.c`.  Change them and recompile.

| Constant | Default | Meaning |
|---|---|---|
| `COLS` | 18 | Number of columns |
| `ROWS` | 14 | Number of rows |
| `MASS` | 1.0 | Mass of each node (kg) |
| `SPRING_K` | 400.0 | Spring constant (N/m) |
| `DAMPING` | 0.6 | Velocity damping coefficient (kg/s) |
| `DT` | 0.0005 | Integration timestep (s) |
| `STEPS_PER_FRAME` | 16 | Physics sub-steps per rendered frame |
| `SPACING` | 60 | Equilibrium grid spacing (px) |
| `MASS_RADIUS` | 7 | Rendered circle radius (px) |

The natural frequency of a single mass on two cardinal springs is
ω = √(2k/m), so with the defaults that is about 28 rad/s.  Reducing `DT`
or increasing `STEPS_PER_FRAME` is advisable if you raise `SPRING_K`
significantly, to keep the integrator stable (Verlet requires ω·Δt ≪ 2).

---

## Physics notes

**Integration.** Velocity Verlet is used in preference to simple Euler
because it is time-reversible and symplectic: energy drifts are bounded
rather than growing monotonically, which matters for long runs.

**Damping.** A linear drag force −γ**v** is applied to each mass at every
step.  This is not physical viscosity but it keeps the simulation from
accumulating floating-point energy indefinitely.  Set `DAMPING` to 0 for
a conservative (energy-preserving) system.

**Diagonal springs.** With only cardinal springs the lattice has zero
resistance to shear — you can deform it into a parallelogram with no
restoring force.  Diagonal springs fix this, making the sheet behave more
like a stiff membrane.  Toggle with **D** to hear and see the difference.

---

## File layout

```
spring_grid.c   — entire simulation (single translation unit)
Makefile        — release and debug targets
README.md       — this file
```
