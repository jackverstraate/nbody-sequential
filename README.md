# N-Body Simulation (C++)

This is a simple N-body simulator that computes gravitational forces between particles and updates their motion over time.

It follows the assignment’s requirements:
- arrays/vectors for mass, position, velocity, force
- O(N^2) force calculation with a softening term
- update velocity and position each time step
- output the simulation state regularly in a TSV format that matches `solar.tsv`

## Build

```bash
make
```

## Run

```bash
./nbody <n_or_file_or_preset> <dt> <steps> <dumpEvery> [output.tsv]
```

- If the first argument is an **integer N**, it creates a random system with N particles.
- If it is a **file path**, it loads the initial state from the **first line** of that TSV file.
- If it is `sem`, it uses a small **Sun/Earth/Moon** preset.

### Examples

Random 100 particles, dt=1, 10000 steps, dump every 10 steps:

```bash
./nbody 100 1 10000 10 out.tsv
```

Load from the provided `solar.tsv` (first line), dt=200, 5000 steps:

```bash
./nbody solar.tsv 200 5000 10 solar_out.tsv
```

Sun/Earth/Moon preset:

```bash
./nbody sem 200 5000 1 sem.tsv
```

## Output format (TSV)

Each output line is one state:

```
N    m x y z vx vy vz fx fy fz   m x y z vx vy vz fx fy fz   ...
```

- First value is the number of particles **N**
- Then for each particle: **mass**, position (**x,y,z**), velocity (**vx,vy,vz**), force (**fx,fy,fz**)
- Values are separated by **tabs**

This should work with the provided plotting script:

```bash
python3 plot.py out.tsv out.pdf 10000
```
