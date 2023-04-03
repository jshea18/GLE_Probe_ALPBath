# LAMMPS Extensions

## Create atoms
This function is a modification of the LAMMPS ```create_atoms``` function which adds the possibility of creating a "single/sphere" object. The "single/sphere" object adds ```N``` atoms to locations on a sphere of radius ```R``` around the point ```x y z```. 

#### Relevant LAMMPS files:
- create_atoms.cpp
- create_atoms.h

#### Sample execution code:
```
create_atoms 2 single/sphere x y z R N units box
```

## Move region
#### Relevant LAMMPS files:
- fix_movereg.cpp
- fix_movereg.h

#### Sample execution code:
```
fix 8 col movereg 1 ${R} solvent file bath_bubble.dat
```

## CDF
#### Relevant LAMMPS files:
- fix_cdf.cpp
- fix_cdf.h

#### Sample execution code:
```
fix 8 col cdf 1000 200 20 200 20 solvent file data_radial_function.dat
```

## Add Activity
This function adds an active propulsion force to particles so that they behave as an active Langevin particle (ALP) without rotational inertia, e.g. an active Brownian particle (ABP) with translational inertia. The equation of motion for such a particle is:

$$m \dot{\mathbf{v}}_n(t) = F_0\mathbf{e}_n(t)-\gamma \mathbf{v_n}(t)+\boldsymbol{\xi}_n(t)$$

where $\mathbf{v}_n(t)$ is the velocity of the particle $n$ at time $t$, $F_0$ is the active force, $\gamma$ is the damping constant, \mathbf{e}_n(t)

#### Relevant LAMMPS files:
- fix_addactivity.cpp
- fix_addactivity.h

#### Sample execution code:
```
fix 4 solvent addactivity style ABP ${Fact} ${D_solvent} RAN2
```

# Analysis code

All analysis code is in the ```correlation.c``` file. This analysis code has many functions including: calculating correlation functions; calculating the memory kernel through different methods such as the inverse Volterra method; calculating the stochastic force

## To compile:
```
g++ correlation.c -o Correlation -O2
```

## To run:

#### Velocity autocorrelation function, force autocorrelation function, and velocity-force correlation function:
Over vectors: 
```
./Correlation colloid.dat 9 90000001 0.001 5000000 10000 1 3 100.0 1 -n -c 3 3 3 vv.dat -c 3 6 3 vf.dat -c 6 6 3 ff.dat
```

Over x:
```
./Correlation colloid.dat 9 90000001 0.001 500000 10000 1 3 100.0 1 -n -c 3 3 1 vv_x.dat -c 3 6 1 vf_x.dat -c 6 6 1 ff_x.dat
```

Over y:
```
./Correlation colloid.dat 9 90000001 0.001 500000 10000 1 3 100.0 1 -n -c 4 4 1 vv_y.dat -c 4 7 1 vf_y.dat -c 7 7 1 ff_y.dat
```

Over z:
```
./Correlation colloid.dat 9 90000001 0.001 500000 10000 1 3 100.0 1 -n -c 5 5 1 vv_z.dat -c 5 8 1 vf_z.dat -c 8 8 1 ff_z.dat
```

#### Memory and stochastic forces:

**Only do this after you have gotten the correlation functions!**

Over x:
```
./Correlation colloid.dat 9 90000001 0.001 500000 10000 1 3 100.0 1 -v vv_x.dat vf_x.dat ff_x.dat volterra_mem_x.dat -s sfx.dat 0
```

Over y:
```
./Correlation colloid.dat 9 90000001 0.001 500000 10000 1 3 100.0 1 -v vv_y.dat vf_y.dat ff_y.dat volterra_mem_y.dat -s sfy.dat 1
```

Over z:
```
./Correlation colloid.dat 9 90000001 0.001 500000 10000 1 3 100.0 1 -v vv_z.dat vf_z.dat ff_z.dat volterra_mem_z.dat -s sfz.dat 2
```

## Meaning of commands for this code:

```colloid.dat```: Inputfile

```9```: Number of columns in File (currently fixed to 9 or 18! (3 dimensions and space, velocity and force; up to two particles))

```90000001```: Length of file

```0.001```: Time step

```500000```: Timestep at which equilibrium is reached

```10000```: Length of the correlation function (how many time steps)

```1```: Coarse-grain Factor (for you always 1)

```3```: Number of dimensions

```100.0```: Mass of particle (mass of two particles must be the same)

```1```: Number of particles (up to 2)

```-n```: Flag for when you do not start in equilibrium

```-c A B C FILENAME```: To calculate correlation functions

```A```: Starting column for observable 1 (Column 1 is 0)

```B```: Starting column for observable 2

```C```: Dimension (number of columns you want to use to calculate the correlation functions)

