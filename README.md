# Engineering a thermostable mutation of Cofilin (1-COF)

The task of the Biomolecular Modeling's Project (A.A. 2026-2027) is to engineer a thermostable mutation of enzyme Cofilin.

Starting from wild type analysis and following the procedure presented in the reference paper[^1] we are tasked to research mutations
that supposedly would make Cofilin protein more thermically stable.

---

<!--toc:start-->
- [Engineering a thermostable mutation of Cofilin (1-COF)](#engineering-a-thermostable-mutation-of-cofilin-1-cof)
  - [TOOLS](#tools)
  - [GUIDE](#guide)
    - [Auxiliary Code](#auxiliary-code)
      - [Export Principal Protein Atoms](#export-principal-protein-atoms)
      - [How To Add Water](#how-to-add-water)
      - [How To Add Ions](#how-to-add-ions)
  - [Notes from the tutorial](#notes-from-the-tutorial)
    - [Electrostatics with MSM](#electrostatics-with-msm)
    - [Electrostatics with PME](#electrostatics-with-pme)
<!--toc:end-->

---


## TOOLS

The tools we used are:

- [VMD](https://www.ks.uiuc.edu/Research/vmd/)
- [NAMD](https://www.ks.uiuc.edu/Research/namd/)

## GUIDE

We present here a step by step guide on what we did.

1. Download the [1-COF](https://www.rcsb.org/structure/1COF) `.pdb` file from the Protein Data Bank;
2. Import `top_all27_prot_lipid.inp` that describes the CHARMM22 force field from a given library;
3. Export the protein atoms and remove noise atoms using the lines in [here](#export-principal-protein-atoms);
4. Write the `.pgn` file to generate the protein structure with the addition of the missing atoms;
5. Run the `.pgn` file by launching `vmd -dispdev text -e cof.pgn` 
6. [Add water to the obtained protein](#how-to-add-water)
7. [Add ions to the obtained water solution](#how-to-add-ions)


### Auxiliary Code
#### Export Principal Protein Atoms

In vmd Plugins --> TkConsole

```
tkcon font 'Courier' 18 # Better readability
set a [atomselect top protein]
$a writepdb cofp.pdb
```

#### How To Add Water

Open the desired protein in vmd: `vmd cof.pdb cof.psf`

In vmd Plugins --> TkConsole

```
package require solvate
solvate cof.psf cof.pdb -t 5 -o cof_wb
```

-t <dist> is the box padding in Armstrong

#### How To Add Ions

In vmd Plugins --> TkConsole

```
package require autoionize
autoionize -psf cof_wb.psf -pdb cof_wb.pdb -sc 0.1 -o cof_ionized
```


## Notes from the tutorial

### Electrostatics with MSM

Multilevel Summation Method (MSM) is a useful method for dealing
with electrostatic interactions in a system that does not have fully pe-
riodic boundary conditions along all three dimensions, in other words,
when employing non-periodic or semi-periodic boundaries. 
Although
MSM also supports fully periodic boundary conditions, the Particle
Mesh Ewald (PME) method discussed in Section 1.5 is more efficient
for smaller periodic systems, such as ubiquitin simulated in a water
box. In the present simulation involving a water sphere, MSM is the
best alternative for electrostatics.

MSM uses particle meshes, similar to PME, arranged as 3-D grids for
evaluating long-range electrostatic interactions. 
Charges are assigned
to nearby grid points, and the interaction of these gridded charges
result in gridded potentials that are interpolated back to the atoms
to determine the long-range atomic forces. 
Unlike PME which uses
a single grid resolution, MSM uses a nesting of grids having coarser
resolution at each successive “level” in real space for resolving the
long-ranged tail of the 1/r Coulomb potential. 
Rather than evaluating the Ewald sum as PME does, MSM uses a continued splitting of the interaction potential and performs interpolation from these
nested grid levels.

- MSM: indicates whether or not the simulation uses the Multilevel
Summation Method; uses values yes and no.

- MSMGridSpacing: The grid spacing between the finest level grid
points determines, in part, the accuracy and efficiency of MSM.
According to an error versus cost analysis, the optimal value is for
the grid spacing to be a bit larger than the average interatomic
spacing; since the performance is extremely sensitive to the grid
spacing, it is best to leave this parameter set to its default value
of 2.5Å.

- MSMxmin, MSMxmax, MSMymin, MSMymax, MSMzmin, MSMzmax: For
non-periodic simulation, MSM needs to know the expected max-
imum and minimum coordinates so as to establish grids that will
contain all atoms throughout the entire simulation. To keep the
atoms from leaving the simulated space, it might be necessary
to use a restraining potential, such as the spherical boundary
restraints used for the water sphere simulation.

Additional MSM parameters that further control accuracy and per-
formance are documented in the [NAMD User’s Guide](https://www.ks.uiuc.edu/Research/namd/3.0.2/ug/).


### Electrostatics with PME

Particle Mesh Ewald (PME) is a useful method for dealing with elec-
trostatic interactions in a system when periodic boundary conditions
are present. The Ewald sum is an efficient way of calculating long
range forces in a periodic system. The particle mesh is a 3-D grid
created in the system over which the system charge is distributed.
From this charge, potentials and forces on atoms in the system are
determined. As a result, your grid size should be chosen such that
it is fine enough to accurately represent the configuration of your
system.

- PME: indicates whether or not the simulation uses the Particle
Mesh Ewald Sum method; uses values yes and no.

- PMEGridSpacing: sets the minimum ratio between the number
of PME grid points along each cellBasisVector and the phys-
ical dimensions. 
Since the grid will replicate the charge distri-
bution in your system, PMEGridSpacing should be chosen to be
large enough so that the grid spacing accurately reproduces your
charge distribution. However, it should not be so large that it
will slow down your simulation to provide unnecessary precision.
Typically, a grid density of slightly more than 1/Å is a good
choice to reproduce charge distribution in biological systems,
where the closest atoms have a bond separation on the order
of 1 Å. 
This corresponds to a PMEGridSpacing of 1.0. NAMD
will then automatically set the PME grid sizes (see below) such
that there is always less than 1.0 Å between grid points and the
sizes contain only small prime factors (e.g. 2, 3, and 5).


Alternatively, one can define the PME grid sizes manually, using
PMEGridSizeX, PMEGridSizeY, and PMEGridSizeZ. 
These set the size
of the PME grid along cellBasisVector1, 2 and 3, respectively
(not the x, y, and z directions as implied). 
For speed in computing
Fast Fourier Transforms, PMEGridSizeX should be chosen so that it
can be factorized by 2, 3, or 5. 
If your cellBasisVector1 = (60,
0, 0), a good choice for PMEGridSizeX might be 64, since 60 Å /
64 = 0.9375 Å and 64 = 26. 
Note that since cellBasisVector is
defined with slightly different values in each direction, the size of the
mesh spacing (in length) will be different in each direction.
Note also that when using the PME method, the command cutoff
dictates the separation between long and short range forces for the
method; it does not simply turn off interactions.

[^1]: <!--TODO: add doi--->
