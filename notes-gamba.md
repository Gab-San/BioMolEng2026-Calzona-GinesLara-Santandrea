Aprire la console:

Plugins->TkConsole

Cambiare font:

`tkcon font 'Courier' 18`

Create a set of atoms:

`set a [atomselect top protein]`

save a into a new pdb file:

`$a writepdb ubqp.pdb`


IMPR section ????

IC section: Internal Coordinates
Distance between 1-2
Angle between 1-2-3
???
Angle between 1-3-4
Distance between 3-4

It is used to guess distance of the missing atom.

Thanks to the topology file, it is possible to assign/generate the 
information about the missing atoms.

## PSF

The psf is storing the complementary information.
In the first section is describing atom by atom describing mass and +

## Opening the final molecule

`vmd <>.pdb <>.psf`

Add the water solution.
`package require solvate`
`solvate ubq.psf ubq.pdb -t 5 -o ubq_wb`

-t adds water in a range of 5 armstrong as a padding to our box.

10 armstrong is a typical cutoff value.

It will generate a new .pdb and a new .psf containing water.

> Ha fatto na cosa incredibile

Modeling -> Add Ions

Adding ions is not deterministic!


## BOH

`set all [atomselect top all]`

measure the size of the system -> `measure minmax $all` ->
`{10.442000389099121 7.947000026702881 -5.327000141143799} {51.11899948120117 50.23099899291992 40.895999908447266}`: minimum and maximum {X,Y,Z} coordinates.

We are interested in the size so max - min = `{41, 42, 46}`.

Calculate the center of the box -> `measure center $all` -> `30.853321075439453 29.084442138671875 17.79519271850586`


`cp ../common/sample.conf .` -> adjust the configuration file
`cp ../common/par_all27_prot_lipid.inp` -> adjust copy parameter file.

The parameter file contains for each bond:
- elastic constant and the reference length;
- elastic constant and angle;
- etc...

## Run the simulation

`namd3 +p2 <>.conf`

## Simulation output

Various

.dcd files are collection of frames of the simulation.

Open with vmd -> `vmd <>.dcd <>.pfs`

In the various tabs alongside graphics one can check:
Periodic and select axes.


Analysis -> RMSD Trajectory Tool -> ALIGN -> Plot -> RMSD.


Why ALIGN? BOH?


Analysis -> Hydrogen Bond: Analyse the energy of the hydrogen bond.

Provide a selection:
if only 1 it will consider bonds between 1-1
if provided 1 and 2 it will consider 1-2

By enabling All hbonds all hydrogen bonds will be evaluated  in detail.

Hbonds are calculated on a geometrical basis. These kind of bonds are empirically defined by geometrical parameters such as:
- distance between acceptor and donor;
- angle between -;
