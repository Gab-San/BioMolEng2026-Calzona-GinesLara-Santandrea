# Engineering a thermostable mutation of Cofilin (1-COF)

The task of the Biomolecular Modeling's Project (A.A. 2026-2027) is to engineer a thermostable mutation of enzyme Cofilin.

Starting from wild type analysis and following the procedure presented in the reference paper[^1] we are tasked to research mutations
that supposedly would make Cofilin protein more thermically stable.

---

<!--toc:start-->
- [Engineering a thermostable mutation of Cofilin (1-COF)](#engineering-a-thermostable-mutation-of-cofilin-1-cof)
  - [TOOLS](#tools)
  - [GUIDE](#guide)
<!--toc:end-->

---


## TOOLS

The tools we used are:

- [VMD](https://www.ks.uiuc.edu/Research/vmd/)
- [NAMD](https://www.ks.uiuc.edu/Research/namd/)

## GUIDE

We present here a step by step guide on what we did.

1. Download the [1-COF](https://www.rcsb.org/structure/1COF) `.pdb` file from the Protein Data Bank;
1. Import `top_all27_prot_lipid.inp` that describes the CHARMM22 force field from a given library;
1. Export the protein atoms and remove noise atoms using the lines in [here](#export-principal-protein-atoms);
1. Write the `.pgn` file to generate the protein structure with the addition of the missing atoms;
1. Run the `.pgn` file by launching `vmd -dispdev text -e cof.pgn` 
1. [Add water to the obtained protein](#how-to-add-water)
1. [Add ions to the obtained water solution](#how-to-add-ions)


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
autoionize -psf cof.psf -pdb cof.pdb -sc 0.1 -o cof_ionized
```


[^1]: <!--TODO: add doi--->
