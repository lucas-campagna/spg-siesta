# spg-siesta

#========================================================================================

## [En]

This implementation makes siesta able to generate input files for BoltzTraP, and study thermoelectricity proprietes.
It is implemented in version siesta-4.1-b3.

## Install

On `path/to/siesta/Obj` folder type:

```console
git clone https://github.com/lucas-campagna/spg-siesta.git
```

Or download it directly from https://github.com/lucas-campagna/spg-siesta clicking in button `Clone or download > Download ZIP`. Save it in `path/to/siesta/Obj` folder as before.

After this, run `./configure` script. It'll modify:

 - `Makefile`, to include the new spg.make in the correct position;
 - `siesta_analysis.F`, to include, in the correct position, the call to the subroutine that makes all the work.

This script will try to compile `spglib` using cmake acording to https://atztogo.github.io/spglib/install.html (accessed 20/09/2018). If it does'nt work for you, we are not responsible, consider the previous link.

## How to use

This extention add up the following pamarmeters to FDF file.

#### ``%block BT.kgrid_Monkhorst_Pack``

Specify the grid of reciprocal space for termoletricity calculations on BoltzTraP.

```
%block BT.kgrid_Monkhorst_Pack
  kx   0   0   0.0
   0  ky   0   0.0
   0   0  kz   0.0
%endblock BT.kgrid_Monkhorst_Pack

```
Where ```kx```, ```ky``` and ```kz``` are integers values.

#### ```BT.Symprec```

Tolerane of symmetry operations. Standard value is 0.5.

#### ```BT.First_Band```, ```BT.Last_Band```

First and last index (Fortran like) of the set of bands to do the calculations. Standard values are 1 and no_u (amount of bands).
