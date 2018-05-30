# spg-siesta

#========================================================================================

## [Pt]

Os arquivos `siesta2bt_espresso.f90`, `siesta2bt_spg_dataset.f90` e `siesta2bt_spg_findsym.f90` possuem a rotina siesta2bt que busca os pontos k pelo método utilizado pelo quantum espresso e os dois últimos por metodos diferentes da biblioteca spglib, sendo que o primeiro faz uso da estrutra
spgdataset e segunda nao.

### Instalação

i) Crie uma pasta no diretório `path/to/siesta/Src/` e copie os arquivos que estão contidos neste repositório, ou, no diretório `Src/`, digite:

git clone https://github.com/lucas-campagna/spg-siesta.git

ii) Copie o arquivo arch-spg.make para o diretorio `Obj/` (ou o que está sendo usado para compilar o siesta).

iii) Inclua o `arch-spg.make` no arquivo `Makefile` do siesta, preste atenção para incluir após a primeira receita, para não correr este risco nós recomendamos que a inclusão seja feita junta ou após a inclusao do arquivo arch.make.

Obs: existem três métodos que podem ser escolhidos como dito anteriormente, para escolher algum copier seu respectivo arquivo para o novo diretório criado no passo i) com nome `siesta2bt.f90`. O método padrão utilizado  aquele que faz o uso das subrotinas `spg_find_primitive` e `spg_get_symmetry` (contido no arquivo `siesta2bt_spg_findsym.f90`).

### Forma de uso

Esta extensão adiciona os seguintes parâmetros os arquivos FDF:

#### ```%block BT.kgrid_Monkhorst_Pack```

Especifíca a grade no espaço recíproco para o cálculo de termoeletricidade com BoltzTraP.

```
%block BT.kgrid_Monkhorst_Pack
  kx   0   0   0.0
   0  ky   0   0.0
   0   0  kz   0.0
%endblock BT.kgrid_Monkhorst_Pack
```

Onde `kx`, `ky` e `kz` são inteiros.

#### ```BT.Symprec```

Tolerância usada nas operaçes de simetria. Valor pardão 0.5.

#### ```BT.First_Band```, ```BT.Last_Band```

Primeiro e último índices (primeira banda indice 1) do conjunto de bandas selecionados. Valores padrões são 1 e no_u (quantidade total de bandas).

#========================================================================================

## [En]

The files `siesta2bt_espresso.f90`, `siesta2bt_spg_dataset.f90` and `siesta2bt_spg_findsym.f90` use the subroutine siesta2bt to search irreducible k-points, each by a method. The first uses the same kind of method used by quantum espresso, the others two use methods provided by spglib, the first one makes use spgdataset structure unlinke the second one that makes use of `spg_find_primitive` and `spg_get_symmetry` subroutines.

## Install

i) Create a folder in the directry `path/to/siesta/Src/` and copie the files that are in this repository, or just enter into your `Src/` directory ant type:

git clone https://github.com/lucas-campagna/spg-siesta.git

ii) Copie the arch-spg.make to your `Obj/` directory (or the used to compile siesta).

iii) Include the `arch-spg.make` file into `Makefile` of siesta, but pay attention to include this after the first recipe, to miss this error we highly recommend you to include this together or after the `arch.make` inclusion.

Obs: there are three methods used here as said before, to choose one copie its respective file to your new directory created in step i) with the name `siesta2bt.f90`. The default method is the one that does use of `spg_find_primitive` and `spg_get_symmetry` subroutines (in `siesta2bt_spg_findsym.f90` file).

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
