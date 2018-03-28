# spg-siesta

#========================================================================================

[Pt]

Os arquivos siesta2bt_espresso.f90, siesta2bt_spg_dataset.f90 e siesta2bt_spg_findsym.f90 possuem a rotina siesta2bt que busca os pontos k
pelo método utilizado pelo quantum espresso e os dois últimos por metodos diferentes da biblioteca spglib, sendo que o primeiro faz uso da estrutra
spgdataset e segunda nao.

- Instalação
Inclua o arch-spg.make no seu arquivo arch.make ou no seu arquivo Makefile, preste atenção para incluir após a primeira receita. Para não correr este risco nós recomendamos que a inclusão seja feita no arquivo arch.make. Em seguida, inclua libsymspg.a em COMP_LIBS e execute make.

#========================================================================================

[En]

The files siesta2bt_espresso.f90, siesta2bt_spg_dataset.f90 and siesta2bt_spg_findsym.f90 use the subroutine siesta2bt to search irreducible k-points, each by a method. The first uses the same kind of method used by quantum espresso, the others two use methods provided by spglib, the first one makes use spgdataset structure unlinke the second one that makes use of find_primitive and get_symmetry subroutines.


- Install:
Include the arch-spg.make file into your arch.make or into the Makefile, but pay attention to include this after the first recipe. To miss this error we highly recommend you to put this in arch.make. Last step include libsymspg.a to COMP_LIBS, then type make.

