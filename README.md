# spg-siesta

#========================================================================================

[Pt]

Os arquivos siesta2bt_espresso.f90, siesta2bt_spg_dataset.f90 e siesta2bt_spg_findsym.f90 possuem a rotina siesta2bt que busca os pontos k pelo método utilizado pelo quantum espresso e os dois últimos por metodos diferentes da biblioteca spglib, sendo que o primeiro faz uso da estrutra
spgdataset e segunda nao.

- Instalação:

Inclua o arch-spg.make no arquivo Makefile do siesta, preste atenção para incluir após a primeira receita. Para não correr este risco nós recomendamos que a inclusão seja feita junta ou após a inclusao do arquivo arch.make. Em seguida, inclua libsymspg.a em COMP_LIBS e execute make.

#========================================================================================

[En]

The files siesta2bt_espresso.f90, siesta2bt_spg_dataset.f90 and siesta2bt_spg_findsym.f90 use the subroutine siesta2bt to search irreducible k-points, each by a method. The first uses the same kind of method used by quantum espresso, the others two use methods provided by spglib, the first one makes use spgdataset structure unlinke the second one that makes use of find_primitive and get_symmetry subroutines.


- Install:

Include the arch-spg.make file into Makefile of siesta, but pay attention to include this after the first recipe. To miss this error we highly recommend you to include this together or after the arch.make inclusion. Last step include libsymspg.a to COMP_LIBS, then type make.

#========================================================================================

[Es]

Los archivos siesta2bt_espresso.f90, siesta2bt_spg_dataset.f90 y siesta2bt_spg_findsym.f90 poseen la rutina siesta2bt que busca los puntos k por el método utilizado por el quantum espresso y los dos últimos por métodos diferentes de la biblioteca spglib, siendo que el primero hace uso de la estructura spgdataset y segunda no.

- Instalación:

Incluya el arch-spg.make en su archivo arch.make o en su archivo de Makefile, preste atención a incluir después de la primera receta. Para no correr este riesgo nosotros recomendamos que la inclusión sea hecha junta o después de la inclusión del archivo arch.make. A continuación, incluya libsymspg.a en COMP_LIBS y ejecute make.
