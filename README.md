# spg-siesta

Os arquivos siesta2bt_espresso.f90 e siesta2bt_spg.f90 possuem a rotina siesta2bt que busca os pontos k
pelo método utilizado pelo quantum espresso e pelo spglib, respectivamente. Copie o arquivo com o método desejado
na pasta Src/ do siesta com o nome siesta2bt.f90 e compile normalmente.

Comando para compilar a partir do diretório Src/:

$ alias scompila='cd ../Obj/;make;mv siesta siesta_bt; cd ../Src/'

$ scompila
