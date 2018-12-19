#!/bin/bash

Src=../../Src
backup_dir=$(grep 'backup_dir=' configure.sh | cut -d '=' -f 2)

echo cp $backup_dir/Makefile ../Makefile
cp $backup_dir/Makefile ../Makefile

echo cp $backup_dir/siesta_analysis.F $Src/siesta_analysis.F 
cp $backup_dir/siesta_analysis.F $Src/siesta_analysis.F 

echo rm $Src/siesta2bt.f90
rm $Src/siesta2bt.f90

echo rm $Src/spglib_f08.f90
rm $Src/spglib_f08.f90

echo rm -rf $backup_dir
rm -rf $backup_dir
