#!/bin/bash
#
# Script to compile Spglib, and replace siesta_analysis.F with
# siesta_analysis_spg.F (this call the new subroutine) with a backup.
#
# It must be run in a path/to/siesta/Obj/spg-siesta like folder.
# 

# ===================================================================
# Useful function
# ===================================================================

function version_check {
  dir=$1
  if [ $(cat ../../version.info) != $(cat $dir/version.info) ]
  then
    echo WARNING: Siesta version little different from this implementation:
    echo Your version:  $(cat ../../version.info)
    echo Implementaion: $(cat $dir/version.info)
    echo You may found erros
    echo Continue? \(y\/n\)
    read a
    test $a != 'y' && exit
  fi 
}

# ===================================================================
# Variables
# ===================================================================

i=1
function i++ { i=$(echo $i+1 | bc); }
Src=../../Src
backup_dir=./backup

# ===================================================================
# Check Siesta version
# ===================================================================

echo 
echo $i\) Checking Siesta version compatibility...; i++
version=$(cat ../../version.info | cut -d '-' -f 2)

version_dir=''

case $version in
  '3.2')
    version_dir='siesta-v3.2'
  ;;
  'v4.0.2')                      
    version_dir='siesta-v4.0.2'
  ;;                             
  '4.1')
    version_dir='siesta-v4.1'
  ;;
  *)
    echo 'Unable version'
    exit
  ;;
esac

version_check $version_dir
version=$(cat $version_dir/version.info)
echo
echo $i\) Using version: $version; i++

# ===================================================================
# Download Spglib
# ===================================================================

if ! [ -d spglib ]
then
  echo 
  echo $i\) Downloading spglib; i++
  git clone https://github.com/atztogo/spglib.git
fi

# ===================================================================
# Compile Spglig
# ===================================================================

echo 
echo $i\) Compiling Spglib; i++
(
export LIBS='-lgomp';
export CFLAGS='-fopenmp'
cd spglib
mkdir -p _build
cd _build
cmake ..
make
)

# ===================================================================
# Prepare files
# ===================================================================

echo 
echo $i\) Preparing files to compilation; i++

echo cp spglib/example/spglib_f08.f90 $Src
cp spglib/example/spglib_f08.f90 $Src

# ===================================================================
# Back up & copy
# ===================================================================

mkdir -p $backup_dir

if ! [ -f $backup_dir/Makefile ]
then
  echo
  echo $i\) Backing up Makefile \(it\'ll be included the spg.make\); i++

  line=$(cat ../Makefile | grep -n 'siesta:' | head -1 | cut -d':' -f 1)
  line=$(echo $line - 1 | bc)

  cp ../Makefile $backup_dir/Makefile

  echo Appending \'include spg-siesta/$version_dir/spg.make\' in line: $line in Makefile

  head -n $line ../Makefile > ../Makefile.tmp
  echo "#Included by spg-siesta/configure " >> ../Makefile.tmp
  echo "include spg-siesta/$version_dir/spg.make" >> ../Makefile.tmp
  tail -n '+'$line ../Makefile >> ../Makefile.tmp
  mv ../Makefile.tmp ../Makefile
fi


if ! [ -f $backup_dir/siesta_analysis.F ]
then
  echo
  echo $i\) Backing up siesta_analysis.F \(it\'ll be added the call to siesta2bt subroutine\); i++
  echo cp $Src/siesta_analysis.F $backup_dir/siesta_analysis.F
  cp $Src/siesta_analysis.F $backup_dir/siesta_analysis.F
  
  echo
  echo $i\) Wrinting files:; i++
  echo cp $version_dir/siesta_analysis.F   $Src/siesta_analysis.F
  cp $version_dir/siesta_analysis.F $Src/siesta_analysis.F
 
  echo cp $version_dir/siesta2bt.f90 $Src/siesta2bt.f90
  cp $version_dir/siesta2bt.f90 $Src/siesta2bt.f90
fi

