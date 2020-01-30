#!/usr/bin/env sh

#Some flags variables
BUILD_FFTW=false
BUILD_FLTK=false
BUILD_RELION=true
BUILD_RELION_GUI=true
N_THREADS=$@

#Some path variables
RELION_HOME=$PWD
PREFIX=/opt/relion-1.2

#External libraries versions
VFFTW=fftw-3.2.2
VFLTK=fltk-1.3.0

# Some other vars
GREEN="\033[32m"
ENDC="\033[0m"



#################### FFTW ###########################
if $BUILD_FFTW; then
  echo -e "$GREEN Compiling $VFFTW ...$ENDC"
  echo -e "See $RELION_HOME/external/fftw_build.log for details"
  cd external
  tar -zxf $VFFTW.tar.gz
  cd $VFFTW
  CXXFLAG="-march=core2 -mfpmath=sse" ./configure --enable-threads --enable-shared prefix=$PREFIX > $RELION_HOME/external/fftw_build.log
  make $N_THREADS >> $RELION_HOME/external/fftw_build.log 
  make install >> $RELION_HOME/external/fftw_build.log 
  cd ../..
fi

#################### FLTK ###########################
if $BUILD_FLTK; then
  echo -e "$GREEN Compiling $VFLTK ...$ENDC"
  echo -e "See $RELION_HOME/external/fltk_build.log for details"
  cd external
  tar -zxf $VFLTK.tar.gz
  cd $VFLTK
  ./configure prefix=$PREFIX > $RELION_HOME/external/fltk_build.log
  make $N_THREADS >> $RELION_HOME/external/fltk_build.log
  make install >> $RELION_HOME/external/fltk_build.log
  cd ../..
fi

#################### RELION ###########################
if $BUILD_RELION; then
  echo -e "$GREEN Compiling relion ...$ENDC"
  echo -e "See $RELION_HOME/relion_build.log for details"
 ./configure prefix=$PREFIX --enable-mpi CPPFLAGS=-I$PREFIX/include  LDFLAGS=-L$PREFIX/lib > $RELION_HOME/relion_build.log
 make $N_THREADS >> $RELION_HOME/relion_build.log
 make install >> $RELION_HOME/relion_build.log
fi

#################### GUI ###########################
if $BUILD_RELION_GUI; then
 cd gui
 echo -e "$GREEN Compiling relion GUI ...$ENDC"
 $PREFIX/bin/fltk-config --compile gui.cpp 
 cp $RELION_HOME/gui/gui $PREFIX/bin/relion 
 cp $RELION_HOME/gui/qsub.csh $PREFIX/bin/qsub.csh
 cd ..
fi

echo "Done!"

