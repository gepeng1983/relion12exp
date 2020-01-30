#!/bin/tcsh
#$ -pe XXXqueueXXX XXXcoresXXX
#$ -e XXXerrfileXXX
#$ -o XXXoutfileXXX
#$ -cwd
#$ -S /bin/tcsh

# Environment
source ~/.cshrc

mpiexec --bynode -n XXXmpinodesXXX  XXXcommandXXX
