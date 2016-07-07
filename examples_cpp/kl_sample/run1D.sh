#!/bin/bash

NARGS=2

usage () 
{
  echo "------------------------------------------------------------------------------------------"
  echo "Usage: `basename $0` {klsampl/klcov} corr_len                                  "
  echo "Options:                                                                                  "
  echo "      (1) klsampl - for this option the run script invokes \"kl_sample.x\" with several   "
  echo "                    command line arguments. To explore the command-line options type:     "
  echo "                                 kl_sample.x -h                                           "
  echo "         (a) - generates multivariate normal (MVN) 1D random field (RF) samples.          "
  echo "               the MVN is based on a covariance matrix                                    "
  echo "                        cov(x1,x2)=exp(-(x1-x2)^2/cl^2)                                   "
  echo "               where cl is the second parameter, \"corr_len\". Each RF sample consist     "
  echo "               of values on a uniform spatial grid over [0,1]. Currently 129 grid points  "
  echo "               are used, but this number can be changed through specific command-line     "
  echo "               arguments for \"kl_sample.x\"                                              "
  echo "         (b) - a new, approximate, covariance matrix is computed from the RF samples. The "
  echo "               Karhunen-Loeve decomposition is then based on this matrix. Several output  "
  echo "               files are generated:                                                       "
  echo "                    mean.dat    - mean of samples generated at (a)                        "
  echo "                    eig.dat     - list of eigenvalues for the KL expansion                "
  echo "                    KLmodes.dat - scaled KL modes (multiplied by the square root of the   "
  echo "                                  corresponding eigenvalue                                "
  echo "                    xi_data.dat - realizations of KL eigenvalues, obtained by projecting  "
  echo "                                  samples on the each KL mode                             "
  echo "      (2) klcov - for this option the workflow is similar to (1), except that KL decomp.  "
  echo "                  is performed using the analytical covariance matrix, and no RF samples  "
  echo "                  are generated. For this case, only the \"eig.dat\" and \"KLmodes.dat\"  "
  echo "                  are output.                                                             "
  echo "------------------------------------------------------------------------------------------"
}
 
if [ $# -ne $NARGS ]
then
  usage
  exit $E_BADARGS
fi

run=$1
clen=$2

declare -a slist=(512 1024 2048 4096 8192 16384 32768 65536 131072)
declare -a slist=(512 8192 131072)
sLen=${#slist[@]}

sigma=5.0
ctype="SqExp"

if [ ${run} == "klsampl" ]
then
  for (( i=0; i<${sLen}; i++ ));
  do
    echo "-----------------------------------------------"
    ./kl_sample.x -p ${slist[$i]} -s ${sigma} -l ${clen}
    rsuff="${clen}_${slist[$i]}"
    paste xgrid.dat samples.dat > samples_${rsuff}.dat
    mv mean.dat     mean_${rsuff}.dat
    mv eig.dat      eig_${rsuff}.dat
    mv KLmodes.dat  KLmodes_${rsuff}.dat
    mv xi_data.dat  xidata_${rsuff}.dat
    mv cov.dat      cov_${rsuff}.dat
    resdir="klsampl_${rsuff}"
    if [ ! -d "${resdir}" ]; then
      mkdir ${resdir}
    fi
    /bin/mv *${rsuff}.dat ${resdir}
  done
fi

if [ ${run} == "klcov" ]
then
  ./kl_sample.x -c ${ctype} -s ${sigma} -l ${clen}
  mv eig.dat     eig_${clen}_${ctype}_anl.dat
  mv KLmodes.dat KLmodes_${clen}_${ctype}_anl.dat
  mv cov.dat     cov_${clen}_${ctype}_anl.dat
  if [ ! -d "klcov_${ctype}_${clen}" ]; then
    mkdir klcov_${ctype}_${clen}
  fi
  mv *${clen}_${ctype}_anl.dat klcov_${ctype}_${clen}/.

fi
