#!/bin/bash

NARGS=2

usage () 
{
  echo "------------------------------------------------------------------------------------------"
  echo "Usage: `basename $0` {klsampl/klsampl_u/klcov/klcov_u} corr_len                           "
  echo "Options:                                                                                  "
  echo "      (1) klsampl - for this option the run script invokes \"kl_sample2D.x\" with several "
  echo "                    command line arguments. To explore the command-line options type:     "
  echo "                                 kl_sample2D.x -h                                         "
  echo "         (a) - generates multivariate normal (MVN) 2D random field (RF) samples.          "
  echo "               the MVN is based on a covariance matrix                                    "
  echo "                        cov(x1,x2)=exp(-(x1-x2)^2/cl^2)                                   "
  echo "               where cl is the second parameter, \"corr_len\". Each RF sample consist     "
  echo "               of values on a structured spatial grid over [0,1]x[0,1]. Currently 65x65   "
  echo "               grid points, are used and the grid is clustered near the boundary.         "
  echo "               Some grid options can be changed from command line and some from inside the"
  echo "               code                                                                       "
  echo "         (b) - a new, approximate, covariance matrix is computed from the RF samples. The "
  echo "               Karhunen-Loeve decomposition is then based on this matrix. Several output  "
  echo "               files are generated:                                                       "
  echo "                    mean.dat    - mean of samples generated at (a)                        "
  echo "                    eig.dat     - list of eigenvalues for the KL expansion                "
  echo "                    KLmodes.dat - scaled KL modes (multiplied by the square root of the   "
  echo "                                  corresponding eigenvalue                                "
  echo "                    xi_data.dat - realizations of KL eigenvalues, obtained by projecting  "
  echo "                                  samples on the each KL mode                             "
  echo "      (2) klsampl_u - for this option the run script invokes \"kl_sample2Du.x\" with several"
  echo "                    command line arguments. To explore the command-line options type:     "
  echo "                                 kl_sample2Du.x -h                                         "
  echo "          - this option is similar to the similar to the one above except it uses and     "
  echo "            an unstructured grid                                                          "
  echo "      (3) klcov - for this option the workflow is similar to (1), except that KL decomp.  "
  echo "                  is performed using the analytical covariance matrix, and no RF samples  "
  echo "                  are generated. For this case, only the \"eig.dat\" and \"KLmodes.dat\"  "
  echo "                  are output.                                                             "
  echo "      (4) klcov_u - similar to (3) except it uses and unstructured grid                   "
  echo "------------------------------------------------------------------------------------------"
}
 
if [ $# -lt $NARGS ]
then
  usage
  exit $E_BADARGS
fi

run=$1
clen=$2

# no. of KL eigenvalues/eigenvectors
nkl=129
if [ $# -ge 3 ]
then
  nkl=$3
fi

# List for number of samples
#declare -a slist=(1024 2048 4096)
declare -a slist=(4096)
sLen=${#slist[@]}

sigma=5.0
ctype="SqExp"

if [ ${run} == "klsampl" ]
then
  for (( i=0; i<${sLen}; i++ ));
  do
    echo "-----------------------------------------------"
    ./kl_sample2D.x -p ${slist[$i]} -s ${sigma} -l ${clen} -e ${nkl}
    rsuff="${clen}_${slist[$i]}"
    paste samples.dat > samples2D_${rsuff}.dat
    mv mean.dat     mean2D_${rsuff}.dat
    mv eig.dat      eig2D_${rsuff}.dat
    mv KLmodes.dat  KLmodes2D_${rsuff}.dat
    mv xi_data.dat  xidata2D_${rsuff}.dat
    mv cov.dat      cov2D_${rsuff}.dat
    resdir="klsampl2D_${rsuff}"
    if [ ! -d "${resdir}" ]; then
      mkdir ${resdir}
    fi
    mv xgrid.dat ygrid.dat xg1d.dat ${resdir}
    /bin/mv *${rsuff}.dat ${resdir}
  done
fi

if [ ${run} == "klcov" ]
then
  ./kl_sample2D.x -c ${ctype} -s ${sigma} -l ${clen}
  mv eig.dat     eig2D_${clen}_${ctype}_anl.dat
  mv KLmodes.dat KLmodes2D_${clen}_${ctype}_anl.dat
  mv cov.dat     cov2D_${clen}_${ctype}_anl.dat
  if [ ! -d "klcov2D_${ctype}_${clen}" ]; then
    mkdir klcov2D_${ctype}_${clen}
  fi
  mv *${clen}_${ctype}_anl.dat klcov2D_${ctype}_${clen}/.

fi

if [ ${run} == "klsampl_u" ]
then
  declare -a slist=(4096)
  sLen=${#slist[@]}
  for (( i=0; i<${sLen}; i++ ));
  do
    echo "-----------------------------------------------"
    ./kl_sample2Du.x -p ${slist[$i]} -s ${sigma} -l ${clen} -e ${nkl}
    rsuff="2Du_${clen}_${slist[$i]}"
    mv samples.dat  samples${rsuff}.dat
    mv mean.dat     mean${rsuff}.dat
    mv eig.dat      eig${rsuff}.dat
    mv KLmodes.dat  KLmodes${rsuff}.dat
    mv xi_data.dat  xidata${rsuff}.dat
    mv cov.dat      cov${rsuff}.dat
    resdir="klsampl${rsuff}"
    if [ ! -d "${resdir}" ]; then
      mkdir ${resdir}
    fi
    /bin/mv *${rsuff}.dat ${resdir}
  done
fi

if [ ${run} == "klcov_u" ]
then
  rsuff="${ctype}_${clen}"
  ./kl_sample2Du.x -c ${ctype} -s ${sigma} -l ${clen} -e ${nkl}
  mv eig.dat       eig2Du_${rsuff}_anl.dat
  mv KLmodes.dat   KLmodes2Du_${rsuff}_anl.dat
  mv cov.dat       cov2Du_${rsuff}_anl.dat
  rdir=klcov2Du_${rsuff}
  if [ ! -d "${rdir}" ]; then
    mkdir ${rdir}
  fi
  mv *${rsuff}_anl.dat ${rdir}/.

fi
