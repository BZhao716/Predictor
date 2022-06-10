#!/bin/bash
pdir=`dirname $0`
name=$1
cov=$2

dn=( uniclust30 )
prefix=${name}-${dn}
echo "prefix is ${prefix}."
path=${path}

hhblits=${pdir}/programs/hhblits/hh-suite-master/hh-suite-master/bin/hhblits
if [ ! -s ${prefix}.aln ];then
	${hhblits} -i ${name}.fasta -d ${pdir}/database/uniclust30_2018_08/uniclust30_2018_08 -oa3m ${prefix}.a3m 
	egrep -v "^>" ${prefix}.a3m|sed 's/[a-z]//g' > ${prefix}.aln
	[ ! -s ${prefix}.aln ] && echo "[error]generating ${prefix}.aln" &&exit
fi
 [ ! -s ${prefix}-${cov}.aln ] && ${pdir}/coverage.sh ${prefix} $cov
