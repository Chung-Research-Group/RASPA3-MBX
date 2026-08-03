#! /bin/bash

struc=$1
T=$2

DIR=$(pwd)

mkdir Result
mkdir Result/mmolg

cd calculation

#while read struc
#do

#while read T
#do

echo "Pressure,uptake,error,Q,Q_err" > ${DIR}/Result/mmolg/output_${struc}_${T}.txt

while read P
do
output_file="${struc}/${T}K/${P}Pa/output/*.txt"
init_flag=$(cat $output_file | grep "Current cycle:" | tail -1 | awk '{print $NF}')
finish_flag=$(cat $output_file | grep "Total simulation time" | wc -l)
if [ "$finish_flag" = 1 ]
then
 q=$(cat $output_file | grep "Abs. loading average" | grep "mol/kg-framework" | awk '{print $4","$6}')
 Q=$(cat $output_file | grep "Enthalpy of adsorption" -A 1 | grep "kJ/mol" | awk '{print $1","$3}')
 echo ${P},${q},${Q} >> ${DIR}/Result/mmolg/output_${struc}_${T}.txt
else
 if [ "$init_flag" = 20000 ]
 then
  q=$(cat $output_file | grep "absolute adsorption:" -A 2 | tail -1 | awk '{print $3}'| tr -d '(')
  echo ${P},${q}, >> ${DIR}/Result/mmolg/output_${struc}_${T}.txt
 else
  q=""
  echo ${P},${q}, >> ${DIR}/Result/mmolg/output_${struc}_${T}.txt
 fi
fi
done < ${DIR}/p_range
#done < ${DIR}/t_range
#done < ${DIR}/list_$1


