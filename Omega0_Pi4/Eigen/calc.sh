#!/bin/bash
for((Ns=2;Ns<=16;Ns++))
do
make clean
sed "s/Ns=2/Ns=${Ns}/g" params_temp.f90  > params.f90
make
./eigen
done