#!/bin/bash

out_file_name='topo_potential'
N_replica_pt=10
grid_step=0.05
grid_max=3.0

rm -f ${out_file_name}

# to simulate imaginary-theta term with multicanonic, choose topo-potetial V(i) = -im_theta*i
im_theta=0.5
lambda=0.5

# print topo-potential to file
for i in $( seq -${grid_max} ${grid_step} ${grid_max} ); do
	awk -v i=${i} 'BEGIN{ printf "%.5lf", i}' >> ${out_file_name}
	for j in $( seq 1 1 ${N_replica_pt} ); do
		awk -v i=${i} -v j=${j} -v t=${lambda} -v n=${N_replica_pt} 'BEGIN{ printf " %.18lf", -t*sqrt(i*i)*(1-(j-1)/(n-1))}' >> ${out_file_name}
	done
	awk 'BEGIN{ printf "\n"}' >> ${out_file_name}
done