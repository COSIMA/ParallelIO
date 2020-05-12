#!/bin/sh -l
# To run:
#    sbatch -A gsd-hpcs --ntasks=4 ./slurm_run_4proc_tests.sh
#SBATCH --time=10      # time limit in minutes                                                                       #SBATCH -o %x.o%j      # Output file name to match the current scheme
set -x

env | grep SLURM

printf 'running PIO tests...\n'

PIO_TESTS='test_intercomm2 test_async_mpi test_spmd test_rearr test_async_simple '\
'test_async_3proc test_async_4proc test_iosystem2_simple test_iosystem2_simple2 '\
'test_iosystem2 test_iosystem3_simple test_iosystem3_simple2 test_iosystem3 test_pioc '\
'test_pioc_unlim test_pioc_putget test_pioc_fill test_darray test_darray_multi '\
'test_darray_multivar test_darray_multivar2 test_darray_multivar3 test_darray_1d '\
'test_darray_3d test_decomp_uneven test_decomps test_darray_async_simple '\
'test_darray_async test_darray_async_many test_darray_2sync test_async_multicomp '\
'test_darray_fill test_darray_vard test_async_1d'

success1=true
success2=true
for TEST in $PIO_TESTS
do
  	 success1=false
    echo "running ${TEST}"
    srun ./${TEST} && success1=true
    if test $success1 = false; then
        break
    fi
done
