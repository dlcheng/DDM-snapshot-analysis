#!/bin/sh

#Lists
declare -a Old_name
declare -a New_name
declare -a Ahf_base

Old_name=( 'Decay-test-18' 'Decay-test-21' 'Decay-test-20' 'Decay-test-22' 'Decay-test-31' 'Decay-test-23' 'Decay-test-24' 'Decay-test-25' 'Decay-test-19' 'Decay-test-27' 'Decay-test-26' 'Decay-test-28' 'Decay-test-29' 'Decay-test-30' 'Decay-test-32-altix' 'Decay-test-33' 'Decay-test-34' 'Decay-test-35' 'Decay-test-36' 'Decay-test-37-altix' 'Decay-test-38-altix' )
New_name=( 'Decay-1' 'Decay-2' 'Decay-3' 'Decay-4' 'Decay-5' 'Decay-6' 'Decay-7' 'Decay-8' 'Decay-9' 'Decay-10' 'Decay-11' 'Decay-12' 'Decay-13' 'Decay-14' 'Decay-15-altix' 'Decay-16' 'Decay-17' 'Decay-18' 'Decay-19' 'Decay-20' 'Decay-21' )
Ahf_base=( 'LTic_decay_test_18' 'LTic_decay_test_21' 'LTic_decay_test_20' 'LTic_decay_test_22' 'LTic_decay_test_31' 'LTic_decay_test_23' 'LTic_decay_test_24' 'LTic_decay_test_25' 'LTic_decay_test_19' 'LTic_decay_test_27' 'LTic_decay_test_26' 'LTic_decay_test_28' 'LTic_decay_test_29' 'LTic_decay_test_30' 'LTic_decay_test_32_altix' 'LTic_decay_test_33' 'LTic_decay_test_34' 'LTic_decay_test_35' 'LTic_decay_test_36' 'LTic_decay_test_37_altix' 'LTic_decay_test_38_altix' )

# total 21 files
for i in $(seq 0 20)
do 
make clean > /dev/null
rm define.h

echo "
/* Only valid for AHF version >=1.0
 * and dark matter only simulations
 * with total number of particles less than 2^32
 */

#define AHF_INPUT_PATH               \"/home/dlcheng/data/dlcheng/DDM/Project-DDM-cosmology/Install/ahf-v1.0-049/Result/Old-cosmo-param/${Old_name[$i]}/z-0/\"
#define GADGET_INPUT_PATH            \"/home/dlcheng/data/dlcheng/DDM/Project-DDM-cosmology/Install/Gadget2/DM-cosmo/Result/Old-cosmo-param/${Old_name[$i]}/\"
#define AHF_FILE_BASE                \"${Ahf_base[$i]}\"
#define GADGET_FILE_BASE             \"snapshot_009\"
#define OUTPUT_PATH                  \"./Results/${New_name[$i]}/z-0/\"

#define RESOLUTED_HALO_PART_NUM       50         /* number of undecayed particle in a halo to be considered as trustworthy */
#define PROFILE_BIN_NUM               30         /* the profile bin number in log scale, physical length */

#define ENVIR_MASS_CRIT               1e8        /* halo with mass less than this value
                                                  * will be classified as enviroment, in unit of M_sun/h */

" > ./define.h

make > /dev/null 2>&1
echo "Analyze ${New_name[$i]}"
time ./Halo_info

done 


