#!/bin/bash

#=============================================================================================
#
#
# Current list of sites: BCI-PANAMA, KORUP-CAMEROON, ZF2-BRAZIL, 
#                        ITURI-DRC, LAMBIR-MALAYSIA, HUAI KHA KHAENG-THAILAND
#
#
# site_lat_c =   9.1543, 5.07389,  -2.60909722,  1.4368,   4.1865,  15.6324
# site_lon_c = 280.1539, 8.85472, 299.7907,     28.5826, 114.017,   99.217
#=============================================================================================

MACH=modex
COMP=ICLM45ED
GITHASH=`git log -n 1 --format=%h`
CASE=BR1x1_${GITHASH}

CROOT=/data/Model_Output/fates_testing/

DIN_LOC_ROOT=/data/Model_Data/cesm_input_datasets
DOMAIN_PATH=${DIN_LOC_ROOT}/share/domains/

WORKDIR=/data/sserbin/Modeling/fates-clm.latest/cime/scripts/
cd $WORKDIR

export CASEROOT=${CROOT}${CASE}
echo "CREATING NEW CASE IN "${CASEROOT}

rm -rf ${CASEROOT}
./create_newcase -case ${CASEROOT} -res 1x1_brazil -compset ${COMP} -mach ${MACH} -compiler gnu -project fates_benchmarking

cd ${CASEROOT}

# Modifying : env_mach_pes.xml
./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val 1
./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val 1
./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val 1
./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val 1
./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val 1
./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val 1
./xmlchange -file env_mach_pes.xml -id NTASKS_ROF -val 1
./xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val 1
./xmlchange -file env_mach_pes.xml -id MAX_TASKS_PER_NODE -val 1
./xmlchange -file env_mach_pes.xml -id TOTALPES -val 1
# Modifying : env_build.xml
./xmlchange -file env_build.xml -id CESMSCRATCHROOT -val ${CASEROOT}
./xmlchange -file env_build.xml -id DEBUG -val FALSE
./xmlchange -file env_build.xml -id EXEROOT -val ${CASEROOT}/bld

# Modifying : env_run.xml
./xmlchange -file env_run.xml -id STOP_N -val 25
./xmlchange -file env_run.xml -id RUN_STARTDATE -val '1500-01-01'
./xmlchange -file env_run.xml -id STOP_OPTION -val nyears
./xmlchange -file env_run.xml -id REST_N -val 1
./xmlchange -file env_run.xml -id DATM_CLMNCEP_YR_START -val 1996
./xmlchange -file env_run.xml -id DATM_CLMNCEP_YR_END -val 2004

./xmlchange -file env_run.xml -id DIN_LOC_ROOT -val ${DIN_LOC_ROOT}
./xmlchange -file env_run.xml -id DIN_LOC_ROOT_CLMFORC -val '$DIN_LOC_ROOT'
./xmlchange -file env_run.xml -id DOUT_S_SAVE_INTERIM_RESTART_FILES -val TRUE
./xmlchange -file env_run.xml -id DOUT_S_SAVE_EVERY_NTH_RESTART_FILE_SET -val 1
./xmlchange -file env_run.xml -id DOUT_S -val TRUE
./xmlchange -file env_run.xml -id RUNDIR -val ${CASEROOT}/run
./xmlchange -file env_run.xml -id BATCHQUERY -val ''
./xmlchange -file env_run.xml -id BATCHSUBMIT -val ''

./xmlchange -file env_run.xml -id PIO_DEBUG_LEVEL -val 0

#./cesm_setup
echo "*** Running case.setup ***"
./case.setup


# Modify run script
#

#hist_fincl2='ED_GPP_GD_SCPF','ED_NPP_LEAF_GD_SCPF',,'ED_NPP_FNRT_GD_SCPF','ED_NPP_BGSW_GD_SCPF','ED_NPP_BGDW_GD_SCPF','ED_NPP_AGSW_GD_SCPF','ED_NPP_AGDW_GD_SCPF','ED_NPP_STOR_GD_SCPF','ED_NPP_SEED_GD_SCPF','ED_DDBH_GD_SCPF','ED_BA_GD_SCPF','ED_NPLANT_GD_SCPF','ED_M1_GD_SCPF','ED_M2_GD_SCPF','ED_M3_GD_SCPF','ED_M4_GD_SCPF','ED_M5_GD_SCPF'
#hist_mfilt             = 480,12
#hist_nhtfrq            = -1,0

cat >> user_nl_clm << \EOF
finidat = ''
hist_mfilt             = 480
hist_nhtfrq            = -1
EOF

# Modify user_nl_datm
cat >> user_nl_datm << EOF
EOF
# building case :
#./${CASE}.build
echo "*** Running case.build ***"
./case.build


cd $WORKDIR
