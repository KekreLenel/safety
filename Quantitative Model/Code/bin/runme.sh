#!/bin/sh
time0=$SECONDS

debug_mode=0        # run fast (0) or debug (1)
make_prm_files=1    # create parameter files (1)
compile_exe=1       # compile executable (1)
solve_model=1       # solve model (1)
create_results=1    # run matlab results script (1) 

# 1. SETUP OF ENVIRONMENT

src_folder="../src/fortran/"

ulimit -s unlimited
export OMP_STACKSIZE=500m

# set variables for for OMP and NAG
export OMP_NUM_THREADS=32
export MKL_NUM_THREADS=1
export NAG_KUSARI_FILE=lic.txt

# set environmental variables for intel fortran compiler
echo  "Set Intel Fortan Compiler vars."
source /opt/intel/oneapi/setvars.sh intel64 || exit 1
echo  "DONE"
echo  ""

# set environmental variables for nag library
nagfdir=/usr/local/NAG/nll6i287nl
echo  "Set NAG vars."
source ${nagfdir}/scripts/nagvars.sh -quiet int64 vendor dynamic || exit 1
echo  "DONE"
echo  ""

if [ $debug_mode -eq 0 ]
then
    # fast compilation 
    fcompile="ifort -qmkl -qopenmp -O3 -xhost ${NAGLIB_FFLAGS} ${NAGLIB_FINCLUDE}"
else 
    # debug compilation 
    fcompile="ifort -qmkl -O0 -qopenmp -init=snan -debug -traceback -check all,noarg_temp_created -nowarn -CB ${NAGLIB_FFLAGS} ${NAGLIB_FINCLUDE}"
fi 

# create output folder  
mkdir ../output
mkdir ../output/tmp

# 2. CREATE PARAMETER FILES

time1=$SECONDS

if [ $make_prm_files -eq 1 ]
then
    echo  "Create parameter files."
    rm ../src/params/*.csv
    matlab -nodisplay -nosplash -nodesktop <../src/params/create_param_files.m > ../output/tmp/param_output.txt 2> ../output/tmp/param_error.txt || exit 1
    echo  "DONE."
fi 
time2=$SECONDS

# 3. COMPILE EXECUTABLE

if [ $compile_exe -eq 1 ]
then
    start_time=$SECONDS
 	echo "Remove previous compilation files."
 	rm main
 	rm *.o
 	rm *.mod
 	rm *_genmod.f90
    echo  "DONE."
 	echo ""

 	echo "Compile executable."

	${fcompile}  ${src_folder}/AuxCodes/base_lib.f90 ${src_folder}/AuxCodes/mod_normal.f90 ${src_folder}/AuxCodes/mod_markov.f90 \
			 ${src_folder}/AuxCodes/mod_smolyak.f90 ${src_folder}/mod_param.f90 ${src_folder}/mod_calc.f90  \
             ${src_folder}/mod_results.f90 ${src_folder}/main.f90 -o main ${NAGLIB_FLINK} || exit 1
 
    echo  "DONE."
	echo ""
fi 
time3=$SECONDS

# 4. SOLVE MODEL ACROSS CALIBRATIONS

# get number of calibrations to run
file="../output/tmp/n_comp.txt"
n_comp=$(cat "$file")
if [ $solve_model -eq 1 ]
then
# loop over calibrations
for i in `seq 1 $n_comp`
do
    # remove previous output files 
	foo="../output/tmp/res_${i}"
	rm -r $foo
	mkdir $foo

 	echo "Run calibration ${i}."
 	rm output.txt
 	./main $i | tee ../output/tmp/output_${i}.txt
    echo  "DONE."
done
fi
time4=$SECONDS

echo "Remove compilation files."
rm *.exe
rm *.o
rm *.mod
rm *_genmod.f90
echo  "DONE."
echo ""

# 5. CREATE RESULTS

if [ $create_results -eq 1 ]
then

    rm -r ../output/figures
    mkdir ../output/figures
    mkdir ../output/figures/addl_figs
    rm -r ../output/tables
    mkdir ../output/tables
    mkdir ../output/tables/num_checks

    cp ../src/params/create_param_files.m ../output/tables 
     
    echo "Run MATLAB results script."
    matlab -nodisplay -nosplash -nodesktop <../src/matlab/main.m > ../output/tmp/results_output.txt 2> ../output/tmp/results_error.txt
    echo "DONE"
    
    source make_pdf.sh

fi

time5=$SECONDS
step1=$(($time1 - $time0))
step2=$(($time2 - $time1))
step3=$(($time3 - $time2))
step4=$(($time4 - $time3))
step5=$(($time5 - $time4))
total=$(($time5 - $time0))
rm ../output/tmp/runtime.txt
printf "RUNTIME \n"                                            >> ../output/tmp/runtime.txt
printf "============================================ \n"       >> ../output/tmp/runtime.txt 
printf "SETUP     $(($step1/60)) min $(($step1%60)) sec \n"    >> ../output/tmp/runtime.txt   
printf "PARAMS    $(($step2/60)) min $(($step2%60)) sec \n"    >> ../output/tmp/runtime.txt 
printf "COMPILE   $(($step3/60)) min $(($step3%60)) sec \n"    >> ../output/tmp/runtime.txt 
printf "RUN       $(($step4/60)) min $(($step4%60)) sec \n"    >> ../output/tmp/runtime.txt 
printf "OUTPUT    $(($step5/60)) min $(($step5%60)) sec \n"    >> ../output/tmp/runtime.txt 
printf "============================================ \n"       >> ../output/tmp/runtime.txt 
printf "TOTAL     $(($total/60)) min $(($total%60)) sec \n"    >> ../output/tmp/runtime.txt 


