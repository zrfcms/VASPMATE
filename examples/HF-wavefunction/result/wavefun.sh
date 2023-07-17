#Add your pbs script head here
#!/bin/bash
#PBS -N panzhaocheng
#PBS -l nodes=2:ppn=24
#PBS -j oe
#PBS -l walltime=240:00:00
#PBS -q normal

source /public/software/profile.d/compiler_intel-composer_xe_2015.2.164.sh
source /public/software/profile.d/mpi_openmpi-1.8.5-intel.sh

cd $PBS_O_WORKDIR
NP=`cat $PBS_NODEFILE|wc -l`

VASPMATE_DIR='/home/zrfcms11/panzhaocheng/VASPMATE/bin' #The path for vaspmate eg. /usr/bin/

for SYS_name in HF
do
cp ./structure/${SYS_name}.vasp ./INPOS

#>>> create POSCAR file
cp INPOS POSCAR

#--------------------------Relaxtion--------------------------
#>>> create INCAR file using VASPMATE and the file incar_rlx and incar_stc are created
${VASPMATE_DIR}/VASPMATE --i rlx
#>>> HSE calculation
#${VASPMATE_DIR}/VASPMATE --i hse stc
cp incar_rlx INCAR
#>>> correct parameters
#${VASPMATE_DIR}/VASPMATE --i_replace ISPIN 2
#${VASPMATE_DIR}/VASPMATE --i_replace ENCUT 400
#>>> create KPOINTS file
${VASPMATE_DIR}/VASPMATE --ka 4000
cp NEWKPT KPOINTS
#>>> create POSAR file
cp INPOS POSCAR
#>>> create POTCAR file
${VASPMATE_DIR}/VASPMATE --pot -PBE
mpirun -np $NP -machinefile $PBS_NODEFILE --mca btl self,sm,openib ~/bin/vasp535 >> log1.vasp 
#--------------------------------------------------------------

#------------------option(self-consistent)---------------------
${VASPMATE_DIR}/VASPMATE --i stc wfn
cp CONTCAR INPOS
cp incar_stc_wfn INCAR
${VASPMATE_DIR}/VASPMATE --i_replace LCHARG T
#${VASPMATE_DIR}/VASPMATE --i_replace ISPIN 2
#${VASPMATE_DIR}/VASPMATE --i_replace ICHARG 2
${VASPMATE_DIR}/VASPMATE --ka 4000
cp NEWKPT KPOINTS
cp INPOS POSCAR
mpirun -np $NP -machinefile $PBS_NODEFILE --mca btl self,sm,openib ~/bin/vasp535 >> log2.vasp 
#--------------------------------------------------------------

${VASPMATE_DIR}/VASPMATE --wfun -k 1 -b 1
#for i in `seq 1 48`
#do
#${VASPMATE_DIR}/VASPMATE --wfun -k 1 -b $i
#done 

done
#rm POSCAR KPOINTS* INCAR* POTCAR
#rm CHG* CONTCAR DOSCAR OSZICAR OUTCAR EIGENVAL PCDAT WAVECAR XDATCAR IBZKPT vasprun.xml
#END
