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

VASPMATE_DIR='/home/zrfcms11/panzhaocheng/VASPMATE/bin'
VASP_CAL='mpirun -np $NP -machinefile $PBS_NODEFILE --mca btl self,sm,openib ~/bin/vasp535'

for SYS_name in SiC 
do
cp ./structure/${SYS_name}.vasp ./INPOS

${VASPMATE_DIR}/VASPMATE --kpt3d 20
#>>> create POSCAR file
cp PRIMPOS POSCAR
cp POSCAR INPOS

#--------------------------Relaxtion--------------------------
#>>> create INCAR file using VASPMATE and the file incar_rlx and incar_stc are created
${VASPMATE_DIR}/VASPMATE --i rlx
${VASPMATE_DIR}/VASPMATE --i stc
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
cp CONTCAR INPOS
cp incar_stc INCAR
${VASPMATE_DIR}/VASPMATE --i_replace LCHARG T
#${VASPMATE_DIR}/VASPMATE --i_replace ISPIN 2
#${VASPMATE_DIR}/VASPMATE --i_replace ICHARG 2
${VASPMATE_DIR}/VASPMATE --ka 4000
cp NEWKPT KPOINTS
cp INPOS POSCAR
mpirun -np $NP -machinefile $PBS_NODEFILE --mca btl self,sm,openib ~/bin/vasp535 >> log2.vasp 
#--------------------------------------------------------------

#>>> get the efermi energy from the self-consistent calculation
${VASPMATE_DIR}/VASPMATE --dos -efermi >> FERMI_LEVEL

#----------------------band calculation--------------------------
#>>> Bandstructure calculation
mkdir band
cd band
cp ../CONTCAR POSCAR
cp ../POTCAR ./
cp ../CHGCAR ./
cp ../NEWKPATH ./
cp ../FERMI_LEVEL ./
#>>> create KPOINTS file
cp NEWKPATH KPOINTS
#>>> HSEcalculaiton
#VASPMATE --kahse 8000 0.05 G
#>>> create INCAR file
${VASPMATE_DIR}/VASPMATE --i pbs
cp incar_pbs INCAR
mpirun -np $NP -machinefile $PBS_NODEFILE --mca btl self,sm,openib ~/bin/vasp535 >> log3.vasp 
#--------------------------------------------------------------

${VASPMATE_DIR}/VASPMATE --band -bg
#>>>select output mode
#${VASPMATE_DIR}/VASPMATE --band -b >> vaspmate.log
#${VASPMATE_DIR}/VASPMATE --band -a 
#${VASPMATE_DIR}/VASPMATE --band -e
#${VASPMATE_DIR}/VASPMATE --band -s 1-4 N
#${VASPMATE_DIR}/VASPMATE --band -m 1-4 N
#${VASPMATE_DIR}/VASPMATE --band -o 1-3 s px py
done
#rm POSCAR KPOINTS* INCAR* POTCAR
#rm CHG* CONTCAR DOSCAR OSZICAR OUTCAR EIGENVAL PCDAT WAVECAR XDATCAR IBZKPT vasprun.xml
#END