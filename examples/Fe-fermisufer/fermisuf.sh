#Add your pbs script head here

VASPMATE_DIR='' #The path for vaspmate eg. /usr/bin/
VASP_RUN=''		#mpirun -np 4 vasp541

for SYS_name in Fe
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
eval $VASP_RUN >> log1.vasp 
#--------------------------------------------------------------

#------------------option(self-consistent)---------------------
${VASPMATE_DIR}/VASPMATE --i VASPMATE --i stc
cp CONTCAR INPOS
cp incar_stc INCAR
${VASPMATE_DIR}/VASPMATE --i_replace LCHARG T
${VASPMATE_DIR}/VASPMATE --i_replace LORBIT 11
#${VASPMATE_DIR}/VASPMATE --i_replace ISPIN 2
#${VASPMATE_DIR}/VASPMATE --i_replace ICHARG 2
${VASPMATE_DIR}/VASPMATE --fskm 64 64 64
#${VASPMATE_DIR}/VASPMATE --fska 8000
#${VASPMATE_DIR}/VASPMATE --fskv 0.05
cp FERMIKPT KPOINTS
cp INPOS POSCAR
eval $VASP_RUN >> log2.vasp 
#--------------------------------------------------------------

${VASPMATE_DIR}/VASPMATE --fsxd
${VASPMATE_DIR}/VASPMATE --fs
#${VASPMATE_DIR}/VASPMATE --fs -b 1 2 -o 1 s
done
#rm POSCAR KPOINTS* INCAR* POTCAR
#rm CHG* CONTCAR DOSCAR OSZICAR OUTCAR EIGENVAL PCDAT WAVECAR XDATCAR IBZKPT vasprun.xml
#END
