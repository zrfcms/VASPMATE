#Add your pbs script head here

VASPMATE_DIR='' #The path for vaspmate eg. /usr/bin/
VASP_RUN=''		#mpirun -np 4 vasp541

for SYS_name in SiC
do
cp ./structure/${SYS_name}.vasp ./INPOS

#>>> create POSCAR file
cp INPOS POSCAR

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
eval $VASP_RUN >> log1.vasp 
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
eval $VASP_RUN >> log2.vasp 
#--------------------------------------------------------------

#----------------------dos calculation--------------------------
#>>> DOS calculation
mkdir dos
cd dos
cp ../CONTCAR POSCAR
cp ../POTCAR ./
cp ../CHGCAR ./
cp ../FERMI_LEVEL ./
#>>> create KPOINTS file
${VASPMATE_DIR}/VASPMATE --ka 4000
cp NEWKPT KPOINTS
#>>> create INCAR file
${VASPMATE_DIR}/VASPMATE --i dos
${VASPMATE_DIR}/VASPMATE --i_replace ICHARG 11
#${VASPMATE_DIR}/VASPMATE --i_replace LORBIT 11
#${VASPMATE_DIR}/VASPMATE --i_replace NEDOS 1000
cp incar_dos INCAR
eval $VASP_RUN >> log3.vasp 
#--------------------------------------------------------------

#>>>select output mode
${VASPMATE_DIR}/VASPMATE  --dos -t >> vaspmate.log
#${VASPMATE_DIR}/VASPMATE --dos -a 
#${VASPMATE_DIR}/VASPMATE --dos -e 
#${VASPMATE_DIR}/VASPMATE --dos -s 1-4 N
#${VASPMATE_DIR}/VASPMATE --dos -m 1-4 N
#${VASPMATE_DIR}/VASPMATE --dos -o 1-3 s px py
#${VASPMATE_DIR}/VASPMATE --dos -bc all
done
#rm POSCAR KPOINTS* INCAR* POTCAR
#rm CHG* CONTCAR DOSCAR OSZICAR OUTCAR EIGENVAL PCDAT WAVECAR XDATCAR IBZKPT vasprun.xml
#END
