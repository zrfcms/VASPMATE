#Add your pbs script head here

VASPMATE_DIR='' #The path for vaspmate eg. /usr/bin/
VASP_RUN=''		#mpirun -np 4 vasp541

for SYS_name in SiC
do
cp ./structure/${SYS_name}.vasp ./INPOS

#--------------------------Relaxtion--------------------------
#>>> create INCAR file using VASPMATE and the file incar_rlx and incar_stc are created
${VASPMATE_DIR}/VASPMATE --i rlx
cp incar_rlx INCAR
#>>> correct parameters
#${VASPMATE_DIR}/VASPMATE --i_replace NSW 100
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
#Bandstructure calculation
cp CONTCAR INPOS
#>>> create KPOINTS file
${VASPMATE_DIR}/VASPMATE --kpt3d 20
#cp NEWKPATH KPOINTS
#>>> HSEcalculaiton
${VASPMATE_DIR}/VASPMATE --kahse 8000 0.05 G
cp NEWKPT KPOINTS
#>>> create INCAR file
${VASPMATE_DIR}/VASPMATE --i stc
cp incar_stc INCAR
#${VASPMATE_DIR}/VASPMATE --i_replace ICHARG 11
#>>> create POSAR file
cp PRIMPOS POSCAR
eval $VASP_RUN >> log2.vasp
#--------------------------------------------------------------

${VASPMATE_DIR}/VASPMATE --dos -efermi >> FERMI_LEVEL

#-------------------hse band calculation-----------------------
#>>> HSE calculation
${VASPMATE_DIR}/VASPMATE --i stc hse
cp incar_stc_hse INCAR
${VASPMATE_DIR}/VASPMATE --i_replace LORBIT 11
eval $VASP_RUN >> log3.vasp 
#--------------------------------------------------------------

#>>>select output mode
${VASPMATE_DIR}/VASPMATE --band -hb >> vaspmate.log
#${VASPMATE_DIR}/VASPMATE --band -ha 
#${VASPMATE_DIR}/VASPMATE --band -he
#${VASPMATE_DIR}/VASPMATE --band -hs 1-4 B
#${VASPMATE_DIR}/VASPMATE --band -ho 1-3 s px py
#${VASPMATE_DIR}/VASPMATE --band -ho 4 all
done
#rm POSCAR KPOINTS* INCAR* POTCAR
#rm CHG* CONTCAR DOSCAR OSZICAR OUTCAR EIGENVAL PCDAT WAVECAR XDATCAR IBZKPT vasprun.xml
#END
