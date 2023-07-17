VASPMATE_DIR='' #The path for vaspmate eg. /usr/bin/
#temperature 300k -T
#pressure 1atm -p
#spin multiplicity -s 1
#ignore the contribution of advection and rotation to the thermodynamics -i 1
#Raising lower frequencies  -l 1
${VASPMATE_DIR}/VASPMATE --thermo -T 300 -p 1 -s 1 -l 1 -i 1
#generate Thermol_Info.dat