NULLSTRING:=
QB = ./QB/QB_dump.cpp ./QB/QB_function.cpp ./QB/QB_neighbor.cpp ./QB/QB_network.cpp ./QB/QB_pbc.cpp ./QB/QB_Readlmc.cpp ./QB/QB_Readlmp.cpp ./QB/QB_Readvasp.cpp ./QB/QB_create.cpp

QSPG = ./QSPG/QSPG_link.cpp ./QSPG/QSPG_func.cpp

SPGLIB =./QSPG/spglib/arithmetic.cpp ./QSPG/spglib/cell.cpp ./QSPG/spglib/debug.cpp ./QSPG/spglib/delaunay.cpp ./QSPG/spglib/hall_symbol.cpp \
 ./QSPG/spglib/kgrid.cpp ./QSPG/spglib/kpoint.cpp ./QSPG/spglib/mathfunc.cpp ./QSPG/spglib/niggli.cpp ./QSPG/spglib/pointgroup.cpp \
./QSPG/spglib/primitive.cpp ./QSPG/spglib/refinement.cpp ./QSPG/spglib/site_symmetry.cpp ./QSPG/spglib/sitesym_database.cpp ./QSPG/spglib/spacegroup.cpp \
./QSPG/spglib/spg_database.cpp ./QSPG/spglib/spglib.cpp ./QSPG/spglib/spin.cpp ./QSPG/spglib/symmetry.cpp ./QSPG/spglib/determination.cpp ./QSPG/spglib/overlap.cpp 

 
SRC = $(QB) $(QSPG) $(SPGLIB) Energy_Minimization_main.cpp
OBJ = $(patsubst  %.cpp, %.o, $(SRC))
OBJ_D = $(patsubst  %.cpp, %.d, $(SRC))

#md: $(OBJ) 
#	g++ -static -o $@ $(OBJ)
all: $(OBJ)

%.o: %.cpp
	g++ -c -o $@ $<

clean:
	rm -f $(OBJ) md
	
