INC := -I./include/ 
CC :=  g++
FFTW_LIB=
FFTW_INCLUDE=
FLAGS := -O2 -std=c++0x -I./spglib  -I./ -I./libmsym -L./lib -I$(FFTW_INCLUDE) -L$(FFTW_LIB)
LFLAGS := -lfftw3 -lspg -lmsym

SOURCE := $(wildcard ./src/*.cpp)

SOURCE_spglib := $(wildcard *.c)
SOURCE_msymlib := $(wildcard *.c)

OBJS := $(patsubst %.cpp, %.o, $(SOURCE))
ALL: libspg msymlib VASPMATE 

VASPMATE: $(OBJS)
	$(CC) -o $@ $(OBJS) $(FLAGS) $(LFLAGS)
	@ echo "Compiling VASPMATE"
	mv VASPMATE ./bin/
	@ echo "Compile successfully!"

./src/%.o: ./src/%.cpp ./include/*.h ./spglib/*.h
	$(CC) -o $@ -c $< \
		$(INC) $(FLAGS) \
		
libspg: spglib/spglib.c
	@ echo "Compiling spglib"
	$(MAKE) -C spglib

msymlib: libmsym/msym.c
	@ echo "Compiling libmsym"
	$(MAKE) -C  libmsym

clean:
	@rm -rf $(OBJS) 

