INC := -I./include/VASPMATE_include/
UNAME_S := $(shell uname -s)
# Check for Linux and Windows (assuming CYGWIN is part of Windows environment)
ifeq ($(UNAME_S),Linux)
    LDFLAGS = -static
else ifeq ($(findstring CYGWIN,$(UNAME_S)),CYGWIN)
    LDFLAGS =
else ifeq ($(findstring MINGW,$(UNAME_S)),MINGW)
    LDFLAGS =
else ifeq ($(UNAME_S),Darwin)
    LDFLAGS =
endif
CC :=  g++
CURRENT_DIR=$(shell pwd)
PREFIX=$(CURRENT_DIR)/fftw
FFTW_VERSION=3.3.10
FFTW_DIR=./fftw-$(FFTW_VERSION)
FFTW_LIB= ./fftw/lib
FFTW_INCLUDE= ./fftw/include
	
#FLAGS := -O2 -std=c++11 -I./spglib  -I./ -I./libmsym -L./lib  -I$(FFTW_INCLUDE) -L$(FFTW_LIB)
FLAGS := -O2 -std=c++0x $(LDFLAGS) -I./spglib -I./sqlite3 -I./ -I./libmsym -L./lib -I$(FFTW_INCLUDE) -L$(FFTW_LIB)
LFLAGS := -lfftw3 -lspg -lmsym -lstdc++ -lsqlite3 -lpthread -ldl

SOURCE := $(wildcard ./src/*.cpp) \
	$(wildcard ./src/VASPMATE_src/*.cpp) \
	$(wildcard ./src/Stochastics_src/*.cpp) \
	$(wildcard ./include/Energy_Minimization/*.cpp) \
	$(wildcard ./include/Energy_Minimization/QB/*.cpp) \
	$(wildcard ./include/Energy_Minimization/QSPG/*.cpp) \
	$(wildcard ./src/Evolutionary_src/*.cpp) \
	$(wildcard ./src/help_src/*.cpp) \
	$(wildcard ./src/MEP_src/*.cpp)
SOURCE_spglib := $(wildcard *.c)
SOURCE_msymlib := $(wildcard *.c)

OBJS := $(patsubst %.cpp, %.o, $(SOURCE))


all: VASPMATE

VASPMATE: check_fftw libspg msymlib sqlit3 $(OBJS)
	$(CC) -o $@ $(OBJS) $(FLAGS) $(LFLAGS)
	@echo "Compiling VASPMATE"
	mv VASPMATE ./bin/
	@echo "Compile successfully!"

%.o: %.cpp
	$(CC) -o $@ -c $< \
		$(INC) $(FLAGS)

#./src/main.o: ./src/main.cpp
#	$(CC) -o $@ -c $< \
#		$(INC) $(FLAGS)

#./src/VASPMATE_src/%.o: ./src/VASPMATE_src/%.cpp ./include/VASPMATE_include/*.h ./spglib/*.h
#	$(CC) -o $@ -c $< \
#		$(INC) $(FLAGS) \

#./src/Stochastics_src/%.o: ./src/Stochastics_src/%.cpp ./include/Stochastics_include/*.h ./spglib/*.h
#	$(CC) -o $@ -c $< \
#		$(INC) $(FLAGS) \

check_fftw:
	@if [ ! -d "$(FFTW_LIB)" ]; then \
	$(MAKE) build_fftw; \
	elif [ ! -d "$(FFTW_INCLUDE)" ]; then \
	$(MAKE) build_fftw; \
	fi

build_fftw: configure_fftw install_fftw compile_fftw

configure_fftw:
	cd $(FFTW_DIR) && chmod +x configure && chmod +x install-sh && ./configure --prefix=$(PREFIX) --disable-dependency-tracking
	@echo "FFTW configured"
	
install_fftw:
	cd $(FFTW_DIR) && make install
	@echo "FFTW installed"

compile_fftw:
	cd $(FFTW_DIR) && make
	@echo "FFTW compiled"
		
libspg: spglib/spglib.c
	@ echo "Compiling spglib"
	$(MAKE) -C spglib

msymlib: libmsym/msym.c
	@ echo "Compiling libmsym"
	$(MAKE) -C  libmsym
	
sqlit3: sqlite3/shell.c sqlite3/sqlite3.c
	@ echo "Compiling sqlite3"
	$(MAKE) -C  sqlite3

clean:
	@rm -rf $(OBJS) 

