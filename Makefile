# NOTES:
# -g	Extra symbolic debugging info for use with the gdb debugger
# -o	Specify output file name
# -c	Compile into object file
# -Wall	Print all warning messages
#
# $@	Name of the target
# $<	First dependency of the target
# $^	All dependencies of the target
#
# Linking occurs when the object files (.o) are put together into an executable (.exe)

INSTALL_DIR=$(HOME)/programs/wham
DLIB_DIR=$(HOME)/programs/dlib/19.4

# Compilers
CXX=g++

# Output binary's name
PROJECT=wham


### !!! Shouldn't have to change anything below this line !!! ###


#############
### Flags ###
#############

# DLIB
DLIB_LIB     := $(DLIB_DIR)/lib
DLIB_INCLUDE := $(DLIB_DIR)/include

# Compiler flags
# - Misc: -Wno-comment -Wno-sign-compare -DPLUMED_MODE
CXXFLAGS += -g -std=c++11 -DCPLUSPLUS -I$(DLIB_INCLUDE) -Wall
# - Optimizations
CXXFLAGS += -O3
# CXXFLAGS += -ffast-math -march=native    # less safe optimizations
# - Debug mode
CXXFLAGS += -DDEBUG

# Check for whether MPI is enabled
ifeq ($(is_mpi_enabled),1)
CXXFLAGS += -DMPI_ENABLED
endif

# Linking
LDFLAGS += -g -L$(DLIB_LIB)
LIBS    += -ldlib -lblas


#############################
### Directories and Files ###
#############################

START_DIR := $(PWD)

SRC_DIR         := src
BUILD_DIR       := build
BUILD_BIN_DIR   := $(BUILD_DIR)/bin
INSTALL_BIN_DIR := $(INSTALL_DIR)/bin

# Get source and target objects
SRC_FILES = $(shell find $(SRC_DIR) -name '*.cpp')
SRC_DIRS  = $(shell find $(SRC_DIR) -type d | sed 's/$(SRC_DIR)/./g' )
OBJECTS   = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRC_FILES))


####################
### Main Program ###
####################

all : buildrepo $(PROJECT)

# Make directories 
buildrepo :
	@{ \
	for dir in $(BUILD_DIR) $(BUILD_BIN_DIR) $(INSTALL_DIR)/bin ; do \
		if [ ! -d $$dir ]; then \
			echo "Making directory $$dir ..." ;\
			mkdir -p $$dir ;\
		fi ;\
	done ;\
	}

# Compile
$(BUILD_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CXX) -o $@ $(CXXFLAGS) -c $<

# Link project
$(PROJECT) : $(OBJECTS)
	$(CXX) -o $(BUILD_BIN_DIR)/$@ $(LDFLAGS) $(OBJECTS) $(LIBS)

.PHONY : install
install :
	@{ \
	if [ ! -d $(INSTALL_BIN_DIR) ]; then \
		echo "Creating $(INSTALL_BIN_DIR) ..." ;\
		mkdir -p $(INSTALL_BIN_DIR) ;\
	fi ;\
	echo "Installing at $(INSTALL_DIR) ..." ;\
	cp $(BUILD_BIN_DIR)/* $(INSTALL_BIN_DIR) ;\
	echo "Done" ;\
	}


###############
### Cleanup ###
###############

# Clean up build directory
.PHONY : clean	
clean :
	rm -rf $(BUILD_DIR)/*

# Clean up install directory
.PHONY : clean_install
clean_install :
	rm -rf $(INSTALL_DIR)/*


########################
### Regression Tests ###
########################

.PHONY: test
test:
	@{ \
	cd test ;\
	./run_tests.sh $(START_DIR)/$(BUILD_BIN_DIR)/$(PROJECT) $(echo_failed_diffs) ;\
	}

