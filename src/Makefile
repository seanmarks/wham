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

# Output program (LINK_TARGET)
PROJECT=wham

# Compiler
CXX= g++


OBJECTS= Bias.o \
         Bins.o \
         main.o \
         OrderParameter.o \
         TimeSeries.o \
         Wham.o \
         WhamDlibWrappers.o 


# DLIB
DLIB_DIR=$(HOME)/programs/dlib/19.4
DLIB_INCL=$(DLIB_DIR)/include
DLIB_LIB=$(DLIB_DIR)/lib

# Compiler flags 
#   Optimizations
CXXFLAGS_OPT= -O3 #-march=native
#   Warnings
CXXFLAGS_WARN= -Wall
CXXFLAGS_DEF= -DCPLUSPLUS
CXXFLAGS= $(CXXFLAGS_DEF) -g -std=c++11 $(CXXFLAGS_OPT) $(CXXFLAGS_WARN) -I $(DLIB_INCL)

# FIXME DEBUG
#CXXFLAGS += -DDEBUG_MODE


# Linking
LFLAGS= -lm -g -I $(DLIB_INCL) -L $(DLIB_LIB) -ldlib -lblas


# Make all
all : $(PROJECT)

# Linking
$(PROJECT) : $(OBJECTS)
	$(CXX) -o ../bin/$@ $^ $(LFLAGS)

# Compiling
%.o : %.cpp	
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -f $(OBJECTS) 
	rm -rf ../bin/$(PROJECT)*
