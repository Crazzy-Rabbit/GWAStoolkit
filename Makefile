#########################################
# GWAStoolkit - Makefile
#########################################

# Compiler
CXX = g++

# Flags

CXXFLAGS = -std=c++11 -O3 -fopenmp -Isrc -DUSE_RMATH
LDFLAGS  = -lz -lRmath -lm

# All source files in current directory

SRC = main.cpp \
        cmds/cmd_rsidImpu.cpp \
        cmds/cmd_convert.cpp \
        cmds/cmd_or2beta.cpp \
        rsidImpu/dbsnp.cpp \
        rsidImpu/gwas.cpp \
        rsidImpu/allele.cpp \
        convert/convert.cpp \
        or2beta/or2beta.cpp \
        utils/args.cpp \
        utils/FormatEngine.cpp \
        utils/gadgets.cpp \
        utils/gwasQC.cpp \
        utils/linereader.cpp \
        utils/log.cpp \
        utils/util.cpp \
        utils/writer.cpp

# Convert .cpp → .o
OBJ = $(SRC:.cpp=.o)

# Output binary
TARGET = rsidImpu

#########################################
all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

#########################################
clean:
	rm -f $(OBJ) $(TARGET)

#########################################
.PHONY: all clean