#########################################
# GWAStoolkit - Makefile
#########################################

CXX = g++
CXXFLAGS = -std=c++11 -O3 -fopenmp -Isrc -DUSE_RMATH
LDFLAGS  = -lz -lm

SRC = \
    src/main.cpp \
    src/cmds/cmd_rsidImpu.cpp \
    src/cmds/cmd_convert.cpp \
    src/cmds/cmd_or2beta.cpp \
    src/cmds/cmd_computeNeff.cpp \
    src/rsidImpu/dbsnp.cpp \
    src/rsidImpu/rsidImpu.cpp \
    src/rsidImpu/allele.cpp \
    src/convert/convert.cpp \
    src/or2beta/or2beta.cpp \
    src/computeNeff/computeNeff.cpp \
    src/utils/args.cpp \
    src/utils/FormatEngine.cpp \
    src/utils/gadgets.cpp \
    src/utils/gwasQC.cpp \
    src/utils/linereader.cpp \
    src/utils/log.cpp \
    src/utils/util.cpp \
    src/utils/writer.cpp \
    src/utils/StatFunc.cpp

OBJ = $(SRC:.cpp=.o)
TARGET = GWAStoolkit

#########################################
all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

#########################################
clean:
	rm -f $(OBJ) $(TARGET)

#########################################
.PHONY: all clean
