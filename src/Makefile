
#Need to install GSL

#Need to install fftw3

#Need to install lapack and lablas for GSL??
#Check to see if I actually need lapack and lablas - doesn't look like I need GSL

#MAY need to install DarkSUSY

# Declaration of variables
FC = gfortran
F_FLAGS = -O3

C = gcc
C_FLAGS = -O3 -w

CC = g++
CC_FLAGS = -O3 -Wall -ansi

#-I../cosmomc/source

#TARGET_DIR = /home/ppxbk2/astroraid/cosmomc/source
#TARGET_DIR = ../test

TARGET_TEST = ../test
TARGET_GENERATE = ..
TARGET_RADON = ..

LIB_GSL = -lm -lstdc++ -lgsl -lgslcblas
PATH_LIB_GSL = /home/bkavanag/gsl-1.9/lib
INC_GSL = /home/bkavanag/gsl-1.9/include

INC_HEADERS = ../include

LIB_FFTW = -lfftw3
PATH_LIB_FFTW = /home/bkavanag/fftw/lib
INC_FFTW = /home/bkavanag/fftw/include

#LIB_LAPACK = -llapack -lblas
PATH_LIB_LAPACK = /home/bkavanag/lapack-3.2

LIB_DS = -ldarksusy -lgfortran -lHB -lFH -lisajet
PATH_LIB_DS = /home/bkavanag/darksusy-5.1.1/lib

LIBS =  -L$(PATH_LIB_GSL) -L$(PATH_LIB_FFTW) -L$(PATH_LIB_LAPACK)
LD_FLAGS = $(LIB_GSL) $(LIB_FFTW) $(LIB_LAPACK)
INCLUDE = -I$(INC_GSL) -I$(INC_HEADERS) -I$(INC_FFTW)



# File names
EXEC = run
SOURCES = $(*.cc)
OBJECTS = $(SOURCES:.cc=.o)

OBJ_FILES =  DMUtils.o Detector_Class.o Event_Class.o Parametrisations.o EventRates.o MatrixUtils.o Distributions.o Astrophysics_Class.o RadonTransform.o

HEADER_FILES = $(OBJ_FILES:.o=.h)
#Main Target

all: GenerateEvents GenerateEvents_IC CalcLikelihood CalcRadon ChainAnalysis ChainAnalysis2 Lib

#NB: Check ordering of library flags...

Lib: libdecodm.a

polyf.o: polyf.c
	$(C) $(C_FLAGS) -c $^ $(LD_FLAGS)

%.o: %.cc
	$(CC) $(CC_FLAGS) $(INCLUDE) $(LIBS) -c $^ $(LD_FLAGS)

Test:
	$(CC) $(CC_FLAGS) $(INCLUDE) $(LIBS) -o $@ $^ test.cc $(LD_FLAGS)

ChainAnalysis: ChainAnalysis.o
	$(CC) $(CC_FLAGS) $(INCLUDE) $(LIBS) -o $@ $^ $(LD_FLAGS)
	cp ChainAnalysis $(TARGET_GENERATE)

ChainAnalysis2: ChainAnalysis2.c CalcLikelihood.cc $(OBJ_FILES) polyf.o
	$(FC) $(F_FLAGS) $(INCLUDE) $(LIBS) -L$(PATH_LIB_DS) -o $@ $^ $(LD_FLAGS) $(LIB_DS)
	cp ChainAnalysis2 $(TARGET_GENERATE)

GenerateEvents: $(OBJ_FILES)
	$(CC) $(CC_FLAGS) $(INCLUDE) $(LIBS) -o $@ $^ GenerateEvents.cc $(LD_FLAGS)
	cp GenerateEvents $(TARGET_GENERATE)
	cp GenerateEvents $(TARGET_TEST)
	rm GenerateEvents

GenerateEvents_IC: $(OBJ_FILES)
	$(CC) $(CC_FLAGS) $(INCLUDE) $(LIBS) -L$(PATH_LIB_DS) -o $@ $^ GenerateEvents_IC.cc $(LD_FLAGS) $(LIB_DS)
	cp GenerateEvents_IC $(TARGET_GENERATE)
	cp GenerateEvents_IC $(TARGET_GENERATE)
	rm GenerateEvents_IC

CalcRadon: CalcRadon.o $(OBJ_FILES)
	$(CC) $(CC_FLAGS) $(INCLUDE) $(LIBS)  -o $@ $^ $(LD_FLAGS)
	cp CalcRadon $(TARGET_RADON)
	rm CalcRadon

CalcLikelihood: CalcLikelihood.cc $(OBJ_FILES)
	$(CC) $(CC_FLAGS) $(INCLUDE) $(LIBS)  -c $< $(LD_FLAGS)

libdecodm.a: CalcLikelihood $(OBJ_FILES)
	ar r ../$@ CalcLikelihood.o $(OBJ_FILES)

clean:
	rm *.o

uninstall:
	rm $(TARGET_GENERATE)/GenerateEvents
	rm $(TARGET_GENERATE)/GenerateEvents_IC
	rm $(TARGET_RADON)/CalcRadon
	rm $(TARGET_TEST)/GenerateEvents
	rm $(TARGET_TEST)/GenerateEvents_IC
	rm ../libdecodm.a



#   $@ name of the target
#   $^ name of all prerequisites with duplicates removed
#   $< name of the first prerequisite
