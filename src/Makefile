
#Need to install GSL

#Need to install fftw3

#Need to install lapack and lablas for GSL??
#Check to see if I actually need lapack and lablas - doesn't look like I need GSL

#MAY need to install DarkSUSY

# Declaration of variables
CC = g++

CC_FLAGS = -O3 -Wall -ansi

#-I../cosmomc/source

#TARGET_DIR = /home/ppxbk2/astroraid/cosmomc/source
#TARGET_DIR = ../test

TARGET_GENERATE = ~/MultiNest_neutrinos
TARGET_RADON = ..

LIB_GSL = -lm -lstdc++ -lgsl -lgslcblas
PATH_LIB_GSL = /home/ppxbk2/gsl-1.9/lib
INC_GSL = /home/ppxbk2/gsl-1.9/include

INC_HEADERS = ../include

LIB_FFTW = -lfftw3
PATH_LIB_FFTW = /home/ppxbk2/fftw/lib
INC_FFTW = /home/ppxbk2/fftw/include

#LIB_LAPACK = -llapack -lblas
PATH_LIB_LAPACK = /home/ppxbk2/lapack-3.2

LIB_DS = -lgfortran -ldarksusy -lHB -lFH -lisajet
PATH_LIB_DS = /home/ppxbk2/DarkSUSY/darksusy-5.1.1/lib

LIBS =  -L$(PATH_LIB_GSL) -L$(PATH_LIB_FFTW) -L$(PATH_LIB_LAPACK)
LD_FLAGS = $(LIB_GSL) $(LIB_FFTW) $(LIB_LAPACK)
INCLUDE = -I$(INC_GSL) -I$(INC_HEADERS) -I$(INC_FFTW)



# File names
EXEC = run
SOURCES = $(*.cc)
OBJECTS = $(SOURCES:.cc=.o)

OBJ_FILES =  DMUtils.o Detector_Class.o Event_Class.o Parametrisations.o EventRates.o MatrixUtils.o RadonTransform.o Distributions.o Astrophysics_Class.o
HEADER_FILES = $(OBJ_FILES:.o=.h)
#Main Target

all: GenerateEvents GenerateEvents_IC CalcLikelihood CalcRadon Lib

#NB: Check ordering of library flags...

Lib: libdecodm.a

%.o: %.cc
	$(CC) $(CC_FLAGS) $(INCLUDE) $(LIBS) -c $^ $(LD_FLAGS)

Test:
	$(CC) $(CC_FLAGS) $(INCLUDE) $(LIBS) -o $@ $^ test.cc $(LD_FLAGS)

GenerateEvents: $(OBJ_FILES)
	$(CC) $(CC_FLAGS) $(INCLUDE) $(LIBS) -o $@ $^ GenerateEvents.cc $(LD_FLAGS)
	cp GenerateEvents $(TARGET_GENERATE)
	rm GenerateEvents

GenerateEvents_IC: $(OBJ_FILES)
	$(CC) $(CC_FLAGS) $(INCLUDE) $(LIBS) -L$(PATH_LIB_DS) -o $@ $^ GenerateEvents_IC.cc $(LD_FLAGS) $(LIB_DS)
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
	rm $(TARGET_GENERATE)GenerateEvents
	rm $(TARGET_RADON)CalcRadon
	rm ../libdecodm.a



#   $@ name of the target
#   $^ name of all prerequisites with duplicates removed
#   $< name of the first prerequisite
