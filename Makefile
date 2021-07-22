# Makefile for VAssemble
#
# Written: April 1st 2012
# Last Update: April 1st 2012
# Author: Eric C. Dykeman

F90 = gfortran
FFLAGS = -O4

# ---------------------------------------------------------------------
# OBJECTS
# ---------------------------------------------------------------------

MAIN = vassembly.o

MODULES = systemvar.o

CLASSES = class_capsid.o

INOUT = readdata.o

MISC = ssareaction.o random.o

OBJECTS = $(MAIN) $(MODULES) $(CLASSES) $(INOUT) $(MISC)

# ---------------------------------------------------------------------
# MAKE COMMANDS
# ---------------------------------------------------------------------

all: vas.x clean

vas.x: $(OBJECTS)
	$(F90) $(FFLAGS) -o VAssembly.x $(OBJECTS)

clean:
	rm -f *.o *.mod

veryclean:
	rm -f *.o *.mod *.x

# ---------------------------------------------------------------------
# COMPILE COMMANDS
# ---------------------------------------------------------------------

# --------- MAIN ------------

vassembly.o: vassembly.f90 $(OBJECTS)
	$(F90) $(FFLAGS) -c vassembly.f90

# -------- MODULES ----------

systemvar.o: systemvar.f90
	$(F90) $(FFLAGS) -c systemvar.f90

# -------- CLASSES ----------

class_capsid.o: class_capsid.f90 $(MODULES)
	$(F90) $(FFLAGS) -c class_capsid.f90

# --------- INOUT -----------

readdata.o: readdata.f90 $(MODULES)
	$(F90) $(FFLAGS) -c readdata.f90
outputconc.o: outputconc.f90 $(MODULES)
	$(F90) $(FFLAGS) -c outputconc.f90

# --------- MISC -----------

ssareaction.o: ssareaction.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ssareaction.f90
random.o: random.f90
	$(F90) $(FFLAGS) -c random.f90
