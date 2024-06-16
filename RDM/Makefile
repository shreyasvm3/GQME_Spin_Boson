bindir = ../bin/
srcdir = ./src
# last edited 06/17/2018  M.S.Church
FC = ifort -O2 -g -traceback
LIB = -llapack -lblas
.SUFFIXES: .f90 .o

#--------------------------------------------------
MOD = params.f90  

SRC = RDM.f90

OBJ = $(MOD) $(SRC)

.f90.o:
	$(FC) -c $<
.f.o:
	$(FC) -c $<

dyn.x  :  $(OBJ)
	$(FC) $(OPTION) -o $@ $(OBJ) $(LIB)

clear   :
	rm *.o
	rm *.x
	rm *.out 
	rm *.mod
	rm fort.*
	rm mqcivr.*
