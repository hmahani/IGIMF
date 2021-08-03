# Start of the makefile
# Defining variables
objects = global.o main.o beginning.o singlemass.o binary.o equivalentmass.o igimf.o mremnant.o IGIMFremnant.o fit.o
#f90comp = ifort
f90comp =gfortran
# Makefile
exemf: $(objects)
	$(f90comp) -o exemf $(objects)
utilmf.mod: global.o global.f90
	 $(f90comp) -c global.f90
global.o: global.f90
	 $(f90comp) -c global.f90
main.o: utilmf.mod main.f90
	 $(f90comp) -c main.f90
beginning.o: utilmf.mod beginning.f90
	 $(f90comp) -c beginning.f90
singlemass.o: utilmf.mod singlemass.f90
	 $(f90comp) -c singlemass.f90
binary.o: utilmf.mod binary.f90
	 $(f90comp) -c binary.f90
equivalentmass.o: utilmf.mod equivalentmass.f90
	 $(f90comp) -c equivalentmass.f90
igimf.o: utilmf.mod igimf.f90
	 $(f90comp) -c igimf.f90
mremnant.o: utilmf.mod mremnant.f90
	 $(f90comp) -c mremnant.f90
IGIMFremnant.o: utilmf.mod IGIMFremnant.f90
	 $(f90comp) -c IGIMFremnant.f90
fit.o: utilmf.mod fit.f90
	 $(f90comp) -c fit.f90

# Cleaning everything
clean:
	 rm utilmf.mod exemf
	 rm $(objects)
	 rm *~ *.dat *.log fort.* 
# End of the makefile
