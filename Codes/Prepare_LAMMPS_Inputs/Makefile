LINKFLAGS_FOR = -std=gnu -pedantic -ffree-form -march=native -O3
COMP_FOR = gfortran
install:
	${COMP_FOR} ${LINKFLAGS_FOR} cif2lammps.f90 -o cif2lammps.exe
all:
	make install
	make execute
	make clean
execute:
	./cif2lammps.exe -i Zn_144_mIm_288_SOD_cluster.cif --cif Zn_144_mIm_288_SOD_cluster.cif --no-charges --search-and-modify-topology --lammps
clean:;         @rm -f *.o *.mod cif2lammps.exe 
