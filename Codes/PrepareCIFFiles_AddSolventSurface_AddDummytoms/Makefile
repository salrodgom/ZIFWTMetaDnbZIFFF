LINKFLAGS_FOR = -std=gnu -pedantic -ffree-form -march=native -O3
COMP_FOR = gfortran
install:
	${COMP_FOR} ${LINKFLAGS_FOR} pdb2cif.f90 -o pdb2cif.exe
all:
	make install
	make execute
	make clean
execute:
	./pdb2cif.exe < Zn_144_mIm_288_SOD_cluster.pdb > pdb2cif.log
	mv p1_fromPDBFile.cif Zn_144_mIm_288_SOD_cluster_Dummies.cif
clean:;         @rm -f *.o *.mod pdb2cif.exe 
