LFLAGS= -larmadillo -lblas -llapack
LFLAGS_UNI = -L/mn/felt/u9/sarahrei/General/Libraries/usr/include/ -lblas -llapack
CFLAGS_UNI =-O3 -I/mn/felt/u9/sarahrei/General/Libraries/usr/include/
CFLAGS =-O3
CC=mpicxx -DMPICH_IGNORE_CXX_SEEK
DEBUG=-Wall -g

all: main dfp main_DMC blocking_analyze blocking Integration density

main_DMC: main_DMC.o Walker.o DMC.o lib.o Hamiltonian.o Radial.o Hermite.o Wavefunction.o Slater.o Jastrow.o  QForce.o  ControlWalkers.o zigrandom.o zignor.o ziggurat.o ini.o normal.o
	$(CC) -o main_DMC main_DMC.o  DMC.o Walker.o Hamiltonian.o Radial.o Hermite.o Wavefunction.o Slater.o Jastrow.o ControlWalkers.o QForce.o lib.o zigrandom.o zignor.o ziggurat.o ini.o normal.o $(LFLAGS_UNI) 


main: main.o qvmc.o lib.o Hamiltonian.o Radial.o Hermite.o Wavefunction.o Slater.o Jastrow.o  QForce.o  zigrandom.o zignor.o ziggurat.o ini.o
	$(CC) -o main main.o  qvmc.o Hamiltonian.o Radial.o Hermite.o Wavefunction.o Slater.o Jastrow.o QForce.o lib.o zigrandom.o zignor.o ziggurat.o ini.o $(LFLAGS_UNI)

dfp: dfp.o qvmc.o lib.o Hamiltonian.o Radial.o Hermite.o Wavefunction.o Slater.o Jastrow.o QForce.o  zigrandom.o zignor.o ziggurat.o  ini.o 
	$(CC) -o dfp dfp.o  qvmc.o Hamiltonian.o Radial.o Hermite.o Wavefunction.o Slater.o Jastrow.o QForce.o lib.o zigrandom.o zignor.o ziggurat.o  ini.o  $(LFLAGS_UNI)
	
blocking_analyze: blocking_analyze.o
	$(CC) -o blocking_analyze blocking_analyze.o $(LFLAGS_UNI)

blocking: blocking.o qvmc.o lib.o Hamiltonian.o Wavefunction.o Radial.o Hermite.o Slater.o Jastrow.o QForce.o  zigrandom.o zignor.o ziggurat.o ini.o
	$(CC) -o blocking blocking.o  qvmc.o Hamiltonian.o Radial.o Hermite.o Wavefunction.o Slater.o Jastrow.o QForce.o lib.o zigrandom.o zignor.o ziggurat.o ini.o $(LFLAGS_UNI)

density: density.o
	$(CC) -o density density.o $(LFLAGS_UNI)

paircorDMC: paircorDMC.o
	$(CC) -o paircorDMC paircorDMC.o $(LFLAGS_UNI)
	
hist3d: hist3d.o
	$(CC) -o hist3d hist3d.o $(LFLAGS_UNI)
	
Integration: Integration.o qvmc.o lib.o Hamiltonian.o Radial.o Hermite.o Wavefunction.o Slater.o Jastrow.o QForce.o  zigrandom.o zignor.o ziggurat.o ini.o
	$(CC) -o Integration Integration.o  qvmc.o Hermite.o Radial.o Hamiltonian.o Wavefunction.o Slater.o Jastrow.o  QForce.o lib.o zigrandom.o zignor.o ziggurat.o ini.o $(LFLAGS_UNI)

main.o: main.cpp 	
	$(CC) -c  $(CFLAGS_UNI)  main.cpp 
	
main_DMC.o: main_DMC.cpp 	
	$(CC) -c  $(CFLAGS_UNI) main_DMC.cpp 

blocking.o: blocking.cpp 	
	$(CC) -c $(CFLAGS_UNI)  blocking.cpp  
	
hist3d.o: hist3d.cpp 	
	$(CC) -c $(CFLAGS_UNI)  hist3d.cpp  

blocking_analyze.o: blocking_analyze.cpp
	$(CC) -c $(CFLAGS_UNI) blocking_analyze.cpp

density.o: density.cpp
	$(CC) -c $(CFLAGS_UNI) density.cpp

paircorDMC.o: paircorDMC.cpp
	$(CC) -c $(CFLAGS_UNI) paircorDMC.cpp

Integration.o: Integration.cpp 	
	$(CC) -c $(CFLAGS_UNI)  Integration.cpp 

dfp.o: dfp.cpp 	
	$(CC) -c $(CFLAGS_UNI)  dfp.cpp     
	  
qvmc.o: qvmc.cpp	
	$(CC) -c $(CFLAGS_UNI) qvmc.cpp 
	
Hermite.o: Hermite.cpp
	$(CC) -c $(CFLAGS_UNI) Hermite.cpp 	
	
Radial.o: Radial.cpp
	$(CC) -c $(CFLAGS_UNI) Radial.cpp 		
	
Hamiltonian.o: Hamiltonian.cpp
	$(CC) -c $(CFLAGS_UNI) Hamiltonian.cpp  
	
Wavefunction.o: Wavefunction.cpp
	$(CC) -c $(CFLAGS_UNI) Wavefunction.cpp  
 
Slater.o: Slater.cpp
	$(CC) -c $(CFLAGS_UNI) Slater.cpp 

Jastrow.o: Jastrow.cpp
	$(CC) -c $(CFLAGS_UNI) Jastrow.cpp
		
QForce.o: QForce.cpp
	$(CC) -c $(CFLAGS_UNI) QForce.cpp 
	
Walker.o: Walker.cpp
	$(CC) -c $(CFLAGS_UNI) Walker.cpp 
	
DMC.o: DMC.cpp
	$(CC) -c $(CFLAGS_UNI) DMC.cpp 
	
ControlWalkers.o: ControlWalkers.cpp
	$(CC) -c $(CFLAGS_UNI) ControlWalkers.cpp 	
	

lib.o: lib.cpp
	$(CC) -c $(CFLAGS_UNI) lib.cpp
	
zignor.o: zignor.c
	$(CC) -c $(CFLAGS_UNI) zignor.c
	
zigrandom.o: zigrandom.c
	$(CC) -c $(CFLAGS_UNI) zigrandom.c
	
ziggurat.o: ziggurat.cpp
	$(CC) -c $(CFLAGS_UNI) ziggurat.cpp
	
ini.o: ini.cpp
	$(CC) -c $(CFLAGS_UNI) ini.cpp
	
normal.o: normal.cpp
	$(CC) -c $(CFLAGS_UNI) normal.cpp
		
clean:
	rm -rf *o main dfp main_DMC blocking_analyze blocking Integration *~
 
