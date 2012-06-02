LFLAGS= -larmadillo -lblas -llapack
LFAGS_UNI = -L/mn/felt/u9/sarahrei/General/Libraries/usr/include/ -lblas -llapack
CFLAGS_UNI =-O3 -I/mn/felt/u9/sarahrei/General/Libraries/usr/include/
CFLAGS =-O3
CC=mpicxx
DEBUG=-Wall -g

all: main dfp main_DMC blocking_analyze blocking Integration density

main_DMC: main_DMC.o Walker.o DMC.o lib.o Hamiltonian.o Radial.o Hermite.o Wavefunction.o Slater.o Jastrow.o  QForce.o  ControlWalkers.o zigrandom.o zignor.o ziggurat.o ini.o normal.o
	$(CC) -o main_DMC main_DMC.o  DMC.o Walker.o Hamiltonian.o Radial.o Hermite.o Wavefunction.o Slater.o Jastrow.o ControlWalkers.o QForce.o lib.o zigrandom.o zignor.o ziggurat.o ini.o normal.o $(LFLAGS) 


main: main.o qvmc.o lib.o Hamiltonian.o Radial.o Hermite.o Wavefunction.o Slater.o Jastrow.o  QForce.o  zigrandom.o zignor.o ziggurat.o ini.o
	$(CC) -o main main.o  qvmc.o Hamiltonian.o Radial.o Hermite.o Wavefunction.o Slater.o Jastrow.o QForce.o lib.o zigrandom.o zignor.o ziggurat.o ini.o $(LFLAGS)

dfp: dfp.o qvmc.o lib.o Hamiltonian.o Radial.o Hermite.o Wavefunction.o Slater.o Jastrow.o QForce.o  zigrandom.o zignor.o ziggurat.o  ini.o 
	$(CC) -o dfp dfp.o  qvmc.o Hamiltonian.o Radial.o Hermite.o Wavefunction.o Slater.o Jastrow.o QForce.o lib.o zigrandom.o zignor.o ziggurat.o  ini.o  $(LFLAGS)
	
blocking_analyze: blocking_analyze.o
	$(CC) -o blocking_analyze blocking_analyze.o $(LFLAGS)

blocking: blocking.o qvmc.o lib.o Hamiltonian.o Wavefunction.o Radial.o Hermite.o Slater.o Jastrow.o QForce.o  zigrandom.o zignor.o ziggurat.o ini.o
	$(CC) -o blocking blocking.o  qvmc.o Hamiltonian.o Radial.o Hermite.o Wavefunction.o Slater.o Jastrow.o QForce.o lib.o zigrandom.o zignor.o ziggurat.o ini.o $(LFLAGS)
	
Integration: Integration.o qvmc.o lib.o Hamiltonian.o Radial.o Hermite.o Wavefunction.o Slater.o Jastrow.o QForce.o  zigrandom.o zignor.o ziggurat.o ini.o
	$(CC) -o Integration Integration.o  qvmc.o Hermite.o Radial.o Hamiltonian.o Wavefunction.o Slater.o Jastrow.o  QForce.o lib.o zigrandom.o zignor.o ziggurat.o ini.o $(LFLAGS)

density: density.o
	$(CC) -o density density.o $(LFLAGS)

main.o: main.cpp 	
	$(CC) -c  $(CFLAGS)  main.cpp 
	
main_DMC.o: main_DMC.cpp 	
	$(CC) -c  $(CFLAGS) main_DMC.cpp 

blocking.o: blocking.cpp 	
	$(CC) -c $(CFLAGS)  blocking.cpp  

blocking_analyze.o: blocking_analyze.cpp
	$(CC) -c $(CFLAGS) blocking_analyze.cpp
	
density.o: density.cpp
	$(CC) -c $(CFLAGS) density.cpp	

Integration.o: Integration.cpp 	
	$(CC) -c $(CFLAGS)  Integration.cpp 

dfp.o: dfp.cpp 	
	$(CC) -c $(CFLAGS)  dfp.cpp     
	  
qvmc.o: qvmc.cpp	
	$(CC) -c $(CFLAGS) qvmc.cpp 
	
Hermite.o: Hermite.cpp
	$(CC) -c $(CFLAGS) Hermite.cpp 	
	
Radial.o: Radial.cpp
	$(CC) -c $(CFLAGS) Radial.cpp 		
	
Hamiltonian.o: Hamiltonian.cpp
	$(CC) -c $(CFLAGS) Hamiltonian.cpp  
	
Wavefunction.o: Wavefunction.cpp
	$(CC) -c $(CFLAGS) Wavefunction.cpp  
 
Slater.o: Slater.cpp
	$(CC) -c $(CFLAGS) Slater.cpp 

Jastrow.o: Jastrow.cpp
	$(CC) -c $(CFLAGS) Jastrow.cpp
		
QForce.o: QForce.cpp
	$(CC) -c $(CFLAGS) QForce.cpp 
	
Walker.o: Walker.cpp
	$(CC) -c $(CFLAGS) Walker.cpp 
	
DMC.o: DMC.cpp
	$(CC) -c $(CFLAGS) DMC.cpp 
	
ControlWalkers.o: ControlWalkers.cpp
	$(CC) -c $(CFLAGS) ControlWalkers.cpp 	
	

lib.o: lib.cpp
	$(CC) -c $(CFLAGS) lib.cpp
	
zignor.o: zignor.c
	$(CC) -c $(CFLAGS) zignor.c
	
zigrandom.o: zigrandom.c
	$(CC) -c $(CFLAGS) zigrandom.c
	
ziggurat.o: ziggurat.cpp
	$(CC) -c $(CFLAGS) ziggurat.cpp
	
ini.o: ini.cpp
	$(CC) -c $(CFLAGS) ini.cpp
	
normal.o: normal.cpp
	$(CC) -c $(CFLAGS) normal.cpp
		
clean:
	rm -rf *o main dfp main_DMC blocking_analyze blocking Integration *~
 
