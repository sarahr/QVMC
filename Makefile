# First line: university account, second line: home laptop

all: main
#all: dfp	
#all: blocking_analyze
#all: blocking
#all: Integration

main: main.o qvmc.o lib.o Hamiltonian.o Radial.o Hermite.o Wavefunction.o Slater.o Jastrow.o  QForce.o  ExpFactor.o zigrandom.o zignor.o ziggurat.o ini.o
	#mpicxx -o main main.o  qvmc.o Hamiltonian.o Radial.o Hermite.o Wavefunction.o Slater.o Jastrow.o ExpFactor.o QForce.o lib.o zigrandom.o zignor.o ziggurat.o  ini.o -L/mn/felt/u9/sarahrei/General/Libraries/usr/include/ -lblas -llapack
	mpicxx -o main main.o  qvmc.o Hamiltonian.o Radial.o Hermite.o Wavefunction.o Slater.o Jastrow.o ExpFactor.o QForce.o lib.o zigrandom.o zignor.o ziggurat.o ini.o -larmadillo -lblas -llapack

dfp: dfp.o qvmc.o lib.o Hamiltonian.o Radial.o Hermite.o Wavefunction.o Slater.o Jastrow.o ExpFactor.o QForce.o  zigrandom.o zignor.o ziggurat.o  ini.o 
	#mpicxx -o dfp dfp.o  qvmc.o Hamiltonian.o Radial.o Hermite.o Wavefunction.o Slater.o Jastrow.o ExpFactor.o QForce.o lib.o zigrandom.o zignor.o ziggurat.o  ini.o  -L/mn/felt/u9/sarahrei/General/Libraries/usr/include/ -lblas -llapack
	mpicxx -o dfp dfp.o  qvmc.o Hamiltonian.o Radial.o Hermite.o Wavefunction.o Slater.o Jastrow.o ExpFactor.o QForce.o lib.o zigrandom.o zignor.o ziggurat.o  ini.o  -larmadillo -lblas -llapack
	
blocking_analyze: blocking_analyze.o
	g++ -o blocking_analyze blocking_analyze.o -L/mn/felt/u9/sarahrei/General/Libraries/usr/include/ -lblas -llapack

blocking: blocking.o qvmc.o lib.o Hamiltonian.o Wavefunction.o Radial.o Hermite.o Slater.o Jastrow.o ExpFactor.o QForce.o  zigrandom.o zignor.o ziggurat.o ini.o
	#mpicxx -o blocking blocking.o  qvmc.o Hamiltonian.o Radial.o Wavefunction.o Slater.o Jastrow.o ExpFactor.o  QForce.o lib.o zigrandom.o zignor.o ziggurat.o ini.o -L/mn/felt/u9/sarahrei/General/Libraries/usr/include/ -lblas -llapack
	mpicxx -o blocking blocking.o  qvmc.o Hamiltonian.o Radial.o Hermite.o Wavefunction.o Slater.o Jastrow.o ExpFactor.o  QForce.o lib.o zigrandom.o zignor.o ziggurat.o ini.o -larmadillo -lblas -llapack
	
Integration: Integration.o qvmc.o lib.o Hamiltonian.o Radial.o Hermite.o Wavefunction.o Slater.o Jastrow.o ExpFactor.o QForce.o  zigrandom.o zignor.o ziggurat.o ini.o
	mpicxx -o Integration Integration.o  qvmc.o Hermite.o Radial.o Hamiltonian.o Wavefunction.o Slater.o Jastrow.o ExpFactor.o QForce.o lib.o zigrandom.o zignor.o ziggurat.o ini.o -L/mn/felt/u9/sarahrei/General/Libraries/usr/include/ -lblas -llapack

main.o: main.cpp 
#	mpicxx -c  -Wall -O3 -I/mn/felt/u9/sarahrei/General/Libraries/usr/include/ main.cpp	
	mpicxx -c  -Wall -O3  main.cpp 

blocking.o: blocking.cpp 
#	mpicxx -c  -Wall -O3 -I/mn/felt/u9/sarahrei/General/Libraries/usr/include/ blocking.cpp	
	mpicxx -c  -Wall -O3  blocking.cpp  

blocking_analyze.o: blocking_analyze.cpp
	g++ -c  -Wall -O3 -I/mn/felt/u9/sarahrei/General/Libraries/usr/include/ blocking_analyze.cpp

Integration.o: Integration.cpp 
#	mpicxx -c  -Wall -O3 -I/mn/felt/u9/sarahrei/General/Libraries/usr/include/ Integration.cpp	
	mpicxx -c  -Wall -O3  Integration.cpp 


dfp.o: dfp.cpp 
#	mpicxx -c  -Wall -O3 -I/mn/felt/u9/sarahrei/General/Libraries/usr/include/ dfp.cpp	
	mpicxx -c  -Wall -O3  dfp.cpp     
	  
qvmc.o: qvmc.cpp	 
#	g++ -c  -Wall -O3 qvmc.cpp -I/mn/felt/u9/sarahrei/General/Libraries/usr/include/	
	g++ -c  -Wall -O3 qvmc.cpp 
	
Hermite.o: Hermite.cpp
#	g++ -c  -Wall -O3 -I/mn/felt/u9/sarahrei/General/Libraries/usr/include/ Hermite.cpp
	g++ -c  -Wall -O3  Hermite.cpp 	
	
Radial.o: Radial.cpp
#	g++ -c  -Wall -O3 -I/mn/felt/u9/sarahrei/General/Libraries/usr/include/ Radial.cpp
	g++ -c  -Wall -O3  Radial.cpp 		
	
Hamiltonian.o: Hamiltonian.cpp
#	g++ -c  -Wall -O3 -I/mn/felt/u9/sarahrei/General/Libraries/usr/include/ Hamiltonian.cpp
	g++ -c  -Wall -O3  Hamiltonian.cpp  
	
Wavefunction.o: Wavefunction.cpp
#	g++ -c  -Wall -O3 -I/mn/felt/u9/sarahrei/General/Libraries/usr/include/ Wavefunction.cpp
	g++ -c  -Wall -O3  Wavefunction.cpp  
 
Slater.o: Slater.cpp
#	g++ -c  -Wall -O3 -I/mn/felt/u9/sarahrei/General/Libraries/usr/include/ Slater.cpp
	g++ -c  -Wall -O3  Slater.cpp 

Jastrow.o: Jastrow.cpp
#	g++ -c  -Wall -O3 -I/mn/felt/u9/sarahrei/General/Libraries/usr/include/ Jastrow.cpp
	g++ -c  -Wall -O3  Jastrow.cpp
	
ExpFactor.o: ExpFactor.cpp
#	g++ -c  -Wall -O3 -I/mn/felt/u9/sarahrei/General/Libraries/usr/include/ ExpFactor.cpp
	g++ -c  -Wall -O3  ExpFactor.cpp
 	
QForce.o: QForce.cpp
#	g++ -c  -Wall -O3 -I/mn/felt/u9/sarahrei/General/Libraries/usr/include/ QForce.cpp
	g++ -c  -Wall -O3  QForce.cpp 

lib.o: lib.cpp
	g++ -c -Wall -O3 lib.cpp
	
zignor.o: zignor.c
	g++ -c -Wall -O3 zignor.c
	
zigrandom.o: zigrandom.c
	g++ -c -Wall -O3 zigrandom.c
	
ziggurat.o: ziggurat.cpp
	g++ -c -Wall -O3 ziggurat.cpp
	
ini.o: ini.cpp
	g++ -c -Wall ini.cpp
		
clean:
	rm -rf *o main Integration blocking blocking_analyze dfp
 
