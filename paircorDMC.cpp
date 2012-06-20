/*
 * Obtain the probability density (Psi^2) of the wavefunction (Psi)
 * Input files contain the radial positions of the particles
 * 
 * Input arguments: 
 * argv[1]: # number of input files
 * argv[2]: # number of particles
 */

#include <cstdlib>
#include <sys/stat.h>
#include</mn/felt/u9/sarahrei/General/Libraries/usr/include/armadillo>
//#include <armadillo>

#define PI_NUM  3.14159265

using namespace std;
using namespace arma;


////////////////////////////////////////////////////////////////////////////////
//                           The main program
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

    int numfiles = atoi(argv[1]);
    int numpart = atoi(argv[2]);
    int n2 = numpart/2;
    
    ofstream ofile;
    ofile.open("paircor.out");
    
    double green_branch;
    int weight, assign, counts = 0;
    int num_bins = 200;
    vec gcor_eq(num_bins);
    vec gcor_uneq(num_bins);
    gcor_eq.zeros();
    gcor_uneq.zeros();
    
    double rmax = 100;
    double boxwidth = rmax/num_bins;
    
    vec paircor(num_bins);
    paircor.zeros();
    double dist;
    int counter = 0;
   
    for (int i = 0; i < numfiles; i++) {

        ostringstream ost;
        ost << "paircorVMC" << i << ".dat";
        cout << i << endl;

        ifstream infile;
        infile.open(ost.str().c_str(), ios::in);


        if (!infile) { // file couldn't be opened
            cerr << "Error: file could not be opened" << endl;
            exit(1);
        }

       while (!infile.eof()) {
         // for(int o = 0; o< 200000; o++){
            counter++;
  
            int col = numpart - 1;
            int count = 1;

            while (col > 0) {
                if (count <= n2) {
                    for (int k = 0; k < n2 - count; k++) {
                        infile >> dist;
                        assign = floor(dist / boxwidth);
                        if (assign < num_bins && counter>1000)
                            gcor_eq(assign) += 1;
                    }
                    
                    for (int k = 0; k < n2; k++) {
                        infile >> dist;
                        assign = floor(dist / boxwidth);
                        if (assign < num_bins && counter>1000)
                            gcor_uneq(assign) += 1;
                    }
                                 
                } // END LOOP count <= n2
                
                else{
                    
                    for (int k = 0; k < numpart - count; k++) {
                        infile >> dist;
                        assign = floor(dist / boxwidth);
                        if (assign < num_bins && counter>1000)
                            gcor_eq(assign) += 1;
                    }
                }
                
                count++;  
                col--;
                  
            }
             infile >> green_branch;

        }
        
        
        infile.close();

    }
    
    
    double histnorm = sum(gcor_eq);
    gcor_eq /= histnorm*boxwidth;
    
    histnorm = sum(gcor_uneq);
    gcor_uneq /= histnorm*boxwidth;
    
       
    for(int i = 0; i< num_bins; i++){
        ofile << (i+0.5)*boxwidth << " " << gcor_eq(i) << " " << gcor_uneq(i) << " " << endl;
    }
    
      
    ofile.close();
    return 0;
}


