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
#include <armadillo>

#define PI_NUM  3.14159265
#define DMC true

using namespace std;
using namespace arma;


////////////////////////////////////////////////////////////////////////////////
//                           The main program
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

    int numfiles = atoi(argv[1]);
    int numpart = atoi(argv[2]);
    
    ofstream ofile;
    ofile.open("density.out");

    int num_bins = 20;
    double r_max = 4;
    double binwidth = r_max/num_bins;
    
    vec pos(numpart);
    vec pos_xy(2*numpart);
    
    
    double green_branch;
    int weight, assign, counts = 0;
    int assign_x, assign_y;
    
    //vec density(num_bins);
    mat density_xy_up(2*num_bins, 2*num_bins);
   mat density_xy_down(2*num_bins, 2*num_bins);
   // density.zeros();
   density_xy_up.zeros();
   density_xy_down.zeros();
    

    for (int i = 0; i < numfiles; i++) {

        ostringstream ost;
        ost << "density" << i << ".dat";

        ifstream infile;
        infile.open(ost.str().c_str(), ios::in);


        if (!infile) { // file couldn't be opened
            cerr << "Error: file could not be opened" << endl;
            exit(1);
        }

        while (!infile.eof()){
            
            for(int k = 0; k< numpart; k++){
               // infile >> pos(k); // RADIAL DENSITY
                infile >> pos_xy(2*k); // XY-density
                infile >> pos_xy(2*k+1);
            }
                infile >> green_branch;
           
               for(int k = 0; k< numpart; k++){   
                    //assign = floor(pos(k)/binwidth); // RADIAL DENSITY
                    //density(assign) += green_branch;
                    assign_x = floor(pos_xy(2*k)/binwidth);
                    assign_y = floor(pos_xy(2*k+1)/binwidth);
                    if(k< numpart/2){
                        density_xy_up(assign_x+num_bins, assign_y+num_bins) += green_branch;
                    }
                    else{
                        density_xy_down(assign_x+num_bins, assign_y+num_bins) += green_branch;
                    }
                    
                    counts += green_branch;
                }
                    
        }
        
        
        infile.close();

    }

    
    /*
    // Normalize radial density
    for(int i = 0; i< num_bins; i++){
        density(i) /= counts*PI_NUM*binwidth*binwidth*(2*i+1);
        ofile << (i+0.5)*binwidth << " " << density(i) << endl;
    }*/

    // Normalize x-y-density
    for(int i = 0; i< 2*num_bins; i++){
        for(int j = 0; j< 2*num_bins; j++){
            density_xy_up(i,j) /= counts*binwidth;
            density_xy_down(i,j) /= counts*binwidth;
             ofile << (i+0.5)*binwidth << " " << (j+0.5)*binwidth << " " << density_xy_up(i,j) << " " << density_xy_down(i,j) << endl;
        }
    }
    
      
    ofile.close();
    return 0;
}

