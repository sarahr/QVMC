/* 
 * File:   hist3d.cpp
 * Author: sa_rei
 *
 * Created on June 15, 2012, 6:48 PM
 */

#include <cstdlib>
#include <armadillo>

using namespace std;
using namespace arma;


////////////////////////////////////////////////////////////////////////////////
//                           The main program
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

    int numfiles = atoi(argv[1]);
    
    ofstream ofile;
    ofile.open("hist3d.out");

    double posx, posy;
    int num_bins = 100;
    double r_max = 80;
    double binwidth = r_max/num_bins;
    
    int assignx, assigny, counts = 0;
    
    double** density = new double*[2*num_bins];
    
    for(int i = 0; i< 2*num_bins; i++)
        density[i] = new double[2*num_bins];
        
    for(int i = 0; i< 2*num_bins; i++)
        for(int j=0; j<2*num_bins; j++)
        density[i][j] = 0;

    for (int i = 0; i < numfiles; i++) {

        ostringstream ost;
        ost << "positionDMC" << i << ".dat";

        ifstream infile;
        infile.open(ost.str().c_str(), ios::in);


        if (!infile) { // file couldn't be opened
            cerr << "Error: file could not be opened" << endl;
            exit(1);
        }

        while (!infile.eof()) {
          //for(int o=0;o<35000; o++){ 
                infile >> posx;
                infile >> posy;

                assignx = floor(posx/binwidth)+ num_bins;
                assigny = floor(posy/binwidth)+ num_bins;
                
                if(assignx < 2*num_bins && assigny < 2*num_bins){
                density[assignx][assigny]+= 1;
                counts += 1;
                }
                }
                          
        infile.close();

    }

    for(int i = 0; i< 2*num_bins; i++){
        for(int j=0; j<2*num_bins; j++){
            density[i][j] /= counts;
               ofile << (i+0.5)*binwidth-r_max << " " << (j+0.5)*binwidth-r_max<< " " << density[i][j] << endl;
        }
    }
    
   
    for(int i = 0; i<2*num_bins; i++)
        delete [] density[i];
    delete [] density;
    
    ofile.close();
    return 0;
}


