/*
 * Performance of blocking with the output data from "blocking"
 * 
 * Input arguments: 
 * argv[1]: # processors used when executing "blocking"
 * argv[2]: # MC cycles used in "blocking"
 */

#include <cstdlib>
#include <sys/stat.h>
#include <armadillo>

using namespace std;
using namespace arma;

///////////////////////////////////////////////////////////////////////////////
//                          Function declarations
///////////////////////////////////////////////////////////////////////////////

void blocking(mat& mc_results, int n, int block_size, vec& res);
double get_mean(int i, int block_size, mat& mc_results);
void get_meanvar(vec& values, int n_blocks, vec& res);

////////////////////////////////////////////////////////////////////////////////
//                           The main program
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

    int numprocs = atoi(argv[1]);
    int n = atoi(argv[2]);
    int local_n;
    ofstream ofile;
    ofile.open("blocking.out");

    // Variables to be adjusted
    int min_blocksize = 2;
    int max_blocksize = 5000;
    int block_samples = 200;
    int block_size;

    int step_length = (max_blocksize - min_blocksize) / block_samples;

    local_n = n / numprocs;

    vec mc_results(n);
    vec res(2); // contains mean and variance

    // Data input
    for (int i = 0; i < numprocs; i++) {

        ostringstream ost;
        ost << "blocks_rank" << i << ".dat";

        ifstream infile;
        infile.open(ost.str().c_str(), ios::in | ios::binary);


        if (!infile) { // file couldn't be opened
            cerr << "Error: file could not be opened" << endl;
            exit(1);
        }

        for (int j = 0; j < local_n; j++) {
            infile >> mc_results(i * local_n + j);
        }

        infile.close();

    }

    // Loop over all block_sizes

    for (int i = 0; i < block_samples; i++) {

        block_size = min_blocksize + i*step_length;
        blocking(mc_results, n, block_size, res);

        ofile << block_size << "\t" << res(0) << "\t" << sqrt(res(1) / 
                (n / block_size)) << endl;

    }


    ofile.close();
    return 0;
}


////////////////////////////////////////////////////////////////////////////////
//                             Functions
////////////////////////////////////////////////////////////////////////////////


/**
 * Compute for all block sizes the vairance in the QVMC algorithm
 * @param mc_results - vector containing the output from program "blocking"
 * @param n - total number of MC cycles
 * @param block_size - size of the blocks
 * @param res - res(0): mean energy, res(1): variance
 */
void blocking(mat& mc_results, int n, int block_size, vec& res){
    
    int n_blocks = n/block_size;
    vec values(n_blocks);
    
    for(int i = 0; i< n_blocks; i++){
        values(i) = get_mean(i, block_size, mc_results);
    }
    
    get_meanvar(values, n_blocks, res);
    
    return;
    
};

/**
 * Compute mean value of block "i"
 * @param i  - index of block
 * @param block_size - size of the block
 * @param mc_results - output data of program "blocking"
 * @return mean value of block "i"
 */
double get_mean(int i, int block_size, mat& mc_results){
    
    double mean_v = 0.0;
    
    for(int j = 0; j< block_size; j++){
        mean_v += mc_results(i*block_size + j);
    }
    
    mean_v /= block_size;
    
    return mean_v;
    
    
}

/**
 * Compute sample variance of the entries of vector "values"
 * @param values - input vector
 * @param n_blocks - size of vector "values"
 * @param res - res(0): mean value, res(1): sample variance
 */
void get_meanvar(vec& values, int n_blocks, vec& res){
    
    res.zeros();
    
    for(int i=0; i< n_blocks; i++){
        res(0) += values(i);
        res(1) += values(i)*values(i);
    }
    
    res(0) /= n_blocks;
    res(1) = res(1)/n_blocks - res(0)*res(0);
    
    return;
       
}


