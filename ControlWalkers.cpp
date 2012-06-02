
#include "ControlWalkers.h"
#include "lib.h"

/**
 * Constructor
 * @param num_part - number of particles
 * @param dim - dimension
 */
ControlWalkers::ControlWalkers(int num_part, int dim) {
    this->numpart = numpart;
    this->dim = dim;

}

ControlWalkers::~ControlWalkers() {
}

/**
 * Copy a walker
 * @param Original - original walker
 * @param Twin - the new copy
 */
void ControlWalkers::copyWalker(Walker *Original, Walker *Twin) {

    Twin->wf_R = Original->wf_R;
    Twin->R_cur = Original->R_cur;
    Twin->R_tr = Original->R_tr;
    Twin->idum = Original->idum;
    Twin->seed = Original->seed;

    Twin->Trial->SlaterPsi->slater = Original->Trial->SlaterPsi->slater;
    Twin->Trial->SlaterPsi->slat_inv = Original->Trial->SlaterPsi->slat_inv;
    Twin->Trial->SlaterPsi->inv_back = Original->Trial->SlaterPsi->inv_back;

    Twin->Trial->JastrowPsi->g_ij = Original->Trial->JastrowPsi->g_ij;
    Twin->Trial->Pos->r = Original->Trial->Pos->r;
    Twin->Trial->Pos->r_int = Original->Trial->Pos->r_int;
    Twin->Trial->Pos_tr->r = Original->Trial->Pos_tr->r;
    Twin->Trial->Pos_tr->r_int = Original->Trial->Pos_tr->r_int;

    Twin->QFo->qf_old = Original->QFo->qf_old;
    Twin->QFo->qf_new = Original->QFo->qf_new;


}
