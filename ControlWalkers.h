/* 
 * File:   ControlWalkers.h
 * Author: sarahrei
 *
 * Created on April 25, 2012, 4:12 PM
 */


#ifndef CONTROLWALKERS_H
#define	CONTROLWALKERS_H

#include "Walker.h"

class ControlWalkers {
private:
    int numpart;
    int dim;

public:
    /**
     * Constructor
     * @param num_part - number of particles
     * @param dim
     */
    ControlWalkers(int num_part, int dim);

    ~ControlWalkers();

    /**
     * Copy a walker
     * @param Original - original walker
     * @param Twin - new copied walker
     */
    void copyWalker(Walker *Original, Walker *Twin);

};

#endif	/* CONTROLWALKERS_H */

