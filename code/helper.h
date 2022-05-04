/*
 * helper.h
 *
 * This file contains functions that can be used everywhere
 * in the simulation.
 * They might make things more easy.
 */

#ifndef __HELPER_H_
#define __HELPER_H_

#include <list>
#include <string>
#include <utility>
#include <vector>

using namespace std;


/*
   This function computes the cartesian product (i.e. everything with everything)
   of a vector of vectors.
*/
list<list<pair<string,float>>>* cartesian_product(
        vector<pair<string,vector<float>*>>& input);

/*
 * This function compares two struct tm objects directly
 * instead of using mktime.
 *
 * Returns 0,  if tm1 == tm2
 * Returns -1, if tm1 <  tm2
 * Returns +1, if tm1 >  tm2
 */
int compare_struct_tm(struct tm* a, struct tm* b);


#endif

