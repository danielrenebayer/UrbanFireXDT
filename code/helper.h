/*
 * helper.h
 *
 * This file contains functions that can be used everywhere
 * in the simulation.
 * They might make things more easy.
 */

#ifndef __HELPER_H_
#define __HELPER_H_

#include <cmath>
#include <list>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace std;


/*
   This function computes the cartesian product (i.e. everything with everything)
   of a vector of vectors.
*/
list<list<pair<string,float>>>* cartesian_product(
        vector<pair<string,shared_ptr<vector<float>>>>& input);

/*
 * This function compares two struct tm objects directly
 * instead of using mktime.
 *
 * Returns 0,  if tm1 == tm2
 * Returns -1, if tm1 <  tm2
 * Returns +1, if tm1 >  tm2
 */
int compare_struct_tm(struct tm* a, struct tm* b);

/**
 * Rounds a float on n decimal places
 */
float round_float_n(float x, unsigned int decimal_places);

/**
 * Rounds a double on 5 decimal places and returns it as a float
 */
float round_float_5(double x);

/**
 * Rounds a float on 5 decimal places
 */
float round_float_5(float x);


#endif

