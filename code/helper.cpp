
#include "helper.h"

#include <list>
#include <string>
#include <utility>
#include <vector>

using namespace std;



/*
 * Internal function, see next function for details.
 * Only used for nice recursive computation.
 */
list<list<pair<string,float>>>* cartesian_product_internal(
        vector<pair<string,vector<float>*>>& input,
        size_t start_index)
{
    // base case
    if (input.size() - start_index == 0)
        return new list<list<pair<string,float>>>();
    if (input.size() - start_index == 1) {
        list<list<pair<string,float>>>* retlist = new list<list<pair<string,float>>>();
        auto& this_only_element = input[input.size()-1];
        string& this_name = this_only_element.first;
        for (float value : *(this_only_element.second)) {
            retlist->emplace_front();
            auto& new_list = retlist->front();
            new_list.emplace_front(this_name, value);
        }
        return retlist;
    }
    // recursive case
    list<list<pair<string,float>>>* new_list = new list<list<pair<string,float>>>();
    // recursive call, starting at the next position
    auto* old_list = cartesian_product_internal(input, start_index+1);
    // now we loop over all values of the current variable and combine
    // the values with the results from the recursive call
    string& currVarName = input[start_index].first;
    for (float value : *(input[start_index].second )) {
        // for every combination in the output
        // add all combinations of the input
        for (auto& old_sublist : *old_list) {
            new_list->push_front(old_sublist); // copy constructor
            auto& new_element = new_list->front();
            new_element.emplace_back(currVarName, value);
        }
    }
    // delete old objects
    delete old_list;
    //
    return new_list;
}


/*
   This function computes the cartesian product (i.e. everything with everything)
   of a vector of vectors.
*/
list<list<pair<string,float>>>* cartesian_product(
        vector<pair<string,vector<float>*>>& input)
{
    return cartesian_product_internal(input, 0);
}



/*
 * This function compares two struct tm objects directly
 * instead of using mktime.
 *
 * Returns 0,  if tm1 == tm2
 * Returns -1, if tm1 <  tm2
 * Returns +1, if tm1 >  tm2
 */
int compare_struct_tm(struct tm* a, struct tm* b) {
    if (a->tm_year > b->tm_year)
        return  1;
    if (a->tm_year < b->tm_year)
        return -1;
    if (a->tm_mon  > b->tm_mon)
        return  1;
    if (a->tm_mon  < b->tm_mon)
        return -1;
    if (a->tm_mday > b->tm_mday)
        return  1;
    if (a->tm_mday < b->tm_mday)
        return -1;
    if (a->tm_hour > b->tm_hour)
        return  1;
    if (a->tm_hour < b->tm_hour)
        return -1;
    if (a->tm_min  > b->tm_min)
        return  1;
    if (a->tm_min  < b->tm_min)
        return -1;
    if (a->tm_sec  > b->tm_sec)
        return  1;
    if (a->tm_sec  < b->tm_sec)
        return -1;
    return 0;
}


