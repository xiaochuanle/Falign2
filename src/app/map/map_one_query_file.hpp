#ifndef __MAP_ONE_QUERY_FILE_HPP
#define __MAP_ONE_QUERY_FILE_HPP

#include "index.hpp"
#include "map_options.hpp"
#include "hbn_outputs.hpp"
#include "query_input.hpp"

void map_one_pore_c_file(HbnMapOptions* options,
    HbnIndex* index,
    QueryReader* queries,
    HbnOutputs* out);

#endif // __MAP_ONE_QUERY_FILE_HPP