#ifndef __INFER_ENZYME_HPP
#define __INFER_ENZYME_HPP

#include "../../corelib/unpacked_seqdb.hpp"
#include "index.hpp"
#include "map_options.hpp"
#include "query_input.hpp"

#include <string>

std::string infer_enzyme_mt(HbnMapOptions* options, HbnIndex* index, QueryReader* queries);

#endif // __INFER_ENZYME_HPP