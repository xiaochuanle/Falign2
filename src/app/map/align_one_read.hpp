#ifndef __ALIGN_ONE_READ_HPP
#define __ALIGN_ONE_READ_HPP

#include "chain_align_list.hpp"
#include "gapped_align.hpp"
#include "map_options.hpp"
#include "pore_c_aux.hpp"
#include "prelim_search.hpp"
#include "query_input.hpp"
#include "trim_overlap_subseq.hpp"

#include <vector>

void
align_one_read(HbnIndex* index,
    PoreCQuery* query,
    HbnMapOptions* options,
    PrelimSearch* prelim_search,
    PoreCGapAlign* gapped_search,
    PoreCAlignChainData* pca_chain_data,
    TrimPcaList* trim_pca_list);

#endif // __ALIGN_ONE_READ_HPP