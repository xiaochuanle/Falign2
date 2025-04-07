#ifndef __SMOOTH_PCA_LIST_HPP
#define __SMOOTH_PCA_LIST_HPP

#include "../../sw/hbn_traceback_aux.h"
#include "chain_align_list.hpp"
#include "pore_c_aux.hpp"
#include "pore_c_traceback.hpp"
#include "query_input.hpp"

#include <string>
#include <vector>

void
smooth_pca_list(std::vector<PoreCAlign>& pca_list,
    EChainType& chain_type,
    std::vector<PoreCAlign>& all_pca_list,
    HbnUnpackedDatabase* updb,
    PoreCQuery* query,
    const int enzyme_size,
    HbnTracebackData* tbck_data);

#endif // __SMOOTH_PCA_LIST_HPP
