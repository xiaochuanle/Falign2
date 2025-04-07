#ifndef __TRIM_OVERLAP_SUBSEQ_H
#define __TRIM_OVERLAP_SUBSEQ_H

#include "chain_align_list.hpp"
#include "gapped_align.hpp"
#include "pore_c_aux.hpp"
#include "pore_c_traceback.hpp"
#include "prelim_search.hpp"
#include "query_input.hpp"
#include "restrict_enzyme_loci_list.hpp"

#include <string>
#include <vector>

void
trim_overlap_subseqs(HbnTracebackData* tbck_data,
    PrelimSearch* prelim_search,
    HbnUnpackedDatabase* updb,
    RestrictEnzymeLociList* reloci_list,
    PoreCQuery* query,
    PoreCAlign* all_pca_a,
    int all_pca_c,
    PoreCAlign* pca_a,
    int pca_c,
    const EChainType chain_type,
    TrimPcaList& trim_pca_list);

#endif // __TRIM_OVERLAP_SUBSEQ_H