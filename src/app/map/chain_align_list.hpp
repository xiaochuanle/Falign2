#ifndef __CHAIN_ALGIN_LIST_HPP
#define __CHAIN_ALGIN_LIST_HPP

#include "pore_c_aux.hpp"

typedef enum {
    ePerfectChain,
    ePseudoPerfectChain,
    eCompleteChain,
    eMaxCovChain
} EChainType;

const char* get_chain_type_name(EChainType type);

EChainType pca_chain_type(const int* vdfa, const int vdfc, PoreCAlign* pca_a, int pca_c);

struct GraphEdge {
    int start;
    int end;
    int weight;
};

struct PoreCAlignChainData {
    std::vector<PoreCAlign> pca_list;
    std::vector<PoreCAlign> chain;

    std::vector<int> dist_list;
    std::vector<GraphEdge> edge_list;
    std::vector<int> pred_list;
};

bool select_pca_chain( PoreCAlign* pca_array, int pca_count, PoreCAlignChainData& data);

#endif // __CHAIN_ALGIN_LIST_HPP