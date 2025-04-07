#ifndef __SIMULATE_PORE_C_FRAG_HPP
#define __SIMULATE_PORE_C_FRAG_HPP

#include "../../corelib/hbn_aux.h"
#include "../../corelib/unpacked_seqdb.hpp"

#include <iostream>
#include <vector>

double gen_one_prob();
int gen_one_positive_int();

struct SimFrag
{
    int chr_id;
    int chunk_id;
    int from;
    int to;
    int strand;
};

inline std::ostream& operator<<(std::ostream& os, const SimFrag& frag)
{
    os << frag.chr_id << '\t' << frag.from << '\t' << frag.to << '\t' << frag.to - frag.from;
    return os;
}

struct SimFragChunk
{
    int chr_id;
    int chunk_id;

    SimFrag* sfa;
    int sfc;
    
    size_t __frag_offset;
};

struct ChrFragInfo
{
    int chr_id;
    
    SimFragChunk* sfca;
    int sfcc;

    SimFrag* sfa;
    int sfc;

    size_t __chunk_offset;
    size_t __frag_offset;
};

struct SimFragLibrary
{
    std::vector<ChrFragInfo>    chr_list;
    std::vector<SimFragChunk>   chunk_list;
    std::vector<SimFrag>        frag_list;

    SimFragLibrary(HbnUnpackedDatabase* updb, const char* enzyme);
};

#endif // __SIMULATE_PORE_C_FRAG_HPPs