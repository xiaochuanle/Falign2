#ifndef __SIMULATE_PORE_C_READ_HPP
#define __SIMULATE_PORE_C_READ_HPP

#include "../map/restrict_enzyme_loci_list.hpp"
#include "simulate_pore_c_frag.hpp"

void gen_simulated_pore_c_reads(SimFragLibrary& frag_library, HbnUnpackedDatabase* updb, 
    const char* enzyme_seq, const char* output);

#endif // __SIMULATE_PORE_C_READ_HPP