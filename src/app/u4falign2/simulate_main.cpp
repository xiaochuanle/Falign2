#include "simulate_main.hpp"

#include "simulate_pore_c_frag.hpp"
#include "simulate_pore_c_read.hpp"

using namespace std;

int simulate_main(int argc, char* argv[])
{
    const char* ref_path = argv[2];
    const char* enzyme = argv[3];
    const char* output = argv[4];
    
    HbnUnpackedDatabase updb(ref_path);
    SimFragLibrary frag_library(&updb, enzyme);

    gen_simulated_pore_c_reads(frag_library, &updb, enzyme, output);

    return 0;
}