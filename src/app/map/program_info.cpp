#include "program_info.hpp"

#include "../../ncbi_blast/str_util/ncbistr.hpp"

using namespace std;

extern "C"
size_t getMemorySizeBytes();
extern "C" size_t getPeakRSS();

HbnProgramInfo::~HbnProgramInfo() {
    gettimeofday(&M_end, NULL);
    double dur = hbn_time_diff(&M_begin, &M_end);
    size_t peak_ram = getPeakRSS();
    string size = NStr::UInt8ToString_DataSize(peak_ram);
    fprintf(stderr, "\n\n");
    fprintf(stderr, "%s Wallclock time: %.2f seconds.\n", M_program.c_str(), dur);
    fprintf(stderr, "%s Peak RAM usage: %s\n", M_program.c_str(), size.c_str());
}