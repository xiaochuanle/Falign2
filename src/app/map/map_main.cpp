#include "../../corelib/unpacked_seqdb.hpp"
#include "../../third-party/htslib-1.19.1/version.h"
#include "index.hpp"
#include "infer_enzyme.hpp"
#include "map_one_query_file.hpp"
#include "map_options.hpp"
#include "program_info.hpp"
#include "query_input.hpp"

#include <sys/utsname.h>
#include <string>
#include <thread>

using namespace std;

extern "C"
size_t getMemorySizeBytes();

void map_program_info(FILE* out)
{
    struct utsname _os_info_buf;
    struct utsname* os_info = nullptr;
    if (uname(&_os_info_buf) == 0) os_info = &_os_info_buf;

    size_t sys_mem_bytes = getMemorySizeBytes();
    string sys_mem = NStr::UInt8ToString_DataSize(sys_mem_bytes);

    int cpu_threads = thread::hardware_concurrency();

    fprintf(out, "\n");
    fprintf(out, "PROGRAM:\n");
    fprintf(out, "  Name:           %s\n", HBN_PACKAGE_NAME);
    fprintf(out, "  Version:        %s\n", HBN_PACKAGE_VERSION);
    fprintf(out, "  htslib:         %s\n", HTS_VERSION_TEXT);
    fprintf(out, "  Description:    Alignment toolkit for Pore-C reads\n");
    fprintf(out, "  Contact:        chenying2016@gmail.com\n");

    fprintf(out, "\n");
    fprintf(out, "SYSTEM:\n");
    if (os_info) {
    fprintf(out, "  Computer:       %s\n", os_info->nodename);
    fprintf(out, "  Name:           %s\n", os_info->sysname);
    fprintf(out, "  Release:        %s\n", os_info->release);
    fprintf(out, "  Version:        %s\n", os_info->version);
    fprintf(out, "  Machine:        %s\n", os_info->machine);
    }
    fprintf(out, "  CPU threads:    %d\n", cpu_threads);
    fprintf(out, "  RAM:            %s\n", sys_mem.c_str());
    fprintf(out, "\n");
}

int map_main(int argc, char* argv[])
{
    vector<string> query_list;
    HbnMapOptions opts;
    if (!opts.parse(argc, argv, query_list)) {
        opts.simple_usage(argc, argv);
        return EXIT_FAILURE;
    }
    HbnProgramInfo program(HBN_PACKAGE_NAME, map_program_info);

    HbnIndex index;
    index.load(opts.reference);
    if (opts.enzyme) index.init_enzyme(opts.enzyme, opts.num_threads);
    HbnOutputs out(argc, argv, index.updb(), opts.output, opts.outfmt);
    QueryReader queries(query_list);
    string last_enzyme;
    while (queries.get_next_query_file()) {
        if (!opts.enzyme) {
            HBN_LOG("Infer enzyme for %s", queries.get_curr_query_fn());
            string enzyme = infer_enzyme_mt(&opts, &index, &queries);
            if (enzyme.empty()) {
                HBN_LOG("WARNING: Could not infer enzyme for %s", queries.get_curr_query_fn());
                HBN_LOG("  We skip this query file and advance to map next query file");
                continue;
            }
            HBN_LOG("Successfully infer enzyme: %s", enzyme.c_str());
            if (enzyme != last_enzyme) {
                index.init_enzyme(enzyme.c_str(), opts.num_threads);
                last_enzyme = enzyme;
            }
        }
        map_one_pore_c_file(&opts, &index, &queries, &out);
    }

    return 0;
}
