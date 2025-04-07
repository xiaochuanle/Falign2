#include "../../corelib/fasta.hpp"

#include <set>
#include <string>

using namespace std;

int extract_chr_main(int argc, char* argv[])
{
    if (argc < 4) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s %s genome chr-list\n", argv[0], argv[1]);
        return 1;
    }
    const char* genome_path = argv[2];
    set<string> chrs;
    for (int i = 3; i < argc; ++i) chrs.insert(argv[i]);

    HbnFastaReader in(genome_path);
    while (1) {
        if (in.ReadOneSeq() == -1) break;
        const string& name = in.name();
        if (chrs.find(name) == chrs.end()) continue;
        
        const char* s = in.name().c_str();
        size_t sl = in.name().size();
        fprintf(stdout, ">");
        hbn_fwrite(s, 1, sl, stdout);
        fprintf(stdout, "\n");

        s = in.sequence().c_str();
        sl = in.sequence().size();
        hbn_fwrite(s, 1, sl, stdout);
        fprintf(stdout, "\n");
    }

    return 0;
}