#include "../../corelib/fasta.hpp"

using namespace std;

int dump_chr_name_main(int argc, char* argv[])
{
    if (argc != 3) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "  %s %s genome\n", argv[0], argv[1]);
        return 1;
    }
    const char* genome_path = argv[2];
    HbnFastaReader in(genome_path);
    bool is_first_chr = true;
    while (1) {
        if (in.ReadOneSeq() == -1) break;
        if (!is_first_chr) fprintf(stdout, " ");
        fprintf(stdout, "%s", in.name().c_str());
        is_first_chr = false;
    }
    fprintf(stdout, "\n");
    return 0;
}
