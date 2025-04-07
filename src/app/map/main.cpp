#include "index.hpp"
#include "map_main.hpp"
#include "program_info.hpp"

using namespace std;

static const char* kIndexCmd = "index";
static const char* kMapCmd = "map";

static void
s_dump_main_usage(int argc, char* argv[])
{
    fprintf(stderr, "\n");
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "  %s subcommand\n", argv[0]);

    fprintf(stderr, "\n");
    fprintf(stderr, "SUBCOMMANDS and DESCRIPTIONS\n");
    fprintf(stderr, "  %s\t\tConstruct an index for a reference genome\n", kIndexCmd);
    fprintf(stderr, "  %s\t\tMap Pore-C reads to the reference genome\n", kMapCmd);

    fprintf(stderr, "\n");
    fprintf(stderr, "DESCRIPTION\n");
    fprintf(stderr, "  Alignment toolkit for Pore-C sequencing reads\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "VERSION\n");
    fprintf(stderr, "  %s\n", HBN_PACKAGE_VERSION);
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        s_dump_main_usage(argc, argv);
        return EXIT_FAILURE;
    }
    if (strcmp(argv[1], kIndexCmd) == 0) {
        return build_index(argc, argv);
    } else if (strcmp(argv[1], kMapCmd) == 0) {
        return map_main(argc, argv);
    } 

    s_dump_main_usage(argc, argv);
    return EXIT_FAILURE;
}