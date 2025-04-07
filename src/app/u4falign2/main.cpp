#include <map>
#include <string>

#include "../../corelib/hbn_aux.h"

using namespace std;

typedef int main_fun_type(int argc, char* argv[]);

int bam_to_frag_bam_main(int argc, char* argv[]);
int ddf_score_stats(int argc, char* argv[]);
int dump_chr_name_main(int argc, char* argv[]);
int eval_sim_frag_bam(int argc, char* argv[]);
int extract_chr_main(int argc, char* argv[]);
int extract_fastq_main(int argc, char* argv[]);
int fix_ngmlr_sam_main(int argc, char* argv[]);
int frag_and_contact_stats_main(int argc, char* argv[]);
int sim_fastq_stats(int argc, char* argv[]);
int simulate_main(int argc, char* argv[]);

class FalignUtilityCommandList
{
public:
    FalignUtilityCommandList() {
        add_one_cmd("bam-to-frag-bam", bam_to_frag_bam_main);
        add_one_cmd("ddf-score-stats", ddf_score_stats);
        add_one_cmd("dump-chr-name", dump_chr_name_main);
        add_one_cmd("eval-sim-frag-bam", eval_sim_frag_bam);
        add_one_cmd("extract-fastq", extract_fastq_main);
        add_one_cmd("extract-chr", extract_chr_main);
        add_one_cmd("fix-ngmlr-sam", fix_ngmlr_sam_main);
        add_one_cmd("frag-bam-stats", frag_and_contact_stats_main);
        add_one_cmd("sim-fastq-stats", sim_fastq_stats);
        add_one_cmd("sim", simulate_main);
    }

    void add_one_cmd(const char* cmdname, main_fun_type* fun) {
        M_cmds[std::string(cmdname)] = fun;
    }

    void dump_cmds(FILE* out = stderr) {
        fprintf(stderr, "\n");
        fprintf(stderr, "COMMANDS:\n");
        for (auto& cmd : M_cmds) {
            fprintf(stderr, "  %s\n", cmd.first.c_str());
        }
    }

    int run_cmd(int argc, char* argv[]) {
        auto pos = M_cmds.find(std::string(argv[1]));
        if (pos == M_cmds.end()) {
            fprintf(stderr, "ERROR: Unrecognised command '%s'\n", argv[1]);
            return 1;
        }
        return (*pos->second)(argc, argv);
    }

private:
    std::map<std::string, main_fun_type*>    M_cmds;
};

int main(int argc, char* argv[])
{
    FalignUtilityCommandList cmd;
    if (argc < 2) {
        cmd.dump_cmds();
        return 1;
    }
    return cmd.run_cmd(argc, argv);
}