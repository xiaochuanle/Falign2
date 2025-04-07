#include "../../corelib/line_reader.hpp"
#include "../../corelib/split_string_by_char.hpp"

using namespace std;

int fix_ngmlr_sam_main(int argc, char* argv[])
{
    if (argc != 4) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s %s sam fixed-sam\n", argv[0], argv[1]);
        return 1;
    }
    const char* input_sam_path = argv[2];
    const char* fixed_sam_path = argv[3];

    hbn_dfopen(out, fixed_sam_path, "w");
    HbnLineReader in(input_sam_path);
    kstring_t line = KS_INITIALIZE;
    string sline;
    vector<pair<const char*, int>> cols;
    size_t n = 0, m = 0;
    while (in.ReadOneLine(&line)) {
        if (line.l == 0) continue;
        if (line.s[0] == '@') {
            hbn_fwrite(line.s, 1, line.l, out);
            fprintf(out, "\n");
            continue;
        }
        ++n;
        sline.assign(line.s, line.s + line.l);
        cols.clear();
        split_string_by_char(sline.c_str(), sline.size(), '\t', cols);
        int mapQ = atoi(cols[4].first);
        if (mapQ < 0 || mapQ > 60) {
            ++m;
            continue;
        }
        fprintf(out, "%s\n", sline.c_str());
    }
    hbn_fclose(out);
    double f = 1.0 * m / n;
    HBN_LOG("Filter %zu/%zu (%g) SAM records", m, n, f);

    return 0;
}
