#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "../../corelib/frag_id.hpp"
#include "../../corelib/line_reader.hpp"
#include "../../corelib/pdqsort.h"
#include "../../htslib/sam.h"
#include "../map/restrict_enzyme_loci_list.hpp"

using namespace std;

static const u16 kSamFlag_unmapped = 4;
static const u16 kSamFlag_reversed = 16;
static const u16 kSamFlag_not_primary = 256;
static const u16 kSamFlag_supplementary_alignment = 2048;

static inline char completement_residue(const char c)
{
    switch (c) {
        case 'a':
        case 'A':
            return 'T';
        case 'c':
        case 'C':
            return 'G';
        case 'g':
        case 'G':
            return 'C';
        case 't':
        case 'T':
            return 'A';
        case 'n':
        case 'N':
            return 'N';
        default:
            return 'A';
    }
}

struct ReadEnzymePosInfo
{
    map<string, int> name2id;
    vector<int> fwd_enzyme_pos_list;
    vector<size_t> fwd_enzyme_pos_offset_list;
    vector<int> fwd_enzyme_pos_cnt_list;

    vector<int> rev_enzyme_pos_list;
    vector<size_t> rev_enzyme_pos_offset_list;
    vector<int> rev_enzyme_pos_cnt_list;

    void get_enzyme_pos_list(int id, int strand, int*& enzyme_pos_list, int& pos_cnt)
    {
        if (strand == FWD) {
            enzyme_pos_list = fwd_enzyme_pos_list.data() + fwd_enzyme_pos_offset_list[id];
            pos_cnt = fwd_enzyme_pos_cnt_list[id];
        } else {
            enzyme_pos_list = rev_enzyme_pos_list.data() + rev_enzyme_pos_offset_list[id];
            pos_cnt = rev_enzyme_pos_cnt_list[id];
        }
    }
};

static bool 
s_get_next_line(HbnLineReader* in, kstring_t* line, string& sline)
{
    if (!in->ReadOneLine(line)) return false;
    sline.assign(line->s, line->s + line->l);
    return true;
}

static void
s_load_read_enzyme_pos_info(const char* fq_path, const char* enzyme_seq, ReadEnzymePosInfo& enzyme_info)
{
    HBN_LOG("Load read enzyme pos from %s", fq_path);
    HbnLineReader in(fq_path);
    kstring_t line = KS_INITIALIZE;
    string name, hdr, fwd_seq, rev_seq, plus, qual;
    vector<int> fwd_enzyme_pos, rev_enzyme_pos;
    int id = 0;
    RestrictEnzyme enzyme;
    RestrictEnzyme_Init(enzyme_seq, &enzyme);

    while (s_get_next_line(&in, &line, hdr)) {
        if (!s_get_next_line(&in, &line, fwd_seq)) HBN_ERR("unexpected end of fastq file");
        if (!s_get_next_line(&in, &line, plus)) HBN_ERR("unexpected end of fastq file");
        if (!s_get_next_line(&in, &line, qual)) HBN_ERR("unexpected end of fastq file");
        hbn_assert(plus[0] == '+');
        hbn_assert(fwd_seq.size() == qual.size());

        int n = hdr.size();
        name.clear();
        for (int i = 1; i < n; ++i) {
            if (isspace(hdr[i])) break;
            name += hdr[i];
        }
        hbn_assert(enzyme_info.name2id.find(name) == enzyme_info.name2id.end());
        enzyme_info.name2id[name] = id;
        ++id;

        fwd_enzyme_pos.clear();
        extract_enzyme_loci_list_for_one_seq(fwd_seq.c_str(), fwd_seq.size(), &enzyme, fwd_enzyme_pos);
        enzyme_info.fwd_enzyme_pos_offset_list.push_back(enzyme_info.fwd_enzyme_pos_list.size());
        enzyme_info.fwd_enzyme_pos_cnt_list.push_back(fwd_enzyme_pos.size());
        enzyme_info.fwd_enzyme_pos_list.insert(enzyme_info.fwd_enzyme_pos_list.end(), fwd_enzyme_pos.begin(), fwd_enzyme_pos.end());

        rev_seq.assign(fwd_seq.rbegin(), fwd_seq.rend());
        for (auto& c : rev_seq) c = completement_residue(c);
        rev_enzyme_pos.clear();
        extract_enzyme_loci_list_for_one_seq(rev_seq.c_str(), rev_seq.size(), &enzyme, rev_enzyme_pos);
        enzyme_info.rev_enzyme_pos_offset_list.push_back(enzyme_info.rev_enzyme_pos_list.size());
        enzyme_info.rev_enzyme_pos_cnt_list.push_back(rev_enzyme_pos.size());
        enzyme_info.rev_enzyme_pos_list.insert(enzyme_info.rev_enzyme_pos_list.end(), rev_enzyme_pos.begin(), rev_enzyme_pos.end());
    }
    HBN_LOG("Done");
    ks_free(&line);
}

static bool 
s_skip_this_bam(bam1_t* bam)
{
    u16 flag = bam->core.flag;
    return (flag & kSamFlag_unmapped)
           ||
           (flag & kSamFlag_not_primary);
}

struct FragBam
{
    int qoff;
    int frag_id;
    bam1_t* bam;
    bam1_t* frag_bam;
};

static bam1_t*
s_read_nect_bam(const char* sam_path, samFile* fp, sam_hdr_t* hdr)
{
	static size_t cnt = 0;
    bam1_t* bam = bam_init1();
    int r = sam_read1(fp, hdr, bam);
    if (r == -1) {
        bam_destroy1(bam);
        return nullptr;
    }
    if (r < 0) {
        bam_destroy1(bam);
	fprintf(stderr, "%s\n", strerror(r));
        HBN_ERR("ERROR at reading BAM record from %s", sam_path);
    }
    //fprintf(stderr, "%zu\t%s\n", cnt++, bam_get_qname(bam));
    return bam;
}

static int
s_resolve_fwd_query_offset(bam1_t* bam)
{
    int qsize = 0;
    int qoff = 0;
    int qend = 0;
    uint32_t* cigar = bam_get_cigar(bam);
    for (int i = 0; i < bam->core.n_cigar; ++i) {
        char op = bam_cigar_opchr(cigar[i]);
        int num = bam_cigar_oplen(cigar[i]);
        if (i == 0 && op == 'S') {
            qoff = num;
            qend = num;
        }
        if (op == 'H' || op == 'S' || op == 'M' || op == '=' || op == 'X' || op == 'I') qsize += num;
        if (op == 'M' || op == '=' || op == 'X' || op == 'I') qend += num;
    }
    u16 flag = bam->core.flag;
    int qdir = (flag & kSamFlag_reversed) ? REV : FWD;
    if (qdir == REV) qoff = qsize - qend;
    return qoff;
}

static inline char s_decode_bam_query_base(int c)
{
    // 1 -> A
    // 2 -> C
    // 4 -> G
    // 8 -> T
    // 15 -> N
    switch (c) {
        case 1:
            return 'A';
        case 2:
            return 'C';
        case 4:
            return 'G';
        case 8:
            return 'T';
        case 15:
            return 'N';
        default:
            HBN_ERR("Illegal BAM base encoded value %s", c);
    }
}

void compute_align_identity_from_bam(bam1_t* bam, int& frag_size, float& pi)
{
    if (bam->core.n_cigar == 0) return;
    frag_size = 0;
    uint32_t* cigar = bam_get_cigar(bam);
    int as_size = 0, match_size = 0;
    for (int i = 0; i < bam->core.n_cigar; ++i) {
        char op = bam_cigar_opchr(cigar[i]);
        int num = bam_cigar_oplen(cigar[i]);
        if (op == 'M' || op == '=') {
            as_size += num;
            match_size += num;
            frag_size += num;
        } else if (op == 'X') {
            as_size += num;
            frag_size += num;
        } else if (op == 'I') {
            if (num < 8) as_size += num;
            frag_size += num;
        } else if (op == 'D') {
            if (num < 8) as_size += num;
        }
    }
    if (as_size < 50) {
	    pi = 0.0;
	    return;
    }
    hbn_assert(as_size > 0);
    pi = 100.0 * match_size / as_size;
}

static bam1_t*
s_bam_to_frag_bam(const int* fwd_enzyme_pos, const int fwd_enzyme_pos_cnt, const int* rev_enzyme_pos, const int rev_enzyme_pos_cnt, RestrictEnzymeLociList* updb_enzyme, 
    sam_hdr_t* hdr, bam1_t* bam, int read_id, int frag_id)
{
    char frag_suffix[HBN_MAX_PATH_LEN];
    frag_id_to_string(read_id, frag_id, bam->core.tid, bam->core.pos, frag_suffix);
    string frag_name = bam_get_qname(bam);
    frag_name += '_';
    frag_name += frag_suffix;

    vector<u32> frag_cigar;
    uint32_t* cigar = bam_get_cigar(bam);
    int qoff = 0, qend = 0, qsize = 0, qH = 0;
    int soff = bam->core.pos, send = soff;
    for (int i = 0; i < bam->core.n_cigar; ++i) {
        char op = bam_cigar_opchr(cigar[i]);
        int num = bam_cigar_oplen(cigar[i]);
        if (i == 0 && op == 'S') {
            qoff = num;
            qend = num;
        }
        if (i == 0 && op == 'H') {
            qH = num;
        }
        if (op == 'M' || op == '=' || op == 'X' || op == 'I') qend += num;
        if (op == 'M' || op == '=' || op == 'X' || op == 'D') send += num;
        if (op != 'H' && op != 'S') frag_cigar.push_back(cigar[i]);
        if (op == 'M' || op == '=' || op == 'X' || op == 'I' || op == 'S' || op == 'H') qsize += num;
    }

    string frag_seq, frag_qv;
    uint8_t* seq = bam_get_seq(bam);
    uint8_t* qual = bam_get_qual(bam);
    for (int i = qoff; i < qend; ++i) {
        int c = bam_seqi(seq, i);
        c = s_decode_bam_query_base(c);
        frag_seq += c;
        frag_qv += qual[i];
    }

    bam1_t* frag_bam = bam_init1();
    u16 flag = bam->core.flag;
    if (flag & kSamFlag_supplementary_alignment) flag -= kSamFlag_supplementary_alignment;
    int r = bam_set1(frag_bam, frag_name.size(), frag_name.c_str(), flag,
        bam->core.tid, bam->core.pos, bam->core.qual, frag_cigar.size(), frag_cigar.data(),
        -1, -1, 0, frag_seq.size(), frag_seq.c_str(), frag_qv.c_str(), 512);
    if (r < 0) HBN_ERR("FAIL at building BAM record");

    int frag_size = 0;
    float pi = 0.0;
    compute_align_identity_from_bam(bam, frag_size, pi);

    char tag[2];
    /// qid
    tag[0] = 'q'; tag[1] = 'i';
    r = bam_aux_append(frag_bam, tag, 'i', sizeof(int), (const uint8_t*)(void*)(&read_id));
    if (r) HBN_ERR("Could not append tag %c%c to BAM record", tag[0], tag[1]);

    /// qdir
    tag[0] = 'q'; tag[1] = 'd';
    int qdir = (bam->core.flag&16) ? REV : FWD;
    r = bam_aux_append(frag_bam, tag, 'i', sizeof(int), (const uint8_t*)(void*)(&qdir));
    if (r) HBN_ERR("Could not append tag %c%c to BAM record", tag[0], tag[1]);

    /// qs
    tag[0] = 'q'; tag[1] = 's';
    int qb = qH + qoff;
    r = bam_aux_append(frag_bam, tag, 'i', sizeof(int), (const uint8_t*)(void*)(&qb));
    if (r) HBN_ERR("Could not append tag %c%c to BAM record", tag[0], tag[1]);

    /// qe
    tag[0] = 'q'; tag[1] = 'e';
    int qe = qend + qH;
    r = bam_aux_append(frag_bam, tag, 'i', sizeof(int), (const uint8_t*)(void*)(&qe));
    if (r) HBN_ERR("Could not append tag %c%c to BAM record", tag[0], tag[1]);

    /// ql
    tag[0] = 'q'; tag[1] = 'l';
    r = bam_aux_append(frag_bam, tag, 'i', sizeof(int), (const uint8_t*)(void*)(&qsize));
    if (r) HBN_ERR("Could not append tag %c%c to BAM record", tag[0], tag[1]);

    /// sid
    tag[0] = 's'; tag[1] = 'i';
    int sid = bam->core.tid;
    r = bam_aux_append(frag_bam, tag, 'i', sizeof(int), (const uint8_t*)(void*)(&sid));
    if (r) HBN_ERR("Could not append tag %c%c to BAM record", tag[0], tag[1]);

    /// ss
    tag[0] = 's'; tag[1] = 's';
    r = bam_aux_append(frag_bam, tag, 'i', sizeof(int), (const uint8_t*)(void*)(&soff));
    if (r) HBN_ERR("Could not append tag %c%c to BAM record", tag[0], tag[1]);    

    /// se
    tag[0] = 's'; tag[1]= 'e';
    r = bam_aux_append(frag_bam, tag, 'i', sizeof(int), (const uint8_t*)(void*)(&send));
    if (r) HBN_ERR("Could not append tag %c%c to BAM record", tag[0], tag[1]);

    /// sl 
    tag[0] = 's'; tag[1] = 'l';
    int sl = sam_hdr_tid2len(hdr, bam->core.tid);
    r = bam_aux_append(frag_bam, tag, 'i', sizeof(int), (const uint8_t*)(void*)(&sl));
    if (r) HBN_ERR("Could not append tag %c%c to BAM record", tag[0], tag[1]);

    const int* enzyme_pos = bam_is_rev(bam) ? fwd_enzyme_pos : rev_enzyme_pos;
    const int enzyme_pos_cnt = bam_is_rev(bam) ? fwd_enzyme_pos_cnt : rev_enzyme_pos_cnt;
    /// vqb
    int vqb = get_enzyme_pos(enzyme_pos, enzyme_pos_cnt, qb);
    tag[0] = 'q'; tag[1] = 'S';
    r = bam_aux_append(frag_bam, tag, 'i', sizeof(int), (const uint8_t*)(void*)(&vqb));
    if (r) HBN_ERR("Could not append tag %c%c to BAM record", tag[0], tag[1]);

    /// vqe
    int vqe = get_enzyme_pos(enzyme_pos, enzyme_pos_cnt, qe);
    tag[0] = 'q'; tag[1] = 'E';
    r = bam_aux_append(frag_bam, tag, 'i', sizeof(int), (const uint8_t*)(void*)(&vqe));
    if (r) HBN_ERR("Could not append tag %c%c to BAM record", tag[0], tag[1]);

	int vid = sid * 2;
    const int* chr_reloci = updb_enzyme->reloci_array + updb_enzyme->seq_reloci_info_array[vid].enzyme_loci_offset;
    const int chr_reloci_cnt = updb_enzyme->seq_reloci_info_array[vid].enzyme_loci_cnt;
    /// vs
    int vsb = get_enzyme_pos(chr_reloci, chr_reloci_cnt, soff);
    tag[0] = 'v'; tag[1] = 'S';
    r = bam_aux_append(frag_bam, tag, 'i', sizeof(int), (const uint8_t*)(void*)(&vsb));
    if (r) HBN_ERR("Could not append tag %c%c to BAM record", tag[0], tag[1]);    

    /// ve
    int vse = get_enzyme_pos(chr_reloci, chr_reloci_cnt, send);
    tag[0] = 'v'; tag[1] = 'E';
    r = bam_aux_append(frag_bam, tag, 'i', sizeof(int), (const uint8_t*)(void*)(&vse));
    if (r) HBN_ERR("Could not append tag %c%c to BAM record", tag[0], tag[1]);    

    return frag_bam;
}

static void
s_build_frag_bam_list(ReadEnzymePosInfo& read_enzyme, RestrictEnzymeLociList* updb_enzyme, 
    sam_hdr_t* hdr, vector<bam1_t*>& bam_list, int read_id, vector<FragBam>& frag_bam_list)
{
    string name = bam_get_qname(bam_list[0]);
    auto pos = read_enzyme.name2id.find(name);
    hbn_assert(pos != read_enzyme.name2id.end());
    int id = pos->second;
    int* fwd_enzyme_pos;
    int* rev_enzyme_pos;
    int fwd_enzyme_pos_cnt, rev_enzyme_pos_cnt;
    read_enzyme.get_enzyme_pos_list(id, FWD, fwd_enzyme_pos, fwd_enzyme_pos_cnt);
    read_enzyme.get_enzyme_pos_list(id, REV, rev_enzyme_pos, rev_enzyme_pos_cnt);

    frag_bam_list.clear();
    FragBam fb;
    for (bam1_t* bam : bam_list) {
        if (s_skip_this_bam(bam)) continue;
        fb.qoff = s_resolve_fwd_query_offset(bam);
        fb.frag_id = -1;
        fb.bam = bam;
        fb.frag_bam = nullptr;
        frag_bam_list.push_back(fb);
    }
    if (frag_bam_list.empty()) return;

    pdqsort(frag_bam_list.begin(), frag_bam_list.end(), [](const FragBam& x, const FragBam& y) { return x.qoff < y.qoff; });
    int N = frag_bam_list.size();
    for (int i = 0; i < N; ++i) frag_bam_list[i].frag_id = i;

    for (int i = 0; i < N; ++i) frag_bam_list[i].frag_bam = s_bam_to_frag_bam(fwd_enzyme_pos, fwd_enzyme_pos_cnt, rev_enzyme_pos, rev_enzyme_pos_cnt, updb_enzyme, hdr, frag_bam_list[i].bam, read_id, i);
}

int bam_to_frag_bam_main(int argc, char* argv[])
{
    if (argc != 7) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s %s reference reads enzyme bam frag-bam\n", argv[0], argv[1]);
        exit(1);
    }
    const char* reference_path = argv[2];
    const char* reads_path = argv[3];
    const char* enzyme_seq = argv[4];
    const char* bam_path = argv[5];
    const char* frag_bam_path = argv[6];

    HbnUnpackedDatabase updb(reference_path);
    RestrictEnzymeLociList* updb_enzyme = RestrictEnzymeLociListNew(&updb, enzyme_seq, 24);
    ReadEnzymePosInfo read_enzyme;
    s_load_read_enzyme_pos_info(reads_path, enzyme_seq, read_enzyme);

    samFile* in = sam_open(bam_path, "rb");
    if (!in) HBN_ERR("Could not open file %s for reading", bam_path);
    hts_set_threads(in, 8);
    sam_hdr_t* hdr = sam_hdr_read(in);
    if (!hdr) HBN_ERR("Could not read BAM header from %s", bam_path);

    samFile* out = sam_open(frag_bam_path, "wb");
    hts_set_threads(out, 8);
    sam_hdr_t* hdrout = sam_hdr_dup(hdr);
    if (sam_hdr_write(out, hdrout)) HBN_ERR("Could not write BAM header to %s", frag_bam_path);

    vector<bam1_t*> bam_list;
    vector<FragBam> frag_bam_list;
    int read_id = 0;
    bam1_t* bam = s_read_nect_bam(bam_path, in, hdr);
    while (bam) {
        bam_list.clear();
        const char* qname = bam_get_qname(bam);
        bam_list.push_back(bam);
        while ((bam = s_read_nect_bam(bam_path, in, hdr))) {
            const char* next_qname = bam_get_qname(bam);
            if (strcmp(qname, next_qname)) break;
            bam_list.push_back(bam);
        }
        s_build_frag_bam_list(read_enzyme, updb_enzyme, hdr, bam_list, read_id++, frag_bam_list);
        for (auto& fb : frag_bam_list) {
            int r = sam_write1(out, hdrout, fb.frag_bam);
            if (r == -1) HBN_ERR("Could not write BAM record to %s", frag_bam_path);
        }
        for (bam1_t* b : bam_list) bam_destroy1(b);
        bam_list.clear();
        for (auto& fb : frag_bam_list) bam_destroy1(fb.frag_bam);
        frag_bam_list.clear();
    }
    sam_close(out);
    sam_hdr_destroy(hdrout);
    sam_close(in);
    sam_hdr_destroy(hdr);

    return 0;
}
