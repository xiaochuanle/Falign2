#include "hbn_outputs.hpp"

#include "../../corelib/frag_id.hpp"
#include "../../corelib/hbn_package_version.h"
#include "../../corelib/pdqsort.h"
#include "../../sw/hbn_traceback_aux.h"

#include <algorithm>
#include <sstream>

using namespace std;

static void
bam_dump_hdr_hd(sam_hdr_t* hdr)
{
    int r = sam_hdr_add_line(hdr, "HD", "VN", SAM_VERSION, "SO", "unknown", "GO", "query", nullptr);
    if (r) HBN_ERR("Fail at adding BAM HD header line");
}

static void
bam_dump_hdr_sq(const char* name, const int size, sam_hdr_t* hdr)
{
    char sizebuf[256];
    snprintf(sizebuf, 256, "%d", size);
    int r = sam_hdr_add_line(hdr, "SQ", "SN", name, "LN", sizebuf, nullptr);
    if (r) HBN_ERR("Fail at adding BAM SQ header line");
}

static void 
bam_dump_hdr_pg(const char* pg_id, const char* pg_name, const char* pg_version, int argc, char* argv[], sam_hdr_t* hdr)
{
    char* cl = stringify_argv(argc, argv);
    int r = sam_hdr_add_line(hdr, "PG", "ID", pg_id, "PN", pg_name, "VN", pg_version, "CL", cl, nullptr);
    free(cl);
    if (r) HBN_ERR("FAIL at adding BAM PG header");
}

static void
extract_pg_name(const char* argv_0, char pg_name[])
{
    int n = strlen(argv_0);
    int s = n;
    while (s) {
        --s;
        if (argv_0[s] == '/') {
            ++s;
            break;
        }
    }

    int i = 0;
    for (; s < n; ++i, ++s) pg_name[i] = argv_0[s];
    pg_name[i] = '\0';
}

static void init_hbn_bam_hdr(int argc, char* argv[], HbnUnpackedDatabase* subjects, sam_hdr_t* hdr)
{
    bam_dump_hdr_hd(hdr);

    int n = subjects->num_targets();
    for (int i = 0; i < n; ++i) {
        const char* name = subjects->target_name(i);
        const int size = subjects->target_size(i);
        bam_dump_hdr_sq(name, size, hdr);
    }

    char pg_name[256];
    extract_pg_name(argv[0], pg_name);
    bam_dump_hdr_pg(pg_name, pg_name, HBN_PACKAGE_VERSION, argc, argv, hdr);
}

static void 
dump_cigar_string(const char* qaln, const char* saln, const int aln_size, ostringstream& os)
{
    const char* q = qaln;
    const char* s = saln;
    int i = 0;
    while (i < aln_size) {
        char type = 'N';
        int cnt = 0;
        if (q[i] == GAP_CHAR) {  // delete from subject (gap in query)
            type = 'D';
            while (i < aln_size && q[i] == GAP_CHAR) {
                hbn_assert(s[i] != GAP_CHAR);
                ++i;
                ++cnt;
            }
        } else if (s[i] == GAP_CHAR) { // insert into subject (gap in subject)
            type = 'I';
            while (i < aln_size && s[i] == GAP_CHAR) {
                hbn_assert(q[i] != GAP_CHAR);
                ++i;
                ++cnt;
            }
        } else if (q[i] == s[i]) {
            type = '=';
            while (i < aln_size && q[i] == s[i]) {
                ++i;
                ++cnt;
            }
        } else { 
            type = 'X';
            while (i < aln_size) {
                bool r = (q[i] == GAP_CHAR) || (s[i] == GAP_CHAR) || (q[i] == s[i]);
                if (r) break;
                ++i;
                ++cnt;
            }
        }
        os << cnt << type;
    }
}

static void
dump_sam_cigar(const int qoff, const int qend, const int qsize,
    const char* qaln, const char* saln, const int aln_size, ostringstream& os)
{
    if (qoff) os << qoff << 'S';
    dump_cigar_string(qaln, saln, aln_size, os);
    if (qend < qsize) os << qsize - qend << 'S';
}

static BOOL 
validate_cigar_and_seq_size(const char* cigar, const int query_size)
{
    const int c_n = strlen(cigar);
    int c_i = 0;
    int q_base = 0;
    while (c_i < c_n) {
        size_t j = c_i + 1;
        while (j < c_n && isdigit(cigar[j])) ++j;
        char op = cigar[j];
        int num = atoi(cigar + c_i);
        switch (op)
        {
        case 'S':
            q_base += num;
            break;
        case 'M':
        case '=':
        case 'X':
            q_base += num;
            break;
        case 'D':
            break;
        case 'I':
            q_base += num;
            break;
        default:
            break;
        }
        c_i = j + 1;
    }    
    if (q_base != query_size) {
        fprintf(stderr, "%s\n", cigar);
        HBN_LOG("cigar and seq_size is inconsistent: %d v.s. %d", q_base, query_size);
        return false;
    }
    return true;
}

static bam1_t*
build_one_read_bam(RestrictEnzymeLociList* reloci_list,
        HbnUnpackedDatabase* subjects,
        const char* query_name,
        const int query_id,
        const char* fwd_query,
        const char* rev_query,
        const char* fwd_qv,
        const char* rev_qv,
        const int query_size,
        PoreCAlign* all_pca_a,
        int all_pca_c,
        TrimPcaList* tpca_list,
        TrimPca* tpca,
        const bool is_primary_map)
{
    PoreCAlign* pca = &tpca->pca;
    uint16_t flag = 0; if (tpca->pca.sdir == REV) flag |= 0x10; if (!is_primary_map) flag |=0x800;
    int qdir = tpca->pca.sdir;
    int sid = tpca->pca.sid;

    const char* qas = tpca_list->align_strings.c_str() + tpca->qas_offset;
    const char* sas = tpca_list->align_strings.c_str() + tpca->sas_offset;
    int as_size = tpca->as_size;

    if (1) {
        const char* Q = (qdir == FWD) ? fwd_query : rev_query;
        const char* S = subjects->target_sequence(sid, FWD);
        validate_aligned_string_0(__FILE__, __FUNCTION__, __LINE__,
            query_id, Q, pca->qoff, pca->qend, qas,
            sid, S, pca->soff, pca->send, sas, as_size, TRUE);
    }

    float pi = pca->pi;
    ostringstream os;
    dump_sam_cigar(pca->qoff, pca->qend, query_size, qas, sas, as_size, os);
    string cigar_s = os.str();
    validate_cigar_and_seq_size(cigar_s.c_str(), query_size);

    uint32_t* cigar_a = nullptr;
    size_t cigar_buffer_size = 0;
    auto num_cigar_op = sam_parse_cigar(cigar_s.c_str(), nullptr, &cigar_a, &cigar_buffer_size);
    if (num_cigar_op < 0) HBN_ERR("FAIL at parsing CIGAR");

    const char* query = nullptr;
    const char* qv = nullptr;
    if (qdir == REV) {
        query = rev_query;
        if (rev_qv) qv = rev_qv;
    } else {
        query = fwd_query;
        if (fwd_qv) qv = fwd_qv;
    }

    bam1_t* bam = bam_init1();
    int r = bam_set1(bam, strlen(query_name), query_name, flag, 
        sid, pca->soff, tpca->pca.map_q, num_cigar_op, cigar_a,
        -1, -1, 0, query_size, query, qv, 512);
    if (r < 0) HBN_ERR("FAIL at building BAM record");
    free(cigar_a);

    char tagname[2];
    /// chaining score
    tagname[0] = 'd'; tagname[1] = 's';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&pca->ddf_score));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// qid
    tagname[0] = 'q'; tagname[1] = 'i';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&query_id));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// qdir
    tagname[0] = 'q'; tagname[1] = 'd';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&qdir));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// qs
    tagname[0] = 'q'; tagname[1] = 's';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&pca->qoff));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// qe
    tagname[0] = 'q'; tagname[1] = 'e';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&pca->qend));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// vqs
    tagname[0] = 'q'; tagname[1] = 'S';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&pca->enzyme_qoff));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// vqe
    tagname[0] = 'q'; tagname[1] = 'E';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&pca->enzyme_qend));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// ql
    tagname[0] = 'q'; tagname[1] = 'l';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&query_size));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// sid
    tagname[0] = 's'; tagname[1] = 'i';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&sid));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// ss
    tagname[0] = 's'; tagname[1] = 's';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&pca->soff));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// se
    tagname[0] = 's'; tagname[1] = 'e';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&pca->send));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// vss
    tagname[0] = 'v'; tagname[1] = 'S';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&pca->enzyme_soff));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// vse
    tagname[0] = 'v'; tagname[1] = 'E';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&pca->enzyme_send));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// sl
    tagname[0] = 's'; tagname[1] = 'l';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&pca->ssize));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    return bam;
}


static bam1_t*
build_one_frag_bam(RestrictEnzymeLociList* reloci_list,
        HbnUnpackedDatabase* subjects,
        const char* query_name,
        const int query_id,
        const char* fwd_query,
        const char* rev_query,
        const char* fwd_qv,
        const char* rev_qv,
        const int query_size,
        PoreCAlign* all_pca_a,
        int all_pca_c,
        TrimPcaList* tpca_list,
        TrimPca* tpca)
{
    PoreCAlign* pca = &tpca->pca;
    int sid = pca->sid;
    int qdir = pca->sdir;
    uint16_t flag = 0; if (qdir == REV) flag |= 0x10;

    const char* qas = tpca_list->align_strings.c_str() + tpca->qas_offset;
    const char* sas = tpca_list->align_strings.c_str() + tpca->sas_offset;
    int as_size = tpca->as_size;

    if (1) {
        const char* Q = (qdir == FWD) ? fwd_query : rev_query;
        const char* S = subjects->target_sequence(sid, FWD);
        validate_aligned_string_0(__FILE__, __FUNCTION__, __LINE__,
            query_id, Q, pca->qoff, pca->qend, qas,
            sid, S, pca->soff, pca->send, sas, as_size, TRUE);
    }

    char frag_prefix[FRAD_ID_SIZE];
    frag_id_to_string(query_id, tpca->frag_id, sid, pca->soff, frag_prefix);
    ostringstream os;
    os << query_name << '_' << frag_prefix;
    string frag_name = os.str();

    float pi = calc_ident_perc(qas, sas, as_size, nullptr, nullptr);
    os.str("");
    dump_cigar_string(qas, sas, as_size, os);
    string cigar_s = os.str();
    validate_cigar_and_seq_size(cigar_s.c_str(), pca->qend - pca->qoff);

    uint32_t* cigar_a = nullptr;
    size_t cigar_buffer_size = 0;
    auto num_cigar_op = sam_parse_cigar(cigar_s.c_str(), nullptr, &cigar_a, &cigar_buffer_size);
    if (num_cigar_op < 0) HBN_ERR("FAIL at parsing CIGAR");

    const char* query = (qdir == FWD) ? fwd_query + pca->qoff : rev_query + pca->qoff;
    const char* qv = nullptr;
    if (qdir == FWD && fwd_qv) qv = fwd_qv + pca->qoff;
    if (qdir == REV && rev_qv) qv = rev_qv + pca->qoff;

    bam1_t* bam = bam_init1();
    int r = bam_set1(bam, frag_name.size(), frag_name.c_str(), flag, 
        sid, pca->soff, pca->map_q, num_cigar_op, cigar_a,
        -1, -1, 0, pca->qend - pca->qoff, query, qv, 512);
    if (r < 0) HBN_ERR("FAIL at building BAM record");
    free(cigar_a);

    char tagname[2];
    /// chaining score
    tagname[0] = 'd'; tagname[1] = 's';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&pca->ddf_score));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// qid
    tagname[0] = 'q'; tagname[1] = 'i';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&query_id));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// qdir
    tagname[0] = 'q'; tagname[1] = 'd';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&qdir));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// qs
    tagname[0] = 'q'; tagname[1] = 's';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&pca->qoff));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// qe
    tagname[0] = 'q'; tagname[1] = 'e';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&pca->qend));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// vqs
    tagname[0] = 'q'; tagname[1] = 'S';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&tpca->pca.enzyme_qoff));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// vqe
    tagname[0] = 'q'; tagname[1] = 'E';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&tpca->pca.enzyme_qend));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// ql
    tagname[0] = 'q'; tagname[1] = 'l';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&query_size));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// sid
    tagname[0] = 's'; tagname[1] = 'i';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&sid));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// ss
    tagname[0] = 's'; tagname[1] = 's';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&pca->soff));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// se
    tagname[0] = 's'; tagname[1] = 'e';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&pca->send));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// vss
    tagname[0] = 'v'; tagname[1] = 'S';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&pca->enzyme_soff));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// vse
    tagname[0] = 'v'; tagname[1] = 'E';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&pca->enzyme_send));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// sl
    tagname[0] = 's'; tagname[1] = 'l';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&pca->ssize));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    return bam;
}

void
dump_one_pca(RestrictEnzymeLociList* reloci_list,
    HbnUnpackedDatabase* subjects,
    const char* query_name,
    const int query_id,
    const int query_size,
    PoreCAlign* all_pca_a,
    int all_pca_c,
    TrimPcaList* tpca_list,
    TrimPca* tpca,
    ostringstream& out)
{
    PoreCAlign* pca = &tpca->pca;
    const char* subject_name = subjects->target_name(pca->sid);
    const int subject_size = subjects->target_size(pca->sid);
    int qdir = pca->sdir;
    int qb = pca->qoff;
    int qe = pca->qend;
    int eqb = pca->enzyme_qoff;
    int eqe = pca->enzyme_qend;
    if (qdir == REV) {
        qb = pca->qsize - pca->qend;
        qe = pca->qsize - pca->qoff;
        eqb = pca->qsize - pca->enzyme_qend;
        eqe = pca->qsize - pca->enzyme_qoff;
    }

    out << query_name;
    out << '\t' << query_size;
    out << '\t' << qb;
    out << '\t' << qe;
    out << '\t' << ((qdir == FWD) ? '+' : '-');

    out << '\t' << subject_name;
    out << '\t' << subject_size;
    out << '\t' << pca->soff;
    out << '\t' << pca->send;

    const char* qas = tpca_list->align_strings.c_str() + tpca->qas_offset;
    const char* sas = tpca_list->align_strings.c_str() + tpca->sas_offset;
    const int as_size = tpca->as_size;
    int mat = 0;
    for (int i = 0; i < as_size; ++i) mat += (qas[i] == sas[i]);

    out << '\t' << mat;
    out << '\t' << as_size;
    out << '\t' << tpca->pca.map_q;

    out << '\t' << "ds:i:" << pca->ddf_score;

    out << '\t' << "qS:i:" << eqb;
    out << '\t' << "qE:i:" << eqe;

    out << '\t' << "vS:i:" << pca->enzyme_soff;
    out << '\t' << "vE:i:" << pca->enzyme_send;

    out << '\t' << "qi:i:" << query_id;

    out << '\n';
}

void HbnOutputs::x_init_sam_hdr(int argc, char* argv[], HbnUnpackedDatabase* reference)
{
    M_sam_hdr = sam_hdr_init();
    init_hbn_bam_hdr(argc, argv, reference, M_sam_hdr);
    int r = sam_hdr_write(M_sam_out, M_sam_hdr);
    if (r) HBN_ERR("FAIL at writing BAM header to file");
}

void HbnOutputs::dump(RestrictEnzymeLociList* reloci_list, 
        HbnUnpackedDatabase* subjects,
        const char* query_name,
        const int query_id,
        const char* fwd_query,
        const char* rev_query,
        const char* fwd_qv,
        const char* rev_qv,
        const int query_size,
        PoreCAlign* all_pca_a,
        int all_pca_c,
        TrimPcaList& pca_list)
{
    TrimPca* a = pca_list.tpca_list.data();
    int c = pca_list.tpca_list.size();
    for (int i = 0; i < c; ++i) a[i].frag_id = i;

    if (M_outfmt == eOutputFmt_BAM || M_outfmt == eOutputFmt_SAM) {
        vector<bam1_t*> bam_list;
        bam1_t* bam;
        int primary_map = -1, max_sc = 0;
        for (int i = 0; i < c; ++i) if (a[i].pca.map_score > max_sc) { max_sc = a[i].pca.map_score; primary_map = i; }
        for (int i = 0; i < c; ++i) {
            bam = build_one_read_bam(reloci_list, subjects, query_name, query_id,
                fwd_query, rev_query, fwd_qv, rev_qv, query_size, all_pca_a, all_pca_c, &pca_list, a + i, (i == primary_map));
            if (!bam) continue;
            bam_list.push_back(bam);
        }
        lock_guard<mutex> lg(M_out_mutex);
        for (auto b : bam_list) {
            int r = sam_write1(M_sam_out, M_sam_hdr, b);
            if (r < 0) HBN_ERR("FAIL at writing BAM record");
            bam_destroy1(b);
        }
    } else if (M_outfmt == eOutputFmt_FragBAM || M_outfmt == eOutputFmt_FragSAM) {
        vector<bam1_t*> bam_list;
        bam1_t* bam;
        for (int i = 0; i < c; ++i) {
            bam = build_one_frag_bam(reloci_list, subjects, query_name, query_id,
                fwd_query, rev_query, fwd_qv, rev_qv, query_size, all_pca_a, all_pca_c, &pca_list, a + i);
            if (!bam) continue;
            bam_list.push_back(bam);
        }
        lock_guard<mutex> lg(M_out_mutex);
        for (auto b : bam_list) {
            int r = sam_write1(M_sam_out, M_sam_hdr, b);
            if (r < 0) HBN_ERR("FAIL at writing BAM record");
            bam_destroy1(b);
        }        
    } else if (M_outfmt == eOutputFmt_PAF) {
        ostringstream os;
        for (int i = 0; i < c; ++i) {
            dump_one_pca(reloci_list, subjects, query_name, query_id, query_size, all_pca_a, all_pca_c, &pca_list, a+i, os);
        }
        lock_guard<mutex> lg(M_out_mutex);
        hbn_fwrite(os.str().c_str(), 1, os.str().size(), M_out);
    }
}