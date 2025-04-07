#include "unpacked_seqdb.hpp"

#include "../ncbi_blast/str_util/ncbistr.hpp"
#include "fasta.hpp"

using namespace std;

HbnUnpackedDatabase::HbnUnpackedDatabase(const char* fn)
{
    HbnFastaReader in(fn);
    TargetInfo target;
    target.id = 0;
    while (in.ReadOneSeq() != -1) {
        if (in.sequence().empty()) continue;

        target._name_offset = M_name_list.size();
        M_name_list.insert(M_name_list.end(), in.name().begin(), in.name().end());
        M_name_list.push_back('\0');

        string& sequence = in.sequence();
        target.seq_size = sequence.size();

        target._seq_offset = M_fwd_seq_list.size();
        M_fwd_seq_list.insert(M_fwd_seq_list.end(), sequence.begin(), sequence.end());
        M_rev_seq_list.insert(M_rev_seq_list.end(), sequence.rbegin(), sequence.rend());

        M_target_list.push_back(target);
        ++target.id;
    }   
    M_num_targets = target.id;

    for (auto& c : M_fwd_seq_list) c = toupper(c);
    for (auto& c : M_rev_seq_list) c = toupper(c);
    for (auto& c : M_rev_seq_list) c = complement_base_table[c];

    for (auto& t : M_target_list) {
        t.name = M_name_list.data() + t._name_offset;
        t.fwd_seq = M_fwd_seq_list.data() + t._seq_offset;
        t.rev_seq = M_rev_seq_list.data() + t._seq_offset;
        M_target_name2id.add_one_name(t.name);
    }
    M_target_name2id.build_name2id_map();

    string size = NStr::UInt8ToString_DataSize(M_fwd_seq_list.size());
    HBN_LOG("Load %d target sequences (%s) from %s", target.id, size.c_str(), fn);
    size_t lb = 0;
    for (auto c : M_fwd_seq_list) if (c >= 'a' && c <= 'z') ++lb;
    fprintf(stderr, "Lower bases: %zu\n", lb);
}

void HbnUnpackedDatabase::dump(FILE* out)
{
    void* s = nullptr;
    size_t cnt;
    size_t size;

    s = (void*)(&M_num_targets);
    cnt = 1;
    size = sizeof(int);
    hbn_fwrite(s, size, cnt, out);

    s = (void*)(M_target_list.data());
    cnt = M_target_list.size();
    size = sizeof(TargetInfo);
    hbn_fwrite(s, size, cnt, out);

    size_t name_list_size = M_name_list.size();
    s = (void*)(&name_list_size);
    cnt = 1;
    size = sizeof(size_t);
    hbn_fwrite(s, size, cnt, out);

    s = (void*)(M_name_list.data());
    cnt = M_name_list.size();
    size = sizeof(char);
    hbn_fwrite(s, size, cnt, out);

    size_t num_bases = M_fwd_seq_list.size();
    s = (void*)(&num_bases);
    cnt = 1;
    size = sizeof(size_t);
    hbn_fwrite(s, size, cnt, out);

    s = (void*)(M_fwd_seq_list.data());
    cnt = M_fwd_seq_list.size();
    size = sizeof(char);
    hbn_fwrite(s, size, cnt, out);
}

void HbnUnpackedDatabase::load(FILE* in)
{
    void* s = nullptr;
    size_t cnt;
    size_t size;

    s = (void*)(&M_num_targets);
    cnt = 1;
    size = sizeof(int);
    hbn_fread(s, size, cnt, in);

    M_target_list.resize(M_num_targets);
    s = M_target_list.data();
    cnt = M_target_list.size();
    size = sizeof(TargetInfo);
    hbn_fread(s, size, cnt, in);

    size_t name_list_size;
    s = (void*)(&name_list_size);
    cnt = 1;
    size = sizeof(size_t);
    hbn_fread(s, size, cnt, in);

    M_name_list.resize(name_list_size);
    s = (void*)M_name_list.data();
    cnt = name_list_size;
    size = sizeof(char);
    hbn_fread(s, size, cnt, in);

    size_t num_bases = 0;
    s = (void*)(&num_bases);
    cnt = 1;
    size = sizeof(size_t);
    hbn_fread(s, size, cnt, in);

    M_fwd_seq_list.resize(num_bases);
    s = (void*)M_fwd_seq_list.data();
    cnt = num_bases;
    size = sizeof(char);
    hbn_fread(s, size, cnt, in);

    M_rev_seq_list.resize(M_fwd_seq_list.size());
    size_t offset = 0;
    for (int i = 0; i < M_num_targets; ++i) {
        const char* seq = M_fwd_seq_list.data() + M_target_list[i]._seq_offset;
        int seq_size = M_target_list[i].seq_size;
        char* rseq = M_rev_seq_list.data() + offset;
        copy(seq, seq + seq_size, rseq);
        reverse(rseq, rseq + seq_size);
        for (int s = 0; s < seq_size; ++s) rseq[s] = complement_base_table[(int)rseq[s]];
        offset += seq_size;
    }
    hbn_assert(offset == M_fwd_seq_list.size());

    M_target_name2id.clear();
    for (int i = 0; i < M_num_targets; ++i) {
        TargetInfo& t = M_target_list[i];
        t.name = M_name_list.data() + t._name_offset;
        t.fwd_seq = M_fwd_seq_list.data() + t._seq_offset;
        t.rev_seq = M_rev_seq_list.data() + t._seq_offset;
        M_target_name2id.add_one_name(t.name);
    }
    M_target_name2id.build_name2id_map();

    string datasize = NStr::UInt8ToString_DataSize(M_fwd_seq_list.size());
    HBN_LOG("Load %d target sequences (%s)", M_num_targets, datasize.c_str());
}
