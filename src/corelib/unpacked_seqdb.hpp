#ifndef __UNPACKED_SEQDB_HPP
#define __UNPACKED_SEQDB_HPP

#include "hbn_aux.h"
#include "seq_name2id_map.hpp"

#include <string>
#include <vector>

class HbnUnpackedDatabase
{
public:

    struct TargetInfo
    {
        int id;
        const char* name;
        const char* fwd_seq;
        const char* rev_seq;
        int seq_size;

        size_t  _name_offset;
        size_t  _seq_offset;
    };

    HbnUnpackedDatabase(const char* fn);

    HbnUnpackedDatabase() {}

    int num_targets() const {
        return M_num_targets;
    }

    const char* target_name(const int id) const {
        x_validate_target_id(id);
        return M_target_list[id].name;
    }

    const char* target_sequence(const int id, const int strand) const {
        x_validate_target_id(id);
        return (strand == FWD) ? M_target_list[id].fwd_seq : M_target_list[id].rev_seq;
    }

    int target_size(const int id) const {
        x_validate_target_id(id);
        return M_target_list[id].seq_size;
    }

    int target_name2id(const char* name) const {
        return M_target_name2id.GetIdFromNameSafe(name);
    }

    const TargetInfo* target_info_list() const {
        return M_target_list.data();
    }

    void dump(FILE* out);
    void load(FILE* in);

private:

    void x_validate_target_id(const int id) const
    {
        if (id < 0 || id >= M_num_targets) {
            HBN_ERR("Target id %d is out of plausible range [%d, %d)", id, 0, M_num_targets);
        }
    }

private:
    int                     M_num_targets;
    std::vector<TargetInfo> M_target_list;
    std::vector<char>       M_name_list;
    std::vector<char>       M_fwd_seq_list;
    std::vector<char>       M_rev_seq_list;
    SeqName2IdMap           M_target_name2id;
};

#endif // __UNPACKED_SEQDB_HPP