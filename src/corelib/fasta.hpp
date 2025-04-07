#ifndef __FASTA_HPP
#define __FASTA_HPP

#include "line_reader.hpp"

#include <string>

class CFastaReader
{
public:
    CFastaReader(const char* path, const bool ignore_ill_formated_seq = true):
        M_LineReader(path),
        M_IgnoreIllFormatedSequence(ignore_ill_formated_seq) {}

    ~CFastaReader() {}

    int ReadOneSeq(std::string& name, std::string& comment, std::string& sequence, std::string& plus, std::string& qual);

    int ReadOneSeq();

    std::string& sequence() {
        return M_Sequence;
    }

    std::string& name() {
        return M_Name;
    }

    std::string& comment() {
        return M_Comment;
    }

    std::string& plus() {
        return M_Plus;
    }

    std::string& qual() {
        return M_Qual;
    }

    const std::string& sequence() const {
        return M_Sequence;
    }

    const std::string& name() const {
        return M_Name;
    }

    const std::string& comment() const {
        return M_Comment;
    }

    const std::string& plus() const {
        return M_Plus;
    }

    const std::string& qual() const {
        return M_Qual;
    }

private:
    bool ParseDefLine(CTempString line, std::string& name, std::string& comment);

    bool CheckDataLine(CTempString line);

    bool ParseDataLine(CTempString line, std::string& sequence);

    void AdvanceToNextSeqHeader();

    bool ParseFastqSeq(std::string& name, std::string& comment, std::string& sequence, std::string& plus, std::string& qual);

private:
    HbnLineReader   M_LineReader;
    std::string     M_Name;
    std::string     M_Comment;
    std::string     M_Sequence;
    std::string     M_Plus;
    std::string     M_Qual;
    bool            M_IgnoreIllFormatedSequence;
};

typedef CFastaReader HbnFastaReader;

#endif // __FASTA_HPP