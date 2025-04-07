#include "fasta.hpp"

#include <algorithm>

using namespace std;

// The FASTA reader uses these heavily, but the standard versions
// aren't inlined on as many configurations as one might hope, and we
// don't necessarily want locale-dependent behavior anyway.

inline bool s_ASCII_IsUpper(unsigned char c)
{
    return c >= 'A'  &&  c <= 'Z';
}

inline bool s_ASCII_IsLower(unsigned char c)
{
    return c >= 'a'  &&  c <= 'z';
}

inline bool s_ASCII_IsAlpha(unsigned char c)
{
    return s_ASCII_IsUpper(c)  ||  s_ASCII_IsLower(c);
}

// The arg *must* be a lowercase letter or this won't work
inline unsigned char s_ASCII_MustBeLowerToUpper(unsigned char c)
{
    return c + ('A' - 'a');
}

inline bool s_ASCII_IsAmbigNuc(unsigned char c)
{
    switch(c) {
    case 'U': case 'u':
    case 'R': case 'r':
    case 'Y': case 'y':
    case 'S': case 's':
    case 'W': case 'w':
    case 'K': case 'k':
    case 'M': case 'm':
    case 'B': case 'b':
    case 'D': case 'd':
    case 'H': case 'h':
    case 'V': case 'v':
    case 'N': case 'n':
        return true;
    default:
        return false;
    }
}

inline static bool s_ASCII_IsUnAmbigNuc(unsigned char c)
{
    switch( c ) {
    case 'A': case 'C': case 'G': case 'T':
    case 'a': case 'c': case 'g': case 't':
        return true;
    default:
        return false;
    }
}

bool CFastaReader::ParseDefLine(CTempString line, std::string& name, std::string& comment)
{
    name.clear();
    comment.clear();

    const char* s = line.data();
    size_t n = line.size();
    hbn_assert(s[0] == '>' || s[0] == '@');
    if (n == 1) {
        HBN_WARN("HbnFastaReader: Defline around (%zu, %s) lacks a proper ID",
            M_LineReader.GetLineNumber(), M_LineReader.GetFileName());
        std::cerr << line << std::endl;
        return FALSE;
    }
    size_t i = 1;
    while (i < n && (!isspace(s[i]))) ++i;
    name.assign(&s[1], i-1);

    while (i < n && isspace(s[i])) ++i;
    if (i < n) comment.assign(&s[i], n-i);

    return true;
}

bool CFastaReader::CheckDataLine(CTempString line)
{
    // make sure the first data line has at least SOME resemblance to
    // actual sequence data
    if (line.empty()) return TRUE;

    size_t good = 0, bad = 0;
    // in case the data has huge sequences all on the first line we do need
    // a cutoff and "70" seems reasonable since it's the default width of
    // CFastaOstream (as of 2017-03-09)
    size_t len_to_check = min<size_t>(line.size(), 70);
    size_t ambig_nuc = 0;
    const char* s = line.data();
    for (size_t pos = 0; pos < len_to_check; ++pos) {
        unsigned char c = s[pos];
        if (s_ASCII_IsAlpha(c) || c == '*') {
            ++good;
            if (s_ASCII_IsAmbigNuc(c)) {
                ++ambig_nuc;
            }
        } else if (c == '-') {
            ++good;
        } else if (isspace(c) || (c >= '0' && c <= '9')) {
            // treat whilespace and digits as neutral
        } else if (c == ';') {
            break; // comment - ignore rest of line
        } else {
            ++bad;
        }
    }
    if (bad >= good / 3 &&
        (len_to_check > 3 || good == 0 || bad > good)) {
        HBN_WARN("HbnFastaReader: Near (%zu, %s), there's a line that doesn't look like"
            " plausible data, but it's not marked as defline or comment.",
            M_LineReader.GetLineNumber(), M_LineReader.GetFileName());
        return FALSE;        
    }
    // warn if more than a certain percentage is ambiguous nucleotides
    const static size_t kWarnPercentAmbiguous = 40; // e.g. "40" means "40%"
    const size_t percent_ambig = (good == 0)?100:((ambig_nuc * 100) / good);
    if( len_to_check > 3 && percent_ambig > kWarnPercentAmbiguous ) {
#if 0
        HBN_WARN( 
            "HbnFastaReader: Start of first data line (%zu, %s) in seq is about "
            "%zu%% ambiguous nucleotides (shouldn't be over %zu%%)",
            HbnLineReaderLineNumber(reader->line_reader),
            reader->filename,
            percent_ambig,
            kWarnPercentAmbiguous);
#endif
    }
    return TRUE;
}

bool CFastaReader::ParseDataLine(CTempString line, std::string& sequence)
{
    if (!CheckDataLine(line)) return false;
    // most lines won't have a comment (';') 
    // so optimize for that case as much as possible
    const size_t s_len = line.size();
    for (size_t pos = 0; pos < s_len; ++pos) {
        const unsigned char c = line[pos];
        if (c == ';') break;
        if (c == '\t' || c == '\v' || c == '\f' || c == '\r' || c == ' ') continue;
        if (nst_nt16_table[c] == 16) {
            char msg[512];
            if (isprint(c)) {
                snprintf(msg, 512, "HbnFastaReader: "
                    "There are invalid nucleotide residue(s) '%c' in input sequence around (%zu, %zu, %s)", 
                    c, pos + 1, M_LineReader.GetLineNumber(), M_LineReader.GetFileName());
            } else {
                snprintf(msg, 512, "HbnFastaReader: "
                    "There are invalid nucleotide residue(s) (ascii code: %d) in input sequence around (%zu, %zu, %s)", 
                    c, pos + 1, M_LineReader.GetLineNumber(), M_LineReader.GetFileName());
            }
            HBN_WARN(msg);
            return false;
        }
        sequence += c;
    }
    return true;
}

void CFastaReader::AdvanceToNextSeqHeader()
{
    while (M_LineReader.ReadOneLine()) {
        CTempString line = *M_LineReader;
        if (line.empty()) continue;
        int c = line[0];
        if (c == '>' || c == '@') {
            M_LineReader.UngetLine();
            break;
        }
    }
}

bool CFastaReader::ParseFastqSeq(std::string& name, std::string& comment, 
    std::string& sequence, std::string& plus, std::string& qual)
{
    /// header
    int find_header = 0;
    while (M_LineReader.ReadOneLine()) {
        CTempString line = *M_LineReader;
        if (line.empty()) continue;
        if (line[0] == '@') {
            find_header = 1;
            break;
        }
    }
    if (!find_header) return TRUE;
    if (!ParseDefLine(*M_LineReader, name, comment)) return false;

    /// sequence
    int find_sequence = 0;
    while (M_LineReader.ReadOneLine()) {
        CTempString line = *M_LineReader;
        if (line.empty()) continue;
        find_sequence = 1;
        break;
    }
    if (!find_sequence) {
        HBN_WARN("HbnFastaReader: sequence data is missing from FASTQ sequence at around (%zu, %s)",
            M_LineReader.GetLineNumber(), M_LineReader.GetFileName());
        return FALSE;
    }
    //cerr << *m_LineReader << endl;
    if (!ParseDataLine(*M_LineReader, sequence)) return false;

    /// '+'
    int find_plus = 0;
    while (M_LineReader.ReadOneLine()) {
        CTempString line = *M_LineReader;
        if (line.empty()) continue;
        find_plus = 1;
        break;
    }
    if (!find_plus) {
        HBN_WARN("HbnFastaReader: plus data line (which shall start with '+') is missing from FASTQ sequence at around (%zu, %s)",
            M_LineReader.GetLineNumber(), M_LineReader.GetFileName());
        return FALSE;        
    }
    CTempString c_plus = *M_LineReader;
    plus.assign(c_plus.data(), c_plus.size());
    if (plus[0] != '+') {
        HBN_WARN("HbnFastaReader: detect an invalid plus data (which shall start with '+') "
            "line from FASTQ sequence at around (%zu, %s).",
            M_LineReader.GetLineNumber(), M_LineReader.GetFileName());
        return FALSE;
    }

    /// qual line
    int find_qual = 0;
    while (M_LineReader.ReadOneLine()) {
        CTempString line = *M_LineReader;
        if (line.empty()) continue;
        find_qual = 1;
        break;
    }
    if (!find_qual) {
        HBN_WARN("HbnFastaReader: "
            "quality scores line is missing from FASTQ sequence at around (%zu, %s)",
            M_LineReader.GetLineNumber(), M_LineReader.GetFileName());
        return FALSE;
    }
    CTempString c_qual = *M_LineReader;
    qual.assign(c_qual.data(), c_qual.size());
    if (qual.size() != sequence.size()) {
        HBN_WARN("HbnFastaReader: "
            "quality line and sequence line do not have the same length at around (%zu, %s).",
            M_LineReader.GetLineNumber(), M_LineReader.GetFileName());
        return FALSE;            
    }    

    return TRUE;
}

#define handle_read_error do { \
    if (M_IgnoreIllFormatedSequence) { \
        HBN_WARN("We skip this sequence and restart scanning from the next one."); \
        AdvanceToNextSeqHeader(); \
    } else { \
        HBN_ERR("We terminate the scanning process and exit with EXIT_FAILURE."); \
    } \
} while(0)

int CFastaReader::ReadOneSeq(std::string& name, std::string& comment, std::string& sequence, std::string& plus, std::string& qual)
{
    if (!M_LineReader.ReadOneLine()) return -1;
    M_LineReader.UngetLine();

    name.clear();
    comment.clear();
    sequence.clear();
    plus.clear();
    qual.clear();

    BOOL need_defline = TRUE;
    while (M_LineReader.ReadOneLine()) {
        CTempString line = *M_LineReader;
        if (line.empty()) continue; //  ignore lines containing only whitespace
        const char c = line[0];
        if (c == '>')  {
            if (need_defline) {
                if (!ParseDefLine(line, name, comment)) {
                    handle_read_error;
                    need_defline = TRUE;
                    continue;
                }
                need_defline = FALSE;
                continue;
            } else {
                M_LineReader.UngetLine();
                // start of the next sequence;
                break;
            }
        } else if (c == '@') {
            M_LineReader.UngetLine();
            if (!ParseFastqSeq(name, comment, sequence, plus, qual)) {
                handle_read_error;
                need_defline = TRUE;
                continue;
            }
            need_defline = FALSE;
            break;
        }

        if (c == '!' || c == '#' || c == ';') {
            // no content, just a comment or blank line
            continue;
        } else if (need_defline) {
            M_LineReader.UngetLine();
            HBN_WARN("HbnFastaReader: Input doesn't start with a defline or comment around (%zu, %s)",
                M_LineReader.GetLineNumber(), M_LineReader.GetFileName());
            handle_read_error;
            need_defline = TRUE;
        }

        if (!ParseDataLine(line, sequence)) {
            handle_read_error;
            need_defline = TRUE;
        }
    }

    if (need_defline) {
        HBN_WARN("HbnFastaReader: Expected defline around (%zu, %s)", M_LineReader.GetLineNumber(), M_LineReader.GetFileName());
        if (M_IgnoreIllFormatedSequence) {
            HBN_WARN("We skipt this error.");
            return FALSE;
        } else {
            HBN_ERR("We terminate the scanning process and exit with EXIT_FAILURE.");
        }
    }

    return sequence.size();
}

int CFastaReader::ReadOneSeq()
{
    return ReadOneSeq(M_Name, M_Comment, M_Sequence, M_Plus, M_Qual);
}