#include "hbn_aux.h"

#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include <zlib.h>
#include <ctype.h>
#include <locale.h>

static const char* hbn_log_level_name[] = {
    [LOG_EMERG]     = "EMERG",
    [LOG_ALERT]     = "ALERT",
    [LOG_CRIT]      = "CRIT",
    [LOG_ERR]       = "ERROR",
    [LOG_WARNING]   = "WARNING",
    [LOG_NOTICE]    = "NOTICE",
    [LOG_INFO]      = "INFO",
    [LOG_DEBUG]     = "DEBUG",
};

void
what_is_time_now(char now[])
{
    time_t ltime;
    time(&ltime);
    ctime_r(&ltime, now);
    size_t n = strlen(now);
    now[n-1] = '\0';
}

static void 
_hbn_dump_message(HBN_LOG_PARAMS_GENERIC, const int level, const char* fmt, va_list args)
{
    char buf[HBN_MAX_PATH_LEN + 1];
    char time_str[256];
    what_is_time_now(time_str);
    vsnprintf(buf, sizeof(buf), fmt, args);
    fprintf(stderr, "[%s %s:%s:%d] <%s> %s\n", time_str, HBN_LOG_ARGS_GENERIC, hbn_log_level_name[level], buf);
}

void
hbn_dump_message(HBN_LOG_PARAMS_GENERIC, const int level, const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    _hbn_dump_message(HBN_LOG_ARGS_GENERIC, level, fmt, args);
    va_end(args);
}

static void
_hbn_dump_message_short(const int level, const char* fmt, va_list args)
{
    char buf[HBN_MAX_PATH_LEN + 1];
    char time_str[256];
    what_is_time_now(time_str);
    vsnprintf(buf, sizeof(buf), fmt, args);
    fprintf(stderr, "[%s NecatCore] <%s> %s\n", time_str, hbn_log_level_name[level], buf);
}

void
hbn_dump_message_short(const int level, const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    _hbn_dump_message_short(level, fmt, args);
    va_end(args);    
}

////////////////////////////////////////////////////////////////////////////////////////////

int safe_gzread(HBN_LOG_PARAMS_GENERIC, gzFile stream, void* buf, unsigned int len)
{
    int ret = gzread(stream, buf, len);
    if (ret < 0) {
        int errnum = 0;
        const char* msg = gzerror(stream, &errnum);
        const char* why = (Z_ERRNO == errnum) ? strerror(errno) : msg;
        HBN_ERR_GENERIC("%s", why);
    }
    return ret;
}

int err_gzread(gzFile stream, void* buf, unsigned int len)
{
    return safe_gzread(HBN_LOG_ARGS_DEFAULT, stream, buf, len);
}

gzFile safe_gzopen(HBN_LOG_PARAMS_GENERIC, const char* path, const char* mode)
{
    gzFile stream;
    if (strcmp(path, "-") == 0) {
        stream = gzdopen(fileno(strstr(mode, "r") ? stdin : stdout), mode);
        /* according to zlib.h, this is the only reason gzdopen can fail */
        if (!stream) HBN_ERR_GENERIC("fail to open file '%s' with mode '%s': out of memory", path, mode);
        return stream;
    }
    if ((stream = gzopen(path, mode)) == 0) {
        const char* why = errno ? strerror(errno) : "out of memory";
        HBN_ERR_GENERIC("fail to open file '%s' with mode '%s': %s", path, mode, why);
    }
    return stream;
}

int safe_gzclose(HBN_LOG_PARAMS_GENERIC, gzFile stream)
{
    int ret = gzclose(stream);
    if (Z_OK != ret) {
        const char* why = (Z_ERRNO == ret) ? strerror(errno) : zError(ret);
        HBN_ERR_GENERIC("%s", why);
    }
    return ret;
}

FILE* safe_fopen(HBN_LOG_PARAMS_GENERIC, const char* path, const char* mode)
{
    FILE* stream = 0;
    if (strcmp(path, "-") == 0) return strstr(mode, "r") ? stdin : stdout;

    if ((stream = fopen(path, mode)) == 0) {
        const char* y = strerror(errno);
        HBN_ERR_GENERIC("fail to open file '%s' with mode '%s': %s", path, mode, y);
    }
    return stream;
}

size_t safe_fwrite(HBN_LOG_PARAMS_GENERIC, const void* buf, size_t size, size_t nmemb, FILE* stream)
{
    size_t ret = fwrite(buf, size, nmemb, stream);
    if (ret != nmemb) {
        HBN_ERR_GENERIC("%s", strerror(errno));
    }
    return ret;
}

size_t safe_fread(HBN_LOG_PARAMS_GENERIC, void* buf, size_t size, size_t nmemb, FILE* stream)
{
    size_t ret = fread(buf, size, nmemb, stream);
    if (ret != nmemb) {
        const char* y = ferror(stream) ? strerror(errno) : "Unexpected end of file";
        HBN_ERR_GENERIC("%s", y);
    }
    return ret;
}

int safe_fclose(HBN_LOG_PARAMS_GENERIC, FILE* stream)
{
    int ret = fclose(stream);
    if (ret != 0) {
        HBN_ERR_GENERIC("%s", strerror(errno));
    }
    return ret;
}

void hbn_exception(const char* expr, HBN_LOG_PARAMS_GENERIC, const char* fmt, ...)
{
    fprintf(stderr, "Assertion Failed At '%s:%s:%d'\n", file, func, line);
    fprintf(stderr, "\tExpression: '%s'\n", expr);
    if (!fmt) return;
    fprintf(stderr, "Context Information:\n");
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fprintf(stderr, "\n");
    abort();
}

off_t hbn_get_file_size(HBN_LOG_PARAMS_GENERIC, const char* path)
{
    struct stat sbuf;
    if (stat(path, &sbuf) == -1) {
        const char* y = strerror(errno);
        HBN_ERR_GENERIC("fail to stat file '%s': %s", path, y);
    }
    return sbuf.st_size;
}

double hbn_time_diff(const struct timeval* begin, const struct timeval* end)
{
    double d = end->tv_sec - begin->tv_sec;
    d += 1.0 * (end->tv_usec - begin->tv_usec) / 1e6;
    return d;
}

u8 nst_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

u8 nst_nt16_table[256] = {
    16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, //15
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, // 31
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 15, 16, 16, // 47
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, // 63
	16, 0, 10, 1, 11, 16, 16, 2, 12, 16, 16, 7, 16, 6, 14, 16,
	16, 16, 4, 9, 3, 16, 13, 8, 16, 5, 16, 16, 16, 16, 16, 16,
	16, 0, 10, 1, 11, 16, 16, 2, 12, 16, 16, 7, 16, 6, 14, 16,
	16, 16, 4, 9, 3, 16, 13, 8, 16, 5, 16, 16, 16, 16, 16, 16,
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16
};

unsigned char complement_base_table[256] = {
	  0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
	 96, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127,
	128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
	144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
	160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
	176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
	192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
	208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
	224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
	240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255
};

void i64_to_string_with_comma_r(i64 n, char str_n[])
{
    static int comma = '\0';
    char rebuf[64];
    char* p = &rebuf[ sizeof(rebuf)-1 ];
    int i = 0;
    if (comma == '\0') {
        struct lconv* lcp = localeconv();
        if (lcp != NULL) {
            if (lcp->thousands_sep != NULL && *lcp->thousands_sep != '\0') {
                comma = *lcp->thousands_sep;
            } else {
                comma = ',';
            }
        }
    }
    int n_is_nonnegative = 1;
    if (n < 0) {
        n_is_nonnegative = 0;
        n *= -1;
    }
    hbn_assert(n >= 0);
    *p = '\0';
    do {
        if (i && i % 3 == 0) *--p = comma;
        *--p = '0' + n % 10;
        n /= 10;
        ++i;
    } while (n);
    if (!n_is_nonnegative) *--p = '-';

    char* q = str_n;
    while (*p != '\0') {
        *q = *p;
        ++q;
        ++p;
    }
    *q = '\0';
}

char* i64_to_string_with_comma(i64 n)
{
    static char n_str[64];
    i64_to_string_with_comma_r(n, n_str);
    return n_str;
}

void u64_to_fixed_width_string_r(u64 n, char n_str[], const int width)
{
    int n_digit = 0;
    size_t m = n;

    do {
        ++n_digit;
        n /= 10;
    } while (n);
    hbn_assert(n_digit <= width, "n = %zu, n_digit = %d, width = %d", m, n_digit, width);

    char* p = n_str;
    int i = n_digit;
    while (i < width) {
        *p = '0';
        ++p;
        ++i;
    }
    sprintf(p, "%zu", m);
}

char*
u64_to_fixed_width_string(u64 n, const int width)
{
    static char n_str[64];
    u64_to_fixed_width_string_r(n, n_str, width);
    return n_str;
}

void
print_digit_with_comma(FILE* out, size_t num)
{
    char buf[64];
    i64_to_string_with_comma_r(num, buf);
    fprintf(out, "%s", buf);
}

void
print_fixed_width_string(FILE* out, const char* s, const int width)
{
    int n = strlen(s);
    fprintf(out, "%s", s);
    for (int i = n; i < width; ++i) fprintf(out, " ");
}

void
print_fixed_width_string_r(const char* s, const int width, char out[])
{
    int n = strlen(s);
    int i = 0;
    for (; i < n; ++i) out[i] = s[i];
    for (; i < width; ++i) out[i] = ' ';
    out[i] = '\0';
}

int
string_is_valid_number(const char * arg)
{
    char state = 1;
    for (char c = *arg; c; c = *(++arg)){
        switch (c) {
            case 'e':
                if (state != 2 && state != 14) return 0;
                state = 9;
                break;
            case '.':
                if (state >> 2) return 0;
                state = 12;
                break;
            case '+':
            case '-':
                if (!(state & 1)) return 0;
                state--;
                break;
            default:
                if (c >= '0' && c <= '9'){
                    if (state >> 2 == 1) return 0;
                    state = "\2\4\012\016"[state >> 2];
                } else {
                    return 0;
                }
                break;
        }
    }
    return state & 2;
}

#ifdef WIN32
#include "windows.h"
#else
#include "unistd.h"
#endif

int hbn_get_cpu_count()
{
#if WIN32
    SYSTEM_INFO sysInfo;
    GetSystemInfo(&sysInfo);
    int allCPUNum_    = sysInfo.dwNumberOfProcessors;
    return allCPUNum_;
#else //linux
    //获取当前系统的所有CPU核数，包含禁用的
    //int allCPUNum_    = sysconf(_SC_NPROCESSORS_CONF);
    //获取当前系统的可用CPU核数
    int enableCPUNum_ = sysconf(_SC_NPROCESSORS_ONLN);
    return enableCPUNum_;
#endif
}

void create_directory(const char* path)
{
    if ((access(path, F_OK) != 0)
        &&
        (mkdir(path, S_IRWXU) != 0)) {
        HBN_ERR("Failed to create directory %s: %s", path, strerror(errno));
    }
}

#ifdef __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
u64 system_ram_bytes()
{
    int mib[2] = { CTL_HW, HW_MEMSIZE };
    u_int namelen = sizeof(mib) / sizeof(mib[0]);
    u64 size;
    size_t len = sizeof(size);
    if (sysctl(mib, namelen, &size, &len, NULL, 0) < 0) {
        perror("sysctl");
    }
    return size;
}
#endif 

#ifdef __linux__
#include <sys/sysinfo.h>
u64 system_ram_bytes()
{
    struct sysinfo info;
    int r = sysinfo(&info);
    if (r != 0) perror("sysinfo");
    return info.totalram;
}
#endif

u64 system_ram_gb()
{
    u64 bytes = system_ram_bytes();
    return bytes >> 30;
}

void
necat_core_dump()
{
    char s_time[256];
    what_is_time_now(s_time);
    fprintf(stderr, "[%s, NecatCore] ", s_time);
}

i64 simple_atoi_n(const char* s, const int sl)
{
    i64 n = 0;
    if (sl == 0) return n;

    int sign = 1;
    int i = 0;
    if (s[i] == '-') {
        sign = -1;
        i = 1;
    }

    for (; i < sl; ++i) {
        int r = (s[i] >= '0' && s[i] <= '9');
        if (!r) break;
        n = n * 10 + s[i] - '0';
    }

    n *= sign;
    return n;
}

double simple_atof_n(const char* s, const int sl)
{
    if (sl == 0) return 0.0;

    int sign = 1;
    int i = 0;
    if (s[i] == '-') {
        sign = -1;
        i = 1;
    }

    i64 x = 0;
    for (; i < sl; ++i) {
        if (s[i] == '.') break;
        int r = (s[i] >= '0' && s[i] <= '9');
        if (!r) break;
        x = x * 10 + s[i] - '0';
    }

    i64 y = 0;
    i64 div = 1;
    for (; i < sl; ++i) {
        if (s[i] == '.') continue;
        int r = (s[i] >= '0' && s[i] <= '9');
        if (!r) break;
        y = y * 10 + s[i] - '0';
        div *= 10;
    }
    double z = 1.0 * y / div;

    z += x;
    if (sign == -1) z = -z;

    return z;
}