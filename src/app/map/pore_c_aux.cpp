#include "pore_c_aux.hpp"

using namespace std;

static const double kSubOvlpFrac = 0.8;//0.8;

size_t gNumRescueMapQ = 0;
void dump_rescue_mapQ()
{
    fprintf(stderr, "rescue mapQ: %zu\n", gNumRescueMapQ);
}

static void
s_set_mapQ_ddf_sc(PoreCInitHit* a, int c)
{
    const int verbose = 0;
    pdqsort(a, a + c, [](const PoreCInitHit& x, const PoreCInitHit& y) { return x.ddf_score > y.ddf_score; });
    for (int i = 0; i < c; ++i) {
        a[i].id = i;
        a[i].parent = i;
        a[i].is_hom = false;
    }   

    for (int i = 1; i < c; ++i) {
	    PoreCInitHit* ai = a + i;
	    int max_ovlp = 0;
	    int max_j = -1;
	    for (int j = i - 1; j >= 0; --j) {
		    PoreCInitHit* aj = a + j;
		    int qb = max(ai->qoff, aj->qoff);
		    int qe = min(ai->qend, aj->qend);
		    int ovlp = qe - qb;
		    if (ovlp > max_ovlp) {
			    max_ovlp = ovlp;
			    max_j = j;
		    }
	    }
	    if (max_j == -1) continue;
	    if (max_ovlp < (ai->qend - ai->qoff) * kSubOvlpFrac) continue;
	    a[i].is_hom = true;
	    a[max_j].is_hom = true;
	    a[i].parent = a[max_j].parent;
    }

    for (int i = 0; i < c; ++i) {
        if (a[i].id != a[i].parent) {
            a[i].ddf_mapQ = 0;
            continue;
        }
        if (!a[i].is_hom) {
            a[i].ddf_mapQ = 60;
            continue;
        }

        if (verbose) cerr << "Set ddf-mapQ for " << a[i] << '\n';
        int sub_ddf_sc = numeric_limits<int>::min();
        int n_ddfsc_sub = 0;
        for (int j = i + 1; j < c; ++j) {
            if (a[j].parent != i) continue;
            if (verbose) cerr << "  **** find hom align " << a[j] << '\n';
            sub_ddf_sc = max(sub_ddf_sc, a[j].ddf_score);
            if (a[j].ddf_score >= a[i].ddf_score) ++n_ddfsc_sub;
        }

        int ddf_mapQ = 0;
        int ddf_score = a[i].ddf_score;
        if (1) {
		    double x = 1.0 * sub_ddf_sc / ddf_score;
		    int mapq = (int)(160.0 * (1.0f - x) * logf(ddf_score));
		    mapq -= (int)(4.343f * logf(n_ddfsc_sub + 1) + .499f);
		    mapq = mapq > 0? mapq : 0;
		    if (ddf_score > sub_ddf_sc && mapq == 0) mapq = 1;

            ddf_mapQ = min(60, mapq);
        }

        if (verbose) fprintf(stderr, "  -------------- ddf-mapQ = %d\n", ddf_mapQ);
        a[i].ddf_mapQ = ddf_mapQ;
    }
}

static void
s_set_mapQ_map_sc(PoreCInitHit* a, int c)
{
    const int verbose = 0;
    pdqsort(a, a + c, [](const PoreCInitHit& x, const PoreCInitHit& y) { return x.map_score > y.map_score; });
    for (int i = 0; i < c; ++i) {
        a[i].id = i;
        a[i].parent = i;
        a[i].is_hom = false;
    }   

    for (int i = 1; i < c; ++i) {
	    PoreCInitHit* ai = a + i;
	    int max_ovlp = 0;
	    int max_j = -1;
	    for (int j = i - 1; j >= 0; --j) {
		    PoreCInitHit* aj = a + j;
		    int qb = max(ai->qoff, aj->qoff);
		    int qe = min(ai->qend, aj->qend);
		    int ovlp = qe - qb;
		    if (ovlp > max_ovlp) {
			    max_ovlp = ovlp;
			    max_j = j;
		    }
	    }
	    if (max_j == -1) continue;
	    if (max_ovlp < (ai->qend - ai->qoff) * kSubOvlpFrac) continue;
	    a[i].is_hom = true;
	    a[max_j].is_hom = true;
	    a[i].parent = a[max_j].parent;
    }

    for (int i = 0; i < c; ++i) {
        if (a[i].id != a[i].parent) {
            a[i].map_mapQ = 0;
            continue;
        }
        if (!a[i].is_hom) {
            a[i].map_mapQ = 60;
            continue;
        }

        if (verbose) cerr << "Set map-mapQ for " << a[i] << '\n';
        int sub_map_sc = numeric_limits<int>::min();
        int n_mapsc_sub = 0;
        for (int j = i + 1; j < c; ++j) {
            if (a[j].parent != i) continue;
            if (verbose) cerr << "  **** find hom align " << a[j] << '\n';
            sub_map_sc = max(sub_map_sc, a[j].map_score);
            if (a[j].map_score >= a[i].map_score) ++n_mapsc_sub;
        }

        int map_mapQ = 0;
        int map_score = a[i].map_score;
        if (1) {
		    double x = 1.0 * sub_map_sc / map_score;
		    double f = 1.0 * a[i].map_score / a[i].ddf_score;
		    double z = 20.0;
		    int mapq = (int)(z * (1.0f - x) * logf(map_score/f));
		    mapq -= (int)(4.343f * logf(n_mapsc_sub + 1) + .499f);
		    mapq = mapq > 0? mapq : 0;
		    if (map_score > sub_map_sc && mapq == 0) mapq = 1;
            map_mapQ = min(60, mapq);
        }

        if (verbose) fprintf(stderr, "  -------------- map-mapQ = %d\n", map_mapQ);
        a[i].map_mapQ = map_mapQ;
    }
}

void set_mapQ_for_init_hits(PoreCInitHit* a, int c)
{
    s_set_mapQ_ddf_sc(a, c);
    s_set_mapQ_map_sc(a, c);
    for (int i = 0; i < c; ++i) a[i].mapQ = max(a[i].ddf_mapQ, a[i].map_mapQ);
}

/////////

void set_mapQ_for_ddf_chains(HbnChainInfo* a, int c)
{
    const int verbose = 0;
    for (int i = 0; i < c; ++i) {
        a[i].id = i;
        a[i].parent = i;
        a[i].is_hom = false;
        a[i].is_valid = true;
        a[i].cnt = 0;
    }   

    for (int i = 1; i < c; ++i) {
	    HbnChainInfo* ai = a + i;
	    int max_ovlp = 0;
	    int max_j = -1;
	    for (int j = i - 1; j >= 0; --j) {
		    HbnChainInfo* aj = a + j;
		    int qb = max(ai->qb, aj->qb);
		    int qe = min(ai->qe, aj->qe);
            if (qe > qb && ai->sid == aj->sid && ai->sdir == aj->sdir) {
                int sb = max(ai->soff, aj->soff);
                int se = min(ai->send, aj->send);
                if (se > sb) {
                    ai->is_valid = false;
                    break;
                }
            }
		    int ovlp = qe - qb;
		    if (ovlp > max_ovlp) {
			    max_ovlp = ovlp;
			    max_j = j;
		    }
	    }
        if (!ai->is_valid) continue;
	    if (max_j == -1) continue;
	    if (max_ovlp < (ai->qe - ai->qb) * kSubOvlpFrac) continue;
	    a[i].is_hom = true;
	    a[max_j].is_hom = true;
	    a[i].parent = a[max_j].parent;
    }

    for (int i = 0; i < c; ++i) {
        if (!a[i].is_valid) continue;
        if (a[i].id != a[i].parent) {
            int p = a[i].parent;
            if (a[p].cnt >= 5) {
                a[i].is_valid = false;
                continue;
            }
            ++a[p].cnt;
            a[i].mapQ = 0;
            continue;
        }
        if (!a[i].is_hom) {
            a[i].mapQ = 60;
            continue;
        }

        if (verbose) cerr << "Set ddf-mapQ for " << a[i] << '\n';
        int sub_ddf_sc = numeric_limits<int>::min();
        int n_ddfsc_sub = 0;
        for (int j = i + 1; j < c; ++j) {
            if (a[j].parent != i) continue;
            if (verbose) cerr << "  **** find hom align " << a[j] << '\n';
            sub_ddf_sc = max(sub_ddf_sc, a[j].sc);
            if (a[j].sc >= a[i].sc) ++n_ddfsc_sub;
        }

        int ddf_mapQ = 0;
        int ddf_score = a[i].sc;
        if (1) {
		    double x = 1.0 * sub_ddf_sc / ddf_score;
		    int mapq = (int)(120.0 * (1.0f - x) * logf(ddf_score));
		    mapq -= (int)(4.343f * logf(n_ddfsc_sub + 1) + .499f);
		    mapq = mapq > 0? mapq : 0;
		    if (ddf_score > sub_ddf_sc && mapq == 0) mapq = 1;

            ddf_mapQ = min(60, mapq);
        }

        if (verbose) fprintf(stderr, "  -------------- ddf-mapQ = %d\n", ddf_mapQ);
        a[i].mapQ = ddf_mapQ;
    }
}
