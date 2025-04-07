#include "map_one_query_file.hpp"

#include "align_one_read.hpp"

using namespace std;

class PoreCMapThreadData
{
public:
    PoreCMapThreadData(HbnMapOptions* options,
        HbnIndex* index,
        QueryReader* queries,
        HbnOutputs* out)
    {
        M_options = options;
        M_index = index;
        M_queries = queries;
        M_out = out;
    }

    int get_next_query(PoreCQuery* query)
    {
        return M_queries->get_next_query(query);
    }

    HbnMapOptions* opts()
    {
        return M_options;
    }

    HbnOutputs* out() {
        return M_out;
    }

    HbnIndex* index() {
        return M_index;
    }

private:
    HbnMapOptions*      M_options;
    HbnIndex*               M_index;
    QueryReader*            M_queries;
    HbnOutputs*             M_out;
};

static void* s_pore_c_map_thread(void* params)
{
    PoreCMapThreadData* data = (PoreCMapThreadData*)(params);
    HbnMapOptions* options = data->opts();
    HbnIndex* index = data->index();
    HbnUnpackedDatabase* updb = index->updb();
    RestrictEnzymeLociList* enzyme_list = index->enzyme_pos_list();
    HbnOutputs* out = data->out();

    PoreCQuery* query = new PoreCQuery(&enzyme_list->enzyme);
    PrelimSearch* prelim_search = new PrelimSearch(options);
    PoreCGapAlign* gapped_search = new PoreCGapAlign;
    PoreCAlignChainData* pca_chain_data = new PoreCAlignChainData;
    TrimPcaList* trim_pca_list = new TrimPcaList;
    int cnt = 0;

    while (data->get_next_query(query)) {
        query->init();
        align_one_read(index, query, options, prelim_search, gapped_search, pca_chain_data, trim_pca_list);
        out->dump(enzyme_list, updb, query->name, query->id, query->fwd_rqs, query->rev_rqs, 
            query->fwd_qv, query->rev_qv, query->size, nullptr, 0, *trim_pca_list);
	++cnt;
	if (cnt == 10) {
		delete query;
		delete prelim_search;
		delete gapped_search;
		delete pca_chain_data;
		delete trim_pca_list;
		query = new PoreCQuery(&enzyme_list->enzyme);
		prelim_search = new PrelimSearch(options);
		gapped_search = new PoreCGapAlign;
		pca_chain_data = new PoreCAlignChainData;
		trim_pca_list = new TrimPcaList;
		cnt = 0;
	}
    }

    delete query;
    delete prelim_search;
    delete gapped_search;
    delete pca_chain_data;
    delete trim_pca_list;

    return nullptr;
}

void map_one_pore_c_file(HbnMapOptions* options,
    HbnIndex* index,
    QueryReader* queries,
    HbnOutputs* out)
{
    PoreCMapThreadData data(options, index, queries, out);
    const int num_threads = options->num_threads;
    pthread_t jobs[num_threads];
    for (int i = 0; i < num_threads; ++i) {
        pthread_create(jobs + i, nullptr, s_pore_c_map_thread, &data);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(jobs[i], nullptr);
    }
    queries->dump_mapped_stats();
    //dump_extended_hits();
    //dump_rescue_mapQ();
}
