ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := falign2
SOURCES  := \
			align_one_read.cpp \
			chain_align_list.cpp \
			gapped_align.cpp \
			hbn_outputs.cpp \
			index_options.cpp \
			index.cpp \
			infer_enzyme.cpp \
			lookup_table.cpp \
			main.cpp \
			map_main.cpp \
			map_one_query_file.cpp \
			map_options.cpp \
			pore_c_aux.cpp \
			pore_c_traceback.cpp \
			prelim_search.cpp \
			program_info.cpp \
			restrict_enzyme_loci_list.cpp \
			smooth_pca_list.cpp \
			trim_overlap_subseq.cpp \
			../../sw/small_object_alloc.cpp \
			../../sw/edlib_wrapper.cpp \
			../../sw/edlib.cpp \
			../../sw/hbn_traceback_aux.c \
			../../sw/kalloc.c \
			../../sw/ksw2_extd2_sse.c \
			../../sw/ksw2_extz2_sse.c \
			../../sw/ksw2_wrapper.cpp \
			../../sw/small_edlib_align.cpp \

SRC_INCDIRS  := .

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lhbn
TGT_PREREQS := libhbn.a

SUBMAKEFILES :=
