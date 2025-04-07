ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := u4falign2
SOURCES  := \
			bam-to-frag-bam.cpp \
			ddf_score.cpp \
			dump_chr_name.cpp \
			eval_sim_frag_bam.cpp \
			extract_chr.cpp \
			extract-fastq.cpp \
			fix_ngmlr_sam.cpp \
			main.cpp \
			sim_fastq_stats.cpp \
			stat_real_frag_bam.cpp \
			simulate_main.cpp \
			simulate_pore_c_frag.cpp \
			simulate_pore_c_read.cpp \
			../map/restrict_enzyme_loci_list.cpp

SRC_INCDIRS  := .

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lhbn
TGT_PREREQS := libhbn.a

SUBMAKEFILES :=