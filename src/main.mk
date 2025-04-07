ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET       := libhbn.a

SOURCES      := \
	./corelib/db_format.c \
	./corelib/fasta.cpp \
	./corelib/frag_id.cpp \
	./corelib/getMemorySize.c \
	./corelib/getRSS.c \
	./corelib/hbn_aux.c \
	./corelib/hbn_package_version.c \
	./corelib/read_list.cpp \
	./corelib/split_string_by_char.cpp \
	./corelib/unpacked_seqdb.cpp \
	./ncbi_blast/c_ncbi_blast_aux.c \
	./ncbi_blast/ncbi_blast_aux.cpp \
	./ncbi_blast/str_util/ncbistr_util.cpp \
	./ncbi_blast/str_util/ncbistr.cpp \
	./ncbi_blast/str_util/str_cmp.cpp \
	./ncbi_blast/str_util/numeric_str_interconv.cpp \
	./ncbi_blast/str_util/str_util.cpp

SRC_INCDIRS  := 

SUBMAKEFILES := \
	./app/map/main.mk \
	./app/u4falign2/main.mk