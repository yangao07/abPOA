#CC          = gcc
EXTRA_FLAGS = -Wno-unused-function -Wno-misleading-indentation
CFLAGS      = -Wall -O3 $(EXTRA_FLAGS)

# for debug
ifneq ($(debug),)
	DFLAGS   =   -D __DEBUG__
endif
# for gdb
ifneq ($(gdb),)
	CFLAGS   = -Wall -g ${DFLAGS} $(EXTRA_FLAGS)
else
	CFLAGS   = -Wall -O3 ${DFLAGS} $(EXTRA_FLAGS)
endif

# for gprof
ifneq ($(pg),)
	PG_FLAG  =   -pg
	CFLAGS  +=   -pg
endif

LIB     = -lm -lz -lpthread
BIN_DIR = ./bin
LIB_DIR = ./lib
INC_DIR = ./include
SRC_DIR = ./src

SOURCE = $(SRC_DIR)/abpoa.c $(SRC_DIR)/abpoa_align.c $(SRC_DIR)/abpoa_graph.c $(SRC_DIR)/simd_abpoa_align.c $(SRC_DIR)/simd_check.c $(SRC_DIR)/utils.c $(SRC_DIR)/abpoa_pog.c
HEADER = $(SRC_DIR)/abpoa.h $(SRC_DIR)/abpoa_align.h $(SRC_DIR)/abpoa_graph.h $(SRC_DIR)/seq.h $(SRC_DIR)/kdq.h $(SRC_DIR)/kseq.h $(SRC_DIR)/simd_instruction.h $(SRC_DIR)/simd_abpoa_align.h $(SRC_DIR)/utils.h
OBJS   = $(SRC_DIR)/abpoa_align.o $(SRC_DIR)/abpoa_graph.o $(SRC_DIR)/simd_abpoa_align.o $(SRC_DIR)/simd_check.o $(SRC_DIR)/utils.o $(SRC_DIR)/abpoa_pog.o

# SIMD label
SIMD_CHECK_D = -D __CHECK_SIMD_MAIN__

FLAG_SSE2     = -msse2
FLAG_SSE41    = -msse4.1
FLAG_AVX2     = -mavx2
FLAG_AVX512F  = -mavx512f
FLAG_AVX512BW = -mavx512bw
SIMD_FLAG     = -march=native

ifneq ($(sse2),)
	SIMD_FLAG=$(FLAG_SSE2)
else ifneq ($(sse41),)
	SIMD_FLAG=$(FLAG_SSE41)
else ifneq ($(avx2),)
	SIMD_FLAG=$(FLAG_AVX2)
else ifneq ($(avx512f),)
	SIMD_FLAG=$(FLAG_AVX512F)
else ifneq ($(avx512bw),)
	SIMD_FLAG=$(FLAG_AVX512BW)
endif

.c.o:
		$(CC) -c $(CFLAGS) $< -o $@

BIN      = $(BIN_DIR)/abPOA
ifneq ($(gdb),)
	BIN  = $(BIN_DIR)/gdb_abPOA
endif
ABPOALIB = $(LIB_DIR)/libabpoa.a
# TODO add example
EXAMPLE  = example


all:       $(BIN) 
abPOA:     $(BIN)
libabpoa:  $(ABPOALIB)
example:   $(EXAMPLE)

$(BIN):$(SRC_DIR)/abpoa.o $(ABPOALIB)
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	$(CC) $(CFLAGS) $< -L$(LIB_DIR) -labpoa $(LIB) -o $@ $(PG_FLAG)

$(EXAMPLE):example.c $(ABPOALIB)
	$(CC) $(CFLAGS) $< -o $@ -I $(INC_DIR) -L $(LIB_DIR) -labpoa $(LIB)

$(ABPOALIB):$(OBJS)
	if [ ! -d $(LIB_DIR) ]; then mkdir $(LIB_DIR); fi
	$(AR) -csr $@ $(OBJS)

$(SRC_DIR)/abpoa.o:$(SRC_DIR)/abpoa.c $(SRC_DIR)/abpoa.h $(SRC_DIR)/abpoa_graph.h $(SRC_DIR)/abpoa_align.h \
                   $(SRC_DIR)/seq.h $(SRC_DIR)/utils.h $(SRC_DIR)/simd_instruction.h
	$(CC) -c $(CFLAGS) $(SIMD_FLAG) $< -o $@

$(SRC_DIR)/simd_check.o:$(SRC_DIR)/simd_check.c $(SRC_DIR)/simd_instruction.h
	$(CC) -c $(CFLAGS) $(SIMD_FLAG) $< -o $@

$(SRC_DIR)/simd_abpoa_align.o:$(SRC_DIR)/simd_abpoa_align.c $(SRC_DIR)/abpoa_graph.h $(SRC_DIR)/abpoa_align.h $(SRC_DIR)/simd_instruction.h $(SRC_DIR)/utils.h
	$(CC) -c $(CFLAGS) $(SIMD_FLAG) $< -o $@

clean:
	rm -f $(SRC_DIR)/*.[oa] $(LIB_DIR)/*.[oa] $(BIN)
