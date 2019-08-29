CC      =	gcc
CFLAGS  =	-Wall -O3 -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function -Wno-misleading-indentation

# for debug
ifneq ($(gdb),)
	#CFLAGS   =	 -Wall -O3 -D __DEBUG__ -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function -Wno-misleading-indentation
	CFLAGS   =	 -g -Wall -D __DEBUG__ -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function -Wno-misleading-indentation
endif

ifneq ($(old),)
	CFLAGS  += -D __OLD__
endif
# for gprof
ifneq ($(pg),)
	PG_FLAG  =   -pg
	CFLAGS  +=   -pg
endif

LIB     =	-lm -lz -lpthread
BIN_DIR =	./bin
LIB_DIR =   ./lib
INC_DIR =   ./include
SRC_DIR =   ./src

SOURCE  =	$(SRC_DIR)/abpoa.c $(SRC_DIR)/abpoa_align.c $(SRC_DIR)/abpoa_graph.c $(SRC_DIR)/simd_abpoa_align.c $(SRC_DIR)/simd_check.c $(SRC_DIR)/utils.c $(SRC_DIR)/abpoa_dot_plot.c $(SRC_DIR)/agglo_hier_clu.c
HEADER  =	$(SRC_DIR)/abpoa.h $(SRC_DIR)/abpoa_align.h $(SRC_DIR)/abpoa_graph.h $(SRC_DIR)/align.h $(SRC_DIR)/kdq.h $(SRC_DIR)/kseq.h $(SRC_DIR)/ksort.h $(SRC_DIR)/simd_instruction.h $(SRC_DIR)/simd_abpoa_align.h $(SRC_DIR)/utils.h
OBJS    =	$(SRC_DIR)/abpoa_align.o $(SRC_DIR)/abpoa_graph.o $(SRC_DIR)/simd_abpoa_align.o $(SRC_DIR)/simd_check.o $(SRC_DIR)/utils.o $(SRC_DIR)/abpoa_dot_plot.o $(SRC_DIR)/agglo_hier_clu.o

# SIMD label
SIMD_CHECK_D	= -D __CHECK_SIMD_MAIN__

SSE41 			= __SSE4_1__
AVX2 			= __AVX2__
AVX512F 		= __AVX512F__
AVX512BW 		= __AVX512BW__

FLAG_SSE2  		= -msse2
FLAG_SSE41      = -msse4.1
FLAG_AVX2       = -mavx2 #TODO -march=native
FLAG_AVX512F    = -mavx512f
FLAG_AVX512BW   = -mavx512bw
SIMD_FLAG       = -msse2

simd_flag := ${shell ./bin/simd_check 2> /dev/null}

ifeq ($(simd_flag), $(AVX512BW))
	SIMD_FLAG = $(FLAG_AVX512BW)
else ifeq ($(simd_flag), $(AVX512F))
	SIMD_FLAG = $(FLAG_AVX512F)
else ifeq ($(simd_flag), $(AVX2))
	SIMD_FLAG = $(FLAG_AVX2)
else ifeq ($(simd_flag), $(SSE41))
	SIMD_FLAG = $(FLAG_SSE41)
endif

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

SIMD_CHECK  	= $(BIN_DIR)/simd_check
BIN     		= $(BIN_DIR)/abPOA
ifneq ($(gdb),)
	BIN = $(BIN_DIR)/gdb_abPOA
endif
ABPOALIB        = $(LIB_DIR)/libabpoa.a
# TODO add example
EXAMPLE         = example


all:		    $(BIN) 
abPOA:     		$(BIN)
gdb_abPOA:		$(BIN)
libabpoa:       $(ABPOALIB)
example:        $(EXAMPLE)

simd_check:$(SIMD_CHECK)
	$(shell ./bin/simd_check > /dev/null)

$(SIMD_CHECK):$(SRC_DIR)/simd_check.c $(SRC_DIR)/simd_instruction.h
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	$(CC) $(SIMD_CHECK_D) $< -o $@

$(BIN):$(SRC_DIR)/abpoa.o $(ABPOALIB)
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	$(CC) $(CFLAGS) $< -L$(LIB_DIR) -labpoa $(LIB) -o $@ $(PG_FLAG)

$(EXAMPLE):example.c $(ABPOALIB)
	$(CC) $(CFLAGS) $< -o $@ -I $(INC_DIR) -L $(LIB_DIR) -labpoa $(LIB)

$(ABPOALIB):$(OBJS)
	if [ ! -d $(LIB_DIR) ]; then mkdir $(LIB_DIR); fi
	$(AR) -csr $@ $(OBJS)

$(SRC_DIR)/abpoa.o:$(SRC_DIR)/abpoa.c $(SRC_DIR)/abpoa.h $(SRC_DIR)/abpoa_graph.h $(SRC_DIR)/agglo_hier_clu.h $(SRC_DIR)/abpoa_align.h \
				   $(SRC_DIR)/align.h $(SRC_DIR)/utils.h $(SRC_DIR)/simd_instruction.h
	$(CC) -c $(CFLAGS) $(SIMD_FLAG) $< -o $@

$(SRC_DIR)/simd_check.o:$(SRC_DIR)/simd_check.c $(SRC_DIR)/simd_instruction.h
	$(CC) -c $(CFLAGS) $(SIMD_FLAG) $< -o $@

$(SRC_DIR)/simd_abpoa_align.o:$(SRC_DIR)/simd_abpoa_align.c $(SRC_DIR)/abpoa_graph.h $(SRC_DIR)/abpoa_align.h $(SRC_DIR)/simd_instruction.h $(SRC_DIR)/utils.h
	$(CC) -c $(CFLAGS) $(SIMD_FLAG) $< -o $@

clean:
	rm -f $(SRC_DIR)/*.[oa] $(LIB_DIR)/*.[oa] $(BIN) #$(SIMD_CHECK)
