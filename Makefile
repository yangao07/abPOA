.PHONY: default
default: abpoa

#CC          = gcc
OS          := $(shell uname)
ARCH        := $(shell arch)
# add -fno-tree-vectorize to avoid certain vectorization errors in O3 optimization
# right now, we are using -O3 for the best performance, and no vectorization errors were found
EXTRA_FLAGS = -Wall -Wno-unused-function -Wno-misleading-indentation -DUSE_SIMDE -DSIMDE_ENABLE_NATIVE_ALIASES# -fno-tree-vectorize

# for debug
ifneq ($(debug),)
	EXTRA_FLAGS  += -D __DEBUG__
endif
ifneq ($(sdebug),)
	EXTRA_FLAGS  += -D __SIMD_DEBUG__
endif

# for gdb
ifneq ($(gdb),)
	OPT_FLAGS = -g
else
	OPT_FLAGS = -O3
endif

CFLAGS        = $(OPT_FLAGS) $(EXTRA_FLAGS)

# for gprof
ifneq ($(pg),)
	PG_FLAG  =   -pg
	CFLAGS  +=   -pg
endif

LIB     = -lm -lz -lpthread
ifneq ($(PREFIX),)
	OUT_PRE_DIR = $(PREFIX)
else
	OUT_PRE_DIR = .
endif

BIN_DIR = $(OUT_PRE_DIR)/bin
LIB_DIR = $(OUT_PRE_DIR)/lib
INC_DIR = ./include
SRC_DIR = ./src

# This is everything that gets bundled into our main library
OBJS = $(addprefix $(SRC_DIR)/, abpoa_align.o abpoa_graph.o abpoa_plot.o abpoa_seed.o abpoa_seq.o abpoa_output.o kalloc.o kstring.o simd_abpoa_align.o simd_check.o utils.o)

# Set default SIMD flags
SIMD_FLAG   = -march=native

# auto-detect some appropriate defaults -- this helps users in the common case of macOS with arm
ifeq ($(ARCH), $(filter $(ARCH), aarch64 arm64))
ifeq ($(OS), Darwin)
	# issues/44 suggests this doesn't work for some users
	SIMD_FLAG = -march=armv8-a+simd -D__AVX2__
else
	SIMD_FLAG = -march=armv8-a+simd -D__AVX2__
endif
endif

# override if user specified
ifneq ($(armv7),) # for ARMv7
	SIMD_FLAG   =  -march=armv7-a -mfpu=neon -D__AVX2__
else
ifneq ($(armv8),) # for ARMv8
ifneq ($(aarch64),) # for Aarch64
	SIMD_FLAG   =  -march=armv8-a+simd -D__AVX2__
else # for Aarch32
	SIMD_FLAG   =  -march=armv8-a+simd -mfpu=auto -D__AVX2__
endif
endif
endif

# some more possible overrides
FLAG_SSE2     = -msse2
FLAG_SSE41    = -msse4.1
FLAG_AVX2     = -mavx2
# FLAG_AVX512F  = -mavx512f
# FLAG_AVX512BW = -mavx512bw

ifneq ($(sse2),)
	SIMD_FLAG=$(FLAG_SSE2)
	py_SIMD_FLAG = SSE2=1
else ifneq ($(sse41),)
	SIMD_FLAG=$(FLAG_SSE41)
	py_SIMD_FLAG = SSE41=1
else ifneq ($(avx2),)
	SIMD_FLAG=$(FLAG_AVX2)
	py_SIMD_FLAG = AVX2=1
#else ifneq ($(avx512f),)
#	SIMD_FLAG=$(FLAG_AVX512F)
#	py_SIMD_FLAG = AVX512f=1
#else ifneq ($(avx512bw),)
#	SIMD_FLAG=$(FLAG_AVX512BW)
#	py_SIMD_FLAG = AVX512BW=1
endif

# SIMD check -- unused
SIMD_CHECK_D = -D __CHECK_SIMD_MAIN__



$(SRC_DIR)/%.o: $(SRC_DIR)/%.c
		$(CC) -c $(CFLAGS) $^ -I$(INC_DIR) -o $@

BIN      = $(BIN_DIR)/abpoa
ifneq ($(gdb),)
	BIN  = $(BIN_DIR)/gdb_abpoa
endif
ABPOALIB = $(LIB_DIR)/libabpoa.a
# TODO add example
EXAMPLE  = example


all:       $(BIN) 
abpoa:     $(BIN)
libabpoa:  $(ABPOALIB)
example:   $(EXAMPLE)

$(BIN):$(SRC_DIR)/abpoa.o $(ABPOALIB)
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	$(CC) $(CFLAGS) $< -I$(INC_DIR) -L$(LIB_DIR) -labpoa $(LIB) -o $@ $(PG_FLAG)

$(EXAMPLE):example.c $(ABPOALIB)
	$(CC) $(CFLAGS) $< -o $@ -I$(INC_DIR) -L$(LIB_DIR) -labpoa $(LIB)

$(ABPOALIB):$(OBJS)
	if [ ! -d $(LIB_DIR) ]; then mkdir $(LIB_DIR); fi
	$(AR) -csr $@ $(OBJS)

$(SRC_DIR)/abpoa.o:$(SRC_DIR)/abpoa.c $(SRC_DIR)/abpoa.h $(SRC_DIR)/abpoa_graph.h $(SRC_DIR)/abpoa_align.h \
                   $(SRC_DIR)/abpoa_seq.h $(SRC_DIR)/utils.h $(SRC_DIR)/simd_instruction.h
	$(CC) -c $(CFLAGS) $(SIMD_FLAG) -I$(INC_DIR) $< -o $@

$(SRC_DIR)/simd_check.o:$(SRC_DIR)/simd_check.c $(SRC_DIR)/simd_instruction.h
	$(CC) -c $(CFLAGS) $(SIMD_FLAG) -I$(INC_DIR) $< -o $@

$(SRC_DIR)/simd_abpoa_align.o:$(SRC_DIR)/simd_abpoa_align.c $(SRC_DIR)/abpoa_graph.h $(SRC_DIR)/abpoa_align.h $(SRC_DIR)/simd_instruction.h $(SRC_DIR)/utils.h
	$(CC) -c $(CFLAGS) $(SIMD_FLAG) -I$(INC_DIR) $< -o $@

install_py: setup.py python/cabpoa.pxd python/pyabpoa.pyx python/README.md
	${py_SIMD_FLAG} python setup.py install
	
sdist: setup.py python/cabpoa.pxd python/pyabpoa.pyx python/README.md
	${py_SIMD_FLAG} python setup.py sdist #bdist_wheel

publish_pypi: clean_py sdist
	twine upload dist/*

clean:
	rm -f $(SRC_DIR)/*.[oa] $(LIB_DIR)/*.[oa] $(BIN)
clean_py:
	rm -rf build/ dist/ pyabpoa.egg-info/ python/pyabpoa.c pyabpoa.cpython-*.so
