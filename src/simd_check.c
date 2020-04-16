#include "simd_instruction.h"
/* source: https://github.com/Mysticial/FeatureDetector */
//  OS Features
int OS_AVX;
int OS_AVX512;

//  Misc.
int HW_MMX;
int HW_x64;
int HW_ABM;
int HW_RDRAND;
int HW_BMI1;
int HW_BMI2;
int HW_ADX;
int HW_PREFETCHWT1;
int HW_MPX;

//  SIMD: 128-bit
int HW_SSE;
int HW_SSE2;
int HW_SSE3;
int HW_SSSE3;
int HW_SSE41;
int HW_SSE42;
int HW_SSE4a;
int HW_AES;
int HW_SHA;

//  SIMD: 256-bit
int HW_AVX;
int HW_XOP;
int HW_FMA3;
int HW_FMA4;
int HW_AVX2;

//  SIMD: 512-bit
int HW_AVX512_F;
int HW_AVX512_PF;
int HW_AVX512_ER;
int HW_AVX512_CD;
int HW_AVX512_VL;
int HW_AVX512_BW;
int HW_AVX512_DQ;
int HW_AVX512_IFMA;
int HW_AVX512_VBMI;

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || defined(_M_IX86)
#if defined(__GNUC__) || defined(__clang__)
#else
#error "No cpuid intrinsic defined for compiler."
#endif
#else
#error "No cpuid intrinsic defined for processor architecture."
#endif

static inline void cpuid(int32_t out[4], int32_t x) {
    __cpuid_count(x, 0, out[0], out[1], out[2], out[3]);
}

static inline uint64_t xgetbv(unsigned int index){
    uint32_t eax, edx;
    __asm__ __volatile__("xgetbv" : "=a"(eax), "=d"(edx) : "c"(index));
    return ((uint64_t)edx << 32) | eax;
}

#define _XCR_XFEATURE_ENABLED_MASK  0

int detect_OS_AVX(){ //  Copied from: http://stackoverflow.com/a/22521619/922184
    int avxSupported = 0;

    int cpuInfo[4];
    cpuid(cpuInfo, 1);

    int osUsesXSAVE_XRSTORE = (cpuInfo[2] & (1 << 27)) != 0;
    int cpuAVXSuport = (cpuInfo[2] & (1 << 28)) != 0;

    if (osUsesXSAVE_XRSTORE && cpuAVXSuport)
    {
        uint64_t xcrFeatureMask = xgetbv(_XCR_XFEATURE_ENABLED_MASK);
        avxSupported = (xcrFeatureMask & 0x6) == 0x6;
    }

    return avxSupported;
}

int detect_OS_AVX512(){
    if (!detect_OS_AVX())
        return 0;

    uint64_t xcrFeatureMask = xgetbv(_XCR_XFEATURE_ENABLED_MASK);
    return (xcrFeatureMask & 0xe6) == 0xe6;
}

void detect_host(){
    //  OS Features
    OS_AVX = detect_OS_AVX();
    OS_AVX512 = detect_OS_AVX512();

    int info[4];
    cpuid(info, 0);
    int nIds = info[0];

    cpuid(info, 0x80000000);
    uint32_t nExIds = info[0];

    //  Detect Features
    if (nIds >= 0x00000001){
        cpuid(info, 0x00000001);
        HW_MMX    = (info[3] & ((int)1 << 23)) != 0;
        HW_SSE    = (info[3] & ((int)1 << 25)) != 0;
        HW_SSE2   = (info[3] & ((int)1 << 26)) != 0;
        HW_SSE3   = (info[2] & ((int)1 <<  0)) != 0;

        HW_SSSE3  = (info[2] & ((int)1 <<  9)) != 0;
        HW_SSE41  = (info[2] & ((int)1 << 19)) != 0;
        HW_SSE42  = (info[2] & ((int)1 << 20)) != 0;
        HW_AES    = (info[2] & ((int)1 << 25)) != 0;

        HW_AVX    = (info[2] & ((int)1 << 28)) != 0;
        HW_FMA3   = (info[2] & ((int)1 << 12)) != 0;

        HW_RDRAND = (info[2] & ((int)1 << 30)) != 0;
    }
    if (nIds >= 0x00000007){
        cpuid(info, 0x00000007);
        HW_AVX2         = (info[1] & ((int)1 <<  5)) != 0;

        HW_BMI1         = (info[1] & ((int)1 <<  3)) != 0;
        HW_BMI2         = (info[1] & ((int)1 <<  8)) != 0;
        HW_ADX          = (info[1] & ((int)1 << 19)) != 0;
        HW_MPX          = (info[1] & ((int)1 << 14)) != 0;
        HW_SHA          = (info[1] & ((int)1 << 29)) != 0;
        HW_PREFETCHWT1  = (info[2] & ((int)1 <<  0)) != 0;

        HW_AVX512_F     = (info[1] & ((int)1 << 16)) != 0;
        HW_AVX512_CD    = (info[1] & ((int)1 << 28)) != 0;
        HW_AVX512_PF    = (info[1] & ((int)1 << 26)) != 0;
        HW_AVX512_ER    = (info[1] & ((int)1 << 27)) != 0;
        HW_AVX512_VL    = (info[1] & ((int)1 << 31)) != 0;
        HW_AVX512_BW    = (info[1] & ((int)1 << 30)) != 0;
        HW_AVX512_DQ    = (info[1] & ((int)1 << 17)) != 0;
        HW_AVX512_IFMA  = (info[1] & ((int)1 << 21)) != 0;
        HW_AVX512_VBMI  = (info[2] & ((int)1 <<  1)) != 0;
    }
    if (nExIds >= 0x80000001){
        cpuid(info, 0x80000001);
        HW_x64   = (info[3] & ((int)1 << 29)) != 0;
        HW_ABM   = (info[2] & ((int)1 <<  5)) != 0;
        HW_SSE4a = (info[2] & ((int)1 <<  6)) != 0;
        HW_FMA4  = (info[2] & ((int)1 << 16)) != 0;
        HW_XOP   = (info[2] & ((int)1 << 11)) != 0;
    }
}


void print_simd_support(void) {
    printf("OS Features:\n");
    printf("    OS AVX      = %d\n", OS_AVX);
    printf("    OS AVX512   = %d\n", OS_AVX512);

    printf("Hardware Features:\n");
    printf("    MMX         = %d\n", HW_MMX);
    printf("    x64         = %d\n", HW_x64);
    printf("    ABM         = %d\n", HW_ABM);
    printf("    RDRAND      = %d\n", HW_RDRAND);
    printf("    BMI1        = %d\n", HW_BMI1);
    printf("    BMI2        = %d\n", HW_BMI2);
    printf("    ADX         = %d\n", HW_ADX);
    printf("    MPX         = %d\n", HW_MPX);
    printf("    PREFETCHWT1 = %d\n", HW_PREFETCHWT1);

    printf("SIMD: 128-bit\n");
    printf("    SSE         = %d\n", HW_SSE);
    printf("    SSE2        = %d\n", HW_SSE2);
    printf("    SSE3        = %d\n", HW_SSE3);
    printf("    SSSE3       = %d\n", HW_SSSE3);
    printf("    SSE4a       = %d\n", HW_SSE4a);
    printf("    SSE4.1      = %d\n", HW_SSE41);
    printf("    SSE4.2      = %d\n", HW_SSE42);
    printf("    AES-NI      = %d\n", HW_AES);
    printf("    SHA         = %d\n", HW_SHA);

    printf("SIMD: 256-bit\n");
    printf("    AVX         = %d\n", HW_AVX);
    printf("    XOP         = %d\n", HW_XOP);
    printf("    FMA3        = %d\n", HW_FMA3);
    printf("    FMA4        = %d\n", HW_FMA4);
    printf("    AVX2        = %d\n", HW_AVX2);

    printf("SIMD: 512-bit\n");
    printf("    AVX512-F    = %d\n", HW_AVX512_F);
    printf("    AVX512-CD   = %d\n", HW_AVX512_CD);
    printf("    AVX512-PF   = %d\n", HW_AVX512_PF);
    printf("    AVX512-ER   = %d\n", HW_AVX512_ER);
    printf("    AVX512-VL   = %d\n", HW_AVX512_VL);
    printf("    AVX512-BW   = %d\n", HW_AVX512_BW);
    printf("    AVX512-DQ   = %d\n", HW_AVX512_DQ);
    printf("    AVX512-IFMA = %d\n", HW_AVX512_IFMA);
    printf("    AVX512-VBMI = %d\n", HW_AVX512_VBMI);

    printf("Summary:\n");
    printf("    Safe to use AVX:     %d\n", HW_AVX && OS_AVX);
    printf("    Safe to use AVX512:  %d\n", HW_AVX512_F && OS_AVX512);
}

int simd_check(void) {
    int simd_flag = 0;
    detect_host();

    if (HW_SSE) simd_flag |= SIMD_SSE;
    if (HW_SSE2) simd_flag |= SIMD_SSE2;
    if (HW_SSE3) simd_flag |= SIMD_SSE3;
    if (HW_SSE41) simd_flag |= SIMD_SSE41;
    if (HW_SSE42) simd_flag |= SIMD_SSE42;
    if (HW_AVX && OS_AVX) simd_flag |= SIMD_AVX;
    if (HW_AVX2 && OS_AVX) simd_flag |= SIMD_AVX2;
    if (HW_AVX512_F && OS_AVX512) simd_flag |= SIMD_AVX512F;
    if (HW_AVX512_BW && OS_AVX512) simd_flag |= SIMD_AVX512BW;

    return simd_flag;
}

#ifdef __CHECK_SIMD_MAIN__
int main(void) {
    char simd_label[6][20] = {"No SIMD", "SSE2 (128 bits)", "SSE4.1 (128 bits)", "AVX2 (256 bits)", "AVX512F (512 bits)", "AVX512BW (512 bits)"};
    int simd_flag = simd_check(), t=0;
    if (simd_flag & SIMD_AVX512BW) printf("__AVX512BW__\n"), t = 5;
    else if (simd_flag & SIMD_AVX512F) printf("__AVX512F__\n"), t = 4;
    else if (simd_flag & SIMD_AVX2) printf("__AVX2__\n"), t = 3;
    else if (simd_flag & SIMD_SSE41) printf("__SSE4_1__\n"), t = 2;
    else if (simd_flag & SIMD_SSE2) printf("__SSE2__\n"), t = 1;
    else printf("NO SIMD\n"), t = 0;

    char msg[100], i;
    fprintf(stderr, "\n");
    sprintf(msg, "==== %s will be used. ====", simd_label[t]);
    for (i = 0; msg[i]; ++i) fprintf(stderr, "="); fprintf(stderr, "\n");
    fprintf(stderr, "%s\n",msg);
    for (i = 0; msg[i]; ++i) fprintf(stderr, "="); fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    return simd_flag;
}
#endif
