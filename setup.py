import os

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension


simde = ['-DUSE_SIMDE', '-DSIMDE_ENABLE_NATIVE_ALIASES']

machine_system = os.popen("uname").readlines()[0].rsplit()[0]
machine_arch = os.popen("uname -m").readlines()[0].rsplit()[0]

if machine_system == "Darwin":
    # note: see https://github.com/pypa/wheel/issues/406
    simd_flag = ['-march=native', '-D__AVX2__', '-mmacosx-version-min=10.9']
    if machine_arch in ["aarch64", "arm64"]:
        os.environ['_PYTHON_HOST_PLATFORM'] = "macosx-10.9-arm64"
        os.environ['ARCHFLAGS'] = "-arch arm64"
    else: # x86_64
        os.environ['_PYTHON_HOST_PLATFORM'] = "macosx-10.9-x86_64"
        os.environ['ARCHFLAGS'] = "-arch x86_64"
else:
    if machine_arch in ["aarch64", "arm64"]:
        simd_flag = ['-march=armv8-a+simd', '-D__AVX2__']
    elif machine_arch in ["aarch32"]:
        simd_flag = ['-march=armv8-a+simd', '-mfpu=auto -D__AVX2__']
    else: # x86_64
        simd_flag=['-march=native']
        if os.getenv('SSE4', False):
            simd_flag=['-msse4.1']
        elif os.getenv('SSE2', False):
            simd_flag=['-msse2']
        elif os.getenv('AVX2', False):
            simd_flag=['-mavx2']
        #elif os.getenv('AVX512F', False):
        #    simd_flag='-mavx512f'
        #elif os.getenv('AVX512BW', False):
        #    simd_flag='-mavx512bw'

src_dir = 'src/'
inc_dir = 'include/'
sources = [
    'abpoa_align.c', 'abpoa_graph.c', 'abpoa_output.c', 'abpoa_plot.c', 'abpoa_seed.c', 'abpoa_seq.c',
    'kalloc.c', 'kstring.c',
    'simd_abpoa_align.c', 'simd_check.c',
    'utils.c']
depends = [
    'abpoa.h', 'abpoa_align.h', 'abpoa_graph.h', 'abpoa_output.h', 'abpoa_seed.h', 'abpoa_seq.h',
    'kalloc.h', 'khash.h', 'kdq.h', 'kseq.h', 'ksort.h', 'kstring.h', 'kvec.h',
    'simd_abpoa_align.h', 'simd_instruction.h',
    'utils.h']

module_src = 'python/pyabpoa.pyx'
module_dep = 'python/cabpoa.pxd'

long_description = open('python/README.md').read()

setup(
    # Information
    name = "pyabpoa",
    description = "pyabpoa: SIMD-based partial order alignment using adaptive band",
    long_description = long_description,
    long_description_content_type="text/markdown",
    version = "1.5.3",
    url = "https://github.com/yangao07/abPOA",
    author = "Yan Gao",
    author_email = "yangao@ds.dfci.harvard.edu",
    license = "MIT",
    keywords = "multiple-sequence-alignment  partial-order-graph-alignment",
    setup_requires=["cython<3"], # see https://github.com/cython/cython/issues/5568
    # Build instructions
    ext_modules = [
        Extension(
            "pyabpoa",
            sources=[module_src] + [src_dir + x for x in sources],
            include_dirs=[inc_dir],
            depends=[module_dep] + [src_dir + x for x in depends],
            libraries=['z', 'm', 'pthread'],
            # extra_compile_args=['-O3', '-Wno-error=declaration-after-statement', '-D __DEBUG__'] + simde + simd_flag
            extra_compile_args=['-O3', '-Wno-misleading-indentation', '-Wno-error=declaration-after-statement'] + simde + simd_flag
    )]
)
