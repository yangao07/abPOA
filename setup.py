try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

cmdclass = {}

try:
    from Cython.Build import build_ext
except ImportError: # without Cython
    module_src = 'python/pyabpoa.c'
else: # with Cython
    module_src = 'python/pyabpoa.pyx'
    cmdclass['build_ext'] = build_ext

import os, platform, sys

sys.path.append('python')

if platform.machine() in ["aarch64", "arm64"]:
    simd_flag='-msse2'
else:
    simd_flag='-march=native'
    if os.getenv('SSE4', False):
        simd_flag='-msse4.1'
    elif os.getenv('SSE2', False):
        simd_flag='-msse2'
    elif os.getenv('AVX2', False):
        simd_flag='-mavx2'

src_dir='src/'
inc_dir='src/'

src=[module_src, src_dir+'abpoa_align.c', src_dir+'abpoa_graph.c', src_dir+'simd_abpoa_align.c', src_dir+'utils.c', src_dir+'simd_check.c', src_dir+'abpoa_pog.c']

long_description = open('python/README.md').read()

setup(
    # Information
    name = "pyabpoa",
    description = "pyabpoa: SIMD-based partial order alignment using adaptive band",
    long_description = long_description,
    long_description_content_type="text/markdown",
    version = "1.0.0a1",
    url = "https://github.com/yangao07/abpoa",
    author = "Yan Gao",
    author_email = "yangao07@hit.edu.cn",
    license = "GLP",
    keywords = "multiple-sequence-alignment  partial-order-graph-alignment",
    # Build instructions
    ext_modules = [Extension("pyabpoa",
                    sources=src,
                    include_dirs=[inc_dir],
                    depends=[src_dir+'abpoa.h', src_dir+'abpoa_align.h', src_dir+'abpoa_graph.h', src_dir+'kdq.h', src_dir+'kseq.h', src_dir+'simd_abpoa_align.h', src_dir+'simd_instruction.h', src_dir+'utils.h', 'python/cabpoa.pxd'],
                    libraries = ['z', 'm', 'pthread'],
                    extra_compile_args=['-O3', '-Wno-error=declaration-after-statement', simd_flag])],
    install_requires=['cython'],
    cmdclass = cmdclass
)
