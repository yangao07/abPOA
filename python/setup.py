try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

cmdclass = {}
long_description = open('README.md').read()

try:
    from Cython.Build import build_ext
except ImportError: # without Cython
    module_src = 'pyabpoa.c'
else: # with Cython
    module_src = 'pyabpoa.pyx'
    cmdclass['build_ext'] = build_ext

import os

simd_flag='-march=native'
if os.getenv('SSE4', False):
    simd_flag='-msse4.1'
elif os.getenv('SSE2', False):
    simd_flag='-msse2'
elif os.getenv('AVX2', False):
    simd_flag='-mavx2'

src_dir='../src/'
inc_dir='../include'

setup(
    # Information
    name = "pyabPOA",
    description = "pyabPOA: SIMD-based Python library for fast partial order alignment using adaptive band",
    long_description = long_description,
    version = "1.0.0",
    url = "https://github.com/yangao07/abpoa",
    author = "Yan Gao",
    author_email = "yangao07@hit.edu.cn",
    license = "GLP",
    keywords = "multiple-sequence-alignment  sequece-to-graph alignment",
    # Build instructions
    ext_modules = [Extension("pyabPOA",
                    sources=[module_src, src_dir+'abpoa_align.c', src_dir+'abpoa_graph.c', src_dir+'simd_abpoa_align.c', src_dir+'utils.c', src_dir+'simd_check.c', src_dir+'abpoa_pog.c'],
                    include_dirs=[inc_dir],
                    depends=[src_dir+'abpoa.h', src_dir+'abpoa_align.h', src_dir+'abpoa_graph.h', src_dir+'kdq.h', src_dir+'kseq.h', src_dir+'simd_abpoa_align.h', src_dir+'simd_instruction.h', src_dir+'utils.h', 'cabpoa.pxd'],
                    libraries = ['z', 'm', 'pthread'],
                    extra_compile_args=['-O3', '-Wno-error=declaration-after-statement', simd_flag])],
    install_requires=['cython'],
    cmdclass = cmdclass
)
