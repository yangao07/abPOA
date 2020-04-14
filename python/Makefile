default: install

simd_flag = 
ifneq ($(AVX2),)
	simd_flag = AVX2=1
endif
ifneq ($(SSE4),)
	simd_flag = SSE4=1
endif
ifneq ($(SSE2),)
	simd_flag = SSE2=1
endif

install: cabpoa.pxd pyabpoa.pyx setup.py README.md
	${simd_flag} python setup.py install

sdist: pyabpoa.c setup.py README.md
	${simd_flag} python setup.py sdist bdist_wheel

publish: clean sdist
	twine upload dist/*

clean:
	rm -rf build dist abpoa.egg-info build pyabPOA.egg-info
	rm -f pyabpoa.c
