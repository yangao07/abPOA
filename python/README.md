# pyabPOA: abPOA python library

`pyabPOA`

## Installation

### Install `pyabPOA` with `pip`

`pyabPOA` can be installed with `pip`:

```
pip install pyabPOA            # first time installation
pip install pyabPOA --upgrade  # update to the latest version
```


### Install `pyabPOA` from source
Alternatively, you can install `pyabPOA` from source:
```
git clone https://github.com/yangao07/abPOA.git
cd abPOA/python
make install 
```
Similarly, you can choose to add `AVX2=1`, `SSE4=1` or `SSE2=1` to the `make` command to manually specify the SIMD instruction you want to use:
```
make install AVX2=1
make install SSE4=1
```
Note that `AVX2` is used if nothing is specified.

### Getting started
After installation, you can run the toy example script to test it:
```
python ./msa_example.py
```

## Usage

