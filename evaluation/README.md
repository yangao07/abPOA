Evaluation of abPOA and SPOA on simulated datasets
==
Two wrap-up programs, [msa_abPOA.c](msa_abPOA.c) and [msa_spoa.cpp](msa_spoa.cpp) were used to make the two libraries allow multiple sets of sequences as the input.

Firstly, make the two libraries:
```
git clone https://github.com/yangao07/abPOA.git
cd abPOA && make # make abPOA library

cd ../
wget https://github.com/rvaser/spoa/releases/download/3.4.0/spoa-v3.4.0.tar.gz
tar -zxvf spoa-v3.4.0.tar.gz && cd spoa-v3.4.0
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make # make spoa library
```
Then, build two wrap-up programs:
```
cd ../abPOA/evaluation
gcc -O3 msa_abPOA.c -I ../include -L ../lib -labpoa -lz -o msa_abPOA
g++ -O3 msa_spoa.cpp -I ../../spoa/include/ -L ../../spoa/build/lib/ -lspoa -o msa_spoa
```

Run the two wrap-up programs:
```
msa_spoa seq_sets.fa -n depth > cons.fa        # SPOA
msa_abPOA seq_sets.fa -n depth -b -1 > cons.fa # abPOA without adaptive banding
msa_abPOA seq_sets.fa -n depth -b 10 > cons.fa # abPOA with adaptive banding
```
