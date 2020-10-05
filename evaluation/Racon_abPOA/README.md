Racon-abPOA
==
To replace SPOA with abPOA in Racon:
```
# download Racon
wget https://github.com/lbcb-sci/racon/releases/download/1.4.13/racon-v1.4.13.tar.gz
tar -zxvf racon-v1.4.13.tar.gz && cd racon-v1.4.13

# download abPOA
git clone https://github.com/yangao07/abPOA.git

# replace SPOA with abPOA:
cp abPOA -r vendor/
cp abPOA/evaluation/Racon_abPOA/Racon_abPOA_CMakeLists.txt CMakeLists.txt
cp abPOA/evaluation/Racon_abPOA/polisher.cpp src/polisher.cpp
cp abPOA/evaluation/Racon_abPOA/polisher.hpp src/polisher.hpp
cp abPOA/evaluation/Racon_abPOA/window.cpp src/window.cpp
cp abPOA/evaluation/Racon_abPOA/window.hpp src/window.hpp

# build Racon-SPOA
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```