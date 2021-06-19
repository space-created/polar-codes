## Build
Build this program with "mbuild" tool (https://www.mathworks.com/help/compiler_sdk/ml_code/mbuild.html)
```bash
mbuild -O -output a.exe optimizer.cpp
```
## Run
```bash
a.exe
```
File format:  
`snr_number` - number of considered signal to noise ratios  
`iteration_number` - number of iterations   
`m` `K` - power of codeword size N, number of information bits  
`v_matrix_size` - number of rows in the constraints matrix  
`constraints matrix`: boolean matrix (v_matrix_size, N = 2^m)  
`SNRb` - signal to noise ration value  
`frozen_bits` - frozen bit boolean vector, 1 means index is frozen otherwise unfrozen