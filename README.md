## Polar subcodes
### Simulation
The "simulation" folder contains implementation of the polar coding and decoding simulation.
Two types of channels can be used to transmit random messages:
1. The channel with additive white gaussian noise (AWGN)
2. Binary erasure channel (BEC)

Decoding is carried out according to the [Tal-Vardi method](https://ieeexplore.ieee.org/document/7055304).
Both the construction of polar codes and the construction of polar subcodes can be used.
EBCH codes are used as the parent code for polar subcodes. More details about polar subcodes: https://arxiv.org/abs/1511.01646.
There is also implementation for the randomized polar subcodes. The algorithm generates random dynamic constraints for the most unreliable frozen channels using 
information indexes.

### Weight spectrum optimizer

The "optimizer" folder contains implementation of the Weight spectrum optimizer.
The sets of dynamic frozen constraints of BCH code are randomly iterated over.
The resulting polar code configuration is used to calculate the weight spectrum
(see https://arxiv.org/pdf/2102.07362.pdf works for short codes:  N <= 128).
Further the Poltirev's upper bound is calculated and it repeats for the next random constraints set.
After a finite number of iterations, a weight spectrum (and the constraints set) corresponding to the smallest value of the Poltyrev boundary will be selected.
It should give the polar subcode configuration with the least probability of decoding error at maximum likelihood and therefore the least probability of the Tal-Vardy list decoding algorithm.
The results and scripts for the plots can be found in the "results" folder.