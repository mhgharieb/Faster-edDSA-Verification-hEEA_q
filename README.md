# Accelerating EdDSA Signature Verification with Faster Scalar Size Halving

This repository contains the source code for the article *"Accelerating EdDSA Signature Verification with Faster Scalar Size Halving."*

## Source Code

The source code is located in the [`src/`](src/) directory, with the following three subdirectories:

- [`src/half_size`](src/half_size/): Contains the implementation of the $\textsf{hEEA\\_approx\\_q}$, devision-based $\textsf{hEEA}$, and $\textsf{hSIZE\\_HGCD}$ algorithms for Curve25519 and Curve448, along with modified versions of the $\textsf{reduce\\_basis}$ algorithm for Curve25519 and Curve448; these are adapted from the [Curve9767](https://github.com/pornin/curve9767/blob/master/src/scalar_amd64.c) implementation.

- [`src/ed25519-donna`](src/ed25519-donna/): Contains the [`ed25519-donna`](https://github.com/floodyberry/ed25519-donna) implementation, as well as the following additional files:
    * [`new_batch_helper.h`](src/ed25519-donna/new_batch_helper.h): Implements several helper functions for performing individual and batch verifications, including the multiplicative inverse function, a function to reformat the output of the `curve25519_hEEA_vartime` function to the format used in `ed25519-donna`, $\textsf{DSM\\_B\\_doublePre}$ of the double-scalar multiplication function with doubled-size precomputation table for $B$, and the two versions, $\textsf{QSM\\_B\\_B'}$ and $\textsf{QSM\\_B}$, of the quadruple-scalar multiplication functions.
    * [`ed25519-donna-open_new.h`](src/ed25519-donna/ed25519-donna-open_new.h): Implements the $\textsf{DSM}$-based individual verification method using $\textsf{DSM\\_B\\_doublePre}$, as well as the two versions of the new $\textsf{QSM}$-based individual verification method, one using $\textsf{QSM\\_B\\_B'}$ and the other using $\textsf{QSM\\_B}$. Additionaly, it implements the new individual verification with $\textsf{QSM\\_B\\_B'}$, but using the optimized version of $\textsf{hSIZE\\_HGCD}$ instead of $\textsf{hEEA\\_approx\\_q}$.
    * [`ed25519-donna-batchverify_new.h`](src/ed25519-donna/ed25519-donna-batchverify_new.h): Implements two versions of the new batch verification method, one using $\textsf{hEEA\\_approx\\_q}$ and the other using the optimized version of $\textsf{hSIZE\\_HGCD}$.

- [`src/inverse25519`](src/inverse25519/): Contains the source codes for three different algorithms to compute the inverse modulo the prime $p = 2^{255}-19$:
    * [`inverse25519/EEA_q`](src/inverse25519/EEA_q/): Contains the implementation of the inverse function using our proposed $\textsf{EEA\\_approx\\_q}$.
    * [`inverse25519/bingcd`](src/inverse25519/bingcd/): Contains the source code of [binGCD](https://github.com/pornin/bingcd).
    * [`inverse25519/inverse25519skylake-20210110`](src/inverse25519/inverse25519skylake-20210110/): Contains the source code of the latest version of [safeGCD](https://gcd.cr.yp.to/software/inverse25519skylake-20210110.tar.gz).


## Benchmarks
### Compilation Dependencies

To compile the code, you will need:

1. Compliers: [Clang](https://clang.llvm.org/) compiler, and [GCC](https://gcc.gnu.org/) compiler
3. Libraries: [GMP](https://gmplib.org/) `libgmp-dev`, and `libssl-dev`

### Compilation
It can be done using the provided `Makefile`. It will generate five executable binaries for running tests and benchmarks:

1. `test_halfSize_ed448`: Tests the correctness of the half-size scalars using $\textsf{hEEA\\_approx\\_q}$, $\textsf{reduce\\_basis}$, devision-based $\textsf{hEEA}$, and $\textsf{hSIZE\\_HGCD}$ for Ed448 and benchmarks over 10,000 random instances of $v$.
2. `test_halfSize_ed25519`: Tests the correctness of the half-size scalars using $\textsf{hEEA\\_approx\\_q}$, $\textsf{reduce\\_basis}$, devision-based  $\textsf{hEEA}$, and $\textsf{hSIZE\\_HGCD}$ for Ed25519 and benchmarks over 10,000 random instances of $v$.
3. `test_singleVerification`: Tests the correctness of the proposed individual verification method and benchmarks over 1,000 random Ed25519 signatures.
4. `test_batchVerification`: Tests the correctness of the proposed batch verification method and benchmarks over 100 random Ed25519 batches, with batch sizes varying from 4 to 128 signatures per batch.
5. `test_inverse25519`: Tests the correctness of the inverse modulo the prime $p = 2^{255}-19$ using  $\textsf{EEA\\_approx\\_q}$, $\textsf{binGCD}$, and $\textsf{safeGCD}$, and benchmarks over 10,000 random instances of $v$.

Execution times are given in clock cycles. Each benchmarking program, except `test_inverse25519`, measures the clock cycle using the `rdtsc` opcode before and after 10 repeated calls to the algorithm being tested with the same input value, then verifies the correctness of the computations before moving to the next input value. This process is repeated for each algorithm being tested. In `test_inverse25519`, clock cycles are measured using the `rdtsc` opcode before and after a sequence of 100 dependent inversions, where each inversion’s output serves as the input for the next.

### Test Environment

- CPU: Intel i7-9750H (Coffee Lake) @ 2.60GHz (TurboBoost disabled)
- Clang 10.0.0
- GCC 9.4.0
- GMP 6.2.0

### Benchmark Outputs
#### `test_halfSize_ed448`
```
Benchmark of half-size-scalars for Ed448:
Number of samples = 10000 
Number of rounds = 10 
───────────────────────────────────────────────────────────────────────
Time (ticks)|   hEEA_div   |    hEEA_q    | Speed up     | Improvement
───────────────────────────────────────────────────────────────────────
Best        | 43327        | 8171         | 5.3025       | 81.14 %
Median      | 56042        | 9899         | 5.6614       | 82.34 %
Average     | 56525.29     | 9899.21      | 5.7101       | 82.49 %
───────────────────────────────────────────────────────────────────────
───────────────────────────────────────────────────────────────────────
Time (ticks)| reduce_basis |    hEEA_q    | Speed up     | Improvement
───────────────────────────────────────────────────────────────────────
Best        | 36814        | 8171         | 4.5054       | 77.80 %
Median      | 43809        | 9899         | 4.4256       | 77.40 %
Average     | 43884.44     | 9899.21      | 4.4331       | 77.44 %
───────────────────────────────────────────────────────────────────────
───────────────────────────────────────────────────────────────────────
Time (ticks)|   GMP_hgcd   |    hEEA_q    | Speed up     | Improvement
───────────────────────────────────────────────────────────────────────
Best        | 7985         | 8171         | 0.9772       | -2.33 %
Median      | 13733        | 9899         | 1.3873       | 27.92 %
Average     | 13760.03     | 9899.21      | 1.3900       | 28.06 %
───────────────────────────────────────────────────────────────────────
Done!
```
#### ``test_halfSize_ed25519``
```
Benchmark of half-size-scalars for Ed25519:
Number of samples = 10000 
Number of rounds = 10 
───────────────────────────────────────────────────────────────────────
Time (ticks)|   hEEA_div   |    hEEA_q    | Speed up     | Improvement
───────────────────────────────────────────────────────────────────────
Best        | 21606        | 2624         | 8.2340       | 87.86 %
Median      | 30483        | 3516         | 8.6698       | 88.47 %
Average     | 30514.30     | 3531.47      | 8.6407       | 88.43 %
───────────────────────────────────────────────────────────────────────
───────────────────────────────────────────────────────────────────────
Time (ticks)| reduce_basis |    hEEA_q    | Speed up     | Improvement
───────────────────────────────────────────────────────────────────────
Best        | 10029        | 2624         | 3.8220       | 73.84 %
Median      | 13762        | 3516         | 3.9141       | 74.45 %
Average     | 13823.46     | 3531.47      | 3.9144       | 74.45 %
───────────────────────────────────────────────────────────────────────
───────────────────────────────────────────────────────────────────────
Time (ticks)|   GMP_hgcd   |    hEEA_q    | Speed up     | Improvement
───────────────────────────────────────────────────────────────────────
Best        | 10991        | 2624         | 4.1886       | 76.13 %
Median      | 17944        | 3516         | 5.1035       | 80.41 %
Average     | 17933.06     | 3531.47      | 5.0781       | 80.31 %
───────────────────────────────────────────────────────────────────────
───────────────────────────────────────────────────────────────────────
Time (ticks)| hgcd_enhance1|    hEEA_q    | Speed up     | Improvement
───────────────────────────────────────────────────────────────────────
Best        | 3110         | 2624         | 1.1852       | 15.63 %
Median      | 3894         | 3516         | 1.1075       | 9.71 %
Average     | 4542.01      | 3531.47      | 1.2862       | 22.25 %
───────────────────────────────────────────────────────────────────────
───────────────────────────────────────────────────────────────────────
Time (ticks)| hgcd_enhance2|    hEEA_q    | Speed up     | Improvement
───────────────────────────────────────────────────────────────────────
Best        | 2118         | 2624         | 0.8072       | -23.89 %
Median      | 2803         | 3516         | 0.7972       | -25.44 %
Average     | 2975.37      | 3531.47      | 0.8425       | -18.69 %
───────────────────────────────────────────────────────────────────────
Done!
```
#### ``test_singleVerification``
```
Benchmark of individual verification:
Number of samples = 10000 
Number of rounds = 20 
───────────────────────────────────────────────────────────────────────────────
Time (ticks)|       DSM_B       |  QSM_B (hEEA_q)  | Speed up     | Improvement
───────────────────────────────────────────────────────────────────────────────
Best        | 152474            | 126659           | 1.2038       | 16.93 %
Median      | 157804            | 132373           | 1.1921       | 16.12 %
Average     | 157752.02         | 132327.36        | 1.1921       | 16.12 %
───────────────────────────────────────────────────────────────────────────────
───────────────────────────────────────────────────────────────────────────────
Time (ticks)|  DSM_B_doublePre  | QSM_B_B' (hEEA_q)| Speed up     | Improvement
───────────────────────────────────────────────────────────────────────────────
Best        | 150756            | 120032           | 1.2560       | 20.38 %
Median      | 156636            | 125075           | 1.2523       | 20.15 %
Average     | 156586.85         | 125044.72        | 1.2522       | 20.14 %
───────────────────────────────────────────────────────────────────────────────
───────────────────────────────────────────────────────────────────────────────
Time (ticks)| QSM_B_B' (hEEA_q) |QSM_B_B' (hgcd_enhance2)|Speed up| Improvement
───────────────────────────────────────────────────────────────────────────────
Best        | 120032            | 119094                 | 1.0079 | 0.78 %
Median      | 125075            | 124372                 | 1.0057 | 0.56 %
Average     | 125044.72         | 124630.93              | 1.0033 | 0.33 %
───────────────────────────────────────────────────────────────────────────────
Done!
```
#### ``test_batchVerification``
```
Benchmark of batch verification:
Number of samples = 100 
Number of rounds = 10 
───────────────────────────────────────Average Time (ticks/verification)──────────────────────────────────────────────
──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
Batch size | Old approach | New using hEEA | Speed up | Improvement || New using GMP_hgcd | Speed up     | Improvement
──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
4          | 131145.15    | 143320.68      | 0.9150   | -9.28 %     || 141546.66          | 0.9876       | 1.25 %
5          | 119775.78    | 123770.02      | 0.9677   | -3.33 %     || 122500.88          | 0.9897       | 1.04 %
6          | 112155.42    | 112148.90      | 1.0001   | 0.01  %     || 111599.26          | 0.9951       | 0.49 %
7          | 108531.84    | 104915.53      | 1.0345   | 3.33  %     || 104428.13          | 0.9954       | 0.47 %
8          | 106204.45    | 102196.24      | 1.0392   | 3.77  %     || 101050.68          | 0.9888       | 1.13 %

.               .               .               .         .               .                    .              .         
.               .               .               .         .               .                    .              .         
.               .               .               .         .               .                    .              .         
.               .               .               .         .               .                    .              .         
124        | 69967.54     | 61723.42       | 1.1336   | 11.78 %     || 60462.49           | 0.9796       | 2.09 %
125        | 70005.60     | 61737.67       | 1.1339   | 11.81 %     || 60473.66           | 0.9795       | 2.09 %
126        | 69982.33     | 61680.68       | 1.1346   | 11.86 %     || 60380.75           | 0.9789       | 2.15 %
127        | 70033.14     | 61570.05       | 1.1375   | 12.08 %     || 60422.10           | 0.9814       | 1.90 %
128        | 70412.98     | 61687.25       | 1.1415   | 12.39 %     || 60383.75           | 0.9789       | 2.16 %
──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
```

#### `test_inverse25519`
```
Benchmark of inverse25519:
Number of samples = 10000 
───────────────────────────────────────────────────────────────────────
Time (ticks)|   binGCD     |     EEA_q    | Speed up     | Improvement
───────────────────────────────────────────────────────────────────────
Best        | 6187         | 4482         | 1.3804       | 27.56 %
Median      | 6204         | 5187         | 1.1961       | 16.39 %
Average     | 6298.85      | 5229.96      | 1.2044       | 16.97 %
───────────────────────────────────────────────────────────────────────
───────────────────────────────────────────────────────────────────────
Time (ticks)|   safeGCD    |     EEA_q    | Speed up     | Improvement
───────────────────────────────────────────────────────────────────────
Best        | 3766         | 4482         | 0.8402       | -19.01 %
Median      | 3886         | 5187         | 0.7492       | -33.48 %
Average     | 3939.35      | 5229.96      | 0.7532       | -32.76 %
───────────────────────────────────────────────────────────────────────
Done!
```
## Citation
To cite this work, please use:
