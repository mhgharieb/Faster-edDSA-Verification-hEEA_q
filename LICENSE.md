MIT License

Copyright (c) 2024 Muhammad ElSheikh (mhgaelsh@uwaterloo.ca). All rights reserved.

This is the source code for the paper "Accelerating EdDSA Signature 
Verification with Faster Scalar Size Halving" by  
    Muhammad ElSheikh,  
    Irem Keskinkurt Paksoy,  
    Murat Cenk,  
    M. Anwar Hasan 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

This license applies to the files:   
           `src/half_size/*`,   
           `src/ed25519-donna/new_batch_helper.h`,  
           `src/ed25519-donna/ed25519-donna-open_new.h`,  
           `src/ed25519-donna/ed25519-donna-batchverify_new.h`,  
           `src/inverse25519/EEA_q/*`,  
           `test/*`
           

## Third‑Party Code

This repository incorporates third‑party code from the following sources. Please refer to the individual license files where applicable:

- **ed25519‑donna**  
  Source: [https://github.com/floodyberry/ed25519-donna](https://github.com/floodyberry/ed25519-donna)  
  License: The code is released into the public domain (or via a CC0‑like dedication).

- **The original code used in `src/half_size/curve448/curve448_reduce_basis_vartime.c` and `src/half_size/curve25519/curve25519_reduce_basis_vartime.c`**  
  Source: [https://github.com/pornin/curve9767/blob/master/src/scalar_amd64.c](https://github.com/pornin/curve9767/blob/master/src/scalar_amd64.c)  
  License: The code is released under the MIT License.
  
  MIT License

  Copyright (c) 2020 Thomas Pornin

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

- **bingcd (`src/inverse25519/bingcd`)**  
  Source: [https://github.com/pornin/bingcd](https://github.com/pornin/bingcd)  
  License: The bingcd component is released under the MIT License. For the full text of its MIT License, please see the file located at `src/inverse25519/bingcd/LICENSE`.
  
  MIT License

  Copyright (c) 2020 Thomas Pornin

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
  
- **Constant‑Time GCD and Modular Inversion Software (`src/inverse25519/inverse25519skylake-20210110`)**  
  Source: [https://gcd.cr.yp.to/software.html](https://gcd.cr.yp.to/software.html)  
  License: The code is released into the public domain (or via a CC0‑like dedication).

## Attribution

- **ed25519‑donna** – by Andrew Moon and collaborators ([ed25519-donna](https://github.com/floodyberry/ed25519-donna))
- **Original code for `*_reduce_basis_vartime` function** – by Thomas Pornin ([pornin/curve9767/blob/master/src/scalar_amd64.c](https://github.com/pornin/curve9767/blob/master/src/scalar_amd64.c))
- **bingcd** – by Thomas Pornin ([pornin/bingcd](https://github.com/pornin/bingcd))
- **Constant‑Time GCD and Modular Inversion Software** – by Daniel J. Bernstein and Bo-Yin Yang ([gcd.cr.yp.to/software.html](https://gcd.cr.yp.to/software.html))



---

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
