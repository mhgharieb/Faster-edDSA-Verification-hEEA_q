b = 256
q = 2**255 - 19
l = 2**252 + 27742317777372353535851937790883648493

def expmod(b,e,m):
  if e == 0: return 1
  t = expmod(b,e/2,m)**2 % m
  if e & 1: t = (t*b) % m
  return t

def inv(x):
  return expmod(x,q-2,q)

d = -121665 * inv(121666)
I = expmod(2,(q-1)/4,q)

def xrecover(y):
  xx = (y*y-1) * inv(d*y*y+1)
  x = expmod(xx,(q+3)/8,q)
  if (x*x - xx) % q != 0: x = (x*I) % q
  if x % 2 != 0: x = q-x
  return x

By = 4 * inv(5)
Bx = xrecover(By)
B = [Bx % q,By % q]

def edwards(P,Q):
  x1 = P[0]
  y1 = P[1]
  x2 = Q[0]
  y2 = Q[1]
  x3 = (x1*y2+x2*y1) * inv(1+d*x1*x2*y1*y2)
  y3 = (y1*y2+x1*x2) * inv(1-d*x1*x2*y1*y2)
  return [x3 % q,y3 % q]

def scalarmult(P,e):
  if e == 0: return [0,1]
  Q = scalarmult(P,e//2)
  Q = edwards(Q,Q)
  if e & 1: Q = edwards(Q,P)
  return Q

def radix255(x):
  x = x % q
  if x + x > q: x -= q
  x = [x,0,0,0,0,0,0,0,0,0]
  bits = [26,25,26,25,26,25,26,25,26,25]
  for i in range(9):
    carry = (x[i] + 2**(bits[i]-1)) / 2**bits[i]
    x[i] -= carry * 2**bits[i]
    x[i + 1] += carry
  result = ""
  for i in range(9):
    result = result+str(x[i])+","
  result = result+str(x[9])
  return result

def radix51(x):
  x = x % q
  y = [0]*5  
  for i in range(5):
    y[i] = x % (1<<51)
    x = x>>51
  result = "{"
  for i in range(4):
    result = result + "{0:#0{1}x}".format(y[i],18) + ","
  result = result + "{0:#0{1}x}".format(y[4],18) + "}"
  return result


# Generate `ge25519_niels_sliding_multiples` in radix51 as in 'ed25519-donna-64bit-tables.h'
Bi = B

print("static const ge25519_niels ge25519_niels_sliding_multiples[32] = {")
for i in range(32):
  line="    {"
  line += radix51(Bi[1]-Bi[0]) + ","
  line += radix51(Bi[1]+Bi[0]) + ","
  if (i!=31):
    line += radix51(2*d*Bi[0]*Bi[1]) + "},"
  else:
    line += radix51(2*d*Bi[0]*Bi[1]) + "}"

  print(line)
  Bi = edwards(B,edwards(B,Bi))

print("};")

Bi = B

print("static const ge25519_niels ge25519_niels_sliding_multiples_2B[64] = {")
for i in range(64):
  line="    {"
  line += radix51(Bi[1]-Bi[0]) + ","
  line += radix51(Bi[1]+Bi[0]) + ","
  if (i!=63):
    line += radix51(2*d*Bi[0]*Bi[1]) + "},"
  else:
    line += radix51(2*d*Bi[0]*Bi[1]) + "}"

  print(line)
  Bi = edwards(B,edwards(B,Bi))

print("};")

# Generate `ge25519_niels_sliding_multiples2` in radix51
# B_par = scalarmult(B, 1<<128)
B_par = scalarmult(B, 1<<126)

print("static const ge25519 ge25519_B_par = {")
print("    " + radix51(B_par[0]) + ",")
print("    " + radix51(B_par[1]) + ",")
print("    " + radix51(1) + "},")
print("    " + radix51(B_par[1]*B_par[0] % q))
print("};")

Bi = B_par
print("static const ge25519_niels ge25519_niels_sliding_multiples2[32] = {")
for i in range(32):
  line="    {"
  line += radix51(Bi[1]-Bi[0]) + ","
  line += radix51(Bi[1]+Bi[0]) + ","
  if (i!=31):
    line += radix51(2*d*Bi[0]*Bi[1]) + "},"
  else:
    line += radix51(2*d*Bi[0]*Bi[1]) + "}"

  print(line)
  Bi = edwards(B_par,edwards(B_par,Bi))

print("};")