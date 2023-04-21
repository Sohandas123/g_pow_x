RUN :
1. make 
2. ./output

OR RUN :
1. gcc g_power_x.c functions.c -o output
2. ./output

FOR Recompilation RUN :
1. make clean



p = p256 = 2^256-2^224+2^192+2^96-1  (NIST suggested P-256)
  = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff (in hex)

This will print an array of four 64-bit long integers. This is the g^x (mod p).
The actual number is :  g^x[0] + g^x[1]*(2^64) + g^x[2]*(2^64)^2 + g^x[3]*(2^64)^3   (in decimal system)
