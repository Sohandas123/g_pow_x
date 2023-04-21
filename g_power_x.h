#ifndef G_POWER_X_DOT_H  
#define G_POWER_X_DOT_H
#include <stdint.h>

void add(long long int a[], long long int b[], long long int c[], int size);  // c = a + b
void mult(long long int a[], long long int b[], long long int m[], int size); // m = a * b
void sub(long long int a[], long long int b[], long long int s[], int size);  // s = a - b  with assumption a > b
void sub_mod_p(long long int a[], long long int b[], long long int s[], int size);  // subtraction in Z_p
int check(long long int x[], long long int y[], int size); // returns 1 when x>y, 2 when x<y , 3 when x=y

// Barrett reduction with ( mod p). It takes input x of double size and returns x_Barret_p = x (mod p) of given size
void Barrett(long long int x[], long long int x_Barrett_p[], int size);

//base conversion 64 bit to 28 bit
void base64to28(uint64_t m64[4], long long int m28[10]);

//base conversion 28 bit to 64 bit
void base28to64(long long int n28[10], uint64_t n64[4]);

//SQUARE AND MULTIPLY. Computes z = y^d (mod p)
void sqm_mod_p(uint64_t y[], uint64_t d[], uint64_t z[]);       


#endif  