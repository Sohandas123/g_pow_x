#include <stdio.h>
#include <stdint.h>
#include "g_power_x.h"


long long int p[10] = {268435455, 268435455, 268435455, 4095, 0, 0, 16777216, 0, 268435455, 15}; // 256-bit prime (NIST P-256)
// p = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
long long int mu[11] = {327680, 3145728, 0, 0, 268435455, 268435439, 268435199, 268431359, 268435455, 1048575, 16777216}; // precomputed "mu" for Barrett reduction
// mu = floor( ((2^28)^20)/p )

/* --------------------------------- Function Definitions ------------------------------------*/

void add(long long int a[], long long int b[], long long int c[], int size) // a and b will be added and the sum will be stored in c
{
    long long int carry = 0;
    for(int i =0; i<size; i++)
        c[i] = 0;
    for (int i = 0; i < size; i++)
    {
        c[i] = a[i] + b[i] + carry;
        carry = (c[i] >> 28) & 1;
        c[i] = c[i] & ((1 << 28) - 1); // (1<<28) - 1 = 2^28 - 1
    }
}


void mult(long long int a[], long long int b[], long long int m[], int size) // a and b will be multiplied and the product will be stored in m
{
    long long int restof28bit;
    for (int i = 0; i < 2 * size; i++)
    {
        m[i] = 0;
    }
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            m[i + j] = m[i + j] + (a[i] * b[j]);
            restof28bit = m[i + j] >> 28;
            m[i + j] = m[i + j] & ((1 << 28) - 1);     // (1<<28) - 1 = 2^28 - 1
            m[i + j + 1] = m[i + j + 1] + restof28bit; // Carry forwarding
        }
    }
}


void sub(long long int a_o[], long long int b_o[], long long int s[], int size) // s = a_o - b_o with the assumption that a_o > b_o , o stands for original
{   long long int a[size];
    long long int b[size];
    int i;
    for (i = 0; i < size; i++)  // initialization
    {
        s[i] = 0;
        a[i] = a_o[i];
        b[i] = b_o[i];
    }

    /* Using 2's complemnet method */
    long long int I[size], b_2s_com[size];
    for (i = 0; i < size; i++)
    {
        I[i] = 0;
        b[i] = (~b[i])&((1<<28)-1);  // 1's complement is done 
    }
    I[0] = 1;  // I = 1
    add(b, I, b_2s_com, size);  // 2's complement of b is done
    add(a, b_2s_com, s, size);  // s = a + b_2s_com
}


// check whether  x > y or x = y or x < y
// return  1 when x>y
//         2 when x<y
//         3 when x=y
int check(long long int x[], long long int y[], int size)
{
    int flag = 3, i;
    for (i = size - 1; i >= 0; i--)
    {
        if (x[i] > y[i])
        {
            flag = 1;
            return flag;
        }
        if (x[i] < y[i])
        {
            flag = 2;
            return flag;
        }
    }
    return flag;
}



void sub_mod_p(long long int a[], long long int b[], long long int s[], int size) // s = a - b  if a>b , s = p - (b-a) if a<b
{
    long long int pr[size];
    for (int i = 0; i < 10; i++)
        pr[i] = p[i];
    if(size>10)
    {
        for (int i = 10; i < size; i++)
        pr[i] = 0;
    }
    int flag = check(a,b,size);  // flag = 1 when a>b , flag = 2 when a<b and flag = 3 when a=b
    if ((flag == 1) || (flag == 3)) // when a >= b
        sub(a, b, s, size);
    if (flag == 2) // when a < b , then subtraction in mod p will be p - (b - a)
    {
        long long int temp[size]; // for storing b - a
        for (int i = 0; i < size; i++)
            temp[i] = 0;
        sub(b, a, temp, size);  // temp = b - a
        sub(pr, temp, s, size); // s = p - temp = p - (b - a) = p + (a - b)
    }
}


void Barrett(long long int x[], long long int x_Barret_p[], int size) // Barrett reduction with ( mod p), size means size of input x[]
{
    long long int q[11]; // since size of "mu" is 28*11 bits
    for (int i = 0; i < 11; i++)
    {
        q[i] = 0; // initialization
    }
    // q = floor(x/{B^(k-1)}) , k = number of 28-bit blocks used to store the value of p, so k = 10, here and B = 2^28, so to divide x by B^9 we will right shift x by 9*28 bits and store it in q
    int i = 0;
    while (i < (size - 9))
    {
        q[i] = x[i + 9];
        i++;
    }                      // q = floor(x/{B^(k-1)}) is done.
    long long int q22[22]; // to store  q*mu , both of them are of 11*28 bits so their product may be of 2*11*28 bit long
    for (int i = 0; i < 22; i++)
    {
        q22[i] = 0; // initialization
    }
    mult(q, mu, q22, 11);  // q22 = q*mu  is done
    long long int q11[11]; // since q11 = floor(q22/{B^(k+1)}) , q22 is of 22*28 bits and k = 10
    for (int i = 0; i < 11; i++)
    {
        q11[i] = 0; // initialization
    }
    int j = 0; // initialization
    while (j < 11)
    {
        q11[j] = q22[j + 11];
        j++;
    } // q11 = floor(q22/{B^(k+1)}) is done
    // r = x - q11p11 (mod B^(k+1)), B = 2^28, k= 10

    long long int q11p11[22] = {0}; // sizeof(q11) = 11 block, sizeof(p) = 10, so we have to extend it to 11 blocks and we say it p11 also have to extend sizeof(x) = 20 to 22 blocks
    long long int p11[11];          // extending block size of p from 10 to 11
    for (int i = 0; i < 10; i++)    // since sizeof(p) is 10 blocks means 10 long long int is used
    {
        p11[i] = p[i];
    }
    p11[10] = 0;
    mult(p11, q11, q11p11, 11);  // q11p11 = p11 * q11

    // r = x - q11p11 (mod B^(k+1)), first we will calculate x(mod B^(k+1)) and q11p11(mod B^(k+1)) then will subtract in  mod B^(k+1)
    long long int x11[11] = {0}; // for storing the value x (mod B^(k+1)), k=10 here
    long long int qp[11] = {0};  // for storing the value q11p11 (mod B^(k+1)), k=10 here
    for (int i = 0; i < 11; i++)
    {
        x11[i] = x[i];
        qp[i] = q11p11[i];
    }
    long long int r11[11] = {0};  // r11 = x11 - qp (mod B^(k+1))
    int flg = check(x11, qp, 11);  // flg = 1 if x11>qp, flg = 2 if x11<qp, flag = 3 if x11=qp
    if((flg == 1) || (flg == 3))  // if x11 >= qp
        sub(x11, qp, r11, 11);
    if(flg == 2)  // if x11 < qp
        {
            long long int tmp[11] = {0};
            sub(qp, x11, tmp, 11);  // temp = qp - x11  is done 
            //to do mod B^(k+1) we will use 2's complement method as follows
            for(int i = 0; i<11; i++)
            {
                tmp[i] = (~tmp[i])&((1<<28)-1);  // 1's complement in 28-bit is done
            }
            long long int I[11] = {0};
            I[0] = 1;  
            add(tmp, I, r11, 11);  // 2's complement is done
        }

    // Checking whether r >= p or not
    int flag = check(r11, p11, 11); // flag = 1 if r11>p11, flag = 2 if r11<p11, flag = 3 if r11=p11

    while ((flag == 1) || (flag == 3)) // do successive subtraction when r >= p
    {
        // printf("\n Entered in Barrett While loop   flag = %d\n", flag);
        long long int temp[11] = {0}; // to copy r11
        for (int i = 0; i < 11; i++)
        {
            temp[i] = r11[i];
            r11[i] = 0;
        }
        sub(temp, p11, r11, 11); // r11 = temp - p11 = r11(old) - p11

        // Checking whether r >= p or not
        flag = check(r11, p11, 11); // flag = 1 when r > p, flag = 2 when r < p, flag = 0 when r =p
    }
    static int count = 0;
    count++;

    for (int i = 0; i < 10; i++)
    {
        x_Barret_p[i] = r11[i];
    }
}
//base conversion 64 bit to 28 bit
void base64to28(uint64_t m64[4], long long int m28[10])
{
    long long int k = ((1<<28) - 1);
    m28[0] = m64[0] & k;
    m28[1] = (m64[0] >> 28) & k;
    m28[2] = (m64[0] >> 56) & k;
    m28[2] = m28[2] | ((m64[1] & ((1<<20) - 1)) << 8 );
    m28[3] = (m64[1] >> 20) & k;
    m28[4] = (m64[1] >> 48) & k;
    m28[4] = m28[4] | ((m64[2] & ((1<<12) - 1)) << 16);
    m28[5] = (m64[2] >> 12) & k;
    m28[6] = (m64[2] >> 40) & k;
    m28[6] = m28[6] | ((m64[3] & ((1<<4) - 1)) << 24);
    m28[7] = (m64[3]>>4) & k;
    m28[8] = (m64[3]>>32) & k;
    m28[9] = m64[3] >> 60;
}

//base conversion 28 bit to 64 bit
void base28to64(long long int n28[10], uint64_t n64[4])
{
    n64[0] = (uint64_t)n28[0] | ((uint64_t)n28[1] << 28 ) | (((uint64_t)n28[2] & ((1<<8)-1)) << 56);
    n64[1] = ((uint64_t)n28[2] >> 8 ) | ((uint64_t)n28[3] << 20) | (((uint64_t)n28[4] & ((1<<16)-1))<< 48);
    n64[2] = ((uint64_t)n28[4] >> 16) | ((uint64_t)n28[5] << 12) | (((uint64_t)n28[6] & ((1<<24)-1)) << 40);
    n64[3] = ((uint64_t)n28[6] >> 24) | ((uint64_t)n28[7] <<  4) | (((uint64_t)n28[8] << 32)) | (((uint64_t)n28[9] & ((1<<4) - 1)) << 60);
}

// Square and multiply ( used for finding inverse in Z*_p)
void sqm_mod_p(uint64_t y64[], uint64_t d64[], uint64_t z64[])
{
    int size = 10;
    long long int y[10] = {0}, d[10] = {0};
    base64to28(y64, y);
    base64to28(d64, d);
    long long int z[10] = {0};  
    z[0] = 1;
    long long int z2[2 * size]; // for storing the value of z^2
    for (int i = 0; i < 2 * size; i++)
    {
        z2[i] = 0; // initialization
    }
    for (int i = size - 1; i >= 0; i--) // since sizeof(d) = 10 block = 10 *28 bit
    {
        for (int j = 27; j >= 0; j--) // since each block is of 28-bit long
        {
            mult(z, z, z2, size); // z2 = z*z, 10 is the size of inputs z's
            for (int j = 0; j < 10; j++)
                z[j] = 0;               // initialization of before storing values gotten from Barrett
            Barrett(z2, z, 2 * size);   // z = z^2 (mod p) is done
            long long int zy[2 * size]; // for storing z*y
            for (int i = 0; i < 2 * size; i++)
                zy[i] = 0; // initialization
            if (((d[i] >> j) & 1) == 1)
            {
                mult(z, y, zy, size); // zy = z*y
                for (int k = 0; k < 10; k++)
                    z[k] = 0;             // initialization before storing values gotten from Barrett
                Barrett(zy, z, 2 * size); // z = zy (mod p) is done
            }
        }
    }
    base28to64(z, z64);
}