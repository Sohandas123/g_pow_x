#include <stdio.h>
#include <stdint.h>
#include "g_power_x.h"

int main()
{
    // square and multiply to calculate  z = g^x (mod p)
    uint64_t g[4] = {14317131134123807145u, 165634889443579316u, 10579839724117548515u, 12689480371216343084u};
    uint64_t x[4] = {18265553712439590882u, 2017884693948405437u, 8064836623372059513u, 16743275605901433557u};
 
    uint64_t z[4] = {0};
    sqm_mod_p(g, x, z); // z = g^x (mod p)

    printf("Each of g^x[i]'s are of 28-bit long except the last g^x[9] which is 4-bit long so that 28*9 + 4 = 256 \n");
    for (int i = 0; i < 4; i++)
    {
       printf(" g^x[%d] = %lu\n", i, z[i]);
    }

    return 0;
}