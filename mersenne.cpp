#include "mersenne.h"

/**
 * Return a Mersenne number
 * Usage: ./mersenne p
 * p: exponent on 2^p-1
 * */
int main(int argc, char const *argv[])
{

    mpz_class m;
    unsigned int p;

    p = (unsigned int)atoi(argv[1]);
    mpz_ui_pow_ui(m.get_mpz_t(), 2, p);
    m = m - 1;

    /*
    if (sizeof(*argv) > 1)
    {
        const char *arg = argv[2];
    }
    */
    gmp_printf("%Zd", m.get_mpz_t());

    /* code */
    return 0;
}
