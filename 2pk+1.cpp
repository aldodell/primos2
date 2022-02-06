#include "2pk+1.h"

int main(int argc, char const *argv[])
{
    mpz_class m, f, k, r;
    unsigned int p, pp;
    p = (unsigned int)atoi(argv[1]);
    pp = 2 * p;

    mpz_ui_pow_ui(m.get_mpz_t(), 2, p);
    m = m - 1;

    f = 2 * p + 1;

    mpz_sqrt(r.get_mpz_t(), m.get_mpz_t());

    while (true)
    {
        if (mpz_divisible_p(m.get_mpz_t(), f.get_mpz_t()) != 0)
        {
            gmp_printf("F: %Zd", f.get_mpz_t());
            break;
        }

        if (f > r)
        {
            gmp_printf("Primo\n");
            break;
        }

        f = f + pp;
    }

    /* code */
    return 0;
}
