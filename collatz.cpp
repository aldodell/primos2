#include "collatz.h"
using namespace std;

vector<mpz_class> collatz(mpz_t n)
{
    vector<mpz_class> vect;
    mpz_t b;
    mpz_init_set(b, n);
    while (mpz_cmp_ui(b, 1) > 0)
    {
        if (mpz_odd_p(b) != 0)
        {
            mpz_mul_ui(b, b, 3); // b = b * 3
            mpz_add_ui(b, b, 1); // b = b + 1
        }
        else
        {
            mpz_div_ui(b, b, 2); // b = b / 2
        }
        vect.emplace_back(b);
    }
    return vect;
}

/**
 * use: collatz n [args]
 * */
int main(int argc, char const *argv[])
{
    mpz_t n;
    mpz_set_str(n, argv[1], 10);
    vector<mpz_class> r = collatz(n);

    if (sizeof(*argv) > 1)
    {
        const char *arg = argv[2];

        //Show steps
        if (strcmp(arg, "-s") == 0)
        {
            printf("%d", (int)r.size());
        }

        //show numbers list
        else if (strcmp(arg, "-l") == 0)
        {
            for (mpz_class p : r)
            {
                gmp_printf("%Zd\n", p.get_mpz_t());
            }
        }
    }

    return 0;
}
