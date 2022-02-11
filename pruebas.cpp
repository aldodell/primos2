#include <stdio.h>
#include <gmpxx.h>
#include <thread>
#include <chrono>
#include <vector>
#include <mutex>
#include <bitset>
#include <cmath>

using namespace std;
using namespace chrono;

int main(int argc, char const *argv[])
{
    int p = 101;
    int k = 1;
    unsigned long cycles = 10000000000;
    mpz_t M;

    mpz_init(M);
    mpz_ui_pow_ui(M, 2, p);
    mpz_sub_ui(M, M, 1);
    int f = 2 * p * k + 1;

    mpz_t D, d;
    mpz_init_set_ui(D, 1000000);
    mpz_init_set_ui(d, 5);

    auto start = high_resolution_clock::now();

    for (int i = 0; i < cycles; i++)
    {
        mpz_divisible_p(D, d);
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    printf("%lld\n", duration.count());

    mpz_t t, t2;
    mpz_init(t);
    mpz_init(t2);

    bool flag0 = true;

    start = high_resolution_clock::now();
    for (int i = 0; i < cycles; i++)
    {
        while (flag0)
        {
            mpz_set(t, d);
            mpz_set(t2, d);

            while (mpz_cmp(t2, D) < 0)
            {
                mpz_set(t, t);
                mpz_mul_2exp(t2, t2, 2);
            }

            if (mpz_cmp(t, D) == 0)
            {
                flag0 = false;
                break;
            }

            mpz_set(D, t);
        }
    }

    stop = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(stop - start);
    printf("%lld\n", duration.count());

    /* code */
    return 0;
}
