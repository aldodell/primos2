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
    int k = 1000000;
    mpz_t M;

    mpz_init(M);
    mpz_ui_pow_ui(M, 2, p);
    mpz_sub_ui(M, M, 1);
    int f = 2 * p * k + 1;
    auto start = high_resolution_clock::now();

    for (int i = 0; i < 1000000000  ; i++)
    {
        mpz_divisible_ui_p(M, f);
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    printf("%lld\n", duration.count());

    /* code */
    return 0;
}
