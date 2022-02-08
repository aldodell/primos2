#include <stdio.h>
#include <gmpxx.h>
#include <thread>
#include <chrono>
#include <vector>
#include <mutex>
#include <bitset>

using namespace std;
using namespace chrono;

mutex mtx;

void factorize(unsigned int p, mpz_class begin, mpz_class final, bool &done, int &threadsCounter)
{

    threadsCounter++; // Increment counter

    mpz_t a, b, alfa, beta, tmp, root, end, offsetBeta;
    unsigned int betaDiff, n;
    betaDiff = p * 2;
    unsigned int p2 = 2 * p;
    unsigned int p1 = p + 1;

    // initial settings
    mpz_init_set(a, begin.get_mpz_t());
    mpz_init_set(end, final.get_mpz_t());
    mpz_init_set_ui(b, 0);
    mpz_init_set_ui(alfa, 0);
    mpz_init_set_ui(beta, 0);
    mpz_init_set_ui(tmp, 0);
    mpz_init_set_ui(offsetBeta, 0);

    // set alfa:
    mpz_ui_pow_ui(alfa, 2, p);  // alfa = 2^p
    mpz_sub_ui(alfa, alfa, 2);  // alfa = alfa - 2
    mpz_div_ui(alfa, alfa, p2); // alfa = alfa / 2p
    mpz_sub(alfa, alfa, a);     // alfa = alfa - a

    // set beta:
    mpz_set_ui(beta, p2);          // beta = 2p
    mpz_add_ui(beta, beta, 1);     // beta = beta +1
    mpz_mul_ui(offsetBeta, a, p2); // offsetBeta = (a-1)*2*p
    mpz_sub_ui(offsetBeta, offsetBeta, p2);
    mpz_add(beta, beta, offsetBeta); // beta = beta + betaOffset

    mtx.lock();
    gmp_printf("\nalfa: %Zd / beta:%Zd /a:%Zd /b:%Zd\n", alfa, beta, a, b);
    mtx.unlock();

    // take time
    auto start = high_resolution_clock::now();
    int z = 0;

    while (!done)
    {
        if (mpz_divisible_p(alfa, beta) != 0)
        {
            done = true;
            mpz_div(b, alfa, beta);

            mtx.lock();
            gmp_printf("\na: %Zd\nb: %Zd\n", a, b);
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<seconds>(stop - start);
            printf("\ntime: %llu\n", duration.count());
            threadsCounter--;
            mtx.unlock();

            break;
        }

        // mpz_sub_ui(alfa, alfa, p2);
        mpz_sub_ui(alfa, alfa, 1);
        mpz_add_ui(beta, beta, betaDiff);
        mpz_add_ui(a, a, 1);

        // Reach square root or beyond final mileston
        if (mpz_cmp(a, end) > 0 || done)
        {
            mtx.lock();
            gmp_printf("\nFinish: init: %Zd / final:%Zd\n", begin.get_mpz_t(), final.get_mpz_t());
            threadsCounter--;
            mtx.unlock();
            break;
        }

        z++;
        if (z == 500000000)
        {

            mtx.lock();
            z = 0;
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(stop - start);
            start = high_resolution_clock::now();
            gmp_printf("\ntime: %llu\na:%Zd", duration.count(), a);
            mtx.unlock();
        }
    }

    // mpz_clears(a, b, alfa, beta, betaDiff, tmp);
}

bool isDivisible(mpz_t a, mpz_t b)
{
    // MULTIPLY B BY TO HOW MANY TIMES FITS ON A
    mpz_t q, r;
    mpz_init_set(q, b);
    mpz_init(r);

    while (mpz_cmp(q, a) < 0)
    {
        mpz_set(r, q);
        mpz_mul_2exp(q, q, 1);
    }

    // subtract a - r;
    mpz_sub(a, a, r);

    // substract a - b
    while (mpz_cmp_ui(a, 0) > 0)
    {
        mpz_sub(a, a, b);
    }

    gmp_printf("%Zd", a);
    return false;
}

bool test(unsigned int p, mpz_t fa, mpz_t fb)
{

    mpz_t r;
    mpz_init_set_ui(r, 0);

    int k = 0; // carry
    int i, j, ij, ij1;
    unsigned int sizea, sizeb;
    sizea = mpz_sizeinbase(fa, 2);
    sizeb = mpz_sizeinbase(fb, 2);

    for (i = 0; i < sizeb; i++)
    {
        if (mpz_tstbit(fb, i))
        {
            for (j = 0; j < sizea; j++)
            {
                ij = i + j;
                ij1 = i + j + 1;
                if (mpz_tstbit(fa, j))
                {
                    if (!mpz_tstbit(r, ij))
                    {
                        mpz_setbit(r, ij);
                    }
                    else
                    {
                        mpz_clrbit(r, ij);
                        k = 1;
                    }
                }

                if (k)
                {

                    if (!mpz_tstbit(r, ij1))
                    {
                        mpz_setbit(r, ij1);
                        k = 0;
                    }
                    else
                    {
                        mpz_clrbit(r, ij1);
                        k = 1;
                    }
                }
            }
        }

        if (!mpz_tstbit(r, i))
            return false;
    }
    gmp_printf("%Zd\n", r);
    bool z = mpz_popcount(r) == p;
    return z;
}

/**
 * @brief Multiplica fa x fb en un array de p bits para evaluar si es un mersenne
 *
 * @param p
 * @param fa
 * @param fb
 * @return true Si fa*fb = 2^p-1 es Mersenne (no necesariamente primo)
 * @return false Si fa*fb 2^p-1 no es Mersenne.
 */
bool test(unsigned int p, unsigned int fa, unsigned int fb)
{
    vector<bool> a;
    vector<bool> b;
    vector<bool> r;

    r.resize(p, 0);

    int k = 0; // carry

    // cast to binary
    int t = fa;
    while (t > 0)
    {
        a.push_back(t & 1);
        t >>= 1;
    }
    t = fb;
    while (t > 0)
    {
        b.push_back(t & 1);
        t >>= 1;
    }
    int i, j;

    for (i = 0; i < b.size(); i++)
    {
        if (b[i])
        {
            for (j = 0; j < a.size(); j++)
            {
                if (a[j])
                {
                    if (!r[i + j])
                    {
                        r[i + j] = 1;
                    }
                    else
                    {
                        r[i + j] = 0;
                        k = 1;
                    }
                }

                if (k)
                {
                    if (!r[i + j + 1])
                    {
                        r[i + j + 1] = 1;
                        k = 0;
                    }
                    else
                    {
                        r[i + j + 1] = 0;
                        k = 1;
                    }
                }
            }
        }
        if (!r[i])
            return false;
    }

    bool w = true;
    for (bool q : r)
    {
        w &= q;
    }

    return w;
}

int main(int argc, char const *argv[])
{

    /*
     test(11, 23, 89);
    mpz_class fa, fb;
    fa = 23;
    fb = 89;
    bool ss = isDivisible(fb.get_mpz_t(), fa.get_mpz_t());

    test(11, fa.get_mpz_t(), fb.get_mpz_t());
    */

    int MAX_THREADS = 4;

    int threadsCounter = 0;
    unsigned int p, i, q;
    p = (unsigned int)atoi(argv[1]);
    q = (unsigned int)atoi(argv[2]);

    MAX_THREADS = q;

    mpz_class begin, final, THRESHOLD, MERSENNE, ROOT;
    bool done = false;

    // Mersenne number.. and OMEGA: root(2^p-2/(2p))
    mpz_ui_pow_ui(MERSENNE.get_mpz_t(), 2, p);
    mpz_sub_ui(MERSENNE.get_mpz_t(), MERSENNE.get_mpz_t(), 1);
    mpz_sqrt(ROOT.get_mpz_t(), MERSENNE.get_mpz_t());

    // split bby threads numbers
    mpz_div_ui(THRESHOLD.get_mpz_t(), ROOT.get_mpz_t(), MAX_THREADS);
    mpz_set_ui(begin.get_mpz_t(), 1);                                     // first "A"
    mpz_add(final.get_mpz_t(), begin.get_mpz_t(), THRESHOLD.get_mpz_t()); // Last "A"

    vector<thread> threads;

    /*
        while (!done)
        {
            while (threadsCounter < MAX_THREADS)
            {
                if (done)
                {
                    break;
                }
    */

    for (int i = 0; i < MAX_THREADS; i++)
    {
        threads.push_back(thread(factorize, p, begin, final, ref(done), ref(threadsCounter)));
        mpz_add_ui(begin.get_mpz_t(), final.get_mpz_t(), 1);
        mpz_add(final.get_mpz_t(), final.get_mpz_t(), THRESHOLD.get_mpz_t());
        /*
        mtx.lock();
        gmp_printf("\nThreads counter: %d\n", threadsCounter);
        mtx.unlock();
        */
    }
    /*
}
}
*/

    for (auto &th : threads)
    {
        if (th.joinable())
        {
            th.join();
        }
    }

    /* code */
    return 0;
}
