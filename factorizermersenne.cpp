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

    mpz_t a, b, alfa, beta, tmp, root, end;
    unsigned int betaDiff, n;
    betaDiff = p * 2;
    unsigned int p2 = 2 * p;
    unsigned int p1 = p + 1;

    //initial settings
    mpz_init_set(a, begin.get_mpz_t());
    mpz_init_set(end, final.get_mpz_t());
    mpz_init_set_ui(b, 0);
    mpz_init_set_ui(alfa, 0);
    mpz_init_set_ui(beta, 0);
    mpz_init_set_ui(tmp, 0);
    mpz_init_set_ui(root, 0);

    //set root
    mpz_ui_pow_ui(root, 2, p);
    mpz_sub_ui(root, root, 1);
    mpz_sqrt(root, root);
    mpz_sub_ui(root, root, 1);
    mpz_div_ui(root, root, p2);

    //If begin > root, exit now!!!
    if (mpz_cmp(begin.get_mpz_t(), root) > 0)
    {
        /*
        mtx.lock();
        gmp_printf("Root reached before start on a thread!\n");
        mtx.unlock();
        */
        return;
    }

    mtx.lock();
    gmp_printf("\ninit: %Zd / final:%Zd\n", begin.get_mpz_t(), final.get_mpz_t());
    mtx.unlock();

    //set alfa:
    mpz_ui_pow_ui(alfa, 2, p);  // alfa = 2^p
    mpz_sub_ui(alfa, alfa, 2);  //alfa = alfa - 2
    mpz_div_ui(alfa, alfa, p2); //alfa = alfa / 2p
    mpz_sub_ui(alfa, alfa, 1);  //alfa = alfa -1

    //set beta:
    mpz_set_ui(beta, p2);      //beta = 2p
    mpz_add_ui(beta, beta, 1); //beta = beta +1

    //take time
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
            mtx.unlock();

            break;
        }

        // mpz_sub_ui(alfa, alfa, p2);
        mpz_sub_ui(alfa, alfa, 1); //Experimental
        mpz_add_ui(beta, beta, betaDiff);
        mpz_add_ui(a, a, 1);

        //Reach square root or beyond final mileston
        if (mpz_cmp(a, root) > 0 || mpz_cmp(a, end) > 0 || done)
        {
            mtx.lock();
            gmp_printf("\nFinish: init: %Zd / final:%Zd\n", begin.get_mpz_t(), final.get_mpz_t());
            threadsCounter--;
            mtx.unlock();
            break;
        }
        z++;
        if (z == 1000000)
        {

            mtx.lock();
            z = 0;
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(stop - start);
            start = high_resolution_clock::now();
            printf("\ntime: %llu\n", duration.count());
            mtx.unlock();
        }
    }

    // mpz_clears(a, b, alfa, beta, betaDiff, tmp);
}

/**
 * @brief Multiplica fa x fb en un array de p bits para evaluar si es un mersenne
 * 
 * @param p 
 * @param fa 
 * @param fb 
 * @return true 
 * @return false 
 */
bool test(unsigned int p, unsigned int fa, unsigned int fb)
{
    vector<bool> a;
    vector<bool> b;
    vector<bool> r;

    r.resize(p, 0);

    int k = 0; //carry

    //cast to binary
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

    test(11, 23, 89);

    int MAX_THREADS = 4;

    int threadsCounter = 0;
    unsigned int p, i, q;
    p = (unsigned int)atoi(argv[1]);
    q = (unsigned int)atoi(argv[1]);

    MAX_THREADS = q;

    mpz_class begin, final, THRESHOLD, OMEGA;
    bool done = false;

    //Mersenne number.. and OMEGA: root(2^p-2/(2p))
    mpz_ui_pow_ui(OMEGA.get_mpz_t(), 2, p);
    mpz_sub_ui(OMEGA.get_mpz_t(), OMEGA.get_mpz_t(), 1);
    mpz_sqrt(OMEGA.get_mpz_t(), OMEGA.get_mpz_t());
    mpz_div_ui(OMEGA.get_mpz_t(), OMEGA.get_mpz_t(), (2 * p));

    mpz_div_ui(THRESHOLD.get_mpz_t(), OMEGA.get_mpz_t(), MAX_THREADS);

    mpz_set_ui(begin.get_mpz_t(), 1);
    mpz_add(final.get_mpz_t(), begin.get_mpz_t(), THRESHOLD.get_mpz_t());

    //factorize(p, begin, final, done, threadsCounter);

    vector<thread> threads;

    while (!done)
    {
        while (threadsCounter < MAX_THREADS)
        {
            if (done)
            {
                break;
            }
            threads.push_back(thread(factorize, p, begin, final, ref(done), ref(threadsCounter)));
            threadsCounter++;
            mpz_add_ui(begin.get_mpz_t(), final.get_mpz_t(), 1);
            mpz_add(final.get_mpz_t(), final.get_mpz_t(), THRESHOLD.get_mpz_t());
        }
    }

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
