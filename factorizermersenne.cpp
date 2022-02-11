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

mutex mtx;

void factorize6(unsigned int p)
{
    mpz_t fa, t, p2, r, ft, M;

    unsigned int da, db, dp, x, p1;
    mpz_init(fa);
    mpz_init(p2);
    mpz_init(t);
    mpz_init(r);
    mpz_init(ft);
    mpz_init(M);

    mpz_set_ui(p2, 2 * p);
    mpz_ui_pow_ui(M, 2, p);
    mpz_sub_ui(M, M, 1);

    //   da = 0;
    //  db = 0;
    p1 = p + 1;

    //  dp = mpz_sizeinbase(fa, 2);

    mpz_set_ui(fa, 1);
    while (true)
    {
    stageA:

        // Construimos el factor A mediante 2p+1
        //  fa = fa + 2p
        mpz_add(fa, fa, p2);
        //   da = mpz_sizeinbase(fa, 2) + 1;

        x = 0; // posicion a evaluar;
        mpz_set(t, fa);
        mpz_set_ui(r, 1);

        while (true)
        {
        stageB:
            x++;

            if (x == p1)
            {
                goto stageA;
            }

            /*
             evaluamos con 0: suponemos que en el factor B, el bit x lo
            seteamos en cero. Luego, no suma nada al factor temporal T
            */
            if (mpz_tstbit(t, x) == 1)
            {
                goto stageB;
            }

            // Evaluamos con 1;
            // cuando el bit x del factor b est'a en 1 sumamos el factor a desplazad
            // una posicion a t

            mpz_mul_2exp(ft, fa, x); //<deszplamiento a la izquierda del factor temporal
            mpz_add(t, t, ft);       // sumamos
            // Si encontramos el bit a la izquierda con 1 entonces agregamos el bit al resultado
            // o factor B
            if (mpz_tstbit(t, x) == 1)
            {
                mpz_setbit(r, x);
                if (mpz_cmp(t, M) == 0)
                {
                    goto finalStage;
                }
                goto stageB;
            }
            else
            {
                goto stageA; // evaluamos otro factor
                // porque ningun valor en el bit x del factor A genera 1 en el factor b en
                // esta posici'on
            }
        }
    finalStage:
        gmp_printf("%Zd\n", fa);
        return;
    }

    /*

        // Construimos el factor B para que quepa en la cantidad de bits
        mpz_set_ui(fb, p);
        db = mpz_sizeinbase(fb, 2);
        mpz_mul_2exp(fb, fb, p + 1 - da - db);
        mpz_add_ui(fb, fb, 1);
        */
}

void factorize5(unsigned int p)
{

    mpz_class fa, b, t, fa2, r, z, root, m;
    int x, i, cycles;
    long long speed, cyclesBoundry;
    unsigned int da, dt; // digitios del factor a y del factor temporal a2

    fa = 1;
    mpz_ui_pow_ui(m.get_mpz_t(), 2, p);
    m = m - 1;

    auto start = high_resolution_clock::now();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cycles = 0;
    cyclesBoundry = 1000000;
    // 36158000001
    while (true)
    {
    begin:
        fa += (2 * p);
        b = 1;
        t = fa; // asumiendo que el multicando termina en 1
        x = 0;
        r = 1;
        // da = (int)mpz_sizeinbase(fa.get_mpz_t(), 2);
        cycles++;
        if (cycles == cyclesBoundry)
        {
            stop = high_resolution_clock::now();
            duration = duration_cast<milliseconds>(stop - start);
            speed = cyclesBoundry / duration.count();
            gmp_printf("Speed: %lld, fa:%Zd\n", speed, fa.get_mpz_t());
            cycles = 0;
            start = high_resolution_clock::now();
        }

        while (true)
        {
            if (mpz_tstbit(t.get_mpz_t(), x) == 0)
            {
                mpz_mul_2exp(fa2.get_mpz_t(), fa.get_mpz_t(), x);
                mpz_add(t.get_mpz_t(), t.get_mpz_t(), fa2.get_mpz_t());
                // dt = (int)mpz_sizeinbase(t.get_mpz_t(), 2);
                if (mpz_cmp(t.get_mpz_t(), m.get_mpz_t()) > 0)
                {

                    goto begin;
                }
                if (mpz_cmp(t.get_mpz_t(), m.get_mpz_t()) == 0)
                {
                    gmp_printf("%Zd\n", fa.get_mpz_t());
                    return;
                }
            }
            x++;
        }
    }
}

void factorize4(unsigned int p)
{

    unsigned long long fa, b, t, fa2, x, r, i, z, root, m;
    fa = 1;
    m = pow(2, p) - 1;
    root = sqrt(m);

    while (true)
    {
    begin:
        fa += (2 * p);
        b = 1;
        t = fa; // asumiendo que el multicando termina en 1
        x = 0;
        r = 1;
        while (true)
        {
            // numero de desplazamientos
            x++;
            b <<= 1;

            // no pas'o la prueba
            if ((t & b) == 0)
            {
                fa2 = fa << x;
                t = t + fa2;
                r = r + (1 << x);
                if (t > m)
                    goto begin;
            }
            z = 0;
            for (i = 0; i < p; i++)
            {

                if ((t & (1 << i)))
                {
                    z++;
                }
            }
            if (z == p)
            {
                printf("%llu\n", fa);
                return;
            }
        }
    }
}

void factorize3(int p)
{
    mpz_t alfa, beta, alfa2, a, b, t, t2, progress, thresold;
    mpz_init(alfa);
    mpz_init(beta);
    mpz_init(alfa2);
    mpz_init(a);
    mpz_init(b);
    mpz_init(t);
    mpz_init(t2);
    mpz_init(progress);
    mpz_init(thresold);

    // init alfa:
    mpz_ui_pow_ui(alfa, 2, p);
    mpz_sub_ui(alfa, alfa, 2);
    mpz_div_ui(alfa, alfa, 2 * p);
    mpz_sub_ui(alfa, alfa, 1);
    mpz_set(alfa2, alfa);

    // init beta:
    int p2 = 2 * p;
    mpz_set_ui(beta, 2 * p + 1);

    // init a
    mpz_set_ui(a, 1);
    int q;

    // setup progress
    mpz_set_ui(thresold, 1000);
    mpz_add_ui(progress, thresold, 0);
    printf("\n[");

    while (true)
    {
        // Evaluate if beta | alfa
        while (true)
        {
            mpz_set(t, beta);
            mpz_set(t2, t);
            while (mpz_cmp(t2, alfa2) < 0)
            {
                mpz_set(t, t2);
                mpz_mul_2exp(t2, t2, 1);
            }

            q = mpz_cmp(alfa2, t);
            mpz_sub(alfa2, alfa2, t);
            if (q < 0)
            {
                break;
            }
            else if (q == 0)
            {
                gmp_printf("]\n\na=%Zd\n", a);
                return;
            }
        }

        mpz_sub_ui(alfa, alfa, 1);
        mpz_set(alfa2, alfa);
        mpz_add_ui(beta, beta, p2);
        mpz_add_ui(a, a, 1);
        if (mpz_cmp(a, progress) == 0)
        {
            // mtx.lock();
            gmp_printf("*");
            mpz_add(progress, progress, thresold);
            // mtx.unlock();
        }
    }
}

/*
q = mpz_cmp(alfa2, beta);
            if (q < 0)
            {
                mpz_set(alfa2, alfa);
                mpz_add_ui(beta, beta, p2);
                mpz_add_ui(a, a, 1);
            }
            */

void factorize2(int p)
{

    mpz_t M, N, f, t, t2, root;
    mpz_init(M);
    mpz_init(N);
    mpz_init(f);
    mpz_init(t);
    mpz_init(t2);

    mpz_set_ui(f, 2 * p + 1);

    // Mersenne
    mpz_ui_pow_ui(M, 2, p);
    mpz_sub_ui(M, M, 1);

    mpz_set(N, M); // N=M
    mpz_set(t, f); // first possible factor 2p+1

    auto start = high_resolution_clock::now();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);

    int q;
    while (true)
    {
        stop = high_resolution_clock::now();
        duration = duration_cast<seconds>(stop - start);
        if (duration.count() > 10)
        {
            mtx.lock();
            gmp_printf("evaluating: %Zd\n", f);
            start = high_resolution_clock::now();
            mtx.unlock();
        }

        mpz_set(t, f);
        mpz_set(t2, f);

        while (mpz_cmp(N, t2) > 0)
        {
            mpz_set(t, t2);
            mpz_mul_2exp(t2, t2, 1);
        }

        mpz_sub(N, N, t);

        q = mpz_cmp(N, f);
        if (q == 0)
        {
            gmp_printf("Factor: %Zd\n", N);
            return;
        }
        else if (q < 0)
        {
            mpz_add_ui(f, f, (2 * p));
            mpz_set(N, M);
        }
    }
}

//   15.729.190.359
// 7.432.339.208.719

void puzzle(unsigned int p)
{
    mpz_t M, root, fa, fb;
    int p2;

    // inits
    mpz_inits(M, fa, fb, root);

    p2 = 2 * p;

    // set Mersenne
    mpz_ui_pow_ui(M, 2, p);

    // root
    mpz_sqrt(root, M);

    // factors
    mpz_set_ui(fa, p2 + 1);
    mpz_set(fb, root);
    mpz_sub_ui(fb, fb, 1);
    mpz_div_ui(fb, fb, p2);
}

void puzzle(unsigned int p, int qq)
{

    vector<bool> fa; // factor a
    vector<bool> fb; // factor b
    vector<bool> r;  // intermediate result;
    int a;           // factor a index
    int b;           // factor b index
    int x;           // factor a index aux
    int y = 0;       // factor b inde aux
    int fan;         // factor a in decimal number
    int fbn;         // factor a in decimal number

    int da; // factor a digits
    int db; // factor b digits
    bool c; // carry
    int t;
    int limitY = 0;

    // Prepare result set:
    r.resize(p, 0);

    // Set decimal numeric factor A
    fan = 1;
    fbn = 0;

    // main loop.
    while (true)
    {

        // Get factor A
        fan += (2 * p);
        // digits a=0
        da = 0;

        // cast to binary factor A
        fa.clear();
        t = fan;
        while (t > 0)
        {
            da++; // counting digits A
            fa.push_back(t & 1);
            t >>= 1;
        }
        db = p - da + 1; // Calculate factor b digits size
        limitY = pow(2, db - 2) - 1;

        // prepare factor b: form 100000...1
        fb.clear();
        fb.push_back(1);
        for (int i = 0; i < db - 2; i++)
        {
            fb.push_back(0);
        }
        fb.push_back(1);

        x = 0;

        y = 2 * p + 1;
        t = 0;
        while (t < (p / 2))
        {
            y += 2 * p;
        }

        while (y < limitY)
        {
            // Across factor a:
            for (a = 0; a < da; a++)
            {

                // if digit on position "a" is 1
                if (fa[a])
                {
                    // across factor b
                    for (b = 0; b < db; b++)
                    {
                        if (fb[b]) // if digit on position "b" on factor b is 1
                        {
                            if (!r[a + b])
                            {
                                r[a + b] = 1;
                                c = 0;
                            }
                            else
                            {
                                r[a + b] = 0;
                                c = 1;
                            }
                        }

                        // Evaluate carry
                        if (c)
                        {
                            if (!r[a + b + 1])
                            {
                                r[a + b + 1] = 1;
                                c = 0;
                            }
                            else
                            {
                                r[a + b + 1] = 0;
                                c = 1;
                            }
                        }
                    }
                }
                // Evaluate if doesn't match with "1" las dgit
                if (!r[x])
                {
                    fb.clear();
                    fb.push_back(1);
                    t = y;
                    while (t > 0)
                    {
                        fb.push_back(t & 1);
                        t >>= 1;
                    }
                    while (fb.size() < db)
                    {
                        fb.push_back(0);
                    }
                    fb.push_back(1);

                    a = -1;
                    b = -1;
                    x = -1;
                    y++;
                    r.clear();
                    r.resize(p, 0);
                }
                x++;
            }

            t = 1;
            for (int i = 0; i < r.size(); i++)
            {
                t = t & r[i];
            }
            if (t)
            {
                goto showResult;
            }
        }
    }

showResult:
    t = 1;
    for (int i = 0; i < db; i++)
    {
        if (fb[i])
        {
            fbn += t;
        }
        t *= 2;
    }
    printf("%d\n\n%d\n", fan, fbn);
}

// 36.793.758.459
// 1.000.000.001
// 398.069.729.532.862
void factorize(unsigned int p, mpz_class begin, mpz_class final, bool &done, int &threadsCounter)
{
    const int OPERATIONS = 1000000000;
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
    gmp_printf("\n*****\n* alfa: %Zd\n* beta:%Zd\n* a:%Zd\n* b:%Zd\n\n", alfa, beta, a, b);
    mtx.unlock();

    // take time
    auto start = high_resolution_clock::now();
    auto stop = high_resolution_clock::now();
    int z = 0;

    while (!done)
    {
        if (mpz_divisible_p(alfa, beta) != 0)
        {
            done = true;
            mpz_div(b, alfa, beta);

            mtx.lock();
            gmp_printf("\na: %Zd\nb: %Zd\n", a, b);
            stop = high_resolution_clock::now();
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
        if (z == OPERATIONS)
        {

            mtx.lock();
            z = 0;
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(stop - start);
            start = high_resolution_clock::now();
            gmp_printf("\ntime: %llu a:%Zd\nSpeed: %d", duration.count(), a, OPERATIONS / duration.count());
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

    // 1.000.000.001
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
    // q = (unsigned int)atoi(argv[2]);

    factorize6(p);

    // puzzle(p);
    return 0;

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

    auto start = high_resolution_clock::now();
    auto stop = high_resolution_clock::now();

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
    stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    printf("Total time: %lld seconds.\n", duration.count());

    /* code */
    return 0;
}
