#include "transforms.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if (_OPENMP >= 200203)
#include <omp.h>
#endif

#include "am_sysdep.h"
//#include "am_types.h"
#include "math_const.h"

#ifndef FFT_UNIT_STRIDE
#define FFT_UNIT_STRIDE 0
#endif
#ifndef FHT_UNIT_STRIDE
#define FHT_UNIT_STRIDE 0
#endif

static void fft_dif_iter(double*, unsigned long);
static void fft_dif_iter_seq(double*, unsigned long);
static void fft_dif_rec(double*, unsigned long, int);

static void ifft_dit_iter(double*, unsigned long);
static void ifft_dit_iter_seq(double*, unsigned long);
static void ifft_dit_rec(double*, unsigned long, int);

static void fht_dif_iter(double*, unsigned long);
static void fht_dif_iter_seq(double*, unsigned long);
static void fht_dif_rec(double*, unsigned long, int);

static void fht_dit_iter(double*, unsigned long);
static void fht_dit_iter_seq(double*, unsigned long);
static void fht_dit_rec(double*, unsigned long, int);

void fft_dif(double* z, unsigned long n)
{
    fft_dif_rec(z, n, 1);
    return;
} 

void ifft_dit(double* z, unsigned long n)
{
    ifft_dit_rec(z, n, 1);
    return;
}

void fht_dif(double* x, unsigned long n)
{
    fht_dif_rec(x, n, 1);
    return;
} 


void fht_dit(double* x, unsigned long n)
{
    fht_dit_rec(x, n, 1);
    return;
} 

void hilbert(double* z, unsigned long n)
{
    double x;
    unsigned long i, n2;

    n2 = n << 1;

    fft_dif(z, n);

    for (i = 6; i < n2; i += 4) {
        z[i] = 0.;
        z[i + 1] = 0.;
    }

    z[0] *= 0.5;
    z[1] *= 0.5;
    if (n > 1) {
        z[2] *= 0.5;
        z[3] *= 0.5;
    }

    ifft_dit(z, n);

    x = 2. / (double)n;
    for (i = 0; i < n2; ++i)
        z[i] *= x;
    return;
}  

static void fft_dif_iter(double* z, unsigned long n)
{
    unsigned long i, n2;

    n2 = n << 1;
    for (i = n; i > 1; i >>= 1) {
        double a, b, c, s, t;
        unsigned long i2, j;
        i2 = i << 1;
        t = TWOPI / (double)i;
        a = sin(0.5 * t);
        a *= 2.0 * a;
        b = sin(t);
        c = 1.0;
        s = 0.0;
        for (j = 0; j < i; j += 2) {
            double tmp;
            unsigned long kr, kmax;
            kmax = n2 + j;
            for (kr = j; kr < kmax; kr += i2) {
                double ur, ui;
                unsigned long ki, mr, mi;
                ki = kr + 1;
                mr = kr + i;
                mi = mr + 1;
                ur = z[kr];
                ui = z[ki];
                z[kr] = ur + z[mr];
                z[ki] = ui + z[mi];
                ur -= z[mr];
                ui -= z[mi];
                z[mr] = ur * c - ui * s;
                z[mi] = ur * s + ui * c;
            }
            tmp = c;
            c -= a * c + b * s;
            s -= a * s - b * tmp;
        }
    }
    return;
}

static void fft_dif_iter_seq(double* z, unsigned long n)
{
    unsigned long i, n2;

    n2 = n << 1;
    for (i = n; i > 1; i >>= 1) {
        double a, b, t;
        unsigned long i2, k;
        i2 = i << 1;
        t = TWOPI / (double)i;
        a = sin(0.5 * t);
        a *= 2.0 * a;
        b = sin(t);
        for (k = 0; k < n2; k += i2) {
            double c, s;
            unsigned long j;
            c = 1.0;
            s = 0.0;
            for (j = 0; j < i; j += 2) {
                double ur, ui, tmp;
                unsigned long kr, ki, mr, mi;
                kr = k + j;
                ki = kr + 1;
                mr = kr + i;
                mi = mr + 1;
                ur = z[kr];
                ui = z[ki];
                z[kr] = ur + z[mr];
                z[ki] = ui + z[mi];
                ur -= z[mr];
                ui -= z[mi];
                z[mr] = ur * c - ui * s;
                z[mi] = ur * s + ui * c;
                tmp = c;
                c -= a * c + b * s;
                s -= a * s - b * tmp;
            }
        }
    }
    return;
}



static void fft_dif_rec(double* z, unsigned long n, int nbranch)
{
    double a, b, c, s, t;
    unsigned long nh, kr;

    if (n == 1)
        return;
    if (n <= (unsigned long)(L1_CACHE_BYTES / (2 * sizeof(double)))) {
        if (FFT_UNIT_STRIDE)
            fft_dif_iter_seq(z, n);
        else
            fft_dif_iter(z, n);
        return;
    }
    t = TWOPI / (double)n;
    a = sin(0.5 * t);
    a *= 2.0 * a;
    b = sin(t);
    c = 1.0;
    s = 0.0;
    for (kr = 0; kr < n; kr += 2) {
        double ur, ui, tmp;
        unsigned long ki, mr, mi;
        ki = kr + 1;
        mr = kr + n;
        mi = mr + 1;
        ur = z[kr];
        ui = z[ki];
        z[kr] = ur + z[mr];
        z[ki] = ui + z[mi];
        ur -= z[mr];
        ui -= z[mi];
        z[mr] = ur * c - ui * s;
        z[mi] = ur * s + ui * c;
        tmp = c;
        c -= a * c + b * s;
        s -= a * s - b * tmp;
    }
    nh = n >> 1;
    nbranch <<= 1;
#if (_OPENMP >= 200203)
#pragma omp parallel sections if (nbranch <= omp_get_max_threads()) num_threads(2)
#endif
    {
#if (_OPENMP >= 200203)
#pragma omp section
#endif
        fft_dif_rec(z, nh, nbranch);
#if (_OPENMP >= 200203)
#pragma omp section
#endif
        fft_dif_rec(z + n, nh, nbranch);
    }
    return;
}  


static void ifft_dit_iter(double* z, unsigned long n)
{
    unsigned long i, n2;

    n2 = n << 1;
    for (i = 2; i <= n; i <<= 1) {
        double a, b, c, s, t;
        unsigned long i2, j;
        i2 = i << 1;
        t = -TWOPI / (double)i;
        a = sin(0.5 * t);
        a *= 2.0 * a;
        b = sin(t);
        c = 1.0;
        s = 0.0;
        for (j = 0; j < i; j += 2) {
            double tmp;
            unsigned long kr, kmax;
            kmax = n2 + j;
            for (kr = j; kr < kmax; kr += i2) {
                double vr, vi;
                unsigned long ki, mr, mi;
                ki = kr + 1;
                mr = kr + i;
                mi = mr + 1;
                vr = z[mr] * c - z[mi] * s;
                vi = z[mr] * s + z[mi] * c;
                z[mr] = z[kr] - vr;
                z[mi] = z[ki] - vi;
                z[kr] += vr;
                z[ki] += vi;
            }
            tmp = c;
            c -= a * c + b * s;
            s -= a * s - b * tmp;
        }
    }
    return;
} 



static void ifft_dit_iter_seq(double* z, unsigned long n)
{
    unsigned long i, n2;

    n2 = n << 1;
    for (i = 2; i <= n; i <<= 1) {
        double a, b, t;
        unsigned long i2, k;
        i2 = i << 1;
        t = -TWOPI / (double)i;
        a = sin(0.5 * t);
        a *= 2.0 * a;
        b = sin(t);
        for (k = 0; k < n2; k += i2) {
            double c, s;
            unsigned long j;
            c = 1.0;
            s = 0.0;
            for (j = 0; j < i; j += 2) {
                double vr, vi, tmp;
                unsigned long kr, ki, mr, mi;
                kr = k + j;
                ki = kr + 1;
                mr = kr + i;
                mi = mr + 1;
                vr = z[mr] * c - z[mi] * s;
                vi = z[mr] * s + z[mi] * c;
                z[mr] = z[kr] - vr;
                z[mi] = z[ki] - vi;
                z[kr] += vr;
                z[ki] += vi;
                tmp = c;
                c -= a * c + b * s;
                s -= a * s - b * tmp;
            }
        }
    }
    return;
}  


static void ifft_dit_rec(double* z, unsigned long n, int nbranch)
{
    double a, b, c, s, t;
    unsigned long  nh, kr;

    if (n == 1)
        return;
    if (n <= (unsigned long)(L1_CACHE_BYTES / (2 * sizeof(double)))) {
        if (FFT_UNIT_STRIDE)
            ifft_dit_iter_seq(z, n);
        else
            ifft_dit_iter(z, n);
        return;
    }
    nh = n >> 1;
    nbranch <<= 1;
#if (_OPENMP >= 200203)
#pragma omp parallel sections if (nbranch <= omp_get_max_threads()) num_threads(2)
#endif
    {
#if (_OPENMP >= 200203)
#pragma omp section
#endif
        ifft_dit_rec(z, nh, nbranch);
#if (_OPENMP >= 200203)
#pragma omp section
#endif
        ifft_dit_rec(z + n, nh, nbranch);
    }
    t = -TWOPI / (double)n;
    a = sin(0.5 * t);
    a *= 2.0 * a;
    b = sin(t);
    c = 1.0;
    s = 0.0;
    for (kr = 0; kr < n; kr += 2) {
        double vr, vi, tmp;
        unsigned long ki, mr, mi;
        ki = kr + 1;
        mr = kr + n;
        mi = mr + 1;
        vr = z[mr] * c - z[mi] * s;
        vi = z[mr] * s + z[mi] * c;
        z[mr] = z[kr] - vr;
        z[mi] = z[ki] - vi;
        z[kr] += vr;
        z[ki] += vi;
        tmp = c;
        c -= a * c + b * s;
        s -= a * s - b * tmp;
    }
    return;
} 

static void fht_dif_iter(double* x, unsigned long n)
{
    unsigned long m;

    for (m = n; m > 1; m >>= 1) {
        double a, b, c, s, t;
        unsigned long i, j, k, mh, mq;
        mh = m >> 1;
        mq = mh >> 1;
        t = PI / (double)mh;
        a = sin(0.5 * t);
        a *= 2.0 * a;
        b = sin(t);
        for (i = 0; i < n; i += m) {
            double* xp;
            xp = x + i;
            for (j = 0, k = mh; j < mh; ++j, ++k) {
                double u, v;
                u = xp[j];
                v = xp[k];
                xp[j] = u + v;
                xp[k] = u - v;
            }
        }
        c = 1.0;
        s = 0.0;
        for (j = 1, k = mh - 1; j < mq; ++j, --k) {
            double tmp;
            double* xj, * xk;
            xj = x + j + mh;
            xk = x + k + mh;
            tmp = c;
            c -= a * c + b * s;
            s -= a * s - b * tmp;
            for (i = 0; i < n; i += m) {
                double u, v;
                u = xj[i];
                v = xk[i];
                xj[i] = u * c + v * s;
                xk[i] = u * s - v * c;
            }
        }
    }
    return;
} 


static void fht_dif_iter_seq(double* x, unsigned long n)
{
    unsigned long m;

    for (m = n; m > 1; m >>= 1) {
        double a, b, t;
        unsigned long i, mh, mq;
        mh = m >> 1;
        mq = mh >> 1;
        t = PI / (double)mh;
        a = sin(0.5 * t);
        a *= 2.0 * a;
        b = sin(t);
        for (i = 0; i < n; i += m) {
            double c, s;
            double* xp;
            unsigned long j, k;
            xp = x + i;
            for (j = 0, k = mh; j < mh; ++j, ++k) {
                double u, v;
                u = xp[j];
                v = xp[k];
                xp[j] = u + v;
                xp[k] = u - v;
            }
            xp += mh;
            c = 1.0;
            s = 0.0;
            for (j = 1, k = mh - 1; j < mq; ++j, --k) {
                double u, v, tmp;
                tmp = c;
                c -= a * c + b * s;
                s -= a * s - b * tmp;
                u = xp[j];
                v = xp[k];
                xp[j] = u * c + v * s;
                xp[k] = u * s - v * c;
            }
        }
    }
    return;
}


static void fht_dif_rec(double* x, unsigned long n, int nbranch)
{
    double a, b, c, s, t;
    unsigned long j, jmax, k, nh, nq;

    if (n == 1)
        return;
    if (n <= (unsigned long)(L1_CACHE_BYTES / sizeof(double))) {
        if (FHT_UNIT_STRIDE)
            fht_dif_iter_seq(x, n);
        else
            fht_dif_iter(x, n);
        return;
    }
    nh = n >> 1;
    nq = nh >> 1;
    t = PI / (double)nh;
    a = sin(0.5 * t);
    a *= 2.0 * a;
    b = sin(t);
    for (j = 0, k = nh; j < nh; ++j, ++k) {
        double u, v;
        u = x[j];
        v = x[k];
        x[j] = u + v;
        x[k] = u - v;
    }
    c = 1.0;
    s = 0.0;
    jmax = nq + nh;
    for (j = nh + 1, k = n - 1; j < jmax; ++j, --k) {
        double u, v, tmp;
        tmp = c;
        c -= a * c + b * s;
        s -= a * s - b * tmp;
        u = x[j];
        v = x[k];
        x[j] = u * c + v * s;
        x[k] = u * s - v * c;
    }
    nbranch <<= 1;
#if (_OPENMP >= 200203)
#pragma omp parallel sections if (nbranch <= omp_get_max_threads()) num_threads(2)
#endif
    {
#if (_OPENMP >= 200203)
#pragma omp section
#endif
        fht_dif_rec(x, nh, nbranch);
#if (_OPENMP >= 200203)
#pragma omp section
#endif
        fht_dif_rec(x + nh, nh, nbranch);
    }
    return;
}


static void fht_dit_iter(double* x, unsigned long n)
{
    unsigned long m;

    for (m = 2; m <= n; m <<= 1) {
        double a, b, c, s, t;
        unsigned long i, j, k, mh, mq;
        mh = m >> 1;
        mq = mh >> 1;
        t = PI / (double)mh;
        a = sin(0.5 * t);
        a *= 2.0 * a;
        b = sin(t);
        c = 1.0;
        s = 0.0;
        for (j = 1, k = mh - 1; j < mq; ++j, --k) {
            double tmp;
            double* xj, * xk;
            xj = x + j + mh;
            xk = x + k + mh;
            tmp = c;
            c -= a * c + b * s;
            s -= a * s - b * tmp;
            for (i = 0; i < n; i += m) {
                double u, v;
                u = xj[i];
                v = xk[i];
                xj[i] = u * c + v * s;
                xk[i] = u * s - v * c;
            }
        }
        for (i = 0; i < n; i += m) {
            double* xp;
            xp = x + i;
            for (j = 0, k = mh; j < mh; ++j, ++k) {
                double u, v;
                u = xp[j];
                v = xp[k];
                xp[j] = u + v;
                xp[k] = u - v;
            }
        }
    }
    return;
}


static void fht_dit_iter_seq(double* x, unsigned long n)
{
    unsigned long m;

    for (m = 2; m <= n; m <<= 1) {
        double a, b, t;
        unsigned long i, mh, mq;
        mh = m >> 1;
        mq = mh >> 1;
        t = PI / (double)mh;
        a = sin(0.5 * t);
        a *= 2.0 * a;
        b = sin(t);
        for (i = 0; i < n; i += m) {
            double c, s;
            double* xp;
            unsigned long j, k;
            xp = x + i + mh;
            c = 1.0;
            s = 0.0;
            for (j = 1, k = mh - 1; j < mq; ++j, --k) {
                double tmp, u, v;
                tmp = c;
                c -= a * c + b * s;
                s -= a * s - b * tmp;
                u = xp[j];
                v = xp[k];
                xp[j] = u * c + v * s;
                xp[k] = u * s - v * c;
            }
            xp -= mh;
            for (j = 0, k = mh; j < mh; ++j, ++k) {
                double u, v;
                u = xp[j];
                v = xp[k];
                xp[j] = u + v;
                xp[k] = u - v;
            }
        }
    }
    return;
}


static void fht_dit_rec(double* x, unsigned long n, int nbranch)
{
    double a, b, c, s, t;
    unsigned long j, jmax, k, nh, nq;

    if (n == 1)
        return;
    if (n <= (unsigned long)(L1_CACHE_BYTES / sizeof(double))) {
        if (FHT_UNIT_STRIDE)
            fht_dit_iter_seq(x, n);
        else
            fht_dit_iter(x, n);
        return;
    }
    nh = n >> 1;
    nq = nh >> 1;
    nbranch <<= 1;
#if (_OPENMP >= 200203)
#pragma omp parallel sections if (nbranch <= omp_get_max_threads()) num_threads(2)
#endif
    {
#if (_OPENMP >= 200203)
#pragma omp section
#endif
        fht_dit_rec(x, nh, nbranch);
#if (_OPENMP >= 200203)
#pragma omp section
#endif
        fht_dit_rec(x + nh, nh, nbranch);
    }
    t = PI / (double)nh;
    a = sin(0.5 * t);
    a *= 2.0 * a;
    b = sin(t);
    jmax = nq + nh;
    c = 1.0;
    s = 0.0;
    for (j = nh + 1, k = n - 1; j < jmax; ++j, --k) {
        double tmp, u, v;
        tmp = c;
        c -= a * c + b * s;
        s -= a * s - b * tmp;
        u = x[j];
        v = x[k];
        x[j] = u * c + v * s;
        x[k] = u * s - v * c;
    }
    for (j = 0, k = nh; j < nh; ++j, ++k) {
        double u, v;
        u = x[j];
        v = x[k];
        x[j] = u + v;
        x[k] = u - v;
    }
    return;
} 


void bitrev_permute(double* z, unsigned long n)
{
    unsigned long i;
    unsigned int ldn = 0;
    unsigned int rshift;

    i = n;
    while (i >>= 1)
        ++ldn;
    rshift = 8 * (unsigned int)sizeof(unsigned long) - ldn;
    for (i = 0; i < n; ++i) {
        unsigned long r;
#if (ULONG_MAX == 0xffffffff) 
        r = ((i & 0x55555555) << 1) | ((i & ~0x55555555) >> 1);
        r = ((r & 0x33333333) << 2) | ((r & ~0x33333333) >> 2);
        r = ((r & 0x0f0f0f0f) << 4) | ((r & ~0x0f0f0f0f) >> 4);
        r = ((r & 0x00ff00ff) << 8) | ((r & ~0x00ff00ff) >> 8);
        r = (r << 16) | (r >> 16);
#elif (ULONG_MAX == 0xffffffffffffffff) 
        r = ((i & 0x5555555555555555) << 1) | ((i & ~0x5555555555555555) >> 1);
        r = ((r & 0x3333333333333333) << 2) | ((r & ~0x3333333333333333) >> 2);
        r = ((r & 0x0f0f0f0f0f0f0f0f) << 4) | ((r & ~0x0f0f0f0f0f0f0f0f) >> 4);
        r = ((r & 0x00ff00ff00ff00ff) << 8) | ((r & ~0x00ff00ff00ff00ff) >> 8);
        r = ((r & 0x0000ffff0000ffff) << 16) |
            ((r & ~0x0000ffff0000ffff) >> 16);
        r = (r << 32) | (r >> 32);
#endif
        r >>= rshift;
        if (r > i) {
            double tmp;
            unsigned long i2;
            i2 = i << 1;
            r <<= 1;
            tmp = z[i2]; z[i2] = z[r]; z[r] = tmp;
            tmp = z[i2 + 1]; z[i2 + 1] = z[r + 1]; z[r + 1] = tmp;
        }
    }
    return;
} 


void bitrev_permute_real(double* x, unsigned long n)
{
    unsigned long i;
    unsigned int ldn = 0;
    unsigned int rshift;

    i = n;
    while (i >>= 1)
        ++ldn;
    rshift = 8 * (unsigned int)sizeof(unsigned long) - ldn;
    for (i = 0; i < n; ++i) {
        unsigned long r;
#if (ULONG_MAX == 0xffffffff) 
        r = ((i & 0x55555555) << 1) | ((i & ~0x55555555) >> 1);
        r = ((r & 0x33333333) << 2) | ((r & ~0x33333333) >> 2);
        r = ((r & 0x0f0f0f0f) << 4) | ((r & ~0x0f0f0f0f) >> 4);
        r = ((r & 0x00ff00ff) << 8) | ((r & ~0x00ff00ff) >> 8);
        r = (r << 16) | (r >> 16);
#elif (ULONG_MAX == 0xffffffffffffffff) 
        r = ((i & 0x5555555555555555) << 1) | ((i & ~0x5555555555555555) >> 1);
        r = ((r & 0x3333333333333333) << 2) | ((r & ~0x3333333333333333) >> 2);
        r = ((r & 0x0f0f0f0f0f0f0f0f) << 4) | ((r & ~0x0f0f0f0f0f0f0f0f) >> 4);
        r = ((r & 0x00ff00ff00ff00ff) << 8) | ((r & ~0x00ff00ff00ff00ff) >> 8);
        r = ((r & 0x0000ffff0000ffff) << 16) |
            ((r & ~0x0000ffff0000ffff) >> 16);
        r = (r << 32) | (r >> 32);
#endif
        r >>= rshift;
        if (r > i) {
            double tmp;
            tmp = x[i]; x[i] = x[r]; x[r] = tmp;
        }
    }
    return;
}  


#if (L1_CACHE_BYTES == 0)
#define MEM_SMALL 0x8000
#else
#define MEM_SMALL L1_CACHE_BYTES
#endif

#if (L2_CACHE_BYTES == 0)
#define MEM_MEDIUM 0x100000
#else
#define MEM_MEDIUM L2_CACHE_BYTES
#endif

#define MEM_LARGE (4 * (MEM_MEDIUM))

#define T_ESTIMATE 0.1
#define T_RUN 0.5     

struct ft_bmark {
    const char* name; 
    void (*func)(double*, unsigned long); 
    size_t mem_size; 
    size_t obj_size; 
} ft_bmarks[] = {
    {"", NULL, 0, 0},
    {"hilbert()", hilbert, MEM_SMALL, 2 * sizeof(double)},
    {"hilbert()", hilbert, MEM_MEDIUM, 2 * sizeof(double)},
    {"hilbert()", hilbert, MEM_LARGE, 2 * sizeof(double)},
    {"", NULL, 0, 0},
    {"fft_dif()", fft_dif, MEM_SMALL, 2 * sizeof(double)},
    {"fft_dif_iter()", fft_dif_iter, MEM_SMALL, 2 * sizeof(double)},
    {"fft_dif_iter_seq()", fft_dif_iter_seq, MEM_SMALL, 2 * sizeof(double)},
    {"ifft_dit()", ifft_dit, MEM_SMALL, 2 * sizeof(double)},
    {"ifft_dit_iter()", ifft_dit_iter, MEM_SMALL, 2 * sizeof(double)},
    {"ifft_dit_iter_seq()", ifft_dit_iter_seq, MEM_SMALL, 2 * sizeof(double)},
    {"bitrev_permute()", bitrev_permute, MEM_SMALL, 2 * sizeof(double)},
    {"", NULL, 0, 0},
    {"fht_dif()", fht_dif, MEM_SMALL, sizeof(double)},
    {"fht_dif_iter()", fht_dif_iter, MEM_SMALL, sizeof(double)},
    {"fht_dif_iter_seq()", fht_dif_iter_seq, MEM_SMALL, sizeof(double)},
    {"fht_dit()", fht_dit, MEM_SMALL, sizeof(double)},
    {"fht_dit_iter()", fht_dit_iter, MEM_SMALL, sizeof(double)},
    {"fht_dit_iter_seq()", fht_dit_iter_seq, MEM_SMALL, sizeof(double)},
    {"bitrev_permute_real()", bitrev_permute_real, MEM_SMALL, sizeof(double)},
    {"", NULL, 0, 0},
    {"fft_dif()", fft_dif, MEM_MEDIUM, 2 * sizeof(double)},
    {"ifft_dit()", ifft_dit, MEM_MEDIUM, 2 * sizeof(double)},
    {"bitrev_permute()", bitrev_permute, MEM_MEDIUM, 2 * sizeof(double)},
    {"fht_dif()", fht_dif, MEM_MEDIUM, sizeof(double)},
    {"fht_dit()", fht_dit, MEM_MEDIUM, sizeof(double)},
    {"bitrev_permute_real()", bitrev_permute_real, MEM_MEDIUM, sizeof(double)},
    {"", NULL, 0, 0},
    {"fft_dif()", fft_dif, MEM_LARGE, 2 * sizeof(double)},
    {"ifft_dit()", ifft_dit, MEM_LARGE, 2 * sizeof(double)},
    {"bitrev_permute()", bitrev_permute, MEM_LARGE, 2 * sizeof(double)},
    {"fht_dif()", fht_dif, MEM_LARGE, sizeof(double)},
    {"fht_dit()", fht_dit, MEM_LARGE, sizeof(double)},
    {"bitrev_permute_real()", bitrev_permute_real, MEM_LARGE, sizeof(double)},
};

void ft_benchmarks(void)
{
    unsigned int i, j;
    for (i = 0; i < (sizeof(ft_bmarks) / sizeof(struct ft_bmark)); ++i) {
        unsigned int ncalls;
        double t, tstart;
        double* z;
        unsigned long nmax, ngrid;
        if (ft_bmarks[i].func == NULL) {
            printf("%s\n", ft_bmarks[i].name);
            continue;
        }
        /*
         * Find ngrid as the largest power of two such that ngrid objects
         * of size obj_size will fit into mem_size.
         */
        nmax = (unsigned long)(ft_bmarks[i].mem_size / ft_bmarks[i].obj_size);
        ngrid = 1;
        while (nmax >>= 1)
            ngrid <<= 1;
        if ((z = (double*)calloc((size_t)ngrid, ft_bmarks[i].obj_size))
            == NULL) {
            fprintf(stderr,
                "calloc() failed in ft_benchmarks() on test %d.\n", i);
            return;
        }
        /*
         * loop to estimate ncalls
         */
        tstart = am_timer(0.0);
        ncalls = 0;
        do {
            ++ncalls;
            ft_bmarks[i].func(z, ngrid);
        } while ((t = am_timer(tstart)) < T_ESTIMATE);
        ncalls *= (unsigned int)(T_RUN / t);
        ncalls = ncalls ? ncalls : 1;
        /*
         * measurement loop
         */
        tstart = am_timer(0.0);
        for (j = 0; j < ncalls; ++j)
            ft_bmarks[i].func(z, ngrid);
        t = am_timer(tstart);
        printf("%22s, n = %ld, %d call%s: %7.0f us per call\n",
            ft_bmarks[i].name,
            (long int)ngrid,
            ncalls,
            ncalls == 1 ? "" : "s",
            1.e6 * t / ncalls);
        free(z);
    }
    return;
} /* ft_benchmarks() */
