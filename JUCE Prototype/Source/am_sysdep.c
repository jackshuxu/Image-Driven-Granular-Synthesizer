
#include <stdlib.h>
#include <time.h>

#if defined (__unix) || defined (__unix__)
#include <unistd.h>
#define UNISTD_INCLUDED
#elif defined _WIN32
#include <windows.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "am_sysdep.h"

void am_sleep(unsigned int seconds)
{
#if defined UNISTD_INCLUDED
    sleep(seconds);
    return;
#elif defined _WIN32
    Sleep((DWORD)(seconds * 1000));
    return;
#else
    clock_t t;

    t = clock() + (clock_t)(seconds * CLOCKS_PER_SEC);
    while (clock() < t)
        ;
    return;
#endif
}   /* am_sleep() */

double am_timer(double tref)
{
    double t;
#ifdef _OPENMP
#if (_OPENMP >= 200203)
    t = omp_get_wtime();
#else
    t = (double)clock() / (double)CLOCKS_PER_SEC;
#endif
#else
    t = (double)clock() / (double)CLOCKS_PER_SEC;
#endif
    return t - tref;
}   /* am_timer() */
