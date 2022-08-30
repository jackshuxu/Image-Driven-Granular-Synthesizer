

#ifndef AM_AM_SYSDEP_H
#define AM_AM_SYSDEP_H

#ifndef L1_CACHE_BYTES
#define L1_CACHE_BYTES  0x8000
#endif

#ifndef L1_CACHE_WAYS
#define L1_CACHE_WAYS   8
#endif

#ifndef L2_CACHE_BYTES
#define L2_CACHE_BYTES  0x100000
#endif

void am_sleep(unsigned int);
double am_timer(double);

#endif /* AM_AM_SYSDEP_H */
