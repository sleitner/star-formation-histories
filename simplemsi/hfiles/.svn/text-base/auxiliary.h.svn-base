#ifndef __AUXILIARY_H__
#define __AUXILIARY_H__

#include <stdlib.h>
#include <stdarg.h>

#define mabs(A) ((A) < 0.0 ? -(A) : (A))
#define nDim		3
#define min(x,y)        ( ((x) < (y)) ? (x): (y) )
#define max(x,y)        ( ((x) > (y)) ? (x): (y) )
#define sign(x,y)       ( (y>=0) ? fabs(x) : -fabs(x) )
#define sqr(x)          ( (x)*(x) ) 
#define cube(x)         ( (x)*(x)*(x) ) 

void init_auxiliary();
void save_auxiliary();

extern unsigned long int rng_seed;

double integrate( double (*f)(double), double a, double b, double epsrel, double epsabs );
double root_finder( double (*f)(double), double a, double b, double epsrel, double epsabs );
int compare_ints( const void *a, const void *b );
int rev_compare_ints( const void *a, const void *b );
int compare_floats( const void *a, const void *b );
int nearest_int( double x );

double my_rand();
double my_rand_gauss(double sigma);

#define my_alloc_bytes(size)     my_alloc_worker(size,__FILE__,__LINE__)
#define my_alloc(type,size)      (type *)my_alloc_worker((size)*sizeof(type),__FILE__,__LINE__)  
#define my_free(ptr) { my_free_worker(ptr,__FILE__,__LINE__); ptr = NULL; }

void *my_alloc_worker(size_t size, const char *file, int line);
void my_free_worker(void *ptr, const char *file, int line);

//-----------------------------------

char *replace(const char *src, const char *from, const char *to);
void my_error( const char *fmt, ... );

void my_debug( const char *fmt, ... );
#define my_assert( x ) if (!(x)) { my_error( "Assertion (%s) failed: %s line %u", #x, __FILE__, __LINE__ ); }

#endif
