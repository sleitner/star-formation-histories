#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "hfiles/auxiliary.h"

unsigned long int rng_seed = 0L;
gsl_rng *my_random_generator;

void init_auxiliary() {
  //char filename[256];
	my_random_generator = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set( my_random_generator, rng_seed );
}

void save_auxiliary() {
}

double gsl_function_wrapper( double x, void *params ) {
        double (*f)(double) = (double (*)(double))params;
        return  f(x);
}

double integrate( double (*f)(double), double a, double b, double epsrel, double epsabs ) {
        double result, error;
        gsl_function F;
        gsl_integration_workspace *w;

	if ( a == b ) {
		return 0.0;
	}

	w = gsl_integration_workspace_alloc(1000);

        F.function = gsl_function_wrapper;
        F.params = f;

	gsl_integration_qag(&F, a, b, epsrel, epsabs, 1000, 6,
                w, &result, &error);
                                                                                                                                                            
        gsl_integration_workspace_free(w);
                                                                                                                                                            
        return result;
}

#define MAX_ITER 1000
                                                                                                                                                            
double root_finder( double (*f)(double), double a, double b, double epsrel, double epsabs ) {
        int status,i;
        double root;
        const gsl_root_fsolver_type *T;
        gsl_root_fsolver *s;
        gsl_function F;
                                                                                                                                                            
        F.function = gsl_function_wrapper;
        F.params = f;
                                                                                                                                                            
        T = gsl_root_fsolver_brent;
        s = gsl_root_fsolver_alloc (T);
        gsl_root_fsolver_set (s, &F, a, b);
                                                                                                                                                            
        for ( i = 0; i < MAX_ITER; i++ ) {
                status = gsl_root_fsolver_iterate (s);
                status = gsl_root_test_interval (gsl_root_fsolver_x_lower(s),
                                                 gsl_root_fsolver_x_upper(s),
                                                 epsrel, epsabs);

                if (status == GSL_SUCCESS) {
                        root = gsl_root_fsolver_root(s);
                        gsl_root_fsolver_free(s);
                        return root;
                }
        }

	my_error("Did not reach root after %u iterations!", MAX_ITER);
	return 0.0;
}

#undef MAX_ITER

int rev_compare_ints( const void *a, const void *b ) {
  return ( *(int *)b - *(int *)a );
}

int compare_ints( const void *a, const void *b ) {
        if ( *(int *)a == -1 ) {
                return 1;
        } else if ( *(int *)b == -1 ) {
                return -1;
        } else {
                return ( *(int *)a - *(int *)b );
        }
}

int compare_floats( const void *a, const void *b ) {
	if ( *(float *)a > *(float *)b ) {
		return 1;
	} else if ( *(float *)a < *(float *)b ) {
		return -1;
	} else {
		return 0;
	}
}

int nearest_int( double x ) {
	int ix;
	double frac;

	ix = (int)x;
	frac = x - (double)ix;

	if ( frac < 0.5 ) {
		return ix;
	} else {
		return ix+1;
	}
}

double my_rand_gauss(double sigma) {
  return gsl_ran_gaussian( my_random_generator,sigma );
}
double my_rand() {
  return gsl_rng_uniform( my_random_generator );
}

void* my_alloc_worker(size_t size, const char *file, int line)
{
  void *ptr;

  if(size > 0)
    {
      ptr = malloc( size );

      if(ptr == NULL)
        {
          my_error( "Failure allocating %d bytes in file %s, line %d", size, file, line );
        }

      return ptr;
    }
  else
    {
      return NULL;
    }
}


void my_free_worker(void *ptr, const char *file, int line)
{
  if(ptr != NULL)
    {
      free(ptr);
    }
}










char* replace2(char* source_str,char* search_str,char* replace_str)
{
  char *ostr, *nstr = NULL, *pdest = "";
  int length, nlen;
  unsigned int nstr_allocated;
  unsigned int ostr_allocated;

  if(!source_str || !search_str || !replace_str){
    printf("Not enough arguments\n");
    return NULL;
  }
  ostr_allocated = sizeof(char) * (strlen(source_str)+1);
  ostr = malloc( sizeof(char) * (strlen(source_str)+1));
  if(!ostr){
    printf("Insufficient memory available\n");
    return NULL;
  }
  strcpy(ostr, source_str);

  while(pdest)
    {
      pdest = strstr( ostr, search_str );
      length = (int)(pdest - ostr);

      if ( pdest != NULL )
	{
	  ostr[length]='\0';
	  nlen = strlen(ostr)+strlen(replace_str)+strlen( strchr(ostr,0)+strlen(search_str) )+1;
	  if( !nstr || /* _msize( nstr ) */ nstr_allocated < sizeof(char) * nlen){
	    nstr_allocated = sizeof(char) * nlen;
	    nstr = malloc( sizeof(char) * nlen );
	  }
	  if(!nstr){
	    printf("Insufficient memory available\n");
	    return NULL;
	  }

	  strcpy(nstr, ostr);
	  strcat(nstr, replace_str);
	  strcat(nstr, strchr(ostr,0)+strlen(search_str));

	  if( /* _msize(ostr) */ ostr_allocated < sizeof(char)*strlen(nstr)+1 ){
	    ostr_allocated = sizeof(char)*strlen(nstr)+1;
	    ostr = malloc(sizeof(char)*strlen(nstr)+1 );
	  }
	  if(!ostr){
	    printf("Insufficient memory available\n");
	    return NULL;
	  }
	  strcpy(ostr, nstr);
	}
    }
  if(nstr)
    free(nstr);
  return ostr;
}






/* http://www.daniweb.com/code/snippet216517.html
* Description:
* Find and replace text within a string.
*
* Parameters:
* src (in) - pointer to source string
* from (in) - pointer to search text
* to (in) - pointer to replacement text
*
* Returns:
* Returns a pointer to dynamically-allocated memory containing string
* with occurences of the text pointed to by 'from' replaced by with the
* text pointed to by 'to'.
*/
char *replace(const char *src, const char *from, const char *to)
{
  size_t size = strlen(src) + 1;
  size_t fromlen = strlen(from);
  size_t tolen = strlen(to);
  char *value = malloc(size);
  char *dst = value;
  if ( value != NULL )
    {
      for ( ;; )
	{
	  const char *match = strstr(src, from);
	  if ( match != NULL )
	    {
	      size_t count = match - src;
	      char *temp; 
	      size += tolen - fromlen;
	      temp = realloc(value, size);
	      if ( temp == NULL )
		{
		  free(value);
		  return NULL;
		}
	      dst = temp + (dst - value);
	      value = temp;
	      memmove(dst, src, count);
	      src += count;
	      dst += count;
	      memmove(dst, to, tolen);
	      src += fromlen;
	      dst += tolen;
	    }
	  else /* No match found. */
	    {
	      /* Copy any remaining part of the string. This includes the null
	       */
	      strcpy(dst, src);
	      break;
	    }
	}
    }
  return value;
}
      
      
      
      
      
      

void my_error( const char *fmt, ... ) {
	char message[256];

	va_list args;

	va_start( args, fmt );
	vsnprintf( message, 256, fmt, args );
	fprintf(stderr, "%s\n", message );
	fflush(stderr);
	va_end( args );
	exit(1);
}

void my_debug( const char *fmt, ... ) {
	char message[256];

	va_list args;

	va_start( args, fmt );
	vsnprintf( message, 256, fmt, args );
	fprintf( stdout, "%s\n", message );
	va_end(args);
}
