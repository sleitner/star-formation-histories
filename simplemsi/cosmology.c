#include "hfiles/auxiliary.h"
#define ASSERT(exp)      my_assert(exp)
#define NEWARR(size)     my_alloc(double,size)
#define DELETE(ptr)      my_free(ptr)

/* ------------------------------------ */

#include <math.h>
#include <string.h>


#ifndef ASSERT
#include <stdio.h>
#define ASSERT(exp) { if(!(exp)) { fprintf(stderr,"Failed assertion %s, line: %d\n",#exp,__LINE__); exit(1); } }
#endif

#ifndef NEWARR
#include <stdlib.h>
#define NEWARR(size)   (double *)malloc((size)*sizeof(double)
#endif

#ifndef DELETE
#include <stdlib.h>
#define DELETE(ptr)    free(ptr)
#endif


#include "hfiles/cosmology.h"


struct CosmologyParameters cosmology_internal_parameters = { 0.28, 0.234, 0.046, 0.72, 0.0, 0.0, 0.7, 0.0, 1 };
const struct CosmologyParameters *cosmology = &cosmology_internal_parameters;


struct CosmologyInternal
{
  int ndex;
  int size;
  double *la;
  double *aUni;
  double *aBox;
  double *tCode;
  double *tPhys;
  double *dPlus;
  double *qPlus;
  double aLow;
  double tCodeOffset;
} cosmology_internal_data = { 200, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 1.0e-2, 0.0 };

#define c cosmology_internal_parameters
#define d cosmology_internal_data


void cosmology_clear_table();
void cosmology_fill_table(double amin, double amax);
void cosmology_fill_table_abox(int istart, int n);


void cosmology_copy(const struct CosmologyParameters *ptr)
{
  c = *ptr;
  cosmology_clear_table();
}


void cosmology_set_OmegaM(double v)
{
  if(v < 1.0e-3) v = 1.0e-3;
  if(fabs(c.OmegaM-v) > 1.0e-5)
    {
      c.OmegaM = v;
      if(c.flat)
	{
	  c.OmegaL = 1.0 - c.OmegaM;
	}
      cosmology_clear_table();
    }
}


void cosmology_set_OmegaL(double v)
{
  if(fabs(c.OmegaL-v) > 1.0e-5)
    {
      c.OmegaL = v;
      if(c.OmegaL<0.0 || c.OmegaL>0.999) c.flat = 0;
      if(c.flat)
	{
	  c.OmegaM = 1.0 - c.OmegaL;
	}
      cosmology_clear_table();
    }
}


void cosmology_set_OmegaB(double v)
{
  if(v < 0.0) v = 0.0;
  if(v > c.OmegaM) v = c.OmegaM;
  if(fabs(c.OmegaB-v) > 1.0e-5)
    {
      c.OmegaB = v;
      cosmology_clear_table();
    }
}


void cosmology_set_h(double v)
{
  if(fabs(c.h-v) > 1.0e-5)
    {
      c.h = v;
      cosmology_clear_table();
    }
}


void cosmology_set_flat(int flat)
{
  if(!flat) flat = 1;
  if(c.flat != flat)
    {
      c.flat = flat;
      if(c.flat && fabs(c.OmegaM+c.OmegaL-1.0)>1.0e-5)
	{
	  c.OmegaL = 1.0 - c.OmegaM;
	  cosmology_clear_table();
	}
    }
}


void cosmology_set_DeltaDC(double v)
{
  if(fabs(v) < 1.0e-3) v = 0.0;
  if(fabs(c.DeltaDC-v) > 1.0e-3)
    {
      c.DeltaDC = v;
      cosmology_clear_table();
    }
}


void cosmology_init()
{
  if(d.size == 0) /* reset only if the state is dirty */
    {
      c.OmegaD = c.OmegaM - c.OmegaB;
      if(c.flat)
	c.OmegaK = 0.0;
      else
	c.OmegaK = 1.0 - (c.OmegaM+c.OmegaL);
      c.OmegaR = 5.837e-5/(c.h*c.h);

      c.Omh2 = c.OmegaM*c.h*c.h;
      c.Obh2 = c.OmegaB*c.h*c.h;

      cosmology_fill_table(d.aLow,1.0);
    }      
}


void cosmology_insure_consistency(double abox, double tcode)
{
  d.tCodeOffset = 0.0;
  d.tCodeOffset = tcode - tCode(inv_aBox(abox));
}


double cosmology_mu(double a)
{
  return sqrt(((a*a*c.OmegaL+c.OmegaK)*a+c.OmegaM)*a+c.OmegaR);
}


double cosmology_dc_factor(double dPlus)
{
  double dc = 1.0 + dPlus*c.DeltaDC;
  return 1.0/pow((dc>0.001)?dc:0.001,1.0/3.0);
}


void cosmology_fill_table_integrate(double a, double y[], double f[])
{
  double mu = cosmology_mu(a);
  double abox = a*cosmology_dc_factor(y[2]);
  
  f[0] = a/(abox*abox*mu);
  f[1] = a/mu;
  f[2] = y[3]/(a*mu);
  f[3] = 1.5*c.OmegaM*y[2]/mu;
}


void cosmology_fill_table_piece(int istart, int n)
{
  int i, j;
  double tPhysUnit = (3.0856775813e17/(365.25*86400))/c.h;  /* 1/H0 in Julian years */

  double x, aeq = c.OmegaR/c.OmegaM;
  double tCodeFac = 1.0/sqrt(aeq);
  double tPhysFac = tPhysUnit*aeq*sqrt(aeq)/sqrt(c.OmegaM);

  double da, a0, y0[4], y1[4];
  double f1[4], f2[4], f3[4], f4[4];


  for(i=istart; i<n; i++)
    {
      d.aUni[i] = pow(10.0,d.la[i]);
    }

  /*
  //  Small a regime, use analytical formulae for matter + radiation model
  */  
  for(i=istart; d.aUni[i]<(d.aLow+1.0e-9) && i<n; i++)
    {
      x = d.aUni[i]/aeq;

      d.tPhys[i] = tPhysFac*2*x*x*(2+sqrt(x+1))/(3*pow(1+sqrt(x+1),2.0));
      d.dPlus[i] = aeq*(x + 2.0/3.0 + (6*sqrt(1+x)+(2+3*x)*log(x)-2*(2+3*x)*log(1+sqrt(1+x)))/(log(64.0)-9));  /* long last term is the decaying mode generated after euality; it is very small for x > 10, I keep ot just for completeness; */
      d.qPlus[i] = d.aUni[i]*cosmology_mu(d.aUni[i])*(1 + ((2+6*x)/(x*sqrt(1+x))+3*log(x)-6*log(1+sqrt(1+x)))/(log(64)-9)); /* this is a*mu*dDPlus/dt/H0 */

      d.aBox[i] = d.aUni[i]*cosmology_dc_factor(d.dPlus[i]);
      d.tCode[i] = 1.0 - tCodeFac*asinh(sqrt(aeq/d.aBox[i]));
    }
  
  /*
  //  Large a regime, solve ODEs
  */
  ASSERT(i > 0);

  tCodeFac = 0.5*sqrt(c.OmegaM);
  tPhysFac = tPhysUnit;

  y1[0] = d.tCode[i-1]/tCodeFac;
  y1[1] = d.tPhys[i-1]/tPhysFac;
  y1[2] = d.dPlus[i-1];
  y1[3] = d.qPlus[i-1];

  for(; i<n; i++)
    {
      a0 = d.aUni[i-1];
      da = d.aUni[i] - a0;

      /*  RK4 integration */
      for(j=0; j<4; j++) y0[j] = y1[j];
      cosmology_fill_table_integrate(a0,y1,f1);

      for(j=0; j<4; j++) y1[j] = y0[j] + 0.5*da*f1[j];
      cosmology_fill_table_integrate(a0+0.5*da,y1,f2);

      for(j=0; j<4; j++) y1[j] = y0[j] + 0.5*da*f2[j];
      cosmology_fill_table_integrate(a0+0.5*da,y1,f3);

      for(j=0; j<4; j++) y1[j] = y0[j] + da*f3[j];
      cosmology_fill_table_integrate(a0+da,y1,f4);

      for(j=0; j<4; j++) y1[j] = y0[j] + da*(f1[j]+2*f2[j]+2*f3[j]+f4[j])/6.0;

      d.tCode[i] = tCodeFac*y1[0];
      d.tPhys[i] = tPhysFac*y1[1];
      d.dPlus[i] = y1[2];
      d.qPlus[i] = y1[3];

      d.aBox[i] = d.aUni[i]*cosmology_dc_factor(d.dPlus[i]);
    }
} 


void cosmology_fill_table(double amin, double amax)
{
  int i, imin, imax, iold;
  double dla = 1.0/d.ndex;
  double lamin, lamax;
  double *old_la = d.la;
  double *old_aUni = d.aUni;
  double *old_aBox = d.aBox;
  double *old_tCode = d.tCode;
  double *old_tPhys = d.tPhys;
  double *old_dPlus = d.dPlus;
  double *old_qPlus = d.qPlus;
  int old_size = d.size;

  if(amin > d.aLow) amin = d.aLow;
  lamin = dla*floor(d.ndex*log10(amin));
  lamax = dla*ceil(d.ndex*log10(amax));

  d.size = 1 + (int)(0.5+d.ndex*(lamax-lamin)); 
  ASSERT(fabs(lamax-lamin-dla*(d.size-1)) < 1.0e-14);

  d.la = NEWARR(d.size);     ASSERT(d.la != NULL);
  d.aUni = NEWARR(d.size);   ASSERT(d.aUni != NULL);
  d.aBox = NEWARR(d.size);   ASSERT(d.aBox != NULL);
  d.tCode = NEWARR(d.size);  ASSERT(d.tCode != NULL);
  d.tPhys = NEWARR(d.size);  ASSERT(d.tPhys != NULL);
  d.dPlus = NEWARR(d.size);  ASSERT(d.dPlus != NULL);
  d.qPlus = NEWARR(d.size);  ASSERT(d.qPlus != NULL);

  /*
  //  New log10(aUni) table
  */
  for(i=0; i<d.size; i++)
    {
      d.la[i] = lamin + dla*i;
    }

  if(old_size == 0)
    {
      /*
      //  Filling the table for the first time
      */
      cosmology_fill_table_piece(0,d.size);
    }
  else
    {
      /*
      //  Find if we need to expand the lower end
      */
      if(lamin < old_la[0])
	{
	  imin = (int)(0.5+d.ndex*(old_la[0]-lamin));
	  ASSERT(fabs(old_la[0]-lamin-dla*imin) < 1.0e-14);
	}
      else imin = 0;

      /*
      //  Find if we need to expand the upper end
      */
      if(lamax > old_la[old_size-1])
	{
	  imax = (int)(0.5+d.ndex*(old_la[old_size-1]-lamin));
	  ASSERT(fabs(old_la[old_size-1]-lamin-dla*imax) < 1.0e-14);
	}
      else imax = d.size - 1;
  
      /*
      //  Re-use the rest
      */
      if(lamin > old_la[0])
	{
	  iold = (int)(0.5+d.ndex*(lamin-old_la[0]));
	  ASSERT(fabs(lamin-old_la[0]-dla*iold) < 1.0e-14);
	}
      else iold = 0;

      memcpy(d.aUni+imin,old_aUni+iold,sizeof(double)*(imax-imin+1));
      memcpy(d.aBox+imin,old_aBox+iold,sizeof(double)*(imax-imin+1));
      memcpy(d.tCode+imin,old_tCode+iold,sizeof(double)*(imax-imin+1));
      memcpy(d.tPhys+imin,old_tPhys+iold,sizeof(double)*(imax-imin+1));
      memcpy(d.dPlus+imin,old_dPlus+iold,sizeof(double)*(imax-imin+1));
      memcpy(d.qPlus+imin,old_qPlus+iold,sizeof(double)*(imax-imin+1));

      DELETE(old_la);
      DELETE(old_aUni);
      DELETE(old_aBox);
      DELETE(old_tCode);
      DELETE(old_tPhys);
      DELETE(old_dPlus);
      DELETE(old_qPlus);

      /*
      //  Fill in additional pieces
      */
      if(imin > 0) cosmology_fill_table_piece(0,imin);
      if(imax < d.size-1) cosmology_fill_table_piece(imax,d.size);
    }
}


void cosmology_clear_table()
{
  if(d.size > 0)
    {
      DELETE(d.la);
      DELETE(d.aUni);
      DELETE(d.aBox);
      DELETE(d.tCode);
      DELETE(d.tPhys);
      DELETE(d.dPlus);
      DELETE(d.qPlus);

      d.size = 0;
      d.la = NULL;
      d.aUni = NULL;
      d.aBox = NULL;
      d.tCode = NULL;
      d.tPhys = NULL;
      d.dPlus = NULL;
      d.qPlus = NULL;
    }
}


void cosmology_check_range(double a)
{
  ASSERT((a > 1.0e-9) && (a < 1.0e9));

  if(d.size == 0) cosmology_init();

  if(a < d.aUni[0])
    {
      cosmology_fill_table(a,d.aUni[d.size-1]);
    }

  if(a > d.aUni[d.size-1])
    {
      cosmology_fill_table(d.aUni[0],a);
    }
}


void cosmology_set_thread_safe_range(double amin, double amax)
{
  cosmology_check_range(amin);
  cosmology_check_range(amax);
}


double cosmology_get_value_from_table(double a, double table[])
{
  int idx = (int)(d.ndex*(log10(a)-d.la[0]));

  ASSERT(idx>=0 && idx<d.size);

  /*
  //  Do it as a function of aUni rather than la to ensure exact inversion
  */
  return table[idx] + (table[idx+1]-table[idx])/(d.aUni[idx+1]-d.aUni[idx])*(a-d.aUni[idx]);
}


int cosmology_find_index(double v, double table[])
{
  int ic, il = 0;
  int ih = d.size - 1;

  if(v < table[0])
    {
      return -1;
    }
  if(v > table[d.size-1])
    {
      return d.size + 1;
    }

  while((ih-il) > 1)
    {
      ic = (il+ih)/2;
      if(v > table[ic]) /* special, not fully optimal form to avoid checking that il < d.size-1 */
	il = ic;
      else
	ih = ic;
    }

  ASSERT(il+1 < d.size);

  return il;
}




/*
//  Direct and inverse functions
*/
#define DEFINE_FUN(name,offset)			\
double name(double a) \
{ \
  cosmology_check_range(a); \
  return cosmology_get_value_from_table(a,d.name) + offset; \
} \
double inv_##name(double v) \
{ \
  int idx; \
  double *table; \
  v -= offset; \
  if(d.size == 0) cosmology_init(); \
  table = d.name; \
  idx = cosmology_find_index(v,table); \
  while(idx < 0) \
    { \
      cosmology_check_range(0.5*d.aUni[0]); \
      table = d.name; \
      idx = cosmology_find_index(v,table); \
    } \
  while(idx > d.size) \
    { \
      cosmology_check_range(2.0*d.aUni[d.size-1]); \
      table = d.name; \
      idx = cosmology_find_index(v,table); \
    } \
  return d.aUni[idx] + (d.aUni[idx+1]-d.aUni[idx])/(table[idx+1]-table[idx])*(v-table[idx]); \
}

DEFINE_FUN(aBox,0.0);
DEFINE_FUN(tCode,d.tCodeOffset);
DEFINE_FUN(tPhys,0.0);
DEFINE_FUN(dPlus,0.0);

