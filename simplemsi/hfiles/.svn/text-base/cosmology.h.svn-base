#ifndef __COSMOLOGY_H__
#define __COSMOLOGY_H__

/*
//  INTERNAL DECLARATIONS
*/

struct CosmologyParameters
{
  double OmegaM;
  double OmegaD;
  double OmegaB;
  double OmegaL;
  double OmegaK;
  double OmegaR;
  double h;
  double DeltaDC;
  int flat;
  double Omh2;
  double Obh2;
};
extern const struct CosmologyParameters *cosmology;


#define COSMOLOGY_DECLARE_PRIMARY_PARAMETER(name) \
void cosmology_set_##name(double value)

#define cosmology_set(name,value)	\
cosmology_set_##name(value)

COSMOLOGY_DECLARE_PRIMARY_PARAMETER(OmegaM);
COSMOLOGY_DECLARE_PRIMARY_PARAMETER(OmegaB);
COSMOLOGY_DECLARE_PRIMARY_PARAMETER(OmegaL);
COSMOLOGY_DECLARE_PRIMARY_PARAMETER(h);
COSMOLOGY_DECLARE_PRIMARY_PARAMETER(DeltaDC);


/*
//  USER INTEFACE:
//
//  To access parameters, use cosmology-><parameter>, for example
//    cosmology->OmegaM for the matter density parameter.
//
//  To set the parameter, call cosmology_set(<parameter>,<value>), for example
//    cosmology_set(OmegaM,0.3);
*/

/*
//  Copy constructor
*/
void cosmology_copy(const struct CosmologyParameters *c);

/*
//  Manual initialization. This does not need to be called, 
//  the initialization is done automatically on the first call
//  to a relevant fucntion.
*/
void cosmology_init();

/*
//  This is for backward compatibility. Because tcode is only defined 
//  up to a constant, an offset may exist between tcode used internally 
//  in the code and the one read from an (old) data file. Call this 
//  function with the values of abox and tcode to set the offset for tcode.
*/
void cosmology_insure_consistency(double abox, double tcode);


// #ifdef COSMOLOGY

/*
//  Set the range of global scale factors for thread-safe
//  calls to direct functions until the argument leaves the range.
*/
void cosmology_set_thread_safe_range(double amin, double amax);

/*
//  Direct functions take the global cosmological scale factor as the argument.
//  These functionsare are thread-safe if called with the argument in the
//  range set by a prior call to cosmology_set_thread_safe_range(...).
//  Calling with the argument outside that range is ok, but breaks
//  thread-safety assurance.
*/
double  aBox(double a);
double tCode(double a);
double tPhys(double a);
double dPlus(double a);

/*
//  Inverse conversions (converting other variables to the global 
//  scale factor aUni). These functions are NOT generally thread-safe.
*/
double inv_aBox(double abox);
double inv_tCode(double tcode);
double inv_tPhys(double tphys);
double inv_dPlus(double dplus);

/*
//  Conversion macros
*/
#define  abox_from_auni(a)   aBox(a)
#define tcode_from_auni(a)  tCode(a)
#define tphys_from_auni(a)  tPhys(a)
#define dplus_from_auni(a)  dPlus(a)

#define auni_from_abox(v)   inv_aBox(v)
#define auni_from_tcode(v)  inv_tCode(v)
#define auni_from_tphys(v)  inv_tPhys(v)
#define auni_from_dplus(v)  inv_dPlus(v)

#define abox_from_tcode(tcode)   aBox(inv_tCode(tcode))
#define tcode_from_abox(abox)    tCode(inv_aBox(abox))

#define tphys_from_abox(abox)    tPhys(inv_aBox(abox))
#define tphys_from_tcode(tcode)  tPhys(inv_tCode(tcode))
#define dplus_from_tcode(tcode)  dPlus(inv_tCode(tcode))

// #endif  /* COSMOLOGY */

#endif  /* __COSMOLOGY_H__ */

