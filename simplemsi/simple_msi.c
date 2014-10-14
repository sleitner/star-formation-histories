#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>

#include "hfiles/defs.h"
#include "hfiles/auxiliary.h"
#include "hfiles/cosmology.h"

typedef struct {
    int idescendent;
    double zform;
    double mu;
} merger_tree_struct;

typedef struct {
    int ntree;
    merger_tree_struct *list;
} merger_list;

#define ngrid_z  30
#define ngrid_lM 16
#define ngrid_mu 99
#define ntrials 1  //ifdef MERGING
double nu[ngrid_z][ngrid_lM][ngrid_mu];
double dfml(double dt, double tij);
//double cumfml(double t);
void pickIMFname(char fileout[]);
double psi_field(double mass, double z, char field_type[3], int iopt, int ierr);
#define Ntable 100
double lin_interpolate(int ind_low, double x, double x_arr[], double y_arr[]  );
double interpolate(double x, double x_arr[], int nels, double y_arr[], int iopt );

const int npossible_progeny=1;  // simple doesn't include merging




/*
 *Mobs_norm=scaling mass normalization
 *t0_Mstar=galaxy mass at initial time
 *t0_sfr=current SF
 *stoch_sigma = the sigma in the stochastic scatter
 *normA0_1sigma=scatter multiplying the observed_median_psi  
 *normA0= factor adjusting psi (mod_in_sigma_normA0-1)*normA0_1sigma at the given t0_Mstar
 *alpha/beta1=scaling slopes
 *zbreak=fit break
 *tau_delta=smoothing timescale, stoch_sigma=0.3
 *psi_field = SFR field
 *phi = galaxy's SFR
 */

double t0_Mstar; //= 5.62341E10; //can be altered by galaxy_name args or -opt
double grand_systemic[ntrials];

const double normA0_1sigma = 0.3;
int mod_in_sigma_normA0;
const double alpha = 3.36; //3.36//2.86 //3.86
const double beta1 = 0.85; //0.85//0.75//0.5,1.00
double psi_beta;
double normA0, t0_Mstar, t0_obs_shift=0;

double tage_from_tphys(double tphys){ return tphys_from_auni(1.0)-tphys-t0_obs_shift; }

double tlb_from_tphys(double tphys){ return tphys_from_auni(1.0)-tphys; }
double tphys_from_tlb(double tlb){ return tphys_from_auni(1.0)-tlb; }
double z_form_Noeske(double Mbar){   return pow(10,-2.7)*pow(Mbar,0.3)-1; }
double tlb_form_Noeske(double Mbar){
    double aform = 1/(1+z_form_Noeske(Mbar));
    double tform = tphys_from_auni( aform );
    return  tlb_from_tphys(tform);
}
double tau_Noeske(double Mbar){ return pow(10,20.7)/Mbar;}
double Mstar_Noeske(double tlb, double Mbar){
    if(tlb<tlb_form_Noeske(Mbar)){
        return Mbar*(1 - exp( -(tlb_form_Noeske(Mbar)-tlb)/tau_Noeske(Mbar) ));  
    }else{
        return 0;
    }
}
double phi_Noeske(double tlb, double Mbar){
    double  phi0;
    if(tlb<tlb_form_Noeske(Mbar)){
        if(c_IRR!=0.0){
            phi0 = Mbar/tau_Noeske(Mbar)/(1-c_IRR);
        }else{
            phi0 = Mbar/tau_Noeske(Mbar)/(1-0.5);
        }
        return phi0*exp(-(tlb_form_Noeske(Mbar)-tlb)/tau_Noeske(Mbar));
    }else{
        return 0;
    }
}

double Mbar_Noeske_from_Mgal(double Mbar){
    double y;
    y = t0_Mstar - Mstar_Noeske(t0_obs_shift, Mbar);
    return y;
}



main(int argc,char **argv)
{
    double t0_sfr, observed_median_psi;
    const double stoch_sigma=0.0;
    char cstoch_sigma[10] ;
    double grandfact ;
    int mod_in_sigma_normA0 = 1; //mod_in_sigma_normA0=0 = -mod_in_sigma_normA0, 1= 0, 2=+mod_in_sigma_normA0      .......... 10-> manual input below
    double userdef_normA0 = 9.9;     //can be altered by galaxy_name args
    double userdef_t0_sfr = 9.9;     //can be altered by galaxy_name args

    double zfp10 = -1;
    double zfp15 = -1;
    double zfp20 = -1;
    double zfp25= -1;
    double zfp50 = -1;
    double mfz1 = -1;
    double tfp10,tfp25,tfp50;
    double tfp15,tfp20;
    double fracM;

  double G0,mltoi,count_time,tot_count_time;
  double Mbar_Noeske;
  double tlb_first, tfirst,dt, tlb;
  double tage, ln_tage;
  int i,j,k, iter;
  FILE *fp,*fpt;
  char ct0_obs_shift[20], cbeta1[20], calpha[20], clt0_Mstar[20],cmod_in_sigma_normA0[20], fileoutform[256],fileoutbeta[256],fileoutsfh[256],file_label[256], cseed[10], cmod_in_sigma_t0_sfr[10], ctau[10], galaxy_name[10], cdummy[10], fileoutmerge[256];
  char field_type[3];
  int igalaxy_name=0;
  int fiarg=1;
  int ierr=0, iopt=-100;
  double psi0;

///////////////////////////// main quantities
#define ntimebins 1000
  static int ilast=ntimebins-1;
  double time_binl[ntimebins], red_binl[ntimebins];
  double Mstar_dot[ntimebins][npossible_progeny];
  double dMLdt[ntimebins][npossible_progeny];
  double massloss_to_t0[ntimebins];  
  double Mstar[ntimebins][npossible_progeny],phi[ntimebins][npossible_progeny]; 
  
///////////////merging stuff
  int idivide[ntimebins][npossible_progeny];
  merger_list *mtree;
  double Mstar_m0[ntimebins] ,phi_m0[ntimebins]; 
  double Mstar_dot_m0[ntimebins],dMLdt_m0[ntimebins];
  int iprog,ipar, nprogeny, imu, inewprog;
  double mu, thalf, prob_merger_dt, Mdot;
  int count_major, count_minor;
  double mtot;

/////////////averaging
  int itrial;
  double cum_Mstar[ntimebins] ,cum_phi[ntimebins]; 
  
///////////////////////////// convolution stuff
  double sum_overflow = 0;
  double sum_last_vespa_bin = 0;
  double tage_last_vespa,tlb_last_vespa;
  double tage_overflow,tlb_overflow;
  
  double Mstar_smoothed[ntimebins], Mstar_dot_smoothed[ntimebins],phi_smoothed_linear[ntimebins];
  double   phi_inp_smoothed1[ntimebins], phi_inp_smoothed5[ntimebins] , phi_inp_smoothed75[ntimebins] ;
  double   CSFH_inp_smoothed1[ntimebins], CSFH_inp_smoothed5[ntimebins] , CSFH_inp_smoothed75[ntimebins], CSFH_m0[ntimebins] ;
  double   Mstar_dot_inp_smoothed1[ntimebins], Mstar_dot_inp_smoothed5[ntimebins] , Mstar_dot_inp_smoothed75[ntimebins] ;
  double   Mstar_inp_smoothed1[ntimebins], Mstar_inp_smoothed5[ntimebins] , Mstar_inp_smoothed75[ntimebins] ;
  double tlb1=0; double tlb2=0; double P_LN; 
  const int resampleCount=1e3;
  double phi_ln;
  int tbin;
  double  resample_bin[resampleCount];  
  double phi_smoothed_ln[resampleCount];
  double resample_binfin; 
  double resample_bin0; 
  double resample_dbin;
  double smoothed_total_mass;

  double dlnt;

  int ifwhm; double sig_err, sig_err10;
  const int nfwhm=3;
  const double fwhm10[]={0.50,0.75,1.00}; 
  double dst, tform, tsample;
/////////////////////////////

  if(argc !=6){
      my_debug("usage:exec <log10(t0_Mstar[Msun])> <shift in lookback time from 0[Gyr]> <fit_to_sfr-m*(ofp[O10pow],ofe[O10powexp],oda[O10interp],kfp,kfe,kda)>");
      my_debug("<extrapolation interp (0=11[fixed at last data],1=01[extrapolated at lowz],2=10[extrapolated at hiz],3=00[extrapolated])>");
      my_debug("<change in SFR-M* normalization (1=default 0/2= change by +-1sigma>");
      my_debug(" -opt can also be used to specify galaxies+SFR directly in the code of simple_derive_mstar.c");
      my_error("exiting");
  }
  t0_Mstar = pow(10.,atof(argv[fiarg++])); 
  t0_obs_shift = atof(argv[fiarg++]); 
  t0_obs_shift *=1e9;
  strcpy( field_type,argv[fiarg++] );
  iopt = atoi(argv[fiarg++]); 
  mod_in_sigma_normA0 = atoi(argv[fiarg++]);
  my_debug("\n options t0_Mstar=%f t0_obs_shift[Gyr]=%f field_type=%s iopt=%d AdjnormA0=%d \n",log10(t0_Mstar),t0_obs_shift/1e9, field_type, iopt,mod_in_sigma_normA0);

  if(t0_Mstar<1e5||t0_Mstar>1e12||mod_in_sigma_normA0>10||mod_in_sigma_normA0<0||t0_obs_shift>3.0e9){my_error("weird options ... double check");}
/*   if(NITER > 1 && c_IRR != 0.00){  */
/*       my_error("dont iterate and use instant recycling!"); */
/*   } */
  if(NITER == 1 && c_IRR == 0.00){
      my_error("if NITER==1 use instant recycling!");
  }
 
  if(argc>fiarg){
      if(strcmp(argv[fiarg],"-opt")!=0){
          if(argc > 1 && mod_in_sigma_normA0 == 10 ){
              strcpy( galaxy_name,argv[fiarg++] );
          }else{
              if(argc > 1){
                  my_error("scatter is off so set mod_in_sigma_normA0=10 and supply a galaxy_name, OR supply options OR let it run with default set params");
              }else{
                  i=sprintf(galaxy_name,"" );
              }
          }
      }
  }else{
      i=sprintf(galaxy_name,"" );
  }

  if(strcmp(galaxy_name,"M33")==0){
    t0_Mstar = 4.5e9;
    userdef_t0_sfr= 0.7;
  }else if(strcmp(galaxy_name,"n2403")== 0){
    t0_Mstar = 14e9;
    userdef_t0_sfr= 1.3;
  }else if(strcmp(galaxy_name,"n2997")== 0){
    t0_Mstar = 160e9;
    userdef_t0_sfr = 5.0;
  }else if(strcmp(galaxy_name,"n4559")== 0){
    t0_Mstar = 6.8e9;
    userdef_t0_sfr= 0.6;
  }else if(strcmp(galaxy_name,"MW")== 0){
    t0_Mstar = 55e9;
    userdef_t0_sfr= 1.45;
  }else if(strcmp(galaxy_name,"M31")== 0){
    t0_Mstar = 100e9;
    userdef_t0_sfr= 1.0;
  }else if(strcmp(galaxy_name,"n891")== 0){
    t0_Mstar = 100e9;
    userdef_t0_sfr= 3.8;
  }else if(strcmp(galaxy_name,"n5746")== 0){
    t0_Mstar = 160e9;
    userdef_t0_sfr= 1.2;
  }else if(strcmp(galaxy_name,"")== 0){
  }else{my_error("bad galaxy name %s",galaxy_name);}
  userdef_normA0=userdef_t0_sfr;

  if(strcmp(galaxy_name,"")!= 0){
   my_debug("arg %s t0_Mstar=%e normA0=%f\n",galaxy_name,t0_Mstar,userdef_t0_sfr);
   igalaxy_name=1;
  }

  if(userdef_normA0==9.90 && mod_in_sigma_normA0==10){my_error("bad userdef_normA0=9.90");}
  //================== done choosing initialization parameters


  //---------------- PICK field_type
  my_assert(   strcmp( field_type, "ofp" ) == 0 || 
               strcmp( field_type, "ofe" ) == 0 ||
               strcmp( field_type, "kfe" ) == 0 ||
               strcmp( field_type, "kfp" ) == 0 ||
               strcmp( field_type, "kda" ) == 0 ||
               strcmp( field_type, "oda" ) == 0 );
  if( psi_field(-1,-1,field_type,-1,0) != 1){my_error("failed in psi init");}
  
  //================== COSMOLOGY
  double auni, abox;
  double age_of_universe;
  init_auxiliary();

  auni=1.0;
  cosmology_set(OmegaM,0.2581);
  cosmology_set(OmegaB,0.0441);
  cosmology_set(OmegaL,.7419);
  cosmology_set(h,0.719);
  cosmology_set(DeltaDC,0.0);
  cosmology_init();
  cosmology_insure_consistency(1.0,0.0);

  age_of_universe = tphys_from_auni(1.0);
  my_debug("age of universe=%e",age_of_universe);
  //==================

  tfirst = tphys_from_auni(0.1);
  tlb_first = tlb_from_tphys( tfirst );
  dt = ( (age_of_universe-t0_obs_shift) - tfirst )/(1.0*ntimebins);
  resample_binfin=log(1e12); 
//  resample_bin0=log((t0_obs_shift+dt*0.01));
  resample_bin0=log(dt*0.1);
  resample_dbin=(resample_binfin-resample_bin0)/(resampleCount-1.0);
  
  tlb_overflow = exp( log( tlb_first ) ) ;
  tage_overflow = exp( log( tage_from_tphys(tfirst) ) ) ;
  
  tlb_last_vespa=10e9;
  tage_last_vespa=tage_from_tphys(tphys_from_tlb(tlb_last_vespa));
  
  my_debug("dt=%e %e %e",dt,(age_of_universe-t0_obs_shift*1e9),tfirst);
  for( i = 0; i < ntimebins ; i++){
      time_binl[i] = tfirst + dt*(i); //first to now
      red_binl[i] = 1/auni_from_tphys(time_binl[i])-1;
  }

  //================== COSMOLOGY

  

  //---------------- PICK t0_sfr and normA0
  psi0 = psi_field(t0_Mstar, red_binl[ilast], field_type, iopt, ierr);

  if(mod_in_sigma_normA0==10){
      normA0 = userdef_normA0/psi0;
  }else{
      normA0 = pow(10,(mod_in_sigma_normA0-1)*normA0_1sigma); 
  }
  t0_sfr = normA0*psi0;   


//-------------------------------------
  
  
  
  //------------------------------------- filename
  
    strcpy(file_label, ""); 

    if(igalaxy_name == 0){ //galaxy parameters
        i=sprintf(clt0_Mstar,"M%2.2f\n",log10(t0_Mstar));
    }else{
        i=sprintf(clt0_Mstar,"%s\n",galaxy_name);
    }
    strncat(file_label,clt0_Mstar,strlen(clt0_Mstar)-1); 
    
    if(mod_in_sigma_normA0!=10){ //star formation rate field systematic changes
        i=sprintf(cmod_in_sigma_normA0,"N%1d\n",mod_in_sigma_normA0 );
    }else{
        i=sprintf(cmod_in_sigma_normA0,"A%2.2f\n",normA0*psi0 );
    }
  strncat(file_label,cmod_in_sigma_normA0,strlen(cmod_in_sigma_normA0)-1);
  
  strncat(file_label,field_type,strlen(field_type)); //psi form
  if(strcmp( field_type, "ofp" ) == 0 || strcmp( field_type, "ofp" ) == 0){
      if(strcmp( field_type, "ofp" ) == 0 ){
          i=sprintf(calpha,"a%1.2f\n",alpha );  strncat(file_label,calpha,strlen(calpha)-1);
      }
      i=sprintf(cbeta1,"b%1.2f\n",beta1 );  strncat(file_label,cbeta1,strlen(cbeta1)-1);
  }
  if(strcmp( field_type, "oda" ) == 0 || strcmp( field_type, "kda" ) == 0){
      i=sprintf(cdummy,"op%d\n",iopt );  strncat(file_label,cdummy,strlen(cdummy)-1);
  }
  i=sprintf(ct0_obs_shift,"dt%1.2f\n",t0_obs_shift/1e9 );  strncat(file_label,ct0_obs_shift,strlen(ct0_obs_shift)-1);

  if(c_IRR!=0.0 && NITER == 1){
      i=sprintf(cdummy,"IRR%.0f\n",c_IRR*100 );  strncat(file_label,cdummy,strlen(cdummy)-1);
  }
    //IMF
  pickIMFname(file_label);


  strncat(file_label,".dat",4);
  //-------------------------------------


  mtree = my_alloc( merger_list, 1);
  mtree->ntree = npossible_progeny;
  mtree->list = my_alloc(merger_tree_struct, mtree->ntree);





  my_assert(npossible_progeny==1);
  for(i=ilast ; i>=0 ; i--){
      cum_phi[i] = 0;
      cum_Mstar[i] = 0;
  }
  for(itrial=0;itrial<ntrials;itrial++){ 
      for(i=ilast ; i>=0 ; i--){
          for(iprog=0;iprog<npossible_progeny;iprog++){
              Mstar[i][iprog] = 0;
              Mstar_dot[i][iprog] = 0;
              phi[i][iprog] = 0;
              dMLdt[i][iprog] = 0;
              idivide[i][iprog] = 0;

              massloss_to_t0[i] = 0;
          }
      }
      for(iprog=0;iprog<npossible_progeny;iprog++){
          mtree->list[iprog].idescendent=0;
          mtree->list[iprog].zform=0;
          mtree->list[iprog].mu=0;
      }
      Mstar[ilast][0]=t0_Mstar;
      Mstar_dot[ilast][0] = t0_sfr;
      nprogeny = 1;
      for(iter=0;iter<NITER;iter++){ 
          for(i=ilast ; i>=0 ; i--){
              iprog = 0;
              while(iprog<nprogeny){
                  phi[i][iprog] = normA0 * psi_field( Mstar[i][iprog] , red_binl[i], field_type, iopt, ierr);
                  
                  if(iter==0){dMLdt[i][iprog] = c_IRR*phi[i][iprog];}  
                  Mstar_dot[i][iprog] =  phi[i][iprog] - dMLdt[i][iprog]; 
                  if( Mstar[i][iprog] + Mstar_dot[i][iprog]*(-dt) < 1e-5*t0_Mstar){ 
                      phi[i][iprog] = 0;
                      Mstar_dot[i][iprog] = 0;
                      Mstar[i][iprog] = 0;
                  }
                  if(i>0){Mstar[i-1][iprog] = Mstar[i][iprog] + Mstar_dot[i][iprog]*(-dt);}
                  iprog++;
              } //iprog
          }//ilast
 
          
          //use Mstar(t) to figure out dMLdt(t):
          //measured phi gives stellar mass born at a given time after all generations of recycling:
          //since phi already includes iterated generations of star formation, there is no need to iterate to get mass lost
          for(i=0;i<ntimebins;i++){
              for(iprog=0;iprog<nprogeny;iprog++){
                  dMLdt[i][iprog] = 0; 
              }
          }
          for(iprog=0;iprog<nprogeny;iprog++){
              ipar = iprog;
              for(i=0;i<ntimebins;i++){//calculate stellar mass lost to time i

                  while( red_binl[i] < mtree->list[ipar].zform ){ // which descendent loses this mass at z[i]?
                      ipar = mtree->list[ipar].idescendent;
                  }

                  G0=0;
                  for(j=0;j<i;j++){
                      if(j==0){ 
                          mltoi = Mstar[0][iprog]/dt * dfml(dt,time_binl[i]-time_binl[0]); //all stars born before start form instantly at start
                      }else {   
                          mltoi = phi[j][iprog] * dfml(dt,time_binl[i]-time_binl[j]);
                      }
                      G0 = G0 + mltoi;
                      if(i==ilast && iter==NITER-1){ massloss_to_t0[j] += G0; }
                  }
                  dMLdt[i][ipar] += G0; 
              }//i
          }//iprog
          my_debug("iter=%d; z=0: %f mdot %f sfr %f dmldt %f ",iter,log10(Mstar[ilast][0]),Mstar_dot[ilast][0],phi[ilast][0],dMLdt[ilast][0]);
          
      }// iter 

      for(i=ilast ; i>=0 ; i--){
          iprog=0;
          for(iprog=0;iprog<nprogeny;iprog++){
              if(itrial == 0){
                  my_assert(iprog==0);
                  phi_m0[i] = phi[i][iprog];
                  Mstar_m0[i] = Mstar[i][iprog];
                  Mstar_dot_m0[i] = Mstar_dot[i][iprog];
                  dMLdt_m0[i] = dMLdt[i][iprog];
                  //                 if(i>ilast-3){my_debug("snl1 %d %f %e",i,phi_m0[i],Mstar_m0[i]);}
              }else{
                  cum_phi[i] += phi[i][iprog]/(ntrials-1.0);
                  cum_Mstar[i] += Mstar[i][iprog]/(ntrials-1.0);
              }
          }
      }
      my_debug("next trial=%d/%d ----------------",itrial+1,ntrials);
  }////////////////////////////itrial
  
/* MAIN LOOP end*/
  
  /******** calculate CSFH_frac *******/
  for(i=0 ; i<ntimebins; i++){
      CSFH_m0[i] = 0;
  }
  double dummy0=0;
  for(i=0 ; i<ntimebins; i++){
      dummy0  += phi_m0[i]*dt;
      CSFH_m0[i]             = dummy0;
  }
  for(i=0 ; i<ntimebins; i++){
      CSFH_m0[i] /= CSFH_m0[ntimebins-1];
  }
  /******** convolve age errors *******/
 

  /******** output *******/
  strcpy(fileoutsfh, "out/avgsfh-"); 
  strncat(fileoutsfh,file_label,strlen(file_label));
  my_debug("/******output: %s *********/ ",fileoutsfh);

  Mbar_Noeske=root_finder( Mbar_Noeske_from_Mgal,1e4,1e16,1e-9,1e-9 );
  if ( (fp=fopen(fileoutsfh,"w")) == NULL ){my_error("couldn't open file");}
  fprintf(fp,"#t0_Mstar=%e (N07Mbar=%e) t0_sfr=%e alpha=%e beta1=%e normA0_1sigma=%e NITER=%d c_IRR=%f",t0_Mstar,Mbar_Noeske,t0_sfr,alpha,beta1,normA0_1sigma, NITER, c_IRR);
  fprintf(fp,"\n");
  fprintf(fp,"#tlb ,red_binl[i],Mstar[i],CSFH_fraction[i],phi[i],dMLdt[i],massloss_to_t0[i]\n",fwhm10[0],fwhm10[1],fwhm10[2] ) ;
  for(i=0 ; i<ntimebins; i++){
      tlb = tlb_from_tphys(time_binl[i]);
      fprintf(fp,"%e %f  %e %f %f   %f %f   \n",
              tlb ,red_binl[i],
              Mstar_m0[i], CSFH_m0[i], phi_m0[i], //5
              dMLdt_m0[i], massloss_to_t0[i] 
//              Mstar_Noeske(tlb,Mbar_Noeske),phi_Noeske(tlb,Mbar_Noeske)  //9
              ) ;  
      
  }
  fclose(fp);



  /******** output formation data *******/
  for(i=0 ; i<ntimebins; i++){
      tlb = tlb_from_tphys(time_binl[i]);
      fracM=Mstar_m0[i]/t0_Mstar;
      if(mfz1  == -1 && red_binl[i] < 1){ mfz1  = fracM; }
      if(zfp10  == -1 && fracM > 0.1    ){ zfp10  = red_binl[i]; tfp10 = tlb_from_tphys(tphys_from_auni(1/(zfp10 + 1.0))); }
      if(zfp15  == -1 && fracM > 0.15   ){ zfp15  = red_binl[i]; tfp15 = tlb_from_tphys(tphys_from_auni(1/(zfp15 + 1.0))); }
      if(zfp20  == -1 && fracM > 0.2    ){ zfp20  = red_binl[i]; tfp20 = tlb_from_tphys(tphys_from_auni(1/(zfp20 + 1.0))); }
      if(zfp25 == -1 && fracM > 0.25   ){ zfp25 = red_binl[i]; tfp25= tlb_from_tphys(tphys_from_auni(1/(zfp25+ 1.0))); }
      if(zfp50  == -1 && fracM > 0.5    ){ zfp50  = red_binl[i]; tfp50 = tlb_from_tphys(tphys_from_auni(1/(zfp50 + 1.0))); }
  }
  //------------------------------------- filename
    strcpy(file_label, ""); 
    if(mod_in_sigma_normA0!=10){ //star formation rate field systematic changes
        i=sprintf(cmod_in_sigma_normA0,"N%1d\n",mod_in_sigma_normA0 );
    }else{
        i=sprintf(cmod_in_sigma_normA0,"A%2.2f\n",normA0*psi0 );
    }
    strncat(file_label,cmod_in_sigma_normA0,strlen(cmod_in_sigma_normA0)-1);
    
    strncat(file_label,field_type,strlen(field_type)); //psi form
    if(strcmp( field_type, "ofp" ) == 0 || strcmp( field_type, "ofp" ) == 0){
        if(strcmp( field_type, "ofp" ) == 0 ){
          i=sprintf(calpha,"a%1.2f\n",alpha );  strncat(file_label,calpha,strlen(calpha)-1);
        }
        i=sprintf(cbeta1,"b%1.2f\n",beta1 );  strncat(file_label,cbeta1,strlen(cbeta1)-1);
    }
    if(strcmp( field_type, "oda" ) == 0 || strcmp( field_type, "kda" ) == 0){
        i=sprintf(cdummy,"op%d\n",iopt );  strncat(file_label,cdummy,strlen(cdummy)-1);
    }
    i=sprintf(ct0_obs_shift,"dt%1.2f\n",t0_obs_shift/1e9 );  strncat(file_label,ct0_obs_shift,strlen(ct0_obs_shift)-1);
    
    if(c_IRR!=0.0 && NITER == 1){
        i=sprintf(cdummy,"IRR%.0f\n",c_IRR*100 );  strncat(file_label,cdummy,strlen(cdummy)-1);
    }
    //IMF
    pickIMFname(file_label);
  strncat(file_label,".dat",4);

  strcpy( fileoutform,"outform/form-");
  strncat(fileoutform,file_label,strlen(file_label));
  if ( (fp=fopen(fileoutform,"a+")) == NULL ){my_error("couldn't open file");}
  fprintf(fp,"%e %e  %f %f %f %f %f  %e %e %e %e %e\n",t0_Mstar,mfz1, zfp10,zfp15,zfp20,zfp25,zfp50, tfp10,tfp15,tfp20,tfp25,tfp50 );
  fclose(fp);
  /******** output formation data *******/

  /******** output beta *******/
  /* tfirst = tphys_from_auni(0.1); */
  /* dt = ( age_of_universe - tfirst )/(1.0*ntimebins); */
  /* for( i = 0; i < ntimebins ; i++){ */
  /*     time_binl[i] = tfirst + dt*(i+1); //first to now */
  /*     red_binl[i] = 1/auni_from_tphys(time_binl[i])-1; */
  /* } */
  /* ilast = ntimebins-1; */
  strcpy(fileoutbeta, "out/betafit-"); 
  strncat(fileoutbeta,field_type,strlen(field_type));
  if(strcmp( field_type, "ofp" ) == 0 || strcmp( field_type, "ofp" ) == 0){
      if(strcmp( field_type, "ofp" ) == 0 ){
          i=sprintf(cdummy,"a%1.2f\n",alpha );  strncat(fileoutbeta,cdummy,strlen(cdummy)-1);
      }
      i=sprintf(cdummy,"b%1.2f\n",beta1 );  strncat(fileoutbeta,cdummy,strlen(cdummy)-1);
  }
  if(strcmp( field_type, "oda" ) == 0 || strcmp( field_type, "kda" ) == 0){
      i=sprintf(cdummy,"op%d\n",iopt );  strncat(fileoutbeta,cdummy,strlen(cdummy)-1);
  }
  strncat(fileoutbeta,".dat",4);
  my_debug("filename: %s",fileoutbeta);

  if ( (fp=fopen(fileoutbeta,"w")) == NULL ){my_error("couldn't open file");}
  double m=1e11, y;
  for(i=ilast ; i>=0 ; i--){
      y=log10( normA0*psi_field(m,red_binl[i], field_type, iopt, ierr)/m )+9;
      fprintf(fp,"%f %f %f\n", red_binl[i] ,y ,psi_beta );
  }
  fclose(fp);
  /******** output beta*******/







}

double dfml(double dt, double tij){
#ifdef IMF_CONROY
  static double c0=0.05; 
  static double T0=3e5;
#elif defined(IMF_JUNGWIERT)
  static double c0=0.055; 
  static double T0=5e6;
#elif defined(IMF_C03)
  static double c0=0.046; 
  static double T0=2.76e5;
#elif defined(IMF_S98)
  static double c0=0.062; 
  static double T0=6.97e6;
#elif defined(IMF_C03STEEP)
  static double c0=0.051; 
  static double T0=1.33e7;
#else
  #error "bad IMF?"
#endif
  double res;

  if(tij-dt > 0){
    res= c0*log((tij+T0)/T0) - c0*log( (tij+T0-dt)/T0  ); // the difference in mass loss rates
  }else if(tij+dt/2 > 0){ //huh?
    res= c0*log((tij+T0)/T0); //this is the last step and all mass loss so  far is fml 
  }else{
    res=0;
  }
  return res; 
}
void pickIMFname(char *fileout){
#ifdef IMF_CONROY
  strncat(fileout,"IMFc",4);
#elif defined(IMF_JUNGWIERT)
  strncat(fileout,"IMFj",4);
#elif defined(IMF_C03)
  strncat(fileout,"IMFc03",6);
#elif defined(IMF_S98)
  strncat(fileout,"IMFs98",6);
#elif defined(IMF_C03STEEP)
  strncat(fileout,"IMFcs3",6);
#else
  #error "bad IMF?"
#endif

}


double psi_field(double mass, double z, char field_type[3], int iopt, int ierr){
    double psi_A0,psi_alpha, psi_zeta,lg_mass;
    double psi_beta_err, psi_alpha_err,psi_A0_err, psi_zeta_err;
    int i;
    static int nels_karim11_z,nels_karim11_mass, nels_oliv10_z;
    const double zbreak=2.0;

    my_assert(iopt>=-1  && iopt<=3);

    if(mass>0){
        lg_mass=log10(mass);
    }else{
        lg_mass=-1;
    }
    
    if( strcmp( field_type, "ofp" ) == 0 || 
        strcmp( field_type, "ofe" ) == 0 || 
        strcmp( field_type, "kfe" ) == 0 || 
        strcmp( field_type, "kfp" ) == 0 ){
        if(iopt==-1){return 1;}
    }
    
        



    
    if(strcmp( field_type, "ofp" ) == 0){//OLIV10_FIT_POW
        my_assert(iopt == 0);
        const double Mobs_norm = 5.62341E10; //pow(10.,10.75);
        psi_A0 = pow(10,-10.30) * pow( t0_Mstar/Mobs_norm,beta1 ) * Mobs_norm; // the empirical relation
        psi_beta=beta1-1;
        if(z < zbreak){
            return psi_A0*pow(mass/t0_Mstar,psi_beta+1)*pow(1.0+z,alpha);
        }else{
            return psi_A0*pow(mass/t0_Mstar,psi_beta+1)*pow(1.0+zbreak,alpha);
        }
    }

    if(strcmp( field_type, "ofe" ) == 0){ //OLIV10_FIT_EXP
        my_assert(iopt == 0);
        psi_beta=beta1-1;
        psi_A0 =    -1.12; 
        psi_alpha=  -3.14; 
        psi_zeta =   4.02; 
        if(z < zbreak ){ 
            return pow(10,psi_A0+11-9.0) * pow( 10,(lg_mass-11)*(psi_beta+1) ) * pow( z+1 , psi_alpha) * exp( psi_zeta * z );
        }else{
            return pow(10,psi_A0+11-9.0) * pow( 10,(lg_mass-11)*(psi_beta+1) ) * pow( zbreak+1 , psi_alpha) * exp( psi_zeta * zbreak );
        }
    }

    if(strcmp( field_type, "kfp" ) == 0){  //KARIM11_MLEFIT_POW
        my_assert(iopt == 0);
        psi_beta = -0.35;
        psi_A0 =   3.23 ;
        psi_alpha= 3.45 ;
        return     psi_A0 * pow(10,(lg_mass-11)*(1+psi_beta)) * pow(z+1,psi_alpha);
    }

    if(strcmp( field_type, "kfe" ) == 0){  //KARIM11_MLEFIT_EXP
        my_assert(iopt == 0);
        psi_beta  = -0.36;    psi_beta_err  = 0.04;
        psi_A0    =    2.17;  psi_A0_err    = 0.13; 
        psi_alpha =  5.61;    psi_alpha_err = 0.08;
        psi_zeta  = -0.99;    psi_zeta_err  = 0.04;

        psi_beta  = psi_beta +   ierr*psi_beta_err;
        psi_A0    = psi_A0 +       ierr*psi_A0_err;
        psi_alpha = psi_alpha + ierr*psi_alpha_err;
        psi_zeta  = psi_zeta +   ierr*psi_zeta_err;
        return     psi_A0 * pow(10,(lg_mass-11)*(1+psi_beta))*pow(z+1,psi_alpha)*exp(psi_zeta*z); //* pow( (z+1)/(karim11_z[zbin]+1), psi_alpha );
    }
    
    

    if(strcmp( field_type, "kda" ) == 0){ //FIELD_KARIM11_DATA_AVG
        static double karim11_z[Ntable],karim11_A0[Ntable],karim11_A0_err[Ntable],karim11_beta[Ntable],karim11_beta_err[Ntable];
        static double karim11_mass[Ntable],karim11_alpha[Ntable],karim11_alpha_err[Ntable];
    
        if(iopt == -1){
            i=1;
            karim11_z[i]=0.3 ; karim11_A0[i]=-1.27; karim11_A0_err[i]=0.03; karim11_beta[i]=-0.44; karim11_beta_err[i]=0.04; i++;
            karim11_z[i]=0.5 ; karim11_A0[i]=-0.90; karim11_A0_err[i]=0.03; karim11_beta[i]=-0.42; karim11_beta_err[i]=0.03; i++;
            karim11_z[i]=0.7 ; karim11_A0[i]=-0.67; karim11_A0_err[i]=0.03; karim11_beta[i]=-0.40; karim11_beta_err[i]=0.03; i++; 
            karim11_z[i]=0.9 ; karim11_A0[i]=-0.48; karim11_A0_err[i]=0.03; karim11_beta[i]=-0.38; karim11_beta_err[i]=0.03; i++; 
            karim11_z[i]=1.1 ; karim11_A0[i]=-0.38; karim11_A0_err[i]=0.03; karim11_beta[i]=-0.46; karim11_beta_err[i]=0.03; i++; 
            karim11_z[i]=1.4 ; karim11_A0[i]=-0.12; karim11_A0_err[i]=0.03; karim11_beta[i]=-0.30; karim11_beta_err[i]=0.03; i++; 
            karim11_z[i]=1.8 ; karim11_A0[i]= 0.10; karim11_A0_err[i]=0.07; karim11_beta[i]=-0.41; karim11_beta_err[i]=0.03; i++; 
            karim11_z[i]=2.25; karim11_A0[i]= 0.22; karim11_A0_err[i]=0.06; karim11_beta[i]=-0.42; karim11_beta_err[i]=0.05; i++; 
            karim11_z[i]=2.75; karim11_A0[i]= 0.43; karim11_A0_err[i]=0.16; karim11_beta[i]=-0.44; karim11_beta_err[i]=0.15; i++;
            nels_karim11_z=i;
            i=1;
            karim11_mass[i]=9.6;  karim11_alpha[i]=3.02; karim11_alpha_err[i]=0.15;i++;
            karim11_mass[i]=10.0; karim11_alpha[i]=3.42; karim11_alpha_err[i]=0.07;i++;
            karim11_mass[i]=10.4; karim11_alpha[i]=3.62; karim11_alpha_err[i]=0.04;i++;
            karim11_mass[i]=10.8; karim11_alpha[i]=3.48; karim11_alpha_err[i]=0.04;i++;
            karim11_mass[i]=11.1; karim11_alpha[i]=3.40; karim11_alpha_err[i]=0.06;i++;
            nels_karim11_mass=i;
            my_debug("initializing COSMOS");
            return 1;
        }       
        
        
        psi_beta=interpolate(z, karim11_z, nels_karim11_z, karim11_beta, iopt );
        psi_A0  =interpolate(z, karim11_z, nels_karim11_z, karim11_A0,   iopt );
    
        return     pow(10,psi_A0+11-9.0) * pow(10,(lg_mass-11)*(1+psi_beta)); //* pow( (z+1)/(karim11_z[zbin]+1), psi_alpha );
/* 0.3  -1.27 0.03 -0.44 0.04 */
/* 0.5  -0.90 0.03 -0.42 0.03 */
/* 0.7  -0.67 0.03  0.40 0.03 */
/* 0.9  -0.48 0.03  0.38 0.03 */
/* 1.1  -0.38 0.03  0.46 0.03 */
/* 1.4  -0.12 0.03  0.30 0.03 */
/* 1.8   0.10 0.07  0.41 0.03 */
/* 2.25  0.22 0.06  0.42 0.05 */
/* 2.75  0.43 0.16  0.44 0.15 */
    }
        

     

    

    if(strcmp( field_type, "oda" ) == 0){ //FIELD_OLIV10_DATA_AVG
        static double oliv10_z[Ntable],oliv10_A0[Ntable],oliv10_A0_err[Ntable],oliv10_beta[Ntable],oliv10_beta_err[Ntable];
        
        if(iopt == -1){
            my_debug("initializing OLIV10");
            i=1;
            oliv10_z[i] = 0.10; oliv10_A0[i]= -1.08; oliv10_A0_err[i]= 0.01; oliv10_beta[i]= 0;     oliv10_beta_err[i]= 0.02;i++;
            oliv10_z[i] = 0.25; oliv10_A0[i]= -0.93; oliv10_A0_err[i]= 0.01; oliv10_beta[i]= 0.01;  oliv10_beta_err[i]= 0.02;i++;
            oliv10_z[i] = 0.35; oliv10_A0[i]= -0.96; oliv10_A0_err[i]= 0.01; oliv10_beta[i]= -0.16; oliv10_beta_err[i]= 0.03;i++;
            oliv10_z[i] = 0.45; oliv10_A0[i]= -0.84; oliv10_A0_err[i]= 0.01; oliv10_beta[i]= -0.12; oliv10_beta_err[i]= 0.03;i++;
            oliv10_z[i] = 0.55; oliv10_A0[i]= -0.75; oliv10_A0_err[i]= 0.01; oliv10_beta[i]= -0.20; oliv10_beta_err[i]= 0.04;i++;
            oliv10_z[i] = 0.70; oliv10_A0[i]= -0.61; oliv10_A0_err[i]= 0.02; oliv10_beta[i]= -0.26; oliv10_beta_err[i]= 0.04;i++;
            oliv10_z[i] = 0.90; oliv10_A0[i]= -0.37; oliv10_A0_err[i]= 0.02; oliv10_beta[i]= -0.30; oliv10_beta_err[i]= 0.05;i++;
            oliv10_z[i] = 1.13; oliv10_A0[i]= -0.17; oliv10_A0_err[i]= 0.04; oliv10_beta[i]= -0.32; oliv10_beta_err[i]= 0.07;i++;
            oliv10_z[i] = 1.38; oliv10_A0[i]=  0.09; oliv10_A0_err[i]= 0.11; oliv10_beta[i]= -0.30; oliv10_beta_err[i]= 0.20;i++;
            oliv10_z[i] = 1.75; oliv10_A0[i]=  0.04; oliv10_A0_err[i]= 0.16; oliv10_beta[i]= 0.17;  oliv10_beta_err[i]= 0.28;i++;
            nels_oliv10_z=i;
            return 1;
        }
        
        psi_beta= interpolate(z, oliv10_z, nels_oliv10_z, oliv10_beta, iopt );
        psi_A0  = interpolate(z, oliv10_z, nels_oliv10_z, oliv10_A0, iopt );
        
        return     pow(10,psi_A0+11-9.0) * pow(10,(lg_mass-11)*(1+psi_beta)); //* pow( (z+1)/(karim11_z[zbin]+1), psi_alpha );
/* 0.1  -1.08 0.01 */
/* 0.25 -0.93 0.01 */
/* 0.35 -0.96 0.01 */
/* 0.45 -0.84 0.01 */
/* 0.55 -0.75 0.01 */
/* 0.7  -0.61 0.02 */
/* 0.9  -0.37 0.02 */
/* 1.13 -0.17 0.04 */
/* 1.33  0.09 0.11 */
    }   


}

double lin_interpolate(int ind_low, double x, double x_arr[], double y_arr[]  ){
    double y,dya,dxa,dx,dy;
    dya = y_arr[ind_low+1] - y_arr[ind_low];
    dxa = x_arr[ind_low+1] - x_arr[ind_low];
    dx = x-x_arr[ind_low];
    dy = dya/dxa*dx;
        return y_arr[ind_low] + dy;
}
double interpolate(double x, double x_arr[], int nels, double y_arr[], int iopt ){
    int bin_low,bin_high, bin;
    
    bin=0;
    while( bin < nels-2 ){ 
        if(x >= x_arr[bin+1]){ //if x is greater than the right side of the bin increase it
            //my_debug("snl %f %f %d\n\n\n\n",x,x_arr[bin],bin);
            bin++; 
        }else{ 
            break;
        }
    }
    bin_low = bin;
    bin_high=bin_low+1;
    if(bin == 0){
        bin = 1;
        bin_low =  bin ;
        bin_high = bin + 1;
    }
    if(bin > nels - 1) {
        bin = nels - 1;
        bin_high = bin; 
        bin_low  = bin - 1;
    }

    if(iopt == 0){ //fixed ends
        if(x > x_arr[nels-1]) {  return y_arr[nels-1]; 
        }else if(x < x_arr[1]) { return y_arr[1]; 
        }else{                   return lin_interpolate(bin_low,x,x_arr,y_arr);}
    }else if(iopt == 1){ //fixed top extrapolated bottom
        if(x > x_arr[nels-1]) {  return y_arr[nels-1]; 
        }else{                   return lin_interpolate(bin_low,x,x_arr,y_arr);}
    }else if(iopt == 2){ //fixed bottom extrapolated top
        if(x < x_arr[1]) {       return y_arr[1]; 
        }else{                   return lin_interpolate(bin_low,x,x_arr,y_arr);}
    }else if(iopt == 3){ //extrapolated 
                                 return lin_interpolate(bin_low,x,x_arr,y_arr);
    }
        
}
    

