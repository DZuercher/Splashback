#ifndef COSMOLOGY_H
#define COSMOLOGY_H

# include "cosmology.h"
# include "gauleg.h"
# include <cmath>
# include <string>
# include <iostream>
# include <fstream>
# include <gsl/gsl_integration.h>
# include <gsl/gsl_sf.h>
# include <gsl/gsl_math.h>
# include <gsl/gsl_roots.h>
# include <gsl/gsl_spline.h>
# include <gsl/gsl_matrix.h>
# include <gsl/gsl_odeiv.h>
# include <gsl/gsl_errno.h>

struct cosmo
{
    double  Om0,Omk,w0,wa,Omb,hval,th,s8,nspec,ximax,cfac;
};

struct gf_par{
    double Omega0,OmegaL,w0,wa;
};

/// march_params object
struct march_params
{
    double fac;
};

class cosmology;

struct qk_params
{
    cosmology *cptr;
    double *k;
    double *z;
    double *xmax;
    int opt;
};

struct projwpk_params
{
    cosmology *cptr;
    double *R;
    double *z;
    double *fkai;
};

/// Object to pass cosmology, cvir, Omega(z), Deltacrit(z)
struct cDel_params
{
    cosmology *cptr;
    double *cvir;
    double *omegaz;
    double *dcz;
    double *Delta;
};

double dTime(double,void*);
double dChi (double,void*);
//double df3 (double,void*);
double dneffint(double, void*);
double dCint(double, void*);
double findmvir(double, void*);
double findksig(double, void*);
double dxi_L(double, void*);
double dPktest_L(double, void*);
double dPktest_NL(double, void*);
double dPktest_zetaNL(double, void*);
double dxi_NL(double, void*);
double E_sq(gf_par&, double&);
double dE_sqdz(gf_par&,double&);
void getall(gf_par&,double&,double&,double&,double&);
double d2lnE_sqdz2(gf_par&,double&);
int gf_func(double, const double[], double [], void*);
int gf_jac(double, const double[], double*, double [], void*);
double findrz(double x, void *params);
class hod;
double findzmax(double x, void *params);
double dwpnl(double x, void* params);
double dwpl(double x, void* params);
double dQk(double,void*);

double dxinlbar(double,void*);
double dxinlbarbar(double,void*);
double dxilbar(double,void*);
double dxilbarbar(double,void*);
double dwpnl_kaiser(double, void*);
double dwpl_kaiser(double, void*);
double dvar_G(double x, void * params);
double dvar_TH(double x, void * params);
double findcmarch(double, void*);
double findcDel(double, void*);
double findcDelp(double, void*);

class cosmology
{
    protected:
    /// Variables
    double Omega0,Omegal,Omegak,w0,wa,Omegab,h,theta,sigma8,ns,d2norm,t0,rho_crit_0,facmtor;
    double xiNLzetamax,cfactor;

    double *x9_16, *w9_16;
    double *x0_2p, *w0_2p;

    int N9_16;
    //double x9_16[N9_16],w9_16[N9_16];

    int N0_2p;
    //double x0_2p[N0_2p],w0_2p[N0_2p];

    bool verbose, mock, takahashicorr, peacockcorr;

    private:
    /// Some constants
    double kmpspMpctoGyr, gee, c, e;

    int Nsigma;
    int Ngf;
    int Npower;
    int Nxi;
    int Nxibar;

    /// Options for various functions
    int opt_mf,opt_b,opt_c,opt_ps_L,opt_ps_NL;

    /// Numerical interpolation units for comoving distance
    bool bool_init_Chi;
    gsl_interp_accel *Chi_acc;
    gsl_spline *Chi_spline;

    /// Eisenstein and Hu 98 Power spectrum variables
    bool bool_initPS;
    double zeq, keq, zd, Rd, Req, sd, ksilk, alphac, betacinv, alphab, betanode, betab;

    /// Smith et al. 2003 Non linear power spectrum variables
    bool bool_initSmith;
    double an,bn,cn,gamman,alphan,betan,mun,nun,f1,f2,f3,ksigma,smithC,smithneff;

    /// Numerical interpolation units for variance
    bool bool_init_varM_TH_spline;
    gsl_interp_accel *varM_TH_num_acc;
    gsl_spline *varM_TH_num_spline;

    /// Numerical interpolation units for growth factor and f3
    bool bool_init_GF;
    gsl_interp_accel *GF_acc;
    gsl_spline *GF_spline;

    /// Numerical interpolation units for linear power spectra
    bool bool_init_PSL0;
    gsl_interp_accel *PSL0_acc;
    gsl_spline *PSL0_spline;
    double PSL0_dlow, PSL0_dhigh;
    double PSL0_xlow,PSL0_ylow,PSL0_xhigh,PSL0_yhigh;
    double kmin, kmax;

    /// Numerical interpolation units for non-linear power spectra
    bool bool_init_PSNL;
    gsl_interp_accel *PSNL_acc;
    gsl_spline *PSNL_spline;
    double PSNL_dlow, PSNL_dhigh;
    double PSNL_xlow,PSNL_ylow,PSNL_xhigh,PSNL_yhigh;

    /// Numerical interpolation units for linear power spectra
    bool bool_init_xiL0;
    gsl_interp_accel *xiL0_acc;
    gsl_spline *xiL0_spline;
    double xiL0_dlow, xiL0_dhigh;
    double xiL0_xlow,xiL0_ylow,xiL0_xhigh,xiL0_yhigh;
    double rmin, rmax;

    /// Numerical interpolation units for nonlinear power spectra
    bool bool_init_xiNL;
    gsl_interp_accel *xiNL_acc;
    gsl_spline *xiNL_spline;
    double xiNL_dlow, xiNL_dhigh;
    double xiNL_xlow,xiNL_ylow,xiNL_xhigh,xiNL_yhigh;
    double zeta_rmax,zetamax;

    /// Numerical interpolation for ukofm
    int Nuk;
    double cmin, cmax, krsmin, krsmax;

    // double uk_krs[Nuk];
    // double uk_c[Nuk];
    // double ukrsc[Nuk][Nuk];
    double *uk_krs, *uk_c, *ukrsc;
    bool bool_inituk;
    gsl_interp_accel *uk_c_acc;
    gsl_interp_accel *uk_krs_acc;

    /// Numerical interpolation units for xiNLbar
    bool bool_init_xiNL_bar;
    gsl_interp_accel *xiNL_bar_acc;
    gsl_spline *xiNL_bar_spline;

    bool bool_init_xiNL_barbar;
    gsl_interp_accel *xiNL_barbar_acc;
    gsl_spline *xiNL_barbar_spline;

    bool bool_init_xiL_bar;
    gsl_interp_accel *xiL_bar_acc;
    gsl_spline *xiL_bar_spline;

    bool bool_init_xiL_barbar;
    gsl_interp_accel *xiL_barbar_acc;
    gsl_spline *xiL_barbar_spline;

    /// f_g(M) swindle
    double fgm_m0;
    double fgm_slp;


    /// Private functions
    void initialize();          //Initializations
    double Chiofz(double);    //Comoving distance h^{-1} Mpc for flat cosmology, Chi for non-flat cosmology
    double Chiofz_num(double);    //Comoving distance h^{-1} Mpc for flat cosmology, Chi for non-flat cosmology
    double growthfactor(double);  //Growth factor dummy

    void init_growthfactor();   // Speed up for growth factor calculations
    void init_Chi();            // Speed up comoving distance calculations
    double initPS_EH();         // Initialize Power spectrum variables

    void fixd2norm();                // Fix the power spectrum normalization


    void init_powerspectra_NL(double); //Speed up for power spectra calculations
    void init_powerspectra_L(); //Speed up for power spectra calculations
    double Delta2_L(double,double); // \Delta^2(k) Power spectrum, k should be in units of h Mpc^{-1} wrapper
    double Delta2_NL(double,double); // \Delta^2(k) Power spectrum, k should be in units of h Mpc^{-1} wrapper
    double Pk_L(double,double);     // \P(k) Power spectrum, k should be in units of h Mpc^{-1} wrapper
    double Pk_NL(double,double);    // Non-linear power spectrum wrapper

    // Correlation function calculations: Tests
    double xi_L(double,double);       // Linear correlation function quadpack
    double xi_NL(double,double);      // Nonlinear correlation function quadpack
    double xi_L_old(double,double);   // Linear correlation function gauleg way
    double xi_NL_old(double,double);  // Nonlinear correlation function gauleg way

    // Global redshift indicator for different splines
    double z_glob;
    double gf_glob;

    //Speed up the correlation calculations
    void init_xi_L(); 
    void init_xi_NL(double);

    double Delta2_EH(double,double); // \Delta^2(k) Power spectrum, k should be in units of h Mpc^{-1}
    double Pk_EH(double,double);     // \P(k) Power spectrum, k should be in units of h Mpc^{-1}
    double TCold_EH(double);         // EH 98 Transfer function for CDM k in Mpc^{-1} 
    double Tb_EH(double);            // EH 98 Transfer function for baryons k in Mpc^{-1}
    double T0master_EH(double,double,double);  // EH 98 Transfer function master for CDM k in Mpc^{-1}
    double TF_EH(double);            // EH 98 Density-wtd Transfer function for CDM+Baryons k in Mpc^{-1}
    double var_TH(double,double);    // \sigma^2(R) R in h^{-1} Mpc assuming a top hat filter in real space
    double var_G(double,double);     // \sigma^2(R) R in h^{-1} Mpc assuming a top hat filter in real space

    double init_Smith(double);           // Initialize the non-linear power spectrum
    double getksigma(double);        // ksigma for Smith et al. 2003
    double getneff(double);          // neff for Smith et al. 2003
    double getC(double);             // C for Smith et al. 2003
    double Delta2NL_S(double,double);    // Non-linear power spectrum Smith et al. 2003
    double PkNL_S(double,double);        // Non-linear power spectrum Smith et al. 2003

    void Delta2_EH_S(double,double,double&,double&); //Both linear and nonlinear power spectrum
    void Pk_EH_S(double,double,double&,double&); //Both linear and nonlinear power spectrum

    void init_varM_TH_spline();          // Initialise the numerical calculation of variance
    double varM_TH(double, double);      // Variance of the density field smoothed with a Tophat filter of mass scale M

    double MF_WA(double,double);         // Warren et al . mass function now obsolete
    double MF_ST(double,double);         // Sheth Tormen mass function now obsolete
    double MF_BH(double,double);         // Bhattacharya mass function experimental
    double MF_TI09(double,double);       // Tinker et al. 2009 mass function also obsolete
    double MF_TI09_350(double,double);   // Tinker et al. 2009 SO(350) mass function, obsolete
    double bias_TWZZ(double, double);    // Tinker et al. 2006 bias function, obsolete

    // Tinker et al. 2010 mass function and bias function variables
    double alpTink;                      // Tinker et al. 2010, mass function normalization
    bool init_Tink;                      // Normalization initialise
    void initTinker(double);             // Initialise the normalization
    double MF_TI10(double,double);       // Tinker et al. 2010 mass function
    double bias_TI10(double, double);    // Tinker et al. 2010 bias
    
    double getzcoll(double);            // Redshift of collapse for a halo of mass M
    // double getmstar();                  // M* defined such that sigma(M*)=1.686
    double getMvir(double, double);     // Get Mvir from M200 and redshift
    double getc200(double, double);     // Get c200
    double c_MAC(double,double);        // Concentration a'la Maccio 
    double munfw(double);              // Function related to mass enclosed for NFW 

    void setrzeta(double);
    double findrzfn(double x, double tgt, double z);

    void ukinit();
    void ukinit2();


    // For Kaiser effect correction
    double xiNL_bar_num(double,double);
    void init_xiNL_bar(double);
    double xiNL_bar(double,double);

    double xiNL_barbar_num(double,double);
    void init_xiNL_barbar(double);
    double xiNL_barbar(double,double);

    double xiL_bar_num(double,double);
    void init_xiL_bar(double);
    double xiL_bar(double,double);

    double xiL_barbar_num(double,double);
    void init_xiL_barbar(double);
    double xiL_barbar(double,double);


    int kbins;
    double hod_kmin, hod_kmax;
    void initQk(double, double[]);
    double *Qk1;
    double *Qk2;
    bool bool_initQk;

    double ukofm(double,double,double); // Fourier transform of NFW profile
    double uskofm(double,double,double,double); // Fourier transform of NFW profile for satellites
    double ukinterp(double,double); // Interpolation routine to find ukofm quickly
    
    double getzeta_rmax();
    double getzetamax();

    // FT(FT(P(k)))
    double Pktest_L(double,double);
    double Pktest_NL(double,double);
    double Pktest_zetaNL(double,double);


    //Friends
    friend double dTime(double,void*);
    friend double dChi (double,void*);
    //friend double df3 (double,void*);
    friend double dneffint(double, void*);
    friend double dCint(double, void*);
    friend double findmvir(double, void*);
    friend double findksig(double, void*);
    friend double dxi_L(double, void*);
    friend double dPktest_L(double, void*);
    friend double dPktest_NL(double, void*);
    friend double dPktest_zetaNL(double, void*);
    friend double dxi_NL(double, void*);
    friend double E_sq(gf_par&, double&);
    friend double dE_sqdz(gf_par&,double&);
    friend void getall(gf_par&,double&,double&,double&,double&);
    friend double d2lnE_sqdz2(gf_par&,double&);
    friend int gf_func(double, const double[], double [], void*);
    friend int gf_jac(double, const double[], double*, double [], void*);
    friend double findrz(double x, void *params);
    friend class hod;
    friend double findzmax(double x, void *params);
    friend double dwpnl(double x, void* params);
    friend double dwpl(double x, void* params);
    friend double dQk(double,void*);
    /// For Kaiser effect 
    friend double dxinlbar(double,void*);
    friend double dxinlbarbar(double,void*);
    friend double dxilbar(double,void*);
    friend double dxilbarbar(double,void*);
    friend double dwpnl_kaiser(double, void*);
    friend double dwpl_kaiser(double, void*);
    friend double dvar_G(double x, void * params);
    friend double dvar_TH(double x, void * params);
    friend double findcmarch(double, void*);
    friend double findcDel(double, void*);
    friend double findcDelp(double, void*);

    public:
    
        // Basic cosmology
        cosmology(); //Constructor
        ~cosmology(); //Destructor
        cosmology(double om0,double omk,double w0,double wa,double omb,double h,double theta,double sigma8,double ns,double ximax,double cfac); //Constructor
        cosmology(cosmo); //Constructor
	void cosmo_free();

        // Different distances
        double Dcofz(double z);    //Comoving distance h^{-1} Mpc
        double Dlofz(double z);     //Luminosity distance h^{-1} Mpc
        double Daofz(double z);     //Angular diameter distance h^{-1} Mpc
        double Daofzlh(double zl, double zh);     //Angular diameter distance h^{-1} Mpc

        // Growth factor and f for redshift space distortions
        double growthfactor_num(double z); // Interpolating routine        
        double dlnDdln1pz(double z);

        // Density parameters as a function of redshift
        double Omega(double z);       //Omega(z) 
        double Omegaw(double z);       //Omegaw(z) 

        // Virial overdensity 
        double Delta_crit(double z);  //Delta_crit(z), Bryan and Norman '98

	void set_optmf(int opt); // Set the mass function option

        // Power spectrum calculations
        double Delta2_L_num(double k,double z);     // Numerical \Delta^2(k) Power spectrum, k should be in units of h Mpc^{-1} wrapper
        double Delta2_NL_num(double k,double z);    // Numerical Non-linear \Delta^2(k),  k should be in units of h Mpc^{-1}

	double xi_L_num(double k,double z);
	double xi_NL_num(double k,double z);

	// Mass and bias function wrappers
        double nofm(double M,double z);
        double bias(double M,double z);
	
	// Variance related functions
	double varM_TH_num(double M,double z);       
	double varM_TH_num_deriv(double M,double z);

	// Group catalog related functions
        // Number density of haloes above mass M
        double Nplus(double M200, double z);
        // Get mass M to obtain a given Number density of haloes above this mass
        double getM(double Nplus,double z);

        //NFW profile related functions
        void modelNFWhalo(double M200,double z,double& Mvir,double& Rvir,double& cvir,double& R200,double& c200); // Radii in physical units
        void modelNFWhalo_com(double M200,double z,double& Mvir,double& Rvir,double& cvir,double& R200,double& c200); // Mvir, Rvir, cvir, R200, c200
        //void modelNFWhalo_com(double,double,double&,double&,double&); //Mvir, Rvir, cvir

        double conc(double Mvir,double z); //Concentration parameter wrapper

        double Eofz(double z);        //Eofz(z) 

	/// Set new z
	void setnew_z(double z);

	/// Access to private variables
	double gets8();
	double getOmb();
	double geth();
	double getns();
    double getxinlzetamax();
    double get_cfac();
    double set_cfac(double cfac);

    /// SDSS survey specific functions
    double getzmax(double xL);
    double getLmin(double z, double L1);

    double Time(double z);        //Time(z) units 1/H0
    double Lookback(double z);    //Lookback time(z) units 1/H0

    double wpnl(double z,double rad,double projmax);
    double wpl(double z,double rad, double projmax);

    double wpnl_kaiser(double z,double rad, double projmax,double fkai);
    double xi_NL_kaiser(double r,double z, double mu,double fkai);

    double wpl_kaiser(double z,double rad, double projmax,double fkai);
    double xi_L_kaiser(double r,double z, double mu,double fkai);

    /// New functionality added renew cosmology
    void renew(cosmo p);
 
    /// Get sound horizon a'la Eisenstein Hu approximation
    double rsound();

    /// Functions to calculate distances given ra and dec
    double get_deltapi(double z1,double z2);
    double get_sinsqang(double x1,double y1,double z1,double x2,double y2,double z2);
    double get_logrp(double x1,double y1,double z1,double x2,double y2,double z2, double Chisq);
    //
    double getmstar();                  // M* defined such that sigma(M*)=1.686
    void pevolve_fixed(double cdel,int opt,double z,double zstart,double&cdelz ,double&fdelz);
    double getcDel(double cvir, double z, double Delta);
    double getRvirfromMvir(double Mvir, double z);
    double getRDelfromMDel(double Mdel, double z, double Del);
    double getcDeltap_from_cDelta(double cDelta, double Delta, double Deltap);

};

/// This is from haloes.cpp
/// Passing cosmology object, Np and z for calculating M from N_{+}
struct np_params
{
    cosmology *cptr;
    double *z;
    double *Np;
};

///Passing cosmology object
struct c_params
{
    cosmology *cptr;
};


///Passing cosmology object and variance
struct coll_params
{
    cosmology *cptr;
    double *sig;
};

///Passing cosmology object, M200 and z
struct mvir_params
{
    cosmology *cptr;
    double *m200;
    double *z;
};

///Passing cosmology object, cvir, Omega(z), Deltacrit(z)
struct c200_params
{
    cosmology *cptr;
    double *cvir;
    double *omegaz;
    double *dcz;
};

///Passing cosmology object, cvir, Omega(z), Deltacrit(z)
struct cDelta_params
{
    double *cDelta;
    double *frac;
};

/// This is for powerspectrum.cpp
/// Passing cosmology object, R and z for variance
struct cvar_params
{
    cosmology *cptr;
    double *R;
    double *z;
    bool *psinit;
};

/// Passing cosmology object, k and z for Pktest_L
struct pk_params
{
    cosmology *cptr;
    double *k;
    double *z;
};

/// Passing cosmology object, r and z for xi_L
struct xi_params
{
    cosmology *cptr;
    double *r;
    double *z;
};

/// Passing cosmology object and z for calculating k_{\sigma}
struct ksig_params
{
    cosmology *cptr;
    double *z;
};

struct projwp_params
{
    cosmology *cptr;
    double *R;
    double *z;
};

struct rz_params
{
    cosmology *cptr;
    double *z;
    double *tgt;
};

///Passing cosmology object
struct z_params
{
    cosmology *cptr;
    double *mag;
};

#endif
