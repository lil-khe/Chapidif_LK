
#include <cstdio> //old stile c IO
#include <iostream> //c++ IO
#include <cmath>   // sin cos etc. in std:: namespace (but we use "using namespace std;")
#include <complex> // allows these functions with complex values
#include <cfloat>  // defines range of variables in std:: namespace 

#include <cstring> // string operations  in std:: namespace 
#include <cstdlib> // various conversions memory management etc
#include <fstream> // i/o with files
#include <cerrno>  // for error messages
#include <cfenv>  // has to do with floating point issues

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/complex.h> // Enables complex number support

//#ifdef _WIN32
//    #include <quadmath.h>
//    #define QUADMATH =1
//#elif __linux__
#if __linux__
    #include <quadmath.h>
    #define QUADMATH =1
#endif  // so no quadmath for MAC and Windows



#  include "Faddeeva.hh"

using namespace std;
using dcomp = complex<double>;

namespace nb = nanobind;
using npArray= nb::ndarray<double, nb::numpy, nb::ndim<1>, nb::c_contig>;
using npArrayComplex=nb::ndarray<std::complex<double>, nb::ndim<1>, nb::c_contig> ; 
//typedef complex<long  double> ldcomp;
//#define MAXOSC  250// max number of oscillators, should match what is happening in calculate.py
//#define MAXGOS  10 // Max number of core levels described by GOS, should match what is happening in calculate.py
//#define MAXBELKACEM 5// maximum number of Belkacem oscillators, should match what is happening in calculate.py
//#define MAXKANEKO 250 // Max number of core levels described by Kaneko, should match what is happening in calculate.py
enum DFSize{
    MAXOSC = 250,// max number of oscillators, should match what is happening in calculate.py
    MAXGOS = 10, // Max number of core levels described by GOS, should match what is happening in calculate.py
    MAXBELKACEM= 5,// maximum number of Belkacem oscillators, should match what is happening in calculate.py
    MAXKANEKO = 250 // Max number of core levels described by Kaneko, should match what is happening in calculate.py
};

const double pi = 3.141592653589793238463;
const double Hartree = 27.211;
const double BohrRadius  = 0.529177;
const double C = 137.036;

double Projectile_Mass;
bool proton;  // true = proton else electron
bool MottCorrection;
bool debugmessageprinted=false;
bool DebugMode=false;
const int NGOSscaling = 80;
const int NTaucscaling =100;
//double Ai[MAXOSC], gammai[MAXOSC], wi[MAXOSC], alphai[MAXOSC],gapi[MAXOSC], Q_Vlasov[MAXOSC], sigmai[MAXOSC], TaucGap[MAXOSC];
double Ai[MAXOSC], gammai[MAXOSC], wi[MAXOSC], alphai[MAXOSC],gapi[MAXOSC], Q_Vlasov[MAXOSC], sigmai[MAXOSC], TaucGap[MAXOSC], Ethi[MAXOSC], deltai[MAXOSC];	///// MODIF 04/08 LK: Add Eth and delta for F edge factor
double AiGOS[MAXGOS ],Edgei_GOS[MAXGOS],Zi_GOS[MAXGOS];
double GOS_ScalingFactor[MAXGOS][NGOSscaling], Tauc_ScalingFactor[MAXOSC][NTaucscaling];
double GOS_Scaling_q[NGOSscaling], Tauc_Scaling_q[NTaucscaling];

double maxEnergyDensityEffect;

double Ai_Belkacem[MAXBELKACEM],wi_Belkacem[MAXBELKACEM],gammai_Belkacem[MAXBELKACEM];
bool  Full_dispersion, Dispersion_relativistic;
bool delayed_dispersion;
int n_i_GOS[MAXGOS ], l_i_GOS[MAXGOS ];

double N_Kaneko[MAXKANEKO ],Edge_ArchubiKaneko[MAXKANEKO],Q_Kaneko[MAXKANEKO],Width_Kaneko[MAXKANEKO];
double w_pl_l[MAXKANEKO],w_pl_0[MAXKANEKO],A_l_Kaneko[MAXKANEKO],A_0_Kaneko[MAXKANEKO],gamma_Kaneko[MAXKANEKO];
int l_Kaneko[MAXKANEKO];

double UnitCellDensity, PL_GOS;
double epsbkg, aTLan;
double w_global, theta_global; 
int I_Os_global;
double quanc8result; double quanc8errest; int quanc8nofun; double quanc8flag;  // not sure why I can not pass addresses here and have ot make the output a global variable
double E_0,p_0,p_detected,v_0, Egerton_rel_cor_factor, gamma_rel, beta_r, b_zero;

double rest_mass_energy;
bool add_Dopplerwidth_to_classical_DF=false;


bool Apply_Mermin_Correction,DirectMethod,  ApplySumRuleToGOS, OriginalKaneko, AddELF; 
bool modelDrude, modelDL, modelMerminLL,modelVlasov, modelTaucLorentz, modelTL_an, modelForouhiBloomer, modelBrendelBormann;
bool modelTauc_Mermin;
bool Orosco_Coimbra_way;
int ExchangeCorrection;
bool UseBornOchkurExchange;	//// MODIF 05/08 LK: Make it possible to use Born-Ochkur factor for exchange
bool UseExchangeForOsc[MAXOSC];	//// MODIF 18/08 LK: Make it possible to use Born-Ochkur factor for exchange (choose for each oscillator)
double EmaxForOsc[MAXOSC];	//// MODIF 18/08 LK: Add a variable to change the Emax of integration
double BindingEnergyForOsc[MAXOSC];	//// MODIF 01/09 LK: Add a variable to compute Emax from the per-oscillator binding energy for individual/ionization type
double Egap;			//// MODIF 19/08 LK: Add a variable to compute Emin from the bandgap energy
double abserr, relerr;
double FirstEnergy ,StepSize_eV,StepSize,Stepsize_qplot ,LastMomentum; 
double q_lower,q_upper;
int NStep, NThetaStep;
double lin_cont_deltaE; // for the integration of the DIIMFP for stopping, also Tauc normalisation
double q2maxFactor, q2surfmaxFactor;
double c_transition, BE_for_exchange;
double theta_max;

double mom_limit_lower[30], mom_limit_upper[30]; 
double FSumRule, BetheSumRule,GOSBetheSumRule,KKSumRule,A_ScalingFactor,w_ScalingFactor;
double theta0,theta1,sigma,PIcoef1,PIcoef2,PIcoef3,surf_ex_factor, fraction_DIIMFP;  // for reels spectrum
//dcomp  Chi_DL(double q, double w, int i);  // function declaration. we generally use function after definition.  this is the exception, hence needs a declaration
dcomp  Chi_DL_double_width(double q, double w_global,int i);
dcomp  oneoverTL(double q, double current_w, double  C, double E0, double alpha, double Egap);
//dcomp Chi_TL_an(double E, double  C, double E0,  double Egap, double a_TL_analytic);


int check_for_exceptions()
 {
    int fe;
    /* testing multiple exceptions: */
    fe = fetestexcept (FE_ALL_EXCEPT);
    if ((fe !=0)&&(fe != 32) )
    {
          if (fe & FE_DIVBYZERO) printf("FE_DIVBYZERO\n");
          if (fe & FE_INVALID)   printf("FE_INVALID\n");
          if (fe & FE_OVERFLOW)  printf("FE_OVERFLOW\n");
          if (fe & FE_UNDERFLOW) printf("FE_UNDERFLOW\n");
          return fe;
  }
 return 0;
}
void my_perror(const char* str)
{
     if (DebugMode)
     {
         if(errno !=0) perror(str);
         check_for_exceptions();
         errno=0;
    }
}

double Poisson(int  n_occurence, double Plambda)
{  double result;
    if (DebugMode) my_perror("before lgamma\n");
    result= exp(n_occurence*log(Plambda)-Plambda-lgamma(n_occurence+1.0));
     if (DebugMode) my_perror("after lgamma\n");
    return result;

}
inline double velocity_from_energy(double E, double mass)  // E in a.u. mass in a.u.
{   
    double gamma = 1.0 + E/(mass*C*C);  //rest mass energy initialised by copyP_to_Vars()
    double v = sqrt(2*E*(1.0+gamma)/(2*gamma*gamma*mass));
    return v;
}
 //if (Full_dispersion)
    //{  
        //double k_f=pow(w*w*3.0/4.0*pi,1.0/3.0);
        //w_at_q_square = w*w + 2.0 * alpha *k_f*k_f /3.0 *Qrecoil +Qrecoil * Qrecoil;  //Lundqvist way, alpha should be 1, 
    //}


inline double recoil_energy( double q) // calculates the energy a free electron gets after it absorbes momentum q
{  //recoil energy after momentum transfer q to a free, stationary electron
    double Q_recoil;
    if((q < 2.0) or !Dispersion_relativistic) 
    {
        Q_recoil = q*q/2.0;
    }
    else
    {
        Q_recoil= sqrt(C * C * q * q + C * C * C * C) - C * C;
    }
    return Q_recoil;
}

double w_at_q(double w,double q,double alpha)
{   double w_q;
    double Qrecoil = recoil_energy(q); 
    if(Full_dispersion)
    {
        double k_f=pow(w*w*3.0/4.0*pi,1.0/3.0);
        w_q= sqrt( w*w + 2.0*alpha*k_f*k_f/3.0*Qrecoil + Qrecoil*Qrecoil);  
    }
    else w_q = w+alpha*Qrecoil;
    return w_q;
}
//// MODIF 04/08 LK: Add the edge factor function
inline double EdgeFactor(double E, double Eth, double delta)
{
    if (delta == 0.0) return 1.0;	// feature off for this oscillator
    return 1.0 / (1.0 + exp(-delta * (E - Eth)));
}
//// END MODIF
inline dcomp reciprocal(dcomp chi)    //works both ways: also gets chi from oneoverchi
{ 
    dcomp eps,oneovereps, oneoverchi;
    eps=chi+dcomp(1.0,0.0);
    oneovereps=dcomp(1.0,0)/eps;
    oneoverchi= oneovereps-dcomp(1.0,0.0);
    return oneoverchi;
}

double  GOSx(int n, int  l, const double Z, double  dE_min, double q, double  w) //results GOS PER ELECTRON
                                                                           // Z = 0 uses dE_min as binding energy
{           
    double Zs, ne, Wl, Ql, Ql2, Ql3, A, bl, dummy, kH, kH2, Qrecoil;
    double c[10]; 
    int j, jmax;
    if (w < dE_min) return 0.0;
    Zs = 0.0; 
    if ((n == 1)&&(Z==1))  Zs = Z;
    else if(n == 1) Zs = Z - 0.3;
    else if (n == 2)  Zs = Z - 4.15;
    else if (n == 3)
    {
        if ((l == 0) || (l == 1))  Zs = Z - 11.25;
        if (l == 2) Zs = Z - 21.15;
    }
    if (Z == 0.0)  Zs = n*sqrt(dE_min / 0.5);

    if (Zs <= 0.0) return 0.0;
    ne = 2.0 * (2.0 * l + 1.0);
 

    //Ql = (q / Zs)*(q / Zs);  original
   
    Qrecoil = recoil_energy(q+0.01);  //subroutine depends on bool Dispersion_relativistic, added 0.001, seems else not stable at 0
  
    Ql = 2.0 * Qrecoil / ( Zs * Zs );
    Wl = w / (0.5*Zs*Zs);
    kH2 = Wl - 1.0 / (n*n);

    if (kH2 > 0.0)
    {  
        kH = sqrt(kH2);
        bl = atan(2.0 * kH / n / (Ql - Wl + 2.0 / (n*n)));
       
        if (bl < 0)  bl = bl + pi;
        double tmp1=0.0;
        if(kH > 0.02) tmp1= exp(-2.0 * pi / kH);// in this way no underflow errors
        A = 16.0 * pow(2.0 / n, 3.0)*exp(-2.0 / kH*bl) / (1.0 - tmp1) / pow((Ql - Wl)*(Ql - Wl) + 4.0 / (n*n)*Ql, 2.0 * n + 1.0);
      
    }
    else
    { 
        kH = sqrt(-kH2);
        double tmp = (Ql - Wl + 2.0 / (n*n) + 2.0 * kH / n) / (Ql - Wl + 2.0 / (n*n) - 2.0 * kH / n);
        bl = -1.0 / kH*log(tmp);
        
        A = 16.0 * pow(2.0 / n, 3.0)*exp(bl) / pow((Ql - Wl)*(Ql - Wl) + 4.0 / (n* n)*Ql, 2.0 * n + 1.0);
    }

    dummy = 0.0;
    jmax = -1;
  
    if (n == 1)
    {
        c[0] = (Ql + Wl / 3);
        jmax = 0;
    }
    if ((n == 2) && (l == 0))
    {
        Ql2 = Ql*Ql;
        c[0] = (19.0 / 60.0 + 8.0 / 15.0 * Ql)*Ql2;
        c[1] = 1.0 / 15.0 * (Ql + 1)*Ql;
        c[2] = 1.0 / 30.0 * (1.0 - 40.0 * Ql)*Ql;
        c[3] = 2.0 / 3.0 * Ql;
        c[4] = 1.0 / 4.0 + 4.0 / 3.0 * Ql;
        c[5] = 1.0 / 3.0;
        jmax = 5;
    }

    if ((n == 2) && (l == 1))
    {
        Ql2 = Ql*Ql;
        c[0] = (17.0 / 20.0 + 4.0 / 5.0 * Ql)*Ql2; // according to Spannish appendix correcting misprint Serra
                                                   //    c[0] = (17.0 / 12.0 + 4.0 / 5.0 * Ql)*Ql2; // according to Sera
        c[1] = (1.0 / 10.0 + 34.0 / 15.0 * Ql)*Ql;
        c[2] = 1.0 / 30.0 * (49.0 + 120.0 * Ql)*Ql;
        c[3] = 1.0 / 6.0 + 2.0 * Ql;
        c[4] = 1.0 / 4.0;
        jmax = 4;
            
    }

    if ((n == 3) && (l == 0))
    {
        Ql2 = Ql*Ql;
        Ql3 = Ql2*Ql;
        c[0] = (528384.0 / 502211745.0 + 561152.0 / 55801305.0 * Ql + 68608.0 / 6200145.0 * Ql2)*Ql3;
        c[1] = (32768.0 / 167403915.0 + 2048.0 / 413343.0 * Ql - 7424.0 / 6200145.0 * Ql2)*Ql2;
        c[2] = (8192.0 / 6200145.0 - 63488.0 / 1240029.0 * Ql - 17408.0 / 98415.0 * Ql2)*Ql2;
        c[3] = (256768.0 / 6200145.0 + 256.0 / 729.0 * Ql)*Ql2;
        c[4] = (11008.0 / 885735.0 + 30976.0 / 98415.0 * Ql + 15744.0 / 10935.0 * Ql2)*Ql;
        c[5] = (2816.0 / 32805.0 - 10912.0 / 10935.0 * Ql)*Ql;
        c[6] = (128.0 / 19683.0 + 1024.0 / 10935.0 * Ql - 64.0 / 27.0 * Ql2);
        c[7] = (208.0 / 2187.0 + 80.0 / 81.0 * Ql);
        c[8] = (32.0 / 81.0 + 4.0 / 3.0 * Ql);
        c[9] = 1.0 / 3.0;
        jmax = 9;
    }


    if ((n == 3) && (l == 1))
    {
        Ql2 = Ql*Ql;
        Ql3 = Ql2*Ql;
        c[0] = (495616.0 / 167403915.0 + 443392.0 / 18600435.0 * Ql + 8192.0 / 413343.0 * Ql2)*Ql3;
        c[1] = (65536.0 / 167403915.0 + 546304.0 / 18600435.0 * Ql + 38912.0 / 413343.0 * Ql2)*Ql2;
        c[2] = (274944.0 / 18600435.0 + 135424.0 / 2066715.0 * Ql + 2048.0 / 6561.0 * Ql2)*Ql2;
        c[3] = (4096.0 / 2657205.0 + 15872.0 / 137781.0 * Ql - 22528.0 / 32805.0 * Ql2)*Ql;
        c[4] = (512.0 / 32805.0 + 2368.0 / 10935.0 * Ql - 512.0 / 243.0 * Ql2)*Ql;
        c[5] = (8992.0 / 32805.0 + 1664.0 / 729.0 * Ql)*Ql;
        c[6] = (224.0 / 6561.0 + 560.0 / 243.0 * Ql + 128.0 / 27.0 * Ql2);
        c[7] = (208.0 / 729.0 + 64.0 / 27.0 * Ql);
        c[8] = 8.0 / 27.0;
        jmax = 8;
    }


    if ((n == 3) && (l == 2))
    {
        Ql2 = Ql*Ql;
        Ql3 = Ql2*Ql;
        c[0] = (253952.0 / 55801305.0 + 904192.0 / 55801305.0 * Ql + 131072.0 / 6200145.0 * Ql2)*Ql3;
        c[1] = (45056.0 / 167403915.0 + 1440256.0 / 18600435.0 * Ql + 149504.0 / 6200145.0 * Ql2)*Ql2;
        c[2] = (140800.0 / 3720087.0 + 4211968.0 / 6200145.0 * Ql + 32768.0 / 98415.0 * Ql2)*Ql2;
        c[3] = (2048.0 / 885735.0 + 2657792.0 / 6200145.0 * Ql + 48128.0 / 32805.0 * Ql2)*Ql;
        c[4] = (75008.0 / 885735.0 + 193472.0 / 98415.0 * Ql + 8192.0 / 3645.0 * Ql2)*Ql;
        c[5] = (256.0 / 59049.0 + 6304.0 / 10935.0 * Ql + 22912.0 / 10935.0 * Ql2);
        c[6] = (736.0 / 19683.0 + 6416.0 / 10935.0 * Ql);
        c[7] = 80.0 / 2187.0;
        jmax = 7;
    }

    for (j = 0; j < (jmax + 1); j++)
    {
        dummy = dummy + c[j] * pow(Wl - Ql, j);
    }

    double result = A*dummy * 2.0 * w / (0.5*(Zs*Zs)*(0.5*(Zs*Zs))) / ne;
    return result;
}


void quanc8(double(*fun)(double), double a, double b, double abserr, double relerr)
/*
estimate the integral of fun(x) from a to b to a user provided tolerance.
an automatic adaptive routine based on the 8-panel newton-cotes rule.

input:
fun     the name of the integrand function subprogram fun(x).
a       the lower limit of integration.
b       the upper limit of integration.(b may be less than a.)
relerr  a relative error tolerance. (should be non-negative)
abserr  an absolute error tolerance. (should be non-negative)

output:
result  an approximation to the integral hopefully satisfying the
least stringent of the two error tolerances.
quanc8errest  an estimate of the magnitude of the actual error.
quanc8nofun   the number of function values used in calculation of result.
quanc8flag    a reliability indicator.  if quanc8flag is zero, then result
probably satisfies the error tolerance.  if quanc8flag is
xxx.yyy , then  xxx = the number of intervals which have
not converged and  0.yyy = the fraction of the interval
left to do when the limit on  quanc8nofun  was approached.

comments:
Alex Godunov (February 2007)
the program is based on a fortran version of program quanc8.f
*/
{
    double w0, w1, w2, w3, w4, area, x0, f0, stone, cor11, temp;
    double qprev,  tolerr;
    double qright[32], f[17], x[17], fsave[9][31], xsave[9][31];
    // double dabs,dmax1;
    int    levmin, levmax, levout, nomax, nofin, lev, nim, i, j;
    int    key;

    //  ***   stage 1 ***   general initialization

    levmin = 1;
    levmax = 30; // was 30
    levout = 6;
    //nomax = 5000;
    nomax = 250000;  // was 5000
    nofin = nomax - 8 * (levmax - levout + 128);
    //  trouble when quanc8nofun reaches nofin

    w0 = 3956.0 / 14175.0;
    w1 = 23552.0 / 14175.0;
    w2 = -3712.0 / 14175.0;
    w3 = 41984.0 / 14175.0;
    w4 = -18160.0 / 14175.0;

    //  initialize running sums to zero.

    quanc8flag = 0.0;
    quanc8result = 0.0;
    cor11 = 0.0;
    quanc8errest = 0.0;
    area = 0.0;
    quanc8nofun = 0;
    if (a == b) return;

    //  ***   stage 2 ***   initialization for first interval

    lev = 0;
    nim = 1;
    x0 = a;
    x[16] = b;
    qprev = 0.0;
    f0 = fun(x0);
    stone = (b - a) / 16.0;
    x[8] = (x0 + x[16]) / 2.0;
    x[4] = (x0 + x[8]) / 2.0;
    x[12] = (x[8] + x[16]) / 2.0;
    x[2] = (x0 + x[4]) / 2.0;
    x[6] = (x[4] + x[8]) / 2.0;
    x[10] = (x[8] + x[12]) / 2.0;
    x[14] = (x[12] + x[16]) / 2.0;
    for (j = 2; j <= 16; j = j + 2)
    {
        f[j] = fun(x[j]);
    }
    quanc8nofun = 9;

    //  ***   stage 3 ***   central calculation

    while (quanc8nofun <= nomax)
    {
        x[1] = (x0 + x[2]) / 2.0;
        f[1] = fun(x[1]);
        for (j = 3; j <= 15; j = j + 2)
        {
            x[j] = (x[j - 1] + x[j + 1]) / 2.0;
            f[j] = fun(x[j]);
        }
        quanc8nofun = quanc8nofun + 8;
        double step = (x[16] - x0) / 16.0;
        double qleft = (w0*(f0 + f[8]) + w1*(f[1] + f[7]) + w2*(f[2] + f[6])
            + w3*(f[3] + f[5]) + w4*f[4]) * step;
        qright[lev + 1] = (w0*(f[8] + f[16]) + w1*(f[9] + f[15]) + w2*(f[10] + f[14])
            + w3*(f[11] + f[13]) + w4*f[12]) * step;
        double qnow = qleft + qright[lev + 1];
        double qdiff = qnow - qprev;
        area = area + qdiff;

        //  ***   stage 4 *** interval convergence test

        double esterr = fabs(qdiff) / 1023.0;
        if (abserr >= relerr*fabs(area))
            tolerr = abserr;
        else
            tolerr = relerr*fabs(area);
        tolerr = tolerr*(step / stone);

        // multiple logic conditions for the convergence test
        //key = 1;
        if (lev < levmin) key = 1;
        else if (lev >= levmax)
            key = 2;
        else if (quanc8nofun > nofin)
            key = 3;
        else if (esterr <= tolerr)
            key = 4;
        else
            key = 1;

        switch (key) {
            // case 1 ********************************* (mark 50)
        case 1:
            //      ***   stage 5   ***   no convergence
            //      locate next interval.
            nim = 2 * nim;
            lev = lev + 1;

            //      store right hand elements for future use.
            for (i = 1; i <= 8; i = i + 1)
            {
                fsave[i][lev] = f[i + 8];
                xsave[i][lev] = x[i + 8];
            }

            //      assemble left hand elements for immediate use.
            qprev = qleft;
            for (i = 1; i <= 8; i = i + 1)
            {
                j = -i;
                f[2 * j + 18] = f[j + 9];
                x[2 * j + 18] = x[j + 9];
            }
            continue;  // go to start of stage 3 "central calculation"
            break;

            // case 2 ********************************* (mark 62)
        case 2:
            quanc8flag = quanc8flag + 1.0;
            break;
            // case 3 ********************************* (mark 60)
        case 3:
            //    ***   stage 6   ***   trouble section
            //    number of function values is about to exceed limit.
            nofin = 2 * nofin;
            levmax = levout;
            quanc8flag = quanc8flag + (b - x0) / (b - a);
            break;
            // case 4 ********************************* (continue mark 70)
        case 4:
            break;
            // default ******************************** (continue mark 70)
        default:
            break;
            // end case section ***********************
        }

        //   ***   stage 7   ***   interval converged
        //   add contributions into running sums.

        quanc8result = quanc8result + qnow;
        quanc8errest = quanc8errest + esterr;
        cor11 = cor11 + qdiff / 1023.0;

        //  locate next interval

        while (nim != 2 * (nim / 2))
        {
            nim = nim / 2;
            lev = lev - 1;
        }
        nim = nim + 1;
        if (lev <= 0) break;  // may exit futher calculation

        //  assemble elements required for the next interval.

        qprev = qright[lev];
        x0 = x[16];
        f0 = f[16];
        for (i = 1; i <= 8; i = i + 1)
        {
            f[2 * i] = fsave[i][lev];
            x[2 * i] = xsave[i][lev];
        }
    }
    //  *** end stage 3 ***   central calculation

    //  ***   stage 8   ***   finalize and return

    quanc8result = quanc8result + cor11;

    //  make sure quanc8errest not less than roundoff level.

    if (quanc8errest == 0.0) return;
    do
    {
        temp = fabs(quanc8result) + quanc8errest;
        quanc8errest = 2.0*quanc8errest;
    } while (temp == fabs(quanc8result));

    return;
}

//# corresponding python code
//def sumg_c(z,u):
    //zplusu=z+u
    //zminusu=z-u
    //log_result = cmath.log((zplusu+1.0)/(zplusu-1.0))
    //plusresult = (1.0 -zplusu*zplusu)*log_result
    //log_result = cmath.log((zminusu+1.0)/(zminusu-1.0))
    //minusresult = (1.0 -zminusu*zminusu)*log_result
    //return(plusresult+minusresult)
    
inline dcomp sumg (dcomp z, dcomp u)
{
    dcomp zplusu=z+u;
    dcomp zminusu=z-u;
    dcomp log_result = log((zplusu + 1.0)/(zplusu - 1.0));
    dcomp plusresult = (1.0 - zplusu * zplusu) * log_result;
    log_result = log((zminusu + 1.0)/(zminusu - 1.0));
    dcomp minusresult = (1.0 - zminusu*zminusu)*log_result;
    return(plusresult+minusresult);
}


//inline dcomp g_c_quad (dcomp A)
//{   __complex128 QA;
    //__real__ QA = A.real();
    //__imag__ QA = A.imag();
    //__complex128 log_result = clogq((QA + 1.0)/(QA - 1.0));
    //__complex128 qout = (1.0 -QA*QA)*(log_result);
    //dcomp out= (dcomp) qout;
    //return out; 
//}

#ifdef QUADMATH
inline dcomp sumg_quad (dcomp z, dcomp u)
{    
    __complex128 zplusu, zminusu, qz,qu;
    __real__ qz = z.real();
    __imag__ qz = z.imag();
    __real__ qu = u.real();
    __imag__ qu = u.imag();
    
     zplusu=qz+qu;
     zminusu=qz-qu;
     __complex128 log_result = clogq((zplusu + 1.0)/(zplusu - 1.0));
     __complex128 plusresult = (1.0 - zplusu * zplusu) * log_result;
    log_result = clogq((zminusu + 1.0)/(zminusu - 1.0));
     __complex128 minusresult = (1.0 - zminusu*zminusu)*log_result;
     dcomp out=(dcomp)(plusresult+minusresult);
    return(out);
}
#endif

dcomp  Chi_DL(double q, double current_w,int i)
{   double gamma;
    double w_q = w_at_q(wi[i], q,alphai[i]);
    if (add_Dopplerwidth_to_classical_DF)
    {
        double k_f=pow(wi[i]*wi[i]*3.0/4.0*pi,1.0/3.0);
        gamma=sqrt(gammai[i]*gammai[i]+ 0.25*q*k_f*q*k_f);
    }    
    else gamma= gammai[i];
    dcomp gamma_c(0.0,gamma+w_q/10000.0);   // make sure not incredibly spikie at large energy loss, which makes it hard for quanc8 to integrate 
    dcomp one_over_chi = Ai[i] * wi[i]*wi[i]/(current_w*current_w - w_q*w_q - gapi[i]*gapi[i]+ current_w*gamma_c );
    dcomp chi = reciprocal(one_over_chi);
    return chi;
}

dcomp  Chi_DL(double A,double q, double current_w, double gamma, double w_p, double alpha,double gap ) //overloaded one for Kaneko etc, when we have not the parameters in an array
{  
    double w_q = w_at_q(w_p, q,alpha);
    dcomp gamma_c(0.0,gamma+w_q/10000.0);   // make sure not incredibly spikie at large energy loss, which makes it hard for quanc8 to integrate 
    dcomp one_over_chi = A * w_p*w_p/(current_w*current_w - w_q*w_q - gap*gap+ current_w*gamma_c );
    dcomp chi = reciprocal(one_over_chi);
    return chi;
}

dcomp  Chi_Lindhard(double q, dcomp omega_c, double  omega0)
{   dcomp u, sumdterms;
    double z;
    double Qrecoil=recoil_energy(q);
    double v_f=pow(omega0*omega0*3.0/4.0*pi,1.0/3.0);
    double sqrt2Q=sqrt(2.0*Qrecoil);
    z =  sqrt2Q/(2*v_f);
    u =  omega_c/ (sqrt2Q *v_f);
   
    if (abs(u) < 500.0*z) // the transition momentum  may need some fine tuning,  MacOs has no quad math.
    {
       sumdterms=sumg(z,u);
    }
    else    
    {   
#ifdef QUADMATH       
        sumdterms=sumg_quad(z,u);
#else
        sumdterms=sumg(z,u);
#endif      
    }
    double prefactor=3.0*omega0*omega0/(2.0*Qrecoil*v_f*v_f);
    dcomp f = 0.5+sumdterms/(8*z);// note d1-d2 differs from Sigmund's f1+if2 (5.155) by 4z
    return (prefactor*f);
}

dcomp Chi_Lindhard_LL(double q, double current_w, int i)
{
    dcomp z1, z2, z3, top, bottom, omega_c, chi, oneoverchi;
    if(!DirectMethod) omega_c=dcomp(current_w, gammai[i]);
    else omega_c = sqrt(current_w*current_w + dcomp(0.0,current_w*gammai[i]));
    dcomp omega_minus=sqrt(omega_c*omega_c-gapi[i]*gapi[i]);
    
    if (q > c_transition*current_w +0.00001)
    {
        z2 =Chi_Lindhard(q,omega_minus,wi[i]);
        if ( Apply_Mermin_Correction)
        {   
            double g_over_w = gammai[i] / current_w;
            z1 = dcomp(1.0, g_over_w);// omega should be unequal 0
            dcomp omega_minus=dcomp(0.0,gapi[i]+1e-10);   //else strange things happen when U=0
            z3 = Chi_Lindhard(q, omega_minus, wi[i]);
            top = z1*z2;
            bottom = dcomp(1.0, 0.0) + dcomp(0, g_over_w)*z2 / z3;
            chi = top / bottom;
        }
        else chi =  z2;
        oneoverchi=reciprocal(chi);  //take A into account for 1/chi
        chi=reciprocal(Ai[i] * oneoverchi);
    }
    else //calculate the equivalent Drude-Lindhard
    {  
        if ( Apply_Mermin_Correction or DirectMethod) chi =  Chi_DL(Ai[i], q, current_w, gammai[i], wi[i], 1.0 ,gapi[i] );
        else  chi =  Chi_DL(Ai[i], q, current_w, 2*gammai[i], wi[i], 1.0,gapi[i] );//plain lindhard, has double the nominal  width
    }
    return (chi);
}

dcomp Chi_Lindhard_LL(double A, double q, double current_w, double w_p, double gamma, double gap)  //overloaded one, for if parameters are not available in arrays
{
    dcomp z1, z2, z3, top, bottom, omega_c, chi, oneoverchi;
    if(!DirectMethod) omega_c=dcomp(current_w, gamma);
    else omega_c = sqrt(current_w*current_w + dcomp(0.0,current_w*gamma));
    dcomp omega_minus=sqrt(omega_c*omega_c-gap*gap);
    
    if (q > c_transition*current_w +0.00001)
    {
        z2 =Chi_Lindhard(q,omega_minus,w_p);
        if ( Apply_Mermin_Correction)
        {   
            double g_over_w = gamma / current_w;
            z1 = dcomp(1.0, g_over_w);// omega should be unequal 0
            dcomp omega_minus=dcomp(0.0,gap+1e-10);   //else strange things happen when U=0
            z3 = Chi_Lindhard(q, omega_minus, w_p);
            top = z1*z2;
            bottom = dcomp(1.0, 0.0) + dcomp(0, g_over_w)*z2 / z3;
            chi = top / bottom;
        }
        else chi =  z2;
        oneoverchi=reciprocal(chi);  //take A into account for 1/chi
        chi=reciprocal(A * oneoverchi);
    }
    else //calculate the equivalent Drude-Lindhard
    {  
        if ( Apply_Mermin_Correction or DirectMethod) chi =  Chi_DL(A, q, current_w, gamma, w_p, 1.0 ,gap );
        else  chi =  Chi_DL(A, q, current_w, 2*gamma, w_p, 1.0,gap);//plain lindhard, has double the nominal  width
    }
    return (chi);
}
dcomp  Chi_Kaneko(double q, dcomp omega_c, double q_mean, double gamma_fudge)
{   
    dcomp Chi,dd,u,z;
    double Qrecoil=recoil_energy(q);
    double sqrt2Q=sqrt(2.0*Qrecoil);
    u=omega_c/ (sqrt2Q* q_mean);
    z = sqrt2Q/ (2.0 * q_mean); 
    dcomp prefactor= 2.0* gamma_fudge * q_mean*q_mean/(q*q*q*dcomp(0,1));  //factor of 2 due to spin, original kaneko means gamma_fudge=0.5
    double relerror=0.0;
    dd= prefactor*( Faddeeva::w(u + z,relerror) -  Faddeeva::w(u - z,relerror));
    Chi =  dd;
    return Chi;
}
dcomp Chi_KanekoDirect (double q, double omega, double gamma, double  alpha, double U,  double gamma_fudge)
{
    dcomp omega_c=sqrt(omega*omega+dcomp(0.0,omega*gamma)-U*U);
    double Q =sqrt(1.0/alpha);
    dcomp Chi = Chi_Kaneko(q, omega_c, Q,  gamma_fudge);
    return Chi;
}

dcomp Chi_Kaneko_Mermin(double q, double omega, double gamma, double  alpha, double U,  double gamma_fudge)
{
    dcomp omega_c,w_c_u,z1, z2, z3, top, bottom;
    double g_over_w;
    double Q =sqrt(1.0/alpha);

    g_over_w = gamma / omega;
    if(DirectMethod) omega_c=sqrt(omega*omega+dcomp(0.0,omega*gamma));
    else  omega_c= dcomp(omega,gamma);
   
    w_c_u=sqrt(omega_c*omega_c-U*U);
 
    z2 = Chi_Kaneko(q, w_c_u, Q,  gamma_fudge) ;
    if(!Apply_Mermin_Correction or  DirectMethod) return z2;
    
    z1 = dcomp(1.0, g_over_w);// omega should be unequal 0 
    z3 = Chi_Kaneko(q, dcomp(0.0, U+1e-10), Q,  gamma_fudge);
    top = z1*z2;
    bottom = dcomp(1.0, 0.0) + dcomp(0, g_over_w)*z2 / z3;
    z1 = top / bottom;   
    return (z1);
}

dcomp diff_chi(double q, double omega, double gamma, double  alpha, double U, double  gamma_fudge)
{  
    dcomp chi_minus,chi_plus;
    double delta_alpha=0.005*alpha;
   
    chi_minus = Chi_Kaneko_Mermin(q, omega, gamma, alpha-delta_alpha,U,gamma_fudge);
    chi_plus  = Chi_Kaneko_Mermin(q, omega, gamma, alpha+delta_alpha,U, gamma_fudge);
 
    return (chi_plus-chi_minus)/(2.0*delta_alpha);
}
dcomp diff2_chi(double q, double omega, double gamma, double  alpha,double U,  double  gamma_fudge)
{  
    dcomp chi_minus,chi_plus;
    double delta_alpha=0.005*alpha;
    chi_minus = diff_chi(q, omega, gamma, alpha-delta_alpha,U,gamma_fudge);
    chi_plus  = diff_chi(q, omega, gamma, alpha+delta_alpha,U,gamma_fudge);
    return (chi_plus-chi_minus)/(2.0*delta_alpha);
}
dcomp diff3_chi(double q, double omega, double gamma, double  alpha, double U, double  gamma_fudge)
{  
    dcomp chi_minus,chi_plus;
    double delta_alpha=0.005*alpha;
    chi_minus = diff2_chi(q, omega, gamma,  alpha-delta_alpha, U, gamma_fudge);
    chi_plus  = diff2_chi(q, omega, gamma, alpha+delta_alpha, U, gamma_fudge);
    return (chi_plus-chi_minus)/(2.0*delta_alpha);
}





dcomp  calculate_chi_AA_LL(double q)  // Archubi-Arista suggested use of  Levine Louie to calculate Kaneko with gap
                                        // revert to plane Kaneko if energy edge  (=U)is 0, calculated at w_global
{
    dcomp chi, sum_chi, sum_oneoverchi,  omega_c;
    double alpha, Volume, Volumefraction,gamma_fudge;
    double ratio,U;

  
 
    sum_chi=dcomp(0.0,0.0);
    sum_oneoverchi=dcomp(0.0,0.0);
    for (int i = 0; i < MAXKANEKO; i++)
    { 
        if ((fabs(N_Kaneko[i]) > 1E-10)) // && (w_global> Edge_ArchubiKaneko[i]) )
        {  
            U =  Edge_ArchubiKaneko[i];
            omega_c = dcomp( w_global, Width_Kaneko[i]);         
            alpha=1.0/(Q_Kaneko[i]*Q_Kaneko[i]);
            gamma_fudge=gamma_Kaneko[i];
            ratio=w_pl_l[i]*w_pl_l[i]* A_l_Kaneko[i]/(A_0_Kaneko[i]*w_pl_0[i]*w_pl_0[i]) ;    
            
            Volume=4.0*pi*N_Kaneko[i]/(w_pl_l[i]*w_pl_l[i]); 
            Volumefraction=Volume* UnitCellDensity; //gamma fudge already incorportaed via w_pl[i]
            if(Volumefraction > 1.0)
            {
                printf("warning, unphysical density Kaneko component  %i volumefraction: %6.4f\n",i, Volumefraction);
            }
            
             if (q < c_transition*w_global+0.00001)
            {
                if ( Apply_Mermin_Correction or DirectMethod )  chi = Chi_DL(Ai[i], q, w_global, Width_Kaneko[i], w_pl_l[i], 1.0,U );
                else chi = Chi_DL(Ai[i],q, w_global, 2* Width_Kaneko[i], w_pl_l[i], 1.0,U ); //plain RPA , has double the nominal  width
            }
            else if( l_Kaneko[i]==0) chi = Chi_Kaneko_Mermin(q, w_global, Width_Kaneko[i], alpha, U, gamma_fudge);
            else if (l_Kaneko[i]==1) chi=  dcomp(1.0,0.0)-ratio*diff_chi(q, w_global, Width_Kaneko[i],   alpha,U, gamma_fudge);
            else if (l_Kaneko[i]==2) chi=  dcomp(1.0,0.0)+ratio*diff2_chi(q, w_global, Width_Kaneko[i],   alpha, U,  gamma_fudge);
            else if (l_Kaneko[i]==3) chi=  dcomp(1.0,0)-ratio*diff3_chi(q, w_global, Width_Kaneko[i],   alpha,U,  gamma_fudge);
            
            if(AddELF)  sum_oneoverchi+=reciprocal(chi*Volumefraction);
            else sum_chi+=Volumefraction*chi;
           
        }
    }
    if(AddELF) sum_chi=reciprocal(sum_oneoverchi);
   
    return sum_chi;
}

double  calculate_Loss_AA_LL(double q)  // Archubi-Arista suggested use of  Levine Louie to calculate Kaneko with gap
                                        // revert to plane Kaneko if energy edge  (=U)is 0, calculated at w_global
                                    
{
    dcomp eps=calculate_chi_AA_LL(q)+dcomp(1.0,0.0);
    dcomp oneovereps= dcomp(1.0,0.0)/eps;
    double loss= -oneovereps.imag();
    return loss;
}                                        

dcomp  Chi_Drude(double q, double current_w,int i)
{   dcomp chi, gamma_c; 
    double w_q;
    w_q = w_at_q(wi[i], q,alphai[i]);
    gamma_c=dcomp(0,gammai[i]+w_q/10000.0);   // make sure not incredibly spikie at large energy loss, which may  make it hard for quanc8 to integrate  
    chi =  -1.0 / (current_w*current_w - w_q*w_q - gapi[i]*gapi[i] + current_w*gamma_c);// U tries to mimick band gap here, just as in LL
    return Ai[i]*chi;
}




dcomp  F_Orosco(dcomp z)
{   double modz2=abs(z)*abs(z);
    return dcomp(0,pi)*Faddeeva::w(z,0.0)+exp(-z*z)*(log(z)+log(-conj(z)/modz2)-dcomp(0,pi));
}

dcomp  ChiBB(double q, double omega, double  gamma, double w, double sigma_Gauss, double U,double  w_p_square)  // Ai[i] corresponds to plasmon frequency square
{
    double w_at_q, Qrecoil;
    dcomp   Chi;
   
    
    Qrecoil = recoil_energy(q);   
                  
    if (Full_dispersion)
    {  
        double k_f=pow(w*w*3.0/4.0*pi,1.0/3.0);
        w_at_q = sqrt( w*w + 2.0 *k_f*k_f /3.0 *Qrecoil +Qrecoil * Qrecoil);  //alpha always taken to be 1 in BB
    }
    else
    {
        w_at_q = w +  Qrecoil;
    }
    gamma=gamma+Qrecoil/1000.0;   // make sure not incredibly spikie at large energy loss, which makes it hard for quanc8 to integrate 
    dcomp omega_c=sqrt(omega*omega+dcomp(0.0,omega*gamma)-U*U);
//// next lines are BB
    if(!Orosco_Coimbra_way)
    {
        dcomp  term1,term2;
        term1=( omega_c-w_at_q)/(sqrt(2.0)*sigma_Gauss);
        term2=( omega_c+w_at_q)/(sqrt(2.0)*sigma_Gauss);

        Chi= dcomp(0.0,1.0)*sqrt(pi)* w_p_square/(sqrt(8)*sigma_Gauss*omega_c)*(Faddeeva::w(term1,0.0)+Faddeeva::w(term2,0.0));
        return Chi;
    }
////end BB 
//now Orosco and Coimbra implementation (eq. 13-15)should be the same as B&B for positive omega
    dcomp z_plus,z_minus;// corr_plus,corr_minus;//, I_plus, I_minus;
    //double alpha_prime= sqrt(omega/2.0)*sqrt(sqrt(omega*omega+gamma*gamma)+w_at_q);
    //double alpha_doubleprime = sqrt(omega/2.0)*sqrt(sqrt(omega*omega+gamma*gamma)-omega)+ 1.0/w_at_q;
    //omega_c=dcomp(alpha_prime,alpha_doubleprime);
    omega_c+=dcomp(0,1.0/w_at_q);
    z_plus= ( omega_c-w_at_q)/(sqrt(2.0)*sigma_Gauss);  //now omega_c has different signs comp. to B&B, not w_at_q!!!
    z_minus= (-omega_c-w_at_q)/(sqrt(2.0)*sigma_Gauss); //now omega_c has different signs, not w_at_q
    //now Orosco and Coimbra implementation (eq. 13-15)should be the same as B&B for positive omega
   // Chi= ( w_p_square/(sqrt(8*pi)*sigma_Gauss*omega_c))*(F_Orosco(z_plus)-F_Orosco(z_minus) ) ;  //multiplied by v_p^2 (Osc. strength A)  in function call
    //now Orosco and Coimbra implementation (eq. 13-13)
    double Chi0= -4*sqrt(pi)*Faddeeva::Dawson(-w_at_q/(sqrt(2.0)*sigma_Gauss));
    double A= w_p_square/(w_at_q*w_at_q);
    Chi=A*(F_Orosco(z_plus)+F_Orosco(z_minus) )/Chi0;
    return Chi;
}


dcomp Chi_BrendelBormann(double q,double current_w,int i)
{   
    double w_q, gamma;
    dcomp   Chi;
    w_q = w_at_q(current_w,q,alphai[i]);
             
    
    gamma=gammai[i]+w_q/1000.0;   // make sure not incredibly spikie at large energy loss, which makes it hard for quanc8 to integrate 
    dcomp omega_c=sqrt(current_w*current_w+dcomp(0.0,current_w*gamma)-gapi[i]*gapi[i]);
//// next lines are BB
    if(!Orosco_Coimbra_way)
    {
        dcomp  term1,term2;
        term1=( omega_c-w_q)/(sqrt(2.0)*sigmai[i]);
        term2=( omega_c+w_q)/(sqrt(2.0)*sigmai[i]);

        Chi= dcomp(0.0,1.0)*sqrt(pi)* Ai[i]/(sqrt(8)*sigmai[i]*omega_c)*(Faddeeva::w(term1,0.0)+Faddeeva::w(term2,0.0));  //Ai[i] corresponds to plasmon energy square
        return Chi;
    }
////end BB 
//now Orosco and Coimbra implementation (eq. 13-15)should be the same as B&B for positive omega
    dcomp z_plus,z_minus;// corr_plus,corr_minus;//, I_plus, I_minus;
    //double alpha_prime= sqrt(omega/2.0)*sqrt(sqrt(omega*omega+gamma*gamma)+w_at_q);
    //double alpha_doubleprime = sqrt(omega/2.0)*sqrt(sqrt(omega*omega+gamma*gamma)-omega)+ 1.0/w_at_q;
    //omega_c=dcomp(alpha_prime,alpha_doubleprime);
    omega_c+=dcomp(0,1.0/w_q);
    z_plus= ( omega_c-w_q)/(sqrt(2.0)*sigmai[i]);  //now omega_c has different signs comp. to B&B, not w_at_q!!!
    z_minus= (-omega_c-w_q)/(sqrt(2.0)*sigmai[i]); //now omega_c has different signs, not w_at_q
    //now Orosco and Coimbra implementation (eq. 13-15)should be the same as B&B for positive omega
   // Chi= ( w_p_square/(sqrt(8*pi)*sigma_Gauss*omega_c))*(F_Orosco(z_plus)-F_Orosco(z_minus) ) ;  //multiplied by v_p^2 (Osc. strength A)  in function call
    //now Orosco and Coimbra implementation (eq. 13-13)
    double Chi0= -4*sqrt(pi)*Faddeeva::Dawson(-w_q/(sqrt(2.0)*sigmai[i]));
    double A= Ai[i]/(w_q*w_q);
    Chi=A*(F_Orosco(z_plus)+F_Orosco(z_minus) )/Chi0;
    return Chi;
}
    
    //python code
    //def BrendelBormann(q,omega,gamma,Amp,U, w0, sigma_gauss): #https://en.wikipedia.org/wiki/Brendel%E2%80%93Bormann_oscillator_model
   
    //w_c_U=  cmath.sqrt(omega*omega-U*U+complex(0,omega*gamma) ) 
    //wq=w0+RecoilEnergy(q) 
 
    //term1=( w_c_U-w0)/(math.sqrt(2.0)*sigma_gauss)
    //term2=( w_c_U+w0)/(math.sqrt(2.0)*sigma_gauss)
    //chi=complex(0,1) * Amp* math.sqrt(math.pi)/(math.sqrt(8)*sigma_gauss*w_c_U)*(special.wofz(term1)+special.wofz(term2))
    //return 1.0 + chi


dcomp Chi_Vlasov(double q, dcomp omega_c, double  Q_Vlasov) 
{   
    double Qrecoil=recoil_energy(q);
    double sqrt2Q=sqrt(2.0*Qrecoil);
    double z =  sqrt2Q/(2*Q_Vlasov);
    dcomp u =  omega_c/ (sqrt2Q *Q_Vlasov);
    dcomp prefactor=Q_Vlasov*Q_Vlasov/(q*q*q*dcomp(0,1));
    dcomp prefactorVlasov=2.0*prefactor*z;

    double relerror=0.0;
    dcomp w_prime_u= -2.0*u* Faddeeva::w(u,relerror)+ dcomp(0.0,2.0)/sqrt(pi);
    dcomp chi= prefactorVlasov * w_prime_u;
    return chi;
}   
    
dcomp Chi_Vlasov_M(double q, double current_w,int i)
{
    dcomp z1, z2, z3, top, bottom, omega_c, chi;
    if(DirectMethod) omega_c=sqrt(current_w*current_w+dcomp(0.0,current_w*gammai[i]));
    else  omega_c=dcomp(current_w, gammai[i]);
    dcomp omega_minus=sqrt(omega_c*omega_c - gapi[i]*gapi[i]);
    if (q > c_transition*current_w+0.00001)
    {
        z2 =Chi_Vlasov(q,omega_minus,Q_Vlasov[i]);
        if ( Apply_Mermin_Correction)
        {   
            double g_over_w = gammai[i] / current_w;
            z1 = dcomp(1.0, g_over_w);// omega should be unequal 0
            dcomp omega_minus=dcomp(0.0,gapi[i]+1e-10);   //else strange things happen when U=0
            z3 = Chi_Vlasov(q, omega_minus, Q_Vlasov[i]);
            top = z1*z2;
            bottom = dcomp(1.0, 0.0) + dcomp(0, g_over_w)*z2 / z3;
            chi = top / bottom;
        }
        else  chi=  z2;
    }
    else //calculate the equivalent Drude-Lindhard
    {   
        double omega0 = sqrt(2*pow(Q_Vlasov[i],3)/sqrt(pi));
        if ( Apply_Mermin_Correction or DirectMethod) chi =  Chi_DL(Ai[i], q, current_w, gammai[i], omega0, 1,gapi[i] );
        else   chi  = Chi_DL(Ai[i], q, current_w, 2*gammai[i], omega0, 1,gapi[i] );//plain lindhard, has double the nominal  width
    }
    
    return chi;
}




dcomp F_TL_an( dcomp alpha, dcomp beta, dcomp gamma, double Egap)
{
    dcomp denom = alpha*(alpha*alpha-beta*beta)*(alpha*alpha-gamma*gamma);
    dcomp term1 = (Egap+alpha)*(Egap+alpha)*log(Egap+alpha);
    dcomp term2 = (Egap-alpha)*(Egap-alpha)*log(Egap-alpha);
    return (term1 - term2)/denom;
}
dcomp Chi_TL_an(double E, double  C, double E0,  double Egap, double a_TL_analytic)
{   dcomp term1,term2,term3;
    dcomp E0c= dcomp(E0,0);
    dcomp b = dcomp(E,a_TL_analytic);
  
    dcomp d=sqrt(E0c*E0c - C*C/4.0)-dcomp(0,C/2.0);
    
    dcomp dprime= conj(d);// is this allowed for an analytical function?
    term1= F_TL_an(b,d,dprime,Egap);
    term2= F_TL_an(d,dprime,b, Egap);
    term3=F_TL_an(dprime,b,d,Egap);
    //dcomp Chi=  C/pi*( F_TL_an(b,d,dprime,Egap)+ F_TL_an(d,dprime,b, Egap)+ F_TL_an(dprime,b,d,Egap)); //Factor A*E0 added in function call
    dcomp Chi=  C/pi*( term1+term2+term3);
    //printf("d %6.4f+%6.4fi, chi  %6.4f+%6.4fi \n",d.real(),d.imag(),Chi.real(),Chi.imag());
   // printf("Energy %6.4f,term1 %6.4f+%6.4fi, term2  %6.4f+%6.4fi term3  %6.4f+%6.4fi, Chi %6.4f+%6.4fi  \n",E,term1.real(),term1.imag(),term2.real(),term2.imag(),term3.real(),term3.imag(),Chi.real(),Chi.imag());
    return Chi;
}
double TaucNormalisation(double A,double  C, double E0, double E0_thisq,  double Egap,  double a_TL_analytic, double q, double A_Mermin, double wp_Mermin)  // can run without copyP_to_Vars
{   int nstep=50000, i;
    double delta =0.02;
    dcomp chi_Mermin;
    double eps2=0.0;// initialized to avoid warning
    double x,TaucFactor,DrudeLimit;
    double stepsize= C/200;
    if (modelTauc_Mermin) 
    {
        DrudeLimit=A_Mermin*wp_Mermin*wp_Mermin/(4.0*pi);
        E0_thisq= wp_Mermin+q*q/2; // approx peak  ridge position
    }
    else 
    { 
        DrudeLimit=A*E0/(4.0*pi);
    }
    
    double sum=0.0; 
    double Elower= E0_thisq-0.5*stepsize;
    double Eupper= E0_thisq+0.5*stepsize;
    if (add_Dopplerwidth_to_classical_DF and !modelTauc_Mermin)
    {
        double k_f=pow( E0* E0*3.0/4.0*pi,1.0/3.0);
       
        C=sqrt(C*C+ 0.25*q*k_f*q*k_f);
    }    
    for ( i=0; i < nstep/2; i++)
    {
        if( Elower > Egap)
        {  
                if (modelTaucLorentz)
                {
                     eps2 = A * E0* C*pow((Elower-Egap),2)/((pow((Elower*Elower-E0_thisq*E0_thisq),2)+C*C*Elower*Elower))/Elower;
                }
                else if (modelTL_an)
                { 
                    dcomp chi_TL= A *E0* Chi_TL_an(Elower, C, E0_thisq, Egap,a_TL_analytic);
                    eps2=chi_TL.imag();
                }
                else  if (modelTauc_Mermin)
                { 
                    chi_Mermin = Chi_Lindhard_LL( A_Mermin, q, Elower, wp_Mermin, C, 0.0);
                    TaucFactor = pow((Elower-Egap)/Elower,2);
                    eps2 = TaucFactor*chi_Mermin.imag();
                }
                sum += Elower* eps2*stepsize*(1.0+delta/2.0); // because the next  step is larger by  1+ delta, effective width this been is (last step+next step)/2
        }
        if (modelTaucLorentz) 
        {
            eps2 = A * E0* C*pow((Eupper-Egap),2)/((pow((Eupper*Eupper-E0_thisq*E0_thisq),2)+C*C*Eupper*Eupper))/Eupper;
        }
        else if (modelTL_an) 
        { 
            dcomp chi_TL= A *E0* Chi_TL_an(Eupper, C, E0_thisq, Egap,a_TL_analytic);
            eps2=chi_TL.imag();
        }
        else if (modelTauc_Mermin)
        { 
            chi_Mermin = Chi_Lindhard_LL( A_Mermin, q, Eupper, wp_Mermin, C, 0.0);
            TaucFactor = pow((Eupper-Egap)/Eupper,2);
            eps2 = TaucFactor*chi_Mermin.imag();
        }
        x=Eupper*eps2*stepsize*(1.0+delta/2.0);
        sum += x;
        stepsize=stepsize*(1.0 + delta);
        Eupper += stepsize;
        Elower -= stepsize;
        
        if(x/sum < 1e-10) break;
       
    } 
    double normvalue= sum/ (2.0 * pi * pi);
  //  printf("normvalue %6.4f, E0 %6.4f, Drude_Limit= %6.4f\n",normvalue, E0, DrudeLimit);
    return normvalue/DrudeLimit; // Normalisation rel. to Drude should approach 1 for large q when Tauc becomes Drude
}

dcomp  ChiTaucLorentz(double E, double  C,  double E0q,  double Egap )  // E00 E0 at q=0 E0q=E0 at q
{    
    double Chi1,Chi2;
    
    double a_ln=(Egap*Egap-E0q*E0q)*E*E+Egap*Egap*C*C-E0q*E0q*(E0q*E0q+3*Egap*Egap);
    double a_atan=(E*E-E0q*E0q)*(E0q*E0q+Egap*Egap)+Egap*Egap*C*C;
    double gamma=sqrt(E0q*E0q-C*C/2);
    double alpha=sqrt(4*E0q*E0q-C*C);
    double zeta4=(E*E-gamma*gamma)*(E*E-gamma*gamma)+alpha*alpha*C*C/4;
    
    if(E> Egap)
    {
        Chi2=C*pow((E-Egap),2)/((pow((E*E-E0q*E0q),2)+C*C*E*E)*E);
       
    }
    else
    {   
        Chi2=0.0;
    }
        
    Chi1  = C/(pi*zeta4)*a_ln/(2*alpha*E0q*E0q)*log((E0q*E0q+Egap*Egap+alpha*Egap)/(E0q*E0q+Egap*Egap-alpha*Egap));
    Chi1 -= a_atan/(pi*zeta4*E0q*E0q)*(pi-atan((2*Egap+alpha)/C)+ atan((-2*Egap+alpha)/C));
    Chi1 += 2/(pi*zeta4*alpha)*Egap*(E*E-gamma*gamma)*(pi+2*atan(2*(gamma*gamma-Egap*Egap)/(alpha*C)));
    Chi1 -= C/(pi*zeta4)*(E*E+Egap*Egap)/E*log(fabs(E-Egap)/(E+Egap));
    Chi1 += 2*C/(pi*zeta4)*Egap*log((fabs(E-Egap)*(E+Egap))/sqrt(pow((E0q*E0q-Egap*Egap),2)+Egap*Egap*C*C));
  //  printf("in Chi TaucLorentz  chi1 %6.5f  chi2 %6.5f\n", Chi1,Chi2);
    return dcomp(Chi1,Chi2);
}
double GetTaucScalingFactor(int Component, double q)
{
    int iq;
    double renorm, before,after;
    for (int i=1; i < NTaucscaling; i++)
    {       
                if (Tauc_Scaling_q[i] > q)
                { 
                    iq=i-1;
                    break;
                }
                else iq=i;
    }  
    if (iq == NTaucscaling-1)  renorm=Tauc_ScalingFactor[Component][0];  // should be close to the Drude limit
    else
    {   
        before = Tauc_ScalingFactor[Component][iq];
        after = Tauc_ScalingFactor[Component][iq + 1];
        double position_in_bin=(q-Tauc_Scaling_q[iq])/(Tauc_Scaling_q[iq+1]-Tauc_Scaling_q[iq]);// in the range of 0 to 1
        double TaucFactor_at_q = before*(1.0 - position_in_bin) + position_in_bin*after;
        renorm=Tauc_ScalingFactor[Component][0]/TaucFactor_at_q;// so renorm is 1 when q=0
    }
    return(renorm);
}


dcomp Chi_TL(double q,double current_w, int i) // also used for TL_an_Sum give local current_w a different name
{          
    dcomp ChiTL;
    double Gamma_this_q;
    // first the shift due to dispersion
    double Qrecoil, currentE0, renorm;
    renorm=GetTaucScalingFactor(i,q);
    Qrecoil = recoil_energy(q);
         
         
    if (add_Dopplerwidth_to_classical_DF)
    {   double a=wi[i];
        if (a< 17/Hartree) a=17/Hartree;  // make it at least the plasmon energy
        double k_f=pow(a* a*3.0/4.0*pi,1.0/3.0); 
        Gamma_this_q=sqrt(gammai[i]*gammai[i]+ 0.25*q*k_f*q*k_f);  //there is a factor of 0.25 herethat I can not justify
    }  
    else
    {
        Gamma_this_q=gammai[i];
    }    
    if (delayed_dispersion)
    {
        double E0_without_gap= sqrt(wi[i]*wi[i]-TaucGap[i]*TaucGap[i]);
        double E0_without_gap_this_q= E0_without_gap+ Qrecoil;
        currentE0=sqrt(E0_without_gap_this_q*E0_without_gap_this_q+TaucGap[i]*TaucGap[i]);
        //Fermi("E0_without_gap %6.4f, E0_without_gap_this_q %6.4f, currentE0  %6.4f, q %6.4f\n",E0_without_gap,E0_without_gap_this_q, currentE0,q);
    }
    else currentE0=wi[i]+Qrecoil;
    if(modelTaucLorentz)ChiTL = Ai[i]*wi[i]*renorm*ChiTaucLorentz(current_w,Gamma_this_q,currentE0, TaucGap[i]);
    if(modelTL_an) ChiTL = Ai[i]*wi[i]* renorm*Chi_TL_an(current_w,Gamma_this_q,currentE0, TaucGap[i], aTLan);       
    
    return ChiTL;
}

double eps2_TaucMermin(double q, double current_w, int i)
    {   
    if(current_w > TaucGap[i])
    {  
       double TaucFactor = pow((current_w-TaucGap[i])/current_w,2);
       double renorm=GetTaucScalingFactor(i,q);
       dcomp chi_Mermin = Chi_Lindhard_LL(q, current_w,i);
       return (renorm*TaucFactor*chi_Mermin.imag()) ;
    }
    else return (0.0);
}

dcomp Chi_TaucMermin(double q, double current_w, int i)
    {   
    dcomp chi_Mermin, chi; 
    
    double TaucFactor = pow((current_w-TaucGap[i])/current_w,2);
    double renorm=GetTaucScalingFactor(i,q);
    if(current_w > TaucGap[i])
    {
        chi_Mermin = Chi_Lindhard_LL(q, current_w,i);
        chi= renorm * TaucFactor*chi_Mermin; // this also scales the real part, which is wrong!
     }   
     
    else chi= dcomp(0,0);
// the real part will deviate from the one calculated here due to fact that for non-neglibeble density effect 
// the real part is affected by the multiplicATION WITH THE TaucFactor. Hence recalculate via KK for non-negligible densities
    double chi1=0.0;
    if(current_w < maxEnergyDensityEffect)
    {   
        int KKSteps=150;
        double current_StepSize=gammai[i]/25;
        double w_below=current_w-current_StepSize;
        double w_above=current_w+current_StepSize;
        for(int j=1; j< KKSteps; j++)
        {   
            if(w_below >TaucGap[i]) chi1 -= current_StepSize* w_below*2/ pi*eps2_TaucMermin(q, w_below, i)/(current_w*current_w - w_below*w_below);// wooton 6.33
            chi1 -= current_StepSize*w_above*2/pi*eps2_TaucMermin(q, w_above, i)/(current_w*current_w - w_above*w_above);
            current_StepSize=current_StepSize*1.025;
            w_below-=current_StepSize;
            w_above+=current_StepSize;
        }
        chi.real(chi1);
    }
    
    return chi;
}

 //if(omega < maxEnergyDensityEffect)
    //{
      
        //double current_StepSize=Edgei_GOS[GOSi]/200;// assume these are core levels  energy resolution required of the order of 0.5% of BE?
        //double w_below=omega-current_StepSize;
        //double w_above=omega+current_StepSize;
        //int KKWidth=500;
        //for(int i=1; i< KKWidth; i++)
        //{
            //if(w_below >StepSize/2)
            //{   
                //double Contribution_below=current_StepSize*W_plasmon_GOS_square*GOSx(n_i_GOS[GOSi], l_i_GOS[GOSi], Zi_GOS[GOSi], Edgei_GOS[GOSi], q, w_below)/(omega*omega- w_below*w_below);
                //chi_real-= Contribution_below;
            //}
            //double Contribution_top=current_StepSize*W_plasmon_GOS_square*GOSx(n_i_GOS[GOSi], l_i_GOS[GOSi], Zi_GOS[GOSi], Edgei_GOS[GOSi], q, w_above)/(omega*omega- w_above*w_above);
            //chi_real-= Contribution_top;
            //current_StepSize=current_StepSize*1.01;
            //w_below-=current_StepSize;
            //w_above+=current_StepSize;
        //}



dcomp Chi_ForouhiBloomer(double q, double current_w, int i)
{   if(q!=0) printf("Forouhi-Bloomer implemented only for q=0\n");
    dcomp sumeps; double n,n_infty,k,Esquare,B,C,B0,C0,Q,Eg ;
    sumeps=dcomp(1.0,0.0);
    k=0.0;
    n_infty=sqrt(epsbkg);
    n=n_infty;
   
    Esquare=current_w*current_w;
    B=wi[i]; 
    C=gammai[i];
    Eg=gapi[i];
    k += Ai[i]*(current_w-Eg)*(current_w-Eg)/(Esquare-B*current_w+C);
    Q=0.5*sqrt(4*C-B*B);
    B0=(Ai[i]/Q)*(-B*B/2+Eg*B-Eg*Eg+C);
    C0=(Ai[i]/Q)*((Eg*Eg+C)*B/2-2*Eg*C);
    n += (B0*current_w+C0)/(Esquare-B*current_w+C);
     
    return dcomp(n*n-k*k- 1.0,2*n*k);
}


dcomp calculate_eps_Belkacem(double q, double omega,double A_i, double w_i, double gamma_i)

{     
    dcomp eps, chi, gamma, fraction1,fraction2;
    double W_plasmon_Belkacem;
    double omega_k,C1,M1,currentPoisson, ratio;//, w_i_dens;
    int n,n_min, n_max;
    
     if(q> c_transition)
    {
         omega_k = recoil_energy(q);
    }
    else  // use a minimum value for q
          //for small q all intensity in n=1 component and then omega_k cancels
    {
        omega_k=c_transition*c_transition/2.0;
    }
    ratio=omega_k/w_i;
    W_plasmon_Belkacem=PL_GOS*sqrt(A_i); //PL_GOS plasmon energy of 1 electron per unit cell. Ai_Belkacem is number of oscillators per unit cell
  //  w_i_dens = sqrt(w_i*w_i+W_plasmon_Belkacem*W_plasmon_Belkacem);//density effect may affect os. frequency  
    n_min=round(omega/w_i)-16;  // was 60
    if (n_min < 1) n_min=1;
    n_max=round(omega/w_i) + 16;  // was 60 divided w_i_dens
    double fraction_covered= Poisson( 0, ratio);
    chi=dcomp(0.0,0.0);
    gamma = dcomp(0.0, gamma_i);//+ omega/10000);
    C1=pow(W_plasmon_Belkacem,2)/(2.0*omega_k);
    for(n=n_min; n< n_max; n++)
    {    
         currentPoisson= Poisson( n, ratio);
         fraction_covered += currentPoisson;
         M1=C1* currentPoisson;
         fraction1= 1.0/(n*w_i - omega - gamma);
         fraction2= 1.0/(n*w_i + omega + gamma);
         chi += M1*(fraction1 + fraction2);
    }
  
    if (DebugMode) printf("Belkacem, fraction covered %6.4f\n",fraction_covered);
    
    eps=dcomp(1.0,0)+chi;
    return eps;
}

 
dcomp calculate_chi_Belkacem(double q)  //calculate at w_global
{   
    dcomp  chi, sum_chi,sum_one_over_eps,one_over_eps;
    sum_chi=dcomp(0.0,0.0);
    sum_one_over_eps=dcomp(1.0,0.0);
    //in version 10 we applied Merminization to Belkacem. removed, does not make sense to me anymore
    for (int i = 0; i < MAXBELKACEM; i++)
    { 
        if (Ai_Belkacem[i] > 1E-10)
        {
             chi=  calculate_eps_Belkacem(q, w_global,Ai_Belkacem[i], wi_Belkacem[i],gammai_Belkacem[i])-dcomp(1.0,0);
             if(AddELF)
             {
                one_over_eps=dcomp(1.0,0.0)/(dcomp(1.0,0.0)+chi);
                sum_one_over_eps+=one_over_eps-dcomp(1.0,0.0);  //add  one over chi's
             }
             else
             {
                 sum_chi+=chi;
             }
        }    
    }
    if(AddELF)
    { 
        sum_chi= dcomp(1.0,0.0)/sum_one_over_eps-dcomp(1.0,0);
    }
    return(sum_chi);

 } 
   
double calculate_loss_Belkacem(double q)  //calculate at w_global
{   dcomp chi_Belk,oneovereps;      
    chi_Belk=calculate_chi_Belkacem( q);
    oneovereps=dcomp(1.0,0.0)/(chi_Belk+dcomp(1.0,0.0));
    return -oneovereps.imag();
}

dcomp calculate_eps_osc(double q)
{
    dcomp chi, oneoverchi;
    dcomp chi_sum = dcomp(0.0,0.0);
    dcomp oneoverchi_sum = dcomp(0.0,0.0);
    for (int i = 0; i < MAXOSC; i++)
    {
        if (fabs(Ai[i]) >1e-90)
        {
        if (modelDrude)              chi = Chi_Drude(q,w_global,i);
        else if(modelDL)             chi = Chi_DL(q,w_global,i);
        else if(modelMerminLL)       chi = Chi_Lindhard_LL(q,w_global,i);
        else if(modelVlasov)         chi = Chi_Vlasov_M(q,w_global,i);
        else if(modelTaucLorentz)    chi = Chi_TL(q,w_global,i);
        else if(modelTL_an)          chi = Chi_TL(q,w_global,i);
        else if(modelTauc_Mermin)    chi = Chi_TaucMermin(q,w_global,i); 
        else if(modelForouhiBloomer) chi = Chi_ForouhiBloomer(q,w_global,i);
        else if(modelBrendelBormann) chi = Chi_BrendelBormann(q,w_global,i);
        if(AddELF)
            {  oneoverchi= reciprocal(chi);
               oneoverchi *= EdgeFactor(w_global, Ethi[i], deltai[i]);	//// MODIF 04/08 LK: Multiply by F for edges
               oneoverchi_sum +=oneoverchi;
            }
            else  chi_sum+=chi;
        }
   }
   if(AddELF) chi_sum=reciprocal(oneoverchi_sum);
   return chi_sum + dcomp(1.0,0.0);
}

dcomp get_eps_component(int i,double q)  // needs checking esp for Drude
{   
    dcomp eps;
   
    if (modelDrude)
    {
        eps = dcomp(1.0,0)+  Chi_Drude(q,w_global,i);
    }
    else if (modelDL)
    {   
        eps =dcomp(1.0,0)+ Chi_DL(q,w_global,i);
    }
    else if (modelMerminLL)
    {
        eps = dcomp(1.0,0)+Chi_Lindhard_LL(q,w_global,i);
    }
    else
    { printf("eps per component not implemented for this model\n");
        
    }
  
  
    return eps;
}


dcomp calculate_chi_GOS_dens(double q, double omega, int GOSi) //// calculates eps1 eps2 from GOS including (optional) density effect
{  
    double chi_real,chi_imag;
    double W_plasmon_GOS_square;
    double  currentscaling;
    int iq;
    if (ApplySumRuleToGOS)
    {
        for (int i=1; i < NGOSscaling; i++)
        {       
                    if (GOS_Scaling_q[i] > q) {iq=i-1;break;}
                    else iq=i;
        }  
        //double q_remainder = fmod(q, GOSScalingStepSize);
        if (iq == NGOSscaling-1)  currentscaling=GOS_ScalingFactor[GOSi][iq];
        else
        {
            double position_in_bin=(q-GOS_Scaling_q[iq])/(GOS_Scaling_q[iq+1]-GOS_Scaling_q[iq]);// in the range of 0 to 1
            double before = GOS_ScalingFactor[GOSi][iq];
            double after  = GOS_ScalingFactor[GOSi][iq + 1];
            currentscaling = before*(1.0 - position_in_bin) + position_in_bin*after;//currentscaling so we obey sum rule
        }
    }
    else currentscaling=1.0;
   

    W_plasmon_GOS_square=PL_GOS*PL_GOS*AiGOS[GOSi]*currentscaling; //PL_GOS plasmon energy of 1 electron per unit cell. Ai is number of GOS electrons per unit cell of comp. GOSi
                                                                   // contains the scaling factor
    chi_real=0.0;
    chi_imag=pi*W_plasmon_GOS_square/(2.0*omega)*GOSx(n_i_GOS[GOSi], l_i_GOS[GOSi], Zi_GOS[GOSi], Edgei_GOS[GOSi], q, omega);

    if(omega < maxEnergyDensityEffect)
    {
        double current_StepSize=Edgei_GOS[GOSi]/200;// assume these are core levels  energy resolution required of the order of 0.5% of BE?
        double w_below=omega-current_StepSize;
        double w_above=omega+current_StepSize;
        int KKWidth=500;
        for(int i=1; i< KKWidth; i++)
        {
            if(w_below >StepSize/2)
            {   
                double Contribution_below=current_StepSize*W_plasmon_GOS_square*GOSx(n_i_GOS[GOSi], l_i_GOS[GOSi], Zi_GOS[GOSi], Edgei_GOS[GOSi], q, w_below)/(omega*omega- w_below*w_below);
                chi_real-= Contribution_below;
            }
            double Contribution_top=current_StepSize*W_plasmon_GOS_square*GOSx(n_i_GOS[GOSi], l_i_GOS[GOSi], Zi_GOS[GOSi], Edgei_GOS[GOSi], q, w_above)/(omega*omega- w_above*w_above);
            chi_real-= Contribution_top;
            current_StepSize=current_StepSize*1.01;
            w_below-=current_StepSize;
            w_above+=current_StepSize;
        }
        
    }
    else //too small to bother
    { 
        chi_real=0.0;
    }
    return dcomp(chi_real,chi_imag);
}

double calculate_loss_GOS_dens(double q)
{ // calculate contribution GOS to loss function including density effect
    dcomp oneovereps, sum_chi, chiGOS;
    double loss=0.0;
    sum_chi=dcomp(0.0,0.0);
        for (int i = 0; i < MAXGOS; i++)
        {
            if (AiGOS[i] > 0.0)
            {
                chiGOS=calculate_chi_GOS_dens(q, w_global,i);
                if(AddELF)
                {
                    oneovereps=dcomp(1.0,0.0)/(chiGOS+dcomp(1.0,0.0));
                    loss-=oneovereps.imag();
                }
                else
                {
                    sum_chi+=chiGOS;
                }
                
            }
        }
        if(!AddELF)
        {
            oneovereps=dcomp(1.0,0.0)/(sum_chi+dcomp(1.0,0.0));
            loss=-oneovereps.imag();
        }
     return loss;  
}

//// MODIF 19/08 LK: Define the inner electron excitation as individual ionization (apply exchange and new integration limit)
double calculate_loss_GOS_dens_MELF(double q)
{ // calculate contribution GOS to loss function including density effect
    double loss=0.0;
        for (int i = 0; i < MAXGOS; i++)
        {
            if (AiGOS[i] > 0.0)
            {
                double EmaxThisGOS = (E_0 + Edgei_GOS[i]) / 2.0;	//// Ionization cutoff
                if (w_global > EmaxThisGOS) continue;
                dcomp chiGOS = calculate_chi_GOS_dens(q, w_global,i);
                dcomp oneovereps = dcomp(1.0,0.0) / (chiGOS+dcomp(1.0,0.0));
                double loss_i = -oneovereps.imag();
                double x = (q * q) / (2.0*E_0);
                loss_i *= (1.0 - x + x*x);
                loss += loss_i;
            }
        }
     return loss;  
}
//// END MODIF

dcomp calculate_chi_allGOS_dens( double q)
{ // calculate eps corresponding to the GOS, incl density effect
    dcomp chiGOS,chiGOS_i, oneoverchi, sum_chi, sum_oneoverchi;
    sum_chi=dcomp(0.0,0.0);
    sum_oneoverchi=dcomp(0.0,0.0);
        for (int i = 0; i < MAXGOS; i++)
        {
            if (AiGOS[i] > 0.0)
            {  if (errno !=0) my_perror("before chi_allGOS an error occured");
              chiGOS_i=calculate_chi_GOS_dens(q, w_global,i);
             
                if(AddELF)
                {
                    oneoverchi=dcomp(1.0,0.0)/(chiGOS_i+dcomp(1.0,0.0))-dcomp(1.0,0.0);
                    sum_oneoverchi+=oneoverchi;
                }
                else
                {
                    sum_chi+=chiGOS_i;
                }
            }
            
        }
        if(AddELF)
        {
            chiGOS=dcomp(1.0,0.0)/(sum_oneoverchi+dcomp(1.0,0.0)) -dcomp(1.0,0.0);
        }
        else
        {
            chiGOS=sum_chi;
        }
     return chiGOS;  
            
}



double  mylossfun(double q)
{
    dcomp eps;
    double result,GOS,Kaneko,Belkacem;
    eps = calculate_eps_osc( q);
    result = eps.imag() / (norm(eps));  //calculate Im [-1/eps],  norm returns sum of squares
    
    GOS=calculate_loss_GOS_dens(q);
    Kaneko= calculate_Loss_AA_LL(q);
    Belkacem=calculate_loss_Belkacem(q);
    result = result + GOS + Kaneko + Belkacem;
   
    return(result/q);  //return epsilon over q, as required for quanc8 integral
}
double  mylossfun_exchange(double q)
{
    dcomp eps;
    double GOS,Kaneko,Belkacem, result, theta;
    double loss_function_direct,loss_function_exchange;
    eps = calculate_eps_osc( q);
    result = eps.imag() / (norm(eps));  //calculate Im [-1/eps],  norm returns sum of squares
    
    GOS=calculate_loss_GOS_dens(q);
    Kaneko= calculate_Loss_AA_LL(q);
    Belkacem=calculate_loss_Belkacem(q);
    loss_function_direct = result + GOS + Kaneko + Belkacem;
    
    
               
    //The detected electron could be the struck electron calculate the ddcs for that process and take into account the interference between both
    // all high-energy stuff, i.e. assuming the struck electron was free and stationary before the collision
    double E_transfer=E_0-w_global;  // energy of the struck electron that is detected, the energy of the projectile (the slow electron)  is (omega - binding energy)  after the collision for these events
    double E_slow = w_global - BE_for_exchange;  // leftover energy of the projectile
    double gamma_slow=1.0+E_slow/rest_mass_energy;
    double v_slow = velocity_from_energy(E_slow,Projectile_Mass);
    double p_slow=gamma_slow*Projectile_Mass*v_slow;
    //p_0^2=p_detected^2+q^2-2  p_detected q cos(theta) cosine rule
    theta = acos((p_0*p_0-p_detected*p_detected-q*q)/(q*p_detected));
    double p_perp= p_detected*sin(theta); // perpendicular momentum of the struck electron that is detected
    if(p_perp < p_slow)
    {   double w_global_direct=w_global;
        w_global = E_transfer;
        double theta_slow=asin(-p_perp/p_slow);// should be compentated by the perpendicular comp. of momentum of the scattered( now slow) electron
        double q_exchange=sqrt(p_slow*p_slow+p_0*p_0- 2*p_0*p_slow*cos(theta_slow));
        eps= calculate_eps_osc( q_exchange);// at w_global
        result = eps.imag() / (norm(eps));  //calculate Im [-1/eps],  norm returns sum of squares
        GOS=calculate_loss_GOS_dens(q_exchange);
        Kaneko= calculate_Loss_AA_LL(q_exchange);
        Belkacem=calculate_loss_Belkacem(q_exchange);
        loss_function_exchange = result + GOS + Kaneko + Belkacem;
        w_global= w_global_direct;//put back original w_global
    }
    else
    {
        loss_function_exchange = 0.0;  
    }    
    result =loss_function_direct  + loss_function_exchange - sqrt(loss_function_direct *loss_function_exchange);
                  
    return(result/q);  //return epsilon over q, as required for quanc8 integral
}

//// MODIF 05/08 LK: Add the Born-Ochkur exchange factor to the loss function
double mylossfun_BornOchkur(double q)
{
    // Born-Ochkur exchange factor: fex = 1 - x + x^2, x = q^2/(2T)
    // T is the projectile kinetic energy, already held in the global E_0
    double x = (q * q) / (2.0 * E_0);
    double fex = 1.0 - x + x * x;
    return mylossfun(q) * fex;
}
//// END MODIF

//// MODIF 19/08 LK: New function to compute ELF as the sum of independant ELF_i with various use of exchange and integration limit
double mylossfun_MELF(double q)
{
    double result = 0.0;
    for (int i = 0; i < MAXOSC; i++)
    {
        if (fabs(Ai[i]) < 1e-90) continue;
        if (w_global > EmaxForOsc[i]) continue;
        dcomp eps_i = get_eps_component(i, q);
        double loss_i = eps_i.imag() / norm(eps_i);
        loss_i *= EdgeFactor(w_global, Ethi[i], deltai[i]);
    
        if (UseExchangeForOsc[i])
        {
            double x = (q * q) / (2.0 * E_0);
            loss_i *= (1.0 - x + x * x);
        }
        result += loss_i;
    }
    double GOS = calculate_loss_GOS_dens_MELF(q);	//// Use the GOS with ionization type excitation
    return (result + GOS) / q;
}    
//// END MODIF


double  fun_neutral(double q)  // have to add Belkacem and Kaneko?
{
//    Just change Z ^ 2 -> (Z - rho(q)) ^ 2
//        for H0 it reads
//            Z = 1
//            rho(q) = 1 / (1.0 + (2 * q) ^ 2) ^ 4
    dcomp eps;
    double result, GOS,rho, factor;
    eps = calculate_eps_osc(q);

    result = eps.imag() / (norm(eps));  //calculate Im [-1/eps],  norm returns sum of squares

    //GOS = calculate_GOS(q);
    GOS=calculate_loss_GOS_dens(q);
    result = result + GOS;
    rho = 1.0 / ((1.0 +  q*q /4.0)*(1.0 +  q*q/4.0));  //see email Pedro March 1 2019
    factor = (1.0 - rho)*(1.0 - rho);
    return(factor*result / q);  //return epsilon over q, as required for quanc8 integral
}



double DSEPfun(double q)
{  //for any incoming direction
    //only for oscillators no contributions from  GOS, Kaneko and Belkacem  for now
    dcomp eps, epsfraction;
    double qs1_plus,qs1_minus, result;
   
    double theta_E; // characteristic angle Egerton EELS (3rd ed. ) eq. 3.28 omega/ gamma E_0
    double theta_scat,tmp;   //calculated here in the small angle approximation
    
    theta_E= (w_global/E_0)*  Egerton_rel_cor_factor;
    tmp=q*q-p_0*p_0*theta_E*theta_E;
    if(tmp < 0.0) 
    {
    // |tmp| values are very small (1e-14) so contrinute nothing
        return 0.0;
    }
    theta_scat=sqrt(tmp)/p_0;
    
    eps = calculate_eps_osc(q);
   
    
    epsfraction = (eps - dcomp(1.0, 0.0))*(eps - dcomp(1.0, 0.0)) / (eps*(eps + dcomp(1.0, 0.0)));
 
    qs1_plus=abs(p_0*theta_scat*cos(theta_global)+p_0* theta_E * sin (theta_global));
    qs1_minus=abs(p_0*theta_scat*cos(theta_global)-p_0 * theta_E * sin (theta_global));
    result = epsfraction.imag() *(qs1_plus+ qs1_minus )/ (q*q*q);
    return (result);
 return 0.0;}

     dcomp calc_total_eps( double omega, double q)
// the addition of the different contributions is an approximation only and not valid at significant densities of more than one type
// problems can be expected if there are  contributions of different models at similar energies
    { 
        dcomp eps,chi_Belk, chi_GOS,chi_Kaneko, totaleps, oneovereps; 
        dcomp oneover_chi_GOS,oneover_chi_Belk,oneover_chi_Kaneko;
        w_global = omega;
        eps=calculate_eps_osc(q);
        chi_Belk=calculate_chi_Belkacem( q);
        chi_GOS=calculate_chi_allGOS_dens(q);
        chi_Kaneko=calculate_chi_AA_LL(q);
        if (AddELF)
        {  
            oneovereps= dcomp(1.0,0.0)/eps; //oscillator contribution
            oneover_chi_Belk   =  dcomp(1.0,0.0)/(chi_Belk   + dcomp(1.0,0.0))-dcomp(1.0,0.0); //add  chi contribution
            oneover_chi_GOS    =  dcomp(1.0,0.0)/(chi_GOS    + dcomp(1.0,0.0))-dcomp(1.0,0.0); //add  chi contribution
            oneover_chi_Kaneko =  dcomp(1.0,0.0)/(chi_Kaneko + dcomp(1.0,0.0))-dcomp(1.0,0.0); //add  chi contribution
            totaleps= dcomp(1.0,0.0)/(oneovereps+oneover_chi_Belk+oneover_chi_GOS+oneover_chi_Kaneko);
        }
        else
        {
            totaleps=eps;//+chi_GOS+chi_Belk+chi_Kaneko;
        }
        if (errno !=0) my_perror("calc_total_eps an error occured");
        return totaleps;
    }
    
    
    double Surface_loss_only(double omega, double q)
    {  
    dcomp eps, epsfraction;
    eps = calc_total_eps(omega, q);
    epsfraction = (eps - dcomp(1.0, 0.0))*(eps - dcomp(1.0, 0.0)) / (eps*(eps + dcomp(1.0, 0.0)));
    return (epsfraction.imag());  
}

int  copyP_to_Vars(double *p,  int modelchoice)
{
    double precision;
    int Param_Offset;
    modelDrude = false; modelDL = false; modelMerminLL = false;  modelVlasov = false;
    modelForouhiBloomer=false;modelTL_an=false; modelTaucLorentz = false; modelTauc_Mermin = false;
    epsbkg=1.0;
    if (modelchoice == 1)
    {
         modelDrude = true;
         epsbkg = p[0]; 
    }
    else if (modelchoice == 2)
    {
        modelDL = true;
    }
    else if (modelchoice == 3)
    {
        modelMerminLL = true; 
    }
     else if (modelchoice == 4)
    {
        modelVlasov = true; 
    }
    else if (modelchoice == 5)
    {  
        modelTaucLorentz = true;
        epsbkg = p[0];  
    }
    else if (modelchoice == 6)
    {   
         modelTL_an = true; 
         aTLan = p[0] / Hartree; 
    }
    else if (modelchoice == 7)
    {  
        modelTauc_Mermin = true; 
        epsbkg = p[0];  
    }
    else if (modelchoice == 8)
    {   
        modelForouhiBloomer = true; 
        epsbkg = p[0];  
    }
     else if (modelchoice == 9)
    {   
        modelBrendelBormann = true; 
        Orosco_Coimbra_way= false;
        epsbkg = p[0];  
    }
    else if (modelchoice == 10)
    {   
        modelBrendelBormann = true; 
        Orosco_Coimbra_way = true;
        epsbkg = p[0];  
    }
    
   

    // now read the components
    for (int i = 0; i < MAXOSC; i++)
    {  
        Ai[i] = p[5* i + 1];
        if(fabs(Ai[i])> 1e-90 )
        {   
            if (modelVlasov)
            {
                Q_Vlasov[i] = p[5 * i + 2];
            }
            else
            {
                wi[i] = p[5 * i + 2] / Hartree;
            }
            if (!modelForouhiBloomer)  gammai[i] = p[5 * i + 3] / Hartree;
            if (modelForouhiBloomer) gammai[i] = p[5 * i + 3] / (Hartree*Hartree);
        
            if (modelDrude or modelBrendelBormann) Ai[i] = Ai[i] / (Hartree*Hartree);// Ai is in eV^2 in DL and BB
            else if  (modelTaucLorentz or modelTL_an or modelTauc_Mermin )  Ai[i] = Ai[i] / (Hartree);// Ai is in eV in TL and TLan
            if(!modelBrendelBormann)
            {
                alphai[i] = p[5 * i + 4];
            }
            else// using the alpha slot for the Gaussian sigma in that case
            { 
                alphai[i]=1.0;
                sigmai[i]=p[5 * i + 4]/Hartree;
            }
            if  (modelTaucLorentz or modelTL_an or modelTauc_Mermin )
            {
                TaucGap[i] = p[5* i + 5] / Hartree;  
                gapi[i]=0.0;  //  gapi in LL model,gapi[]=0 gives Mermin
            }
            else
            {
                 gapi[i] = p[5* i + 5] / Hartree;
            }

            if (modelTauc_Mermin)
            {   

               double A_Drude= wi[i]*Ai[i];
               wi[i]= sqrt(A_Drude+wi[i]*wi[i]);
               Ai[i]=A_Drude/(wi[i]*wi[i]);
               printf(" copyP Mermin param wi  %6.4f  Ai  %6.4f\n",wi[i]*Hartree,Ai[i]);
            }
        }
        
    }
   
    Param_Offset = 5* MAXOSC;
    for (int i = 0; i < MAXGOS; i++)
    {   
        AiGOS[i] = p[4 *  i+ Param_Offset + 1];
        Zi_GOS[i] = p[4 * i + Param_Offset + 2];
        int nl = (int)round(p[4 * i + Param_Offset + 3]);
        if (nl == 10)
        {
            n_i_GOS[i ] = 1;
            l_i_GOS[i ] = 0;
        }
        else if (nl == 20)
        {
            n_i_GOS[i] = 2;
            l_i_GOS[i] = 0;
        }
        else if (nl == 21)
        {
            n_i_GOS[i ] = 2;
            l_i_GOS[i ] = 1;
        }
        else if (nl == 30)
        {
            n_i_GOS[i] = 3;
            l_i_GOS[i ] = 0;
        }
        else if (nl == 31)
        {
            n_i_GOS[i] = 3;
            l_i_GOS[i] = 1;
        }
        else if (nl == 32)
        {
            n_i_GOS[i] = 3;
            l_i_GOS[i] = 2;
        }
        else
        {
            n_i_GOS[i] = 0;  //this should never occur
            l_i_GOS[i] = 0;
        }
        Edgei_GOS[i] = p[4 * i + Param_Offset + 4] / Hartree;

    }
    Param_Offset += 4 * MAXGOS;
     for (int i = 0; i < MAXBELKACEM; i++)
    {
        Ai_Belkacem[i]      = p[4 * i + Param_Offset + 1];
        wi_Belkacem[i]      = p[4 * i + Param_Offset + 2]/Hartree;
        gammai_Belkacem[i]  = p[4 * i + Param_Offset + 3]/Hartree;
    }
    Param_Offset += 4 *  MAXBELKACEM;
    

    for (int i = 0; i < MAXKANEKO; i++)
    {
        N_Kaneko[i]        = p[6 * i + Param_Offset + 1];
        Q_Kaneko[i]           = p[6 * i + Param_Offset + 2];
        Width_Kaneko[i]       = p[6 * i + Param_Offset + 3]/Hartree;
        Edge_ArchubiKaneko[i] = p[6 * i + Param_Offset + 4]/Hartree;
        l_Kaneko[i]= (int)round(p[6 * i + Param_Offset + 5]);
        gamma_Kaneko[i]       = p[6 * i + Param_Offset + 6];
        int l=l_Kaneko[i];
        //double fraction_of_shell_filled=N_Kaneko[i]/(2.0*(2*l+1));
        
        double Q=Q_Kaneko[i];
        double alpha = 1/(Q*Q);
        
        if((fabs(N_Kaneko[i]) > 1E-10)&&(!OriginalKaneko) )
            {  
            int  doublefact;
            if(l==0) doublefact=1;
            else if (l==1) doublefact=3;
            else if (l==2) doublefact=5*3;
            else  if (l==3) doublefact=7*5*3;
            else
            {
                printf(" lvalue %i not implemented for Kaneko \n",l);
                doublefact=0;  // to avoid warning used uninitialized
            }
           
            A_l_Kaneko[i]=pow(2.0,l)* pow(alpha,(2.0*l+3.0)/2.0)/(doublefact* pow(pi,1.5));
            A_0_Kaneko[i]=pow(alpha/pi,1.5);

           w_pl_l[i]=sqrt(gamma_Kaneko[i]*2*(2*l+1)*doublefact*pow(Q,3)*exp(l)/(pow(2.0*l,l)*sqrt(pi)) );// does not depend on occupation level, assumes all levels full
           w_pl_0[i]=sqrt(gamma_Kaneko[i]*2*pow(Q,3)/sqrt(pi));
        }
        else if((fabs(N_Kaneko[i]) > 1E-10)&&(OriginalKaneko) )
        {   Q_Kaneko[i]= Q_Kaneko[i]*pow(N_Kaneko[i],1.0/3.0);// variable Q_kaneko now contains Kaneko's  q_mean
            w_pl_0[i]=sqrt(pow(Q_Kaneko[i],3)/(sqrt(pi)) );  // no spin maximum occupation = 1
            w_pl_l[i]=w_pl_0[i]; // implicitly assumes l=0
            l_Kaneko[i]=0;   // implicitly assumes l=0
            gamma_Kaneko[i]=0.5* gamma_Kaneko[i];
         
           // Q_Kaneko[i]= Q_Kaneko[i]*pow(N_Kaneko[i],1.0/3.0);// variable Q_kaneko now contains Kaneko's  q_mean
           // printf("original Kaneko, plasmon energy %6.3f q_mean %6.2f\n",w_pl_l[i]*Hartree, Q_Kaneko[i] );
        }
    }

        

    Param_Offset += 6*MAXKANEKO;
    
    //// MODIF 04/08 LK: add and shift the reading of the parameters
    for (int i = 0; i < MAXOSC; i++)
    {
        Ethi[i] = p[2*i + Param_Offset + 1] / Hartree;
        deltai[i] = p[2*i + Param_Offset + 2] * Hartree;	// delta has units 1/energy
    }
    Param_Offset += 2*MAXOSC;
    //// END MODIF
    
    E_0= p[Param_Offset + 1]/Hartree;

    
    precision = p[Param_Offset + 2];
    abserr = 1e-55;
    //relerr = 1.0e-5;
    relerr = 2.0e-4 / precision;

    q2maxFactor = p[Param_Offset + 3];
    q2surfmaxFactor = p[Param_Offset + 4];
    NStep = (int)p[Param_Offset + 5];
    FirstEnergy = p[Param_Offset + 6] / Hartree;
    StepSize_eV = p[Param_Offset + 7];
    StepSize = StepSize_eV / Hartree;
    
       
    
    UnitCellDensity = p[Param_Offset + 8];// note this index is out of order, late addition  per ^3
    UnitCellDensity = UnitCellDensity*BohrRadius*BohrRadius*BohrRadius;// now per a.u.^3
    PL_GOS = sqrt(4.0*pi*UnitCellDensity);
    lin_cont_deltaE=p[Param_Offset + 9];
   
    q_lower=p[Param_Offset + 10];
    q_upper=p[Param_Offset + 11];
    LastMomentum=p[Param_Offset + 12];
    Stepsize_qplot=p[Param_Offset + 13];
    c_transition=p[Param_Offset + 14];  //where small q limit is reached (Mermin, Belkacem rpa-like)
    ExchangeCorrection=p[Param_Offset + 15];
    if(p[Param_Offset + 16]==0.0)
    {
        OriginalKaneko=false;
    }
    else
    {
        OriginalKaneko=true;
    }
    if(p[Param_Offset + 17]==0.0)
    {
        AddELF=false;
    }
    else
    {
        AddELF=true;
    }
    if(p[Param_Offset + 18]==0.0)
    {
        ApplySumRuleToGOS=false;
    }
    else
    {
        ApplySumRuleToGOS=true;
    }
  
    if(p[Param_Offset + 19]==0.0)
    {
        Apply_Mermin_Correction=false;
        DirectMethod=false;
    }
    else if (p[Param_Offset + 19]==1.0)
    {
        Apply_Mermin_Correction=true;
        DirectMethod=false;
    }
    else
    {   Apply_Mermin_Correction=false;
        DirectMethod=true;
    }
    if(p[Param_Offset + 20]==0.0)
    {
        Full_dispersion=false;
    }
    else
    {
        Full_dispersion=true;
    }
    maxEnergyDensityEffect = p[Param_Offset + 21]/Hartree;
   
    // for REELS spectrum 
    sigma=p[Param_Offset + 22]/2.355;
    PIcoef1= p[Param_Offset + 23];
    PIcoef2= p[Param_Offset + 24];
    PIcoef3= p[Param_Offset + 25];
    theta0= p[Param_Offset + 26]*pi/180.0;  //now in rad
    theta1= p[Param_Offset + 27]*pi/180.0;
    surf_ex_factor= p[Param_Offset + 28];
    //p[Param_Offset + 29] not used at the moment

    fraction_DIIMFP= p[Param_Offset + 30];
    if(p[Param_Offset + 31]==1.0)
    { 
        proton=true;
        Projectile_Mass=1836.15;
        ExchangeCorrection=0;  // make sure  we are not bothered by exchange correction for protons
    }
    else
    {   
        proton=false;  //electron
        Projectile_Mass=1.0;
    }
    rest_mass_energy=Projectile_Mass*C*C;

    if(p[Param_Offset + 32]==1.0)
    {
        Dispersion_relativistic= true;
    }
    else
    {
         Dispersion_relativistic= false;
    }    
    theta_max= p[Param_Offset + 33]/1000.0;  // now in rad
    NThetaStep= p[Param_Offset + 34];
    MottCorrection = false;
    if(p[Param_Offset + 35]==1.0 && proton) MottCorrection = true;
    
    Egerton_rel_cor_factor = (E_0+Projectile_Mass*C*C) / (E_0 +2.0*Projectile_Mass*C*C);// Egerton 3rd edition page 427 for theta_E
 
    v_0 = velocity_from_energy(E_0,Projectile_Mass);
    gamma_rel=1.0+E_0/rest_mass_energy;
    beta_r= v_0/C;
    b_zero=pow(1.0-sqrt(1.0-beta_r*beta_r),2);  // for the Moller cross section
   
    p_0=gamma_rel*Projectile_Mass*v_0;
 
    BE_for_exchange=p[Param_Offset + 36];
   // printf("BE_for_exchange%6.4f\n",BE_for_exchange);
    BE_for_exchange=BE_for_exchange/Hartree;
    int ExchangeCorrectionMode= (int)p[Param_Offset + 37];
    if(ExchangeCorrection >0)  ExchangeCorrection+=ExchangeCorrectionMode;//ExchangeCorrection=0 no exchange cor.
                                                                          //ExchangeCorrection=1 Ashley method
                                                                          //ExchangeCorrection=2 SBethe method
                                                                      
    delayed_dispersion= (bool) p[Param_Offset + 38]; //if 1 then dispersion calculated from energy of the level is  assumed shifted as in the LL model    
    add_Dopplerwidth_to_classical_DF= (bool) p[Param_Offset + 39];
    
    UseBornOchkurExchange = (bool) p[Param_Offset + 41];	//// MODIF 05/08 LK: Read the new parameter
    //// MODIF 05/08 LK: Born-Ochkur exchange factor is derived from a non-relativistic formalism. If beta > 0.3 (arbitrary) then it prints a warning for each energy where this condition is raised.
    if (UseBornOchkurExchange && beta_r > 0.3)
    {
    fprintf(stderr,
        "WARNING: Born-Ochkur exchange factor used at beta=%.3f (E_0=%.2f keV) -- "
        "Born-Ochkur is a non-relativistic approximation, check its reliability.\n",
        beta_r, E_0 * Hartree / 1000.0);
    }
    //// END MODIF
    
    //// MODIF 19/08 LK: Read the bandgap and the valence binding energy and determine the max energy for each oscillator
    Egap = p[Param_Offset + 50 + 2*MAXOSC] / Hartree;
    for (int i = 0; i < MAXOSC; i++)
    {
        UseExchangeForOsc[i] = (bool) p[Param_Offset + 50 + i];
        int excitation_type = (int) p[Param_Offset + 50 + MAXOSC + i];
        BindingEnergyForOsc[i] = p[Param_Offset + 50 + 2*MAXOSC + 1 + i] / Hartree;   //// MODIF LK 01/09: new per-osc slot
        if (excitation_type == 1) EmaxForOsc[i] = E_0;
        else if (excitation_type == 2) EmaxForOsc[i] = (E_0 + BindingEnergyForOsc[i]) / 2.0;
        else EmaxForOsc[i] = E_0;
    }
    ////END MODIF

    DebugMode= p[Param_Offset + 40];
    if(DebugMode )printf("debug mode is on\n");
    return 0;//end copyP_to_Vars
}

    

double DIIMFP_at_omega(double omega)
{
    w_global = omega;
    double q2used, Qrecoil_max;
    
    bool neutral=false;
    double E_1=E_0-omega;
    if(E_1 <=0.0) return 0.0;
 
    double v_1 = velocity_from_energy(E_1,Projectile_Mass); // projectile velocity, used in DIIMFP and DSEP calculation
    double p_1=(1.0+E_1/rest_mass_energy)*Projectile_Mass*v_1;
   
    double q1 = p_0 - p_1; //integration boundaries
    double q2 = p_0 + p_1;
    
    Qrecoil_max= recoil_energy(q2);
   

    if ((Qrecoil_max > 2*omega)and !modelVlasov) //calculate the momentum of an electron with energy omega, i.e. the maximum transferred momentum, factor 2 beacuse struck electron not stationary
    {  
         double totalE=omega+C*C;
         q2 =sqrt(pow(totalE,2)-pow(C,4))/C+5.0;// the additional amount 5 is because electrons are not stationary, so better go out a bit further
    }
    q2used=q2;

    if(q1 <  q_lower) q1=q_lower;  //this is for when we want to calculate partial diimfp's
    if(q2used <   q_lower) q2used=q_lower;
    if(q1 >  q_upper) q1=q_upper;
    if(q2used > q_upper)  q2used=q_upper;
    if (errno !=0) my_perror("DIIMFP_at_omega: before quanc8 an error occured");
    if (neutral && (Projectile_Mass > 1.0) )
        quanc8(fun_neutral, q1, q2used, abserr, relerr);
    else if (UseBornOchkurExchange)					//// MODIF 05/08 LK: If use of the Born-Ochkur exchange factor
        quanc8(mylossfun_BornOchkur, q1, q2used, abserr, relerr);	//// MODIF 05/08 LK: Then integrate over the loss function WITH this factor
    else if (ExchangeCorrection == 1)
        quanc8(mylossfun_exchange, q1,  q2used, abserr, relerr); 
    else
        quanc8(mylossfun, q1,  q2used, abserr, relerr); 

    double DIIMFP_au=2.0* quanc8result /(pi*v_0*v_0);
    if  (ExchangeCorrection == 2)  
    { 
        double MoellerFactor= 1.0+ pow((omega/(E_0-omega+BE_for_exchange)),2);
        MoellerFactor += -(1.0-b_zero)*omega/(E_0-omega+BE_for_exchange)+b_zero*omega*omega/(E_0*E_0);
        DIIMFP_au= DIIMFP_au*MoellerFactor;
    } 
    if(DebugMode)
    {
        if(quanc8flag > 0.0)
        {
            printf("omega %6.3f, quanc8flag %9.6f q1 %6.4g q2 %6.4g\n",omega*Hartree,quanc8flag,q1,q2used);
        } 
        if (errno !=0) printf("omega %6.4f\n",omega);
        if (errno !=0) my_perror("DIIMFP_at_omega: an error occured");
    }   
   // if(DIIMFP_au < 1e-200) DIIMFP_au=1e-200;
    return DIIMFP_au;
}

//// MODIF 19/08 LK:
double DIIMFP_MELF_at_omega(double omega)
{
    w_global = omega;
    double q2used, Qrecoil_max;
    
    double E_1=E_0-omega;
    if(E_1 <=0.0) return 0.0;
 
    double v_1 = velocity_from_energy(E_1,Projectile_Mass); // projectile velocity, used in DIIMFP and DSEP calculation
    double p_1=(1.0+E_1/rest_mass_energy)*Projectile_Mass*v_1;
   
    double q1 = p_0 - p_1; //integration boundaries
    double q2 = p_0 + p_1;
    
    Qrecoil_max= recoil_energy(q2);
   

    if ((Qrecoil_max > 2*omega)and !modelVlasov) //calculate the momentum of an electron with energy omega, i.e. the maximum transferred momentum, factor 2 beacuse struck electron not stationary
    {  
         double totalE=omega+C*C;
         q2 =sqrt(pow(totalE,2)-pow(C,4))/C+5.0;// the additional amount 5 is because electrons are not stationary, so better go out a bit further
    }
    q2used=q2;

    if(q1 <  q_lower) q1=q_lower;  //this is for when we want to calculate partial diimfp's
    if(q2used <   q_lower) q2used=q_lower;
    if(q1 >  q_upper) q1=q_upper;
    if(q2used > q_upper)  q2used=q_upper;
    if (errno !=0) my_perror("DIIMFP_at_omega: before quanc8 an error occured");

    quanc8(mylossfun_MELF, q1,  q2used, abserr, relerr); 

    double DIIMFP_au=2.0* quanc8result /(pi*v_0*v_0);

    if(DebugMode)
    {
        if(quanc8flag > 0.0)
        {
            printf("omega %6.3f, quanc8flag %9.6f q1 %6.4g q2 %6.4g\n",omega*Hartree,quanc8flag,q1,q2used);
        } 
        if (errno !=0) printf("omega %6.4f\n",omega);
        if (errno !=0) my_perror("DIIMFP_at_omega: an error occured");
    }   
   // if(DIIMFP_au < 1e-200) DIIMFP_au=1e-200;
    return DIIMFP_au;
}
//// END MODIF
 
double  calc_DSEP(double *DSEPresult, double BeamE,  double theta )
{  
    double  omega, q2used,sep,current_stepsize; 
    double prefactor;
    theta_global=theta;

    omega=FirstEnergy;
    double LastEnergy=FirstEnergy+NStep*StepSize;
    prefactor=2/(pi*v_0*v_0*cos(theta));
    current_stepsize=StepSize;
    if(LastEnergy > 0.9*BeamE)LastEnergy=0.9*BeamE;
    sep=0.0;
    for(int i = 0; omega<LastEnergy; i++)
    {   
        double E_1=E_0-omega;
        w_global = omega;
        //double rel_cor_factor_1 =  (1 + E_1 / (2.0*Projectile_Mass*C*C)) / pow((1 + E_1 / (Projectile_Mass*C*C)),2);    

        //double v_1 = sqrt(2.0*  rel_cor_factor_1*E_1/Projectile_Mass); // projectile velocity, used in DIIMFP and DSEP calculation
        //double gamma_relativistic_1=1.0+E_1/rest_mass_energy;
        //double p_1=gamma_relativistic_1*Projectile_Mass*v_1;
        double v_1 = velocity_from_energy(E_1,Projectile_Mass); // projectile velocity, used in DIIMFP and DSEP calculation
        double p_1=(1.0+E_1/rest_mass_energy)*Projectile_Mass*v_1;
        double q1 = p_0 -p_1; 
        double q2 = p_0 + p_1;
      
        if (q2 > q1* q2maxFactor)
            q2used = q1* q2maxFactor;
        else
            q2used = q2;
        quanc8(DSEPfun, q1, q2used, abserr, relerr);
        DSEPresult[i] = prefactor* quanc8result/ Hartree;  // now in eV-1
        
        sep=sep+DSEPresult[i]*current_stepsize*Hartree;
        omega=omega+current_stepsize;
                
    }//end loop over energy losses
    return sep;
}
double DDCS_excl_retardation(double omega, double theta)
{   
 
    dcomp eps, oneovereps;  
    double ddcs_exchange;
    double  E_detected= E_0-omega; // energy of the fast electron
    double v_detected = velocity_from_energy(E_detected,Projectile_Mass);
    double gamma_relativistic_detected=1.0+E_detected/rest_mass_energy;
    p_detected=gamma_relativistic_detected*Projectile_Mass*v_detected;
    

    double q=sqrt(p_detected*p_detected+p_0*p_0- 2*p_0*p_detected*cos(theta));
    if(q/p_0 <  1e-6)// small angle limit, maybe previous q not accurate due to small difference large numbers (seems not to be an issue)
    {
        double q_along = (p_0-p_detected*cos(theta));
        double q_perp= p_detected*sin(theta);
        q=sqrt(q_perp*q_perp+q_along*q_along);
    }
    eps= calc_total_eps( omega, q);
    oneovereps=dcomp(1.0,0.0)/ eps;
    double ddcs=-oneovereps.imag()*p_0*p_0/(q*q);
    if (ExchangeCorrection == 1)
    // we may have detected a struck electron after events with energy transfer E_0-omega
    // these collisions are indistinguishable from detection of fast electrons after energy transfer omega,
    //  then the scattered electron (not detected)  is the slow one, and travels at 90-theta_detector   Ashley method
            {   
                if(omega > 0.5*(E_0-BE_for_exchange))
                {
                    ddcs=0.0;
                }
                else
                {   //The detected electron could be the struck electron calculate the ddcs for that process and take into account the interference between both
                    // all high-energy stuff, i.e. assuming the struck electron was free and stationary before the collision
                    double E_transfer=E_0-omega;  // energy of the struck electron that is detected, the energy of the projectile (the slow electron)  is (omega - binding energy)  after the collision for these events
                    double E_slow = omega- BE_for_exchange;  // leftover energy of the projectile
                    double gamma_slow=1.0+E_slow/rest_mass_energy;
                    double v_slow = velocity_from_energy(E_slow,Projectile_Mass);
                    double p_slow=gamma_slow*Projectile_Mass*v_slow;
                    
                    double p_perp= p_detected*sin(theta); // perpendicular momentum of the struck electron that is detected
                    if(p_perp < p_slow)
                    {
                        double theta_slow=asin(-p_perp/p_slow);// should be compentated by the perpendicular comp. of momentum of the scattered( now slow) electron
                        double q_exchange=sqrt(p_slow*p_slow+p_0*p_0- 2*p_0*p_slow*cos(theta_slow));
                        dcomp eps_exchange = calc_total_eps( E_transfer, q_exchange);
                        dcomp oneovereps_exchange= dcomp(1.0,0.0)/ eps_exchange;
                        ddcs_exchange = -oneovereps_exchange.imag()*p_0*p_0/(q_exchange*q_exchange);
                    }
                    else
                    {
                        ddcs_exchange=0.0;  
                    }    
                    ddcs =ddcs + ddcs_exchange - sqrt(ddcs*ddcs_exchange);
                }    
            }
            else if (ExchangeCorrection == 2)
            {
                if(omega > 0.5*(E_0-BE_for_exchange))
                 {
                    ddcs=0.0;
                }
                else
                {
                    double MoellerFactor= 1.0+ pow((omega/(E_0-omega+BE_for_exchange)),2);
                    MoellerFactor += -(1.0-b_zero)*omega/(E_0-omega+BE_for_exchange)+b_zero*omega*omega/(E_0*E_0);
                    ddcs= ddcs*MoellerFactor;
                }
                
                
            }
           
    return ddcs;  // note  prefactor = 1.0/(pi * pi * v_0 * v_0 * UnitCellDensity) is added at function call;
}

double DDCS_incl_retardation(double omega, double theta)
{   
    double  v_over_c_square= v_0*v_0/(C*C);
    dcomp eps, oneovereps;  
    double  E_1= E_0-omega;
   
    double v_1 = velocity_from_energy(E_1,Projectile_Mass); // projectile velocity, used in DIIMFP and DSEP calculation
    double p_1=(1.0+E_1/rest_mass_energy)*Projectile_Mass*v_1;
    double theta_e= (p_0-p_1)/p_0; 

    double q=sqrt(p_1*p_1+p_0*p_0- 2*p_0*p_1*cos(theta));
    if(q/p_0 < 1e-6)  // small angle limit, maybe previous q not accurate due to small difference large numbers (seems not to be an issue)
    {
        
        double q_along = (p_0-p_1*cos(theta));
        double q_perp= p_1*sin(theta);
        double q_small_angle=sqrt(q_perp*q_perp+q_along*q_along);
        q=q_small_angle;
    }
    eps= calc_total_eps( omega, q);
    oneovereps=dcomp(1.0,0.0)/ eps;
    double a = v_over_c_square*eps.real()-1.0;
    double b=eps.imag()*v_over_c_square;
    double d=theta*theta-theta_e*theta_e * a ;
    double e= theta_e*theta_e*eps.imag()*v_over_c_square;
    double top= theta*theta+theta_e*theta_e*(a*a+b*b);
    double bottom=d*d+e*e;
    double first= -oneovereps.imag();  //Egerton 3.70
    double ddcs=first*top/bottom;
    return ddcs;  // note  prefactor = 1.0/(pi * pi * v_0 * v_0 * UnitCellDensity) is added at function call;
}


  int Eps1Eps2(npArray ParameterArray, npArray omegaArray, npArrayComplex EpsArray, double q, int modelchoice)
    {   
        size_t NPoints = omegaArray.shape(0);
        double* omega= omegaArray.data();
        dcomp* Eps = EpsArray.data();
        copyP_to_Vars(ParameterArray.data(), modelchoice);

        for (int i = 0; i < NPoints; i++)
        {  
            Eps[i]=calc_total_eps( omega[i]/Hartree, q);
        }
        my_perror("eps1eps2 an error occured");
        return 0;
      
    }
    
    dcomp single_Eps1Eps2(npArray ParameterArray, double q, double my_omega, int modelchoice)
    {
        dcomp eps;
        copyP_to_Vars(ParameterArray.data(), modelchoice);
        eps=calc_total_eps( my_omega/Hartree, q);
        my_perror("eps1eps2 an error occured");
    
        return eps;
    }
    //epslib.Eps1Eps2_q(self.ParArray,self.x_axis, eps, self.Energy_qplot, self.DFChoice)
    int  Eps1Eps2_q(npArray ParameterArray, npArray q_Array, npArrayComplex EpsArray, double omega, int modelchoice)
    {
        size_t NPoints = q_Array.shape(0);
        double* q= q_Array.data();
        dcomp* Eps = EpsArray.data();
       
        copyP_to_Vars(ParameterArray.data(), modelchoice);  
        NStep=LastMomentum/Stepsize_qplot;
        for (int i = 0; i < NPoints; i++)
        {
            Eps[i]=calc_total_eps( omega/Hartree, q[i]);
        }
        if (errno !=0) my_perror(" Eps1Eps2_q an error occured");
        return 0;
    }
    
    int Loss_wide(npArray ParameterArray,   npArray np_w_array,  npArray np_ELF,  double q, int modelchoice,  npArray np_ResultArray)
    {
        //also calculates the constants for the simple DL formulae
        dcomp eps;
        double sumr1 = 0.0;
        double sumr2 = 0.0;
        double sumr3 = 0.0;
        double sumr4 = 0.0;
        double sumr5 = 0.0;
        double sumr6 = 0.0;
        double sumBethe=0.0;
        double delta= 0.005;
        double *w_array = np_w_array.data();
        double *ResultArray= np_ResultArray.data();
        double *ELF=np_ELF.data();
        copyP_to_Vars(ParameterArray.data(), modelchoice);  
        double omega = -0.125/Hartree;  // working in atomic units
        double step = 0.25/Hartree;
        double average_step =step;
       
        for (int i = 0; i < 1300; i++)
        {   
            omega += step;
            eps=calc_total_eps( omega, q);  
            ELF[i]= eps.imag() / (norm(eps));  //calculate Im [-1/eps],  norm returns sum of squares
            w_array[i] = omega*Hartree;  // output in eV
            sumr1 += 2.0 / pi * ELF[i] *  log(omega) * average_step; 
            sumr2 += 2.0 / pi * ELF[i] * average_step;
            sumr3 += 2.0 / pi * ELF[i] * omega * log(omega) * average_step; 
            sumr4 += 2.0 / pi * ELF[i] * omega *average_step;
            sumr5 += 2.0 / pi * ELF[i] * omega * omega * log(omega) * average_step; 
            sumr6 += 2.0 / pi * ELF[i] * omega * omega * average_step;
            sumBethe += 1.0/ (2.0*pi*pi)*ELF[i] * omega * average_step;
            average_step=step;
            step = step * (1+delta);
            average_step= (average_step+step)/2.0;  // average of step before and after
        }
        ResultArray[0] = sumr2 / exp(sumr1 / sumr2);  // self.C0
        ResultArray[1] = sumr4 / pow(exp(sumr3 / sumr4),2);//self.C1
        ResultArray[2] = sumr6 / pow(exp(sumr5 / sumr6),3); // self.C2
        ResultArray[3] = Hartree * exp(sumr1 / sumr2);  //self.I0 in eV
        ResultArray[4] = Hartree * exp(sumr3 / sumr4);  //self.MIE in eV
        ResultArray[5] = Hartree * exp(sumr4 / sumr5);  //self.Istraggling in eV
        ResultArray[6] =  sumBethe;
        if (errno !=0) my_perror("Loss_wide an error occured");
        return 0;
    }
    
    
    
    int  SurfLossFunc(npArray ParameterArray,  npArray np_SurfLoss, double q, int modelchoice)
//gos contributions are not included
    {   double *SurfLoss=np_SurfLoss.data();
        copyP_to_Vars(ParameterArray.data(), modelchoice);  ;  
     
        // we evaluate dielectric function, there decided to add either Chi or elf.
        for (int i = 0; i < NStep; i++)
        {
            double omega = FirstEnergy + StepSize*i;
            SurfLoss[i]=Surface_loss_only(omega, q);
        }

        if (errno !=0) my_perror("SurfLossFunc: an error occured");
        return 0;
    }
    

    double Lossfunction_inclGOS_atE(npArray ParameterArray,  double omega, int modelchoice)
    {
        double  GOS,KanekoLoss, Loss_at_E;
        dcomp eps,oneovereps;
       
        copyP_to_Vars(ParameterArray.data(), modelchoice);  
 
        w_global = omega/Hartree;
        eps=calculate_eps_osc(0.01);
        oneovereps=dcomp(1.0,0)/eps;
        Loss_at_E = -oneovereps.imag();
    
        GOS=calculate_loss_GOS_dens(0.01);
   
        Loss_at_E +=GOS;
        KanekoLoss=calculate_Loss_AA_LL(0.01);
        Loss_at_E +=KanekoLoss;
    //    if (errno !=0) my_perror(" lossfunction_inclGOS_atE an error occured");
        return(Loss_at_E);
    }
    
    int  lossfunction_inclGOS(double *p, double *lossfunction, double q, int modelchoice)
    {
        
        dcomp eps,oneovereps, chi_Belk;
        copyP_to_Vars(p, modelchoice);  
 
        for (int i = 0; i < NStep; i++)
        {   
            double omega = FirstEnergy + StepSize*i;
            w_global = omega;
            eps=calculate_eps_osc(q);
            oneovereps=dcomp(1.0,0.0)/eps;
            lossfunction[i] =-oneovereps.imag();

            lossfunction[i] +=calculate_loss_GOS_dens(q);

            chi_Belk=calculate_chi_Belkacem( q);
            oneovereps=dcomp(1.0,0.0)/(chi_Belk+dcomp(1.0,0.0));
            lossfunction[i] -=oneovereps.imag();
            lossfunction[i] +=calculate_Loss_AA_LL(q);
            
        }
       // if (errno !=0) my_perror(" lossfunction_inclGOS an error occured");
        return 0;
    }


    int  DIIMFP(npArray ParameterArray, npArray npDIIMFPArray, int modelchoice, npArray npResultArray)
    {
       //diimfp calculated for all energy losses with first fixed (for DIIMFP distribution) and then  variable stepsize, (for stopping straggling and IMFP calculation).
       double* DIIMFPresult= npDIIMFPArray.data();
       double* ResultArray=  npResultArray.data();
       double inv_lambda=1E-99;  // so we never get divide by 0
       double stopping=0.0;
       double straggling=0.0;
       copyP_to_Vars(ParameterArray.data(), modelchoice);
       double omega=0.5*StepSize;
       double Current_StepSize=StepSize;
       double Effective_Stepsize=Current_StepSize;
       double DIIMFP_au;
       double omega_max= 2*beta_r*beta_r*gamma_rel*gamma_rel*C*C;

       double MottFactor= (1- beta_r*beta_r*omega/omega_max); // correction factor SalvatPRA eq.16, not sure about this one, esp. for electrons
       
       bool KeepGoing=true;
       int i=0;
       while (KeepGoing)
       {    
            DIIMFP_au= DIIMFP_at_omega(omega);
           
 
            if(MottCorrection)DIIMFP_au *= MottFactor; // correction factor Salvat PRA  106 032809 eq.16
            inv_lambda     += DIIMFP_au * Effective_Stepsize;
            stopping       += DIIMFP_au * Effective_Stepsize * omega;
            straggling     += DIIMFP_au * Effective_Stepsize * omega * omega;
            if(i < NStep)
            {
                DIIMFPresult[i] = DIIMFP_au/ (Hartree*BohrRadius);
                i++;
            }
            else
            {   Effective_Stepsize =Current_StepSize;
                Current_StepSize=StepSize+lin_cont_deltaE*omega;
                Effective_Stepsize =(Effective_Stepsize+Current_StepSize)/2.0;  //average of previous and this stepsize
            }
            omega=omega+Current_StepSize;
            
            if (omega < omega_max) MottFactor= (1- beta_r*beta_r*omega/omega_max);// keep Mott factor the same for omega > omega_max
         

            if((ExchangeCorrection  != 0) && (omega > 0.5*(E_0+BE_for_exchange))) KeepGoing=false; 
            if(omega > 1.5* omega_max+2.0) KeepGoing=false; 
        }
        if (errno !=0) my_perror("halfway _for_stopping: an error occured");
        ResultArray[0]=(1.0/(inv_lambda+1.0e-99))*BohrRadius;  //IMFP now in Angstrom;
        ResultArray[1]=Hartree* (stopping)/BohrRadius;  //stopping now in eV/Angstrom
        ResultArray[2]= Hartree*Hartree* (straggling)/BohrRadius;  //straggling  now in eV^2/Angstrom
        if(DebugMode)
        {
            if (errno !=0) my_perror("DIIMFP_for_stopping at end: an error occured");
        }
        return 0;
    }
    
//// MODIF 19/08 LK:
    int DIIMFP_MELF(npArray ParameterArray, npArray npDIIMFPArray, int modelchoice, npArray npResultArray)
    {
       //diimfp calculated for all energy losses with first fixed (for DIIMFP distribution) and then  variable stepsize, (for stopping straggling and IMFP calculation).
       double* DIIMFPresult= npDIIMFPArray.data();
       double* ResultArray=  npResultArray.data();
       double inv_lambda=1E-99;  // so we never get divide by 0
       double stopping=0.0;
       double straggling=0.0;
       copyP_to_Vars(ParameterArray.data(), modelchoice);
       double omega= Egap + 0.5*StepSize;		//// !!! IMPORTANT !!! Egap is the new min energy for integration !!! was: double omega=0.5*StepSize;
       double Current_StepSize=StepSize;
       double Effective_Stepsize=Current_StepSize;
       double DIIMFP_au;
       double omega_max= 2*beta_r*beta_r*gamma_rel*gamma_rel*C*C;

       double MottFactor= (1- beta_r*beta_r*omega/omega_max); // correction factor SalvatPRA eq.16, not sure about this one, esp. for electrons
       
       bool KeepGoing=true;
       int i=0;
       while (KeepGoing)
       {    
            DIIMFP_au= DIIMFP_MELF_at_omega(omega);
           
 
            if(MottCorrection)DIIMFP_au *= MottFactor; // correction factor Salvat PRA  106 032809 eq.16
            inv_lambda     += DIIMFP_au * Effective_Stepsize;
            stopping       += DIIMFP_au * Effective_Stepsize * omega;
            straggling     += DIIMFP_au * Effective_Stepsize * omega * omega;
            if(i < NStep)
            {
                DIIMFPresult[i] = DIIMFP_au/ (Hartree*BohrRadius);
                i++;
            }
            else
            {   Effective_Stepsize =Current_StepSize;
                Current_StepSize=StepSize+lin_cont_deltaE*omega;
                Effective_Stepsize =(Effective_Stepsize+Current_StepSize)/2.0;  //average of previous and this stepsize
            }
            omega=omega+Current_StepSize;
            
            if (omega < omega_max) MottFactor= (1- beta_r*beta_r*omega/omega_max);// keep Mott factor the same for omega > omega_max
            if(omega > 1.5* omega_max+2.0) KeepGoing=false; 
        }
        if (errno !=0) my_perror("halfway _for_stopping: an error occured");
        ResultArray[0]=(1.0/(inv_lambda+1.0e-99))*BohrRadius;  //IMFP now in Angstrom;
        ResultArray[1]=Hartree* (stopping)/BohrRadius;  //stopping now in eV/Angstrom
        ResultArray[2]= Hartree*Hartree* (straggling)/BohrRadius;  //straggling  now in eV^2/Angstrom
        if(DebugMode)
        {
            if (errno !=0) my_perror("DIIMFP_for_stopping at end: an error occured");
        }
        return 0;
    }
//// END MODIF

 int  DIIMFP_variable_step(npArray ParameterArray, const npArray npomega_values,const npArray npstep_values, npArray npDIIMFPArray, int modelchoice, npArray npResultArray)
    {
       //diimfp calculated for all energy losses with first fixed (for DIIMFP distribution) and then  variable stepsize, (for stopping straggling and IMFP calculation).
       double *DIIMFPresult= npDIIMFPArray.data();
       double *omega_values= npomega_values.data();// should be in Hartrees
       double *step_values= npstep_values.data();// should be in Hartrees
       double *ResultArray=  npResultArray.data();
       double inv_lambda=1E-99;  // so we never get divide by 0
       double stopping=0.0;
       double straggling=0.0;
       double last_step, next_step;
       size_t omega_size = npomega_values.shape(0);
      
       copyP_to_Vars(ParameterArray.data(), modelchoice);
      
     //  double Effective_Stepsize = omega_values[1] - omega_values[0];
       double DIIMFP_au;
       double omega_max= 2*beta_r*beta_r*gamma_rel*gamma_rel*C*C;
       double MottFactor=1.0;
       bool KeepGoing=true;
       int i=0;
       double step=step_values[0];
       double omega=omega_values[0];
       while (KeepGoing)
       {    
            DIIMFP_au= DIIMFP_at_omega(omega);
           
            if(MottCorrection)
            {  
                if (omega < omega_max) MottFactor= (1- beta_r*beta_r*omega/omega_max);// keep Mott factor the same for omega > omega_max
                DIIMFP_au *= MottFactor; // correction factor Salvat PRA  106 032809 eq.16 
               
            }

            inv_lambda     += DIIMFP_au * step;
            stopping       += DIIMFP_au * step * omega;
            straggling     += DIIMFP_au * step * omega * omega;
            if(i < omega_size)  DIIMFPresult[i] = DIIMFP_au/ (Hartree*BohrRadius);
            i++;
            if(i<  omega_size)
            {  
                step =step_values[i];;  //average of previous and this stepsize
                omega=omega_values[i];
                last_step=step;
            }
            else
            { 
                 next_step = 1.1*last_step;
                 step = (next_step+last_step)/2.0;
                 omega += next_step;
                 last_step=next_step;
            }

            if((ExchangeCorrection  != 0) && (omega > 0.5*(E_0+BE_for_exchange))) KeepGoing=false; 
            if(omega > 1.5* omega_max+2.0) KeepGoing=false; 
        }
        ResultArray[0]=(1.0/(inv_lambda+1.0e-99))*BohrRadius;  //IMFP now in Angstrom;
        ResultArray[1]=Hartree* (stopping)/BohrRadius;  //stopping now in eV/Angstrom
        ResultArray[2]= Hartree*Hartree* (straggling)/BohrRadius;  //straggling  now in eV^2/Angstrom
        if(DebugMode)
        {
            if (errno !=0) my_perror("DIIMFP_for_stopping at end: an error occured");
        }
        return 0;
    }    
    
	double  DSEP(npArray ParameterArray, npArray np_DSEPresult,  int modelchoice)
	{
        //projectile=1: electron  projectile = 2: proton
        double *DSEPresult = np_DSEPresult.data();
		double sep;
 		copyP_to_Vars(ParameterArray.data(), modelchoice);
		sep=calc_DSEP(DSEPresult,E_0,theta0);
        if (errno !=0) my_perror("DSEP: an error occured");
        return sep;
	}
    
    void eps_Scaling_init(npArray ParameterArray, int modelchoice)  
    {   double GOSScalingStepSize=0.03; // used in eps_Scaling_init
        double TaucScalingStepSize=0.025; // used in eps_Scaling_init
        copyP_to_Vars(ParameterArray.data(), modelchoice);
        double currentFSum, currentE,current_q,current_q_stepsize,stepsize_au, E0_this_q, extra_width_due_to_disperion;
        double delta=0.02;
        
        for (int i = 0; i < MAXGOS; i++)
        {
            if (AiGOS[i] > 0)
            {   
                current_q=0.0;
                current_q_stepsize=GOSScalingStepSize;
                for(int iq=0;iq< NGOSscaling; iq++)
                {   
                    stepsize_au=Edgei_GOS[i]/5000;  // energy   step equal to 1/500 of the edge energy
                    
                    GOS_Scaling_q[iq]= current_q;
                    currentE=Edgei_GOS[i]+0.5*stepsize_au;
                    currentFSum = 0.0;
                    stepsize_au=Edgei_GOS[i]/5000;  // energy   step equal to 1/500 of the edge energy
                    while(currentE < (200*Edgei_GOS[i]+current_q*current_q/2.0))  // so we integrate up to 200 times the  edge energy+q^2/2
                    {
                        currentFSum += GOSx(n_i_GOS[i], l_i_GOS[i], Zi_GOS[i], Edgei_GOS[i], current_q, currentE)*stepsize_au*(1.0+delta/2.0);  // stepsize is average of previous+next step
                        stepsize_au=stepsize_au*(1.0+delta);
                        currentE=currentE+stepsize_au;
                    }
                    GOS_ScalingFactor[i][iq] = 1.0 / currentFSum;
                    if( iq % 5 == 0)printf(" i  %i,current_q %6.3f,scalingfactor, currentFSum %6.3f gos sc fac %6.4f\n",iq, current_q, currentFSum,GOS_ScalingFactor[i][iq] );
                    current_q +=current_q_stepsize;
                    current_q_stepsize *=1.05;
               }
            }        
        }
        if (errno !=0) my_perror("an error occured");
        if (modelTaucLorentz or modelTL_an or modelTauc_Mermin)
        {   
            for (int i = 0; i < MAXOSC; i++)
            {
                if (abs(Ai[i]) > 1e-50)
                {   current_q_stepsize=TaucScalingStepSize;
                    current_q=0.0;
                    for(int iq=0;iq<NTaucscaling; iq++)
                    {   
                        Tauc_Scaling_q[iq]= current_q;
                        double Qrecoil =recoil_energy(current_q);
                        if(delayed_dispersion)
                        {
                            double E0_without_gap= sqrt(wi[i]*wi[i]-TaucGap[i]*TaucGap[i]);
                            double E0_without_gap_this_q= E0_without_gap+alphai[i]*Qrecoil;
                            E0_this_q=sqrt(E0_without_gap_this_q*E0_without_gap_this_q+TaucGap[i]*TaucGap[i]);
                        }
                        else
                        { 
                             E0_this_q=wi[i]+alphai[i]*Qrecoil;
                        }
           //             extra_width_due_to_disperion=5*sqrt(E0_this_q-wi[i]);
                       extra_width_due_to_disperion=0.0;
                       
                        Tauc_ScalingFactor[i][iq] =  TaucNormalisation(Ai[i], gammai[i]+extra_width_due_to_disperion, wi[i],E0_this_q,TaucGap[i],aTLan,current_q, Ai[i],wi[i]);
                       
                        if( iq % 5 == 0)printf(" i  %i,current_q %6.3f, Tauc scalingfactor  %6.5g\n",iq, current_q, Tauc_ScalingFactor[i][iq] );
                        current_q +=current_q_stepsize;
                        current_q_stepsize *=1.05;
                   }
                }        
            }
        }
        return;
    }
    
    double TaucSumRule(double A, double  C, double E0,  double Egap)  // this does NOT call copyP_to_vars, called from calculate to get number of electrons per U.C.
    {     
        modelTaucLorentz = true; modelTL_an = false;   modelTauc_Mermin = false;
        double correctionfactor= TaucNormalisation(A, C, E0,E0,  Egap,0.0, 0.0,0.0,0.0);
        printf("TaucLorentz sum rule %6.4f\n",  correctionfactor);
        return correctionfactor;
    }
    double TL_an_SumRule(double A, double  C, double E0,  double Egap, double a_TL_analytical)  // this does NOT call copyP_to_vars,  called from calculate to get number of electrons per U.C.
    {   
        modelTauc_Mermin= false;  modelTL_an = true;  modelTaucLorentz = false;
        double  correctionfactor = TaucNormalisation(A, C, E0,E0,  Egap, a_TL_analytical, 0.0,0.0,0.0);
        printf("tl an sum rule %6.4f\n",  correctionfactor);
        return  correctionfactor;
    }   
    double TaucMermin_SumRule(double A, double  C, double E0,  double Egap)  // this does NOT call copyP_to_vars,  called from calculate to get number of electrons per U.C.
    {  
        double A_Drude= E0*A;
        double w_p= sqrt(A_Drude+E0*E0);
        double A_Mermin=A_Drude/(w_p*w_p);
        modelTauc_Mermin= true;  modelTL_an = false;  modelTaucLorentz = false;
        double  correctionfactor = TaucNormalisation(0.0, C, 0.0,0.0,  Egap, 0.0,0.0, A_Mermin,w_p);
        printf("tauc Mermin sum  rule %6.4f\n",  correctionfactor);
        return  correctionfactor;
    }   

     int  Kramers_Kronig_eps1_from_eps2(double FirstE, double DeltaE,dcomp eps_last,  int length,  npArray eps2Array, npArray eps1Array)  // does not depend on  copyP_to_Vars
     {  
        double omega_prime;
        double* eps2= eps2Array.data();
        double* eps1= eps1Array.data();  // calculated using kk from eps2
        printf("epslast %6.4f , %6.4f\n" , eps_last.real(),eps_last.imag());
        printf("firstE %6.4f  deltaE %6.4f\n",FirstE, DeltaE);
        printf("eps2[4]  %6.4f, eps2[13] %6.4f, eps2[22] %6.4f\n",eps2[4],eps2[13],eps2[22]);
        printf("eps1[4]  %6.4f, eps1[13] %6.4f, eps1[22] %6.4f\n",eps1[4],eps1[13],eps1[22]);
        for (int i = 0; i < 100; i++) printf("i %i eps2[i]  %6.4f\n",i,eps2[i]); 
        for (int i = 0; i < length; i++)
        {
           double KKsum = 0.0;
           double omega = i * DeltaE + FirstE;  
           for (int k = 0; k < length; k++)
           {
                omega_prime = k * DeltaE + FirstE;
                if (i != k)
                     KKsum += (omega_prime*(eps2[k]-eps_last.imag()) / (omega_prime * omega_prime - omega * omega)) * DeltaE;// Wooten eq. 6.52 
           }                                                                                                               
           eps1[i]  =2.0* KKsum / pi + eps_last.real(); //the factor 2 because we go from 0 to large rather than -large to large
        }
        if (errno !=0) my_perror("Kramers_Kronig_eps1_from_eps2: an error occured");
        return 0;
     }  
     
    int  Kramers_Kronig_eps2_from_eps1(double FirstE, double DeltaE,dcomp eps_last,  int length, const npArray eps1Array, npArray eps2Array) // does not depend on  copyP_to_Vars
     {  
        double omega_prime;
        double* eps1= eps1Array.data();
        double* eps2= eps2Array.data();  // calculated using kk from eps1
        for (int i = 0; i < length; i++)
        {
           double KKsum = 0.0;
           double omega= i * DeltaE + FirstE;  
           for (int k = 0; k < length; k++)
           {
                omega_prime = k * DeltaE + FirstE;
                if (i != k)
                     KKsum += omega * (eps1[k]- eps_last.real()) / (omega_prime*omega_prime - omega*omega) * DeltaE;// Wooten eq. 6.52 
           }                                                                                                               
           eps2[i]  = - 2.0* KKsum / pi+eps_last.imag(); //the factor 2 because we go from 0 to large rather than -large to large
        }
        if (errno !=0) my_perror("Kramers_Kronig_eps2_from_eps1: an error occured");
        return 0;
     }  
     
      
double   stopping_longitudinal(npArray ParameterArray,npArray np_DDCS,  int modelchoice) // stopping from DCS, no retardation
    {   
        double *DDCS= np_DDCS.data();
        double  omega, stopping, DDCS_au;
        int ienergy,itheta;
        copyP_to_Vars( ParameterArray.data(),  modelchoice); 
        double crosssection=0.0;
        double theta_step=theta_max/NThetaStep;
        double prefactor = 1.0/(pi * pi * v_0 * v_0 * UnitCellDensity);
        stopping=0.0;
        for(itheta=0; itheta < NThetaStep; itheta++)
        {
            double theta=(0.5+itheta)*theta_step;
           
            for(ienergy =0;ienergy < NStep; ienergy++)
            {   
                omega= FirstEnergy +ienergy*StepSize;
                DDCS_au= prefactor * DDCS_excl_retardation(omega,theta);
                crosssection+= 2 * pi * sin(theta) * theta_step * DDCS_au *StepSize; 
                stopping += 2 * pi * sin(theta) * theta_step * DDCS_au * omega * StepSize; 
                DDCS[itheta*NStep+ienergy] = DDCS_au * (BohrRadius*BohrRadius)/Hartree;// convert from a.u.^2/Hartree to Angstrom^2/eV    ;
            }
        }
        crosssection /=(BohrRadius*BohrRadius);
        printf(" cross section %6.4e\n",crosssection);
        stopping *= UnitCellDensity;// now in a.u. energy/ a.u. distance
        stopping *= Hartree / BohrRadius ;  // now in eV/Angstrom
        if (errno !=0) my_perror("DDCS_longitudinal: an error occured"); 
        
        return stopping;
    }
    
double   stopping_total(npArray ParameterArray,npArray np_DDCS,   int modelchoice)  // including Cerenkov, i.e. retardation effects
    {   double *DDCS= np_DDCS.data();
        double  omega, stopping, DDCS_au;
        int ienergy,itheta;
        copyP_to_Vars( ParameterArray.data(),  modelchoice); 
        if ( ExchangeCorrection > 0) printf("echange correction not implemented for the case with retardation\n");
        stopping=0.0;
        double theta_step=theta_max/NThetaStep;
        double prefactor = 1.0/(pi * pi * v_0 * v_0 * UnitCellDensity);
        for(itheta=0; itheta < NThetaStep; itheta++)
        {
            double theta=(0.5+itheta)*theta_step;
            for(ienergy =0;ienergy < NStep; ienergy++)
            {
                omega= FirstEnergy +ienergy*StepSize;
                DDCS_au= prefactor * DDCS_incl_retardation(omega,theta);
                stopping += 2 * pi * sin(theta) * theta_step * DDCS_au * omega * StepSize; 
                DDCS[itheta*NStep+ienergy] = DDCS_au * (BohrRadius*BohrRadius)/Hartree;// convert from a.u.^2/Hartree to Angstrom^2/eV    ;// convert from a.u.^2/Hartree to Angstrom^2/eV    ;
            }
        }
        stopping *=  UnitCellDensity;// now in a.u. energy/ a.u. distance
        stopping *= Hartree / BohrRadius ;  // now in eV/Angstrom
        if (errno !=0) my_perror("DDCS_total: an error occured");    
        return stopping;
    }    
    
int  DDCS_at_theta(npArray ParameterArray, npArray np_DDCS, npArray  np_DDCS_incl_ret,  int modelchoice, double theta)  
    {   double  result_au, result_au_incl_ret;
        double *DDCS = np_DDCS.data();
        double *DDCS_incl_ret = np_DDCS_incl_ret.data();
        
        copyP_to_Vars( ParameterArray.data(),  modelchoice); 
        theta=theta/1000.0; // now in rad  
        double prefactor = 1.0/(pi * pi * v_0 * v_0 * UnitCellDensity);
       
      
        for(int ienergy =0;ienergy < NStep; ienergy++)
            {                   
                double omega = FirstEnergy +ienergy*StepSize;
                result_au = prefactor *  DDCS_excl_retardation(omega, theta);  //Egerton 3.32
                DDCS[ienergy]=result_au*(BohrRadius*BohrRadius)/Hartree;// convert from a.u.^2/Hartree to Angstrom^2/eV   
                result_au_incl_ret = prefactor *  DDCS_incl_retardation(omega, theta);  //Egerton 3.32
                DDCS_incl_ret[ienergy]=result_au_incl_ret*(BohrRadius*BohrRadius)/Hartree;// convert from a.u.^2/Hartree to Angstrom^2/eV   
            }    
        if (errno !=0) my_perror("DDCS_at_theta: an error occured");
        return 0;    
    } 
    
int  DDCS_at_omega(npArray ParameterArray, npArray np_xaxis, npArray np_DDCS,  npArray np_DDCS_incl_ret,  int modelchoice,  double omega)
    {   
        double *xaxis = np_xaxis.data();
        double *DDCS = np_DDCS.data();
        double *DDCS_incl_ret = np_DDCS_incl_ret.data();
        
        double  result_au, result_au_incl_ret;
        int itheta;
        copyP_to_Vars( ParameterArray.data(),  modelchoice); 
    
        double prefactor = 1.0/(pi * pi * v_0 * v_0 * UnitCellDensity);
        omega=omega/Hartree;
        for(itheta=0; itheta < NThetaStep; itheta++)
        {
            double theta=xaxis[itheta]/1000.0;   // in rad          
            
            result_au = prefactor *  DDCS_excl_retardation(omega, theta);  //Egerton 3.32
            DDCS[itheta]=result_au*(BohrRadius*BohrRadius)/Hartree;// convert from a.u.^2/Hartree to Angstrom^2/eV   
            result_au_incl_ret = prefactor *  DDCS_incl_retardation(omega, theta);  //Egerton 3.32
            DDCS_incl_ret[itheta]=result_au_incl_ret*(BohrRadius*BohrRadius)/Hartree;// convert from a.u.^2/Hartree to Angstrom^2/eV   
        }
        if (errno !=0) my_perror("DDCS_at_omega: an error occured");
        return 0;
    }
    double DCS(npArray ParameterArray, npArray np_xaxis, npArray np_DCS, npArray np_DCS_incl_ret, int modelchoice)
    {
        int itheta;
        double *xaxis = np_xaxis.data();
        double *DCS = np_DCS.data();
        double *DCS_incl_ret = np_DCS_incl_ret.data();
        copyP_to_Vars( ParameterArray.data(),  modelchoice);  
        double prefactor = 1.0/(pi * pi * v_0 * v_0 * UnitCellDensity);
        double stopping= 0.0;
        double tot_cross=0.0;

        double firststep=xaxis[0]-xaxis[1];
        double secondstep=xaxis[1]-xaxis[2];
        double incrementfactor=firststep/secondstep;
        double omega_max=FirstEnergy +NStep*StepSize;
        if (omega_max > E_0) omega_max = E_0;
        
        double thetastep=firststep*(1.0+incrementfactor)/2;  // as the previous step was larger
        for(itheta=0; itheta < NThetaStep; itheta++)
        {
             DCS[itheta]=0;
             DCS_incl_ret[itheta]=0;
            double sintheta= sin(xaxis[itheta]);
            double omega=FirstEnergy;

            double CurrentEStepSize=StepSize;
            while (omega < omega_max)  
            {    
                double theta=xaxis[itheta];   
                double result_au = prefactor *  DDCS_excl_retardation(omega, theta);
                double increment= result_au * CurrentEStepSize;
                DCS[itheta] += increment;
                tot_cross+=increment*thetastep*2*pi*sintheta;
                stopping+=increment*thetastep*2*pi*sintheta*omega;
                
                double result_incl_au = prefactor *  DDCS_incl_retardation(omega, theta);
                
                double increment_incl= result_incl_au * CurrentEStepSize;
                DCS_incl_ret[itheta] += increment_incl;
                CurrentEStepSize *= 1.025;
                omega+=CurrentEStepSize;
            }
            thetastep/=incrementfactor;
        }
 
        if (errno !=0) my_perror("DCS: an error occured");

        
        // Now put things in eV and Angstrom units
        for(itheta=0; itheta < NThetaStep; itheta++)
        {   DCS[itheta]=DCS[itheta]*BohrRadius*BohrRadius;
            DCS_incl_ret[itheta]=DCS_incl_ret[itheta]*BohrRadius*BohrRadius;
        }
        double lambda=BohrRadius/(tot_cross*UnitCellDensity);
        printf(" cs  %6.4e  lambda %6.4e stopping %6.4e\n",tot_cross*(BohrRadius*BohrRadius),lambda, 
            stopping*UnitCellDensity*Hartree /BohrRadius);
        return tot_cross*BohrRadius*BohrRadius;
    }
    
    
     int  calc_REELS(npArray ParameterArray,   double StartReels, int N_REELS_Step, int modelchoice, npArray np_REELSresult, 
         int EELS, double EELSPathLength)
//works in eV
	{
        
	double lambda; 
	//projectile=1: electron  projectile = 2: proton
	double GE,arg;
    double sprob[3];
    double *REELSresult=np_REELSresult.data();
	copyP_to_Vars( ParameterArray.data(),  modelchoice);  
    double  *DIIMFPresult = new double [NStep];
    double  *DSEPresultIn = new double [NStep];
    double  *DSEPresultOut = new double [NStep];
    double  *NormDSEPresult = new double [NStep];
    double  *normDIIMFP = new double[NStep];
    double  *vSignalconv = new double[2*NStep];

      // Dimensions of the 3D array    
    int MaxPartInt=10;
    double *PartialIntensity= new double [MaxPartInt];
    int MaxSurfEx=3;

   
    double*** matrix2 = new double**[N_REELS_Step];
 
    for (int i = 0; i < N_REELS_Step; i++) {
 
        // Allocate memory blocks for
        // rows of each 2D array
        matrix2[i] = new double*[MaxPartInt];
 
        for (int j = 0; j < MaxPartInt; j++) {
 
            // Allocate memory blocks for
            // columns of each 2D array
            matrix2[i][j] = new double[MaxSurfEx];
        }
    }
    FirstEnergy=0.5* StepSize;  // This makes sure FirstEnergy is sensible, not what in user interface
    double inv_lambda  =1.0e-99;
    for (int i = 0; i < NStep; i++)
		{
			double omega    = FirstEnergy + StepSize*i;
            double DIIMFP_au= DIIMFP_at_omega(omega);
            DIIMFPresult[i] = DIIMFP_au/ (Hartree*BohrRadius);// now eV/angstrom
            inv_lambda     +=  DIIMFPresult[i]*StepSize_eV;// now in angstrom^{-1}
        }
    lambda=1.0/inv_lambda;
    if(EELS == 0)
    {   
        for (int j = 0; j < MaxPartInt; j++)
            PartialIntensity[j] = 1.0 + PIcoef1*j + PIcoef2*j*j + PIcoef3*j*j*j; 
    }   
    else
    {   
        double Lambda_Poisson=EELSPathLength/lambda;
        double zero_loss=Poisson(0, Lambda_Poisson);
        for (int j = 0; j < MaxPartInt; j++)
        {   
            PartialIntensity[j] =Poisson(j, Lambda_Poisson)/zero_loss;  // so elastic peak intensity = 1 
        }   
    } 
    for (int i=0; i< NStep; i++)
    {
        normDIIMFP[i]=DIIMFPresult[i]*lambda*fraction_DIIMFP;
    }
    double SurfExProbin  = calc_DSEP(DSEPresultIn,E_0, theta0);
    double SurfExProbout = calc_DSEP(DSEPresultOut,E_0, theta1);
    if(SurfExProbin > 1.0) SurfExProbin=1.0;
    if(SurfExProbout> 1.0) SurfExProbout=1.0;

    if(surf_ex_factor > 0.0 && SurfExProbin!=0.0 && SurfExProbout !=0.0)
    {
        for (int i=0; i< NStep; i++)
        { // take the average of the two similar DSEPs
            NormDSEPresult[i]=0.5*StepSize*Hartree*(DSEPresultIn[i]/SurfExProbin+DSEPresultOut[i]/SurfExProbout);
        }
        SurfExProbin=SurfExProbin*surf_ex_factor;
        SurfExProbout=SurfExProbout*surf_ex_factor;
        sprob[0]= (1.0-SurfExProbin)*(1.0-SurfExProbout);
        sprob[1]= (1.0-SurfExProbin)*SurfExProbout+(SurfExProbin)*(1.0-SurfExProbout);
        sprob[2]= SurfExProbin*SurfExProbout;
    }
     else
     {
        for (int i=0; i< NStep; i++)
        { // take the average of the two similar DSEPs
            NormDSEPresult[i]=0.0;
        }
        sprob[0]=1.0; sprob[1]=0.0; sprob[2]=0.0;
    }
   

    for (int ipartint = 0; ipartint <MaxPartInt; ipartint++)
	{
		if (ipartint == 0)  // fill up the first element with the spectrometer response,  a Gaussian
		{
			for (int i = 0; i < N_REELS_Step; i++)
			{
				GE = StartReels + StepSize_eV*i;

				//GE = StepSize_eV*(i - 20);
				arg = GE*GE / (2.0*sigma*sigma);
				if (arg < 40)
				{
					matrix2[i][0][0] = exp(-arg); // take maximum to be 1, norm exp elastic peak to 1 for easy comp
				}
				else
				{
					matrix2[i][0][0] = 0.0;
				}
			}
			
		}
		else  //  fill up subsequent elements with  a Gaussian convoluted i times with NDIIMFP
		{
			for (int j = 0; j < 2 * NStep; j++) vSignalconv[j] = 0.0;
			for (int i = 0; i < NStep; i++)
				for (int j = 0; j < NStep; j++){
					vSignalconv[j + i] = vSignalconv[j + i] + matrix2[i][ipartint - 1][0] * normDIIMFP[j] * StepSize_eV;
				}

			// put in matrix2[i][ipartint][0]
			for (int j = 0; j < N_REELS_Step; j++) matrix2[j][ipartint][0] = vSignalconv[j];
		}
// add convolution with normalised DSEP to next columns    

        for (int j = 0; j < 2 * NStep; j++) vSignalconv[j] = 0.0;
        for (int i = 0; i < NStep; i++){
            
            for (int j = 0; j < NStep; j++){
                vSignalconv[j + i] = vSignalconv[j + i] + matrix2[i][ipartint][0] * NormDSEPresult[j];
            }
        }
        // put in matrix2[i][ipartint][1]
        for (int j = 0; j < N_REELS_Step; j++) matrix2[j][ipartint][1] = vSignalconv[j];


        for (int j = 0; j < 2 * NStep; j++) vSignalconv[j] = 0.0;
        for (int i = 0; i < NStep; i++){
            for (int j = 0; j < NStep; j++){
                vSignalconv[j + i] = vSignalconv[j + i] + matrix2[i][ipartint][1] * NormDSEPresult[j];
            }
        }
        // put in matrix2[i][ipartint][2]
        for (int j = 0; j < N_REELS_Step; j++) matrix2[j][ipartint][2] = vSignalconv[j];
    }
    fstream myfile;
    myfile.open("loss_dist.txt",fstream::out);
    myfile << "Energy:\t No surf loss, \t surf loss\t two surf loss" <<std::endl; 
    
    for (int j = 0; j < N_REELS_Step; j++)
	{
        double no_surfplasmon,onesurfplasmon,twosurfplasmon; 
        no_surfplasmon=0.0;onesurfplasmon=0.0;twosurfplasmon=0.0; 
		for (int jj = 0; jj < MaxPartInt; jj++)
        {
            no_surfplasmon+=PartialIntensity[jj] * matrix2[j][jj][0]; // for convenience we keep max height elastic peak at 1
            onesurfplasmon+=(sprob[1]/sprob[0])*PartialIntensity[jj] * matrix2[j][jj][1];
            twosurfplasmon+= (sprob[2]/sprob[0])*PartialIntensity[jj] * matrix2[j][jj][2];
        }
        REELSresult[j] =  no_surfplasmon+onesurfplasmon+twosurfplasmon; 
        myfile << StartReels+j*StepSize_eV << "\t"<< no_surfplasmon << "\t"<<onesurfplasmon<< "\t"<<twosurfplasmon<<std::endl;  
    }    
    myfile.close();
    myfile.open("PartInt.txt",fstream::out);
    myfile << "Part.Int.:\t"; 
    for (int j=0; j< MaxPartInt;j++)// Prints row of x
    {        
        myfile << j << "\t";  
    }

    myfile<< std::endl;

    for (int j=0; j< N_REELS_Step;j++) //This variable is for each row below the x 
    {        
        myfile << StartReels+j*StepSize_eV << "\t";

        for (int jj=0; jj<MaxPartInt;jj++)
        {                      
            myfile <<matrix2[j][jj][0] << "\t";
        }
        myfile<<std::endl;
    }
    myfile.close();

    delete[] NormDSEPresult;
    delete[] DSEPresultIn;
    delete[] DSEPresultOut;
    delete[] DIIMFPresult;
    delete[] vSignalconv;
    delete[] normDIIMFP;
    for (int i = 0; i < N_REELS_Step; i++) {
        for (int j = 0; j < MaxPartInt; j++) {
            delete[] matrix2[i][j];
        }
        delete[] matrix2[i];
    }
    delete[] matrix2;
    if (errno !=0) my_perror("calc_REELS: an error occured");
    return 0;
    }
    
void cumulative_trapezoid(const npArray y_array,const  npArray x_array, npArray result_array, double begin)
{   size_t length = x_array.size();// all arrays need to have the same length
    printf("integrating length array  %i\n", length);
    double *x=x_array.data();
    double *y=y_array.data();
    double *result=result_array.data();
    

    result[0]=y[0]*(x[0]-begin);// because it is the cumulative sum starting as begin, not starting at x[0]. This assumes y=constant for x< x[0]
    for (int i=1; i< length-1 ;i++)
    { double contribution= 0.5*(y[i]+y[i-1])*(x[i]-x[i-1]);
       result[i]=result[i-1]+0.5*(y[i]+y[i-1])*(x[i]-x[i-1]);
    }
   
    return;
};

void test(){
    printf("nanobind seems correctly installed!\n");
        std::cout << "The C library was compiled by: ";

#if defined(__clang__)
    std::cout << "Clang (version " << __clang_version__ << ")" << std::endl;
#elif defined(__GNUC__)
    std::cout << "GCC (version " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__ << ")" << std::endl;
   // does not work printf("GCC (version  %s . %s .%s )\n", __GNUC__ , __GNUC_MINOR__ , __GNUC_PATCHLEVEL__ );
#elif defined(_MSC_VER)
#elif defined(_MSC_VER)
    std::cout << "Microsoft Visual C++ (version " << _MSC_VER << ")" << std::endl;
#elif defined(__VERSION__)
    // A generic fallback for other compilers that define __VERSION__ as a string
    std::cout << "Unknown/Other Compiler (version string: " << __VERSION__ << ")" << std::endl;
#else
    std::cout << "An unknown or highly unusual compiler." << std::endl;
#endif
}
NB_MODULE(epslib, m) {
    m.def("test", &test);
    m.def( "Eps1Eps2", &Eps1Eps2,"calculate real and imag part eps at fixed q");
    m.def( "single_Eps1Eps2", &single_Eps1Eps2,"calculate real and imag part eps for one omega, q combination");
    m.def( "Eps1Eps2_q", &Eps1Eps2_q,"calculate real and imag part eps at fixed omega");
    m.def("eps_Scaling_init",&eps_Scaling_init, "calculate q dependent normalisation factors");  
    m.def("Kramers_Kronig_eps1_from_eps2", &Kramers_Kronig_eps1_from_eps2, "calculate eps1 from eps2 using Kramers Kronig");
    m.def("Kramers_Kronig_eps2_from_eps1", &Kramers_Kronig_eps2_from_eps1, "calculate eps2 from eps1 using Kramers Kronig");
    m.def("DIIMFP", &DIIMFP,"calc differential inverse mean free path");
    m.def("DIIMFP_MELF", &DIIMFP_MELF,"calc differential inverse mean free path as the sum of independent per-oscillator ELF contributions (MELF)");
    m.def("Loss_wide",&Loss_wide,"loss array over wide range, non-linear grid");
    m.def("DDCS_at_omega",&DDCS_at_omega,"double diff cross section at energy omega");
    m.def("DDCS_at_theta",&DDCS_at_theta, "double diff cross section at angle theta");
    m.def("DCS",&DCS, "diff cross section");
    m.def("stopping_total",&stopping_total, "total stopping (incl cerenkov from dcs");
    m.def("stopping_longitudinal", &stopping_longitudinal, "stopping (without cerenkov) from dcs");
    m.def("SurfLossFunc", &SurfLossFunc, "calculate surface loss function");
    m.def("DSEP", &DSEP, "calculate differential surface excitation probability");
    m.def("calc_REELS",&calc_REELS,"calculate an (r)eels spectrum");
    m.def("TaucMermin_SumRule", &TaucMermin_SumRule, "calculate sum rule for Tauc-truncated mermin DF");
    m.def("TL_an_SumRule", &TL_an_SumRule, "calculate sum rule for Tauc_analytic DF"); 
    m.def("TaucSumRule", &TaucSumRule, "calculate sum rule for Tauc DF");
    m.def("DIIMFP_variable_step",&DIIMFP_variable_step,"calculate DIIMFP with, if required, variable step size");
    m.def("Lossfunction_inclGOS_atE", &Lossfunction_inclGOS_atE, "calculate loss function at single omega value");
   // m.def("cumulative_trapezoid",cumulative_trapezoid,"modifiied cumulative trapezoid integration");
}
