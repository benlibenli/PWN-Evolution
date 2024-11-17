// Develpoed by Ben Li in March 2023
// The evolution of SNR, RS and PWN is mainly based on the model by Joseph D. Gelfand, et al. (2009)
#include "arb.h"
#include "arf.h"
#include "arb_hypgeom.h"
#include<iostream>
#include<fstream>
#include <iomanip>
#include<cmath>
using namespace std;
#define pi 3.1415926
#define yr 31536000.0
#define pc 3.08e18
#define c 3e10
#define mH 1.674e-24
#define M_sun 1.98847e33

// pulsar name: J1834-0845

// Spin-down Luminosity
double n0 = 3.0;
double n = (n0+1.0)/(n0-1.0);
double P = 2.48;
double Pdot = 7.96e-12;
double Edot = 2.1e34;
double tauc = P/2.0/Pdot/yr;
double P0 = 0.1;
double Tage = 2.0*tauc/(n0-1.0)*(1.0-pow(P0/P, n0-1.0));
double tau0 = 2.0*tauc/(n0-1.0) - Tage;
double Ls0 = Edot*pow(P/P0, n0+1.0);  // Ls0 = Edot*pow(1.0+Tage/tau0, (n0+1.0)/(n0-1.0));


// Magnetic luminosity
double Rns = 1.2e6;
double B_int = 1e16;
double EB_int = pow(B_int, 2.0)*pow(Rns, 3.0)/6.0;
double t0 = tau0;  // onset of energy injection
double alpha = n;
double Lm0 = (alpha-1.0)*EB_int/(t0*yr);

// Magnetism parameter
double eta_B = 0.1;
double sigma = eta_B/(1.0-eta_B);


// SNR parameters
double E_SN = 1e51;
double M_ej = 11.3*M_sun;
double V0 = sqrt(10.0/3.0*E_SN/M_ej);
double alphao = 6.0/5.0;
double coeff[5] = {1.0, 0.2-0.044898*n, 0.12-0.045438*n+0.004434*n*n, 0.088-0.043397*n+0.007204*n*n-0.000383*pow(n, 3.0), 0.0704-0.041023*n+0.008888*n*n-0.000802*pow(n, 3.0)+0.000024*pow(n, 4.0)};
double rho_ISM = 0.5*mH;
double R_ch = pow(M_ej, 1.0/3.0)*pow(rho_ISM, -1.0/3.0);  // characteristic length
double t_ch = pow(E_SN, -1.0/2.0)*pow(M_ej, 5.0/6.0)*pow(rho_ISM, -1.0/3.0) / yr;  // characteristic time
double t_ST = 0.52*t_ch;  // time when SNR enters Sedov phase
double t_core = 0.25*t_ch;  // time when reverse shock enters the constant density core
double t_r;  // time when reverse shock hits the pwn
double t_s;  // time when the pressure inside the PWN meets the Sedov solution 


int main(){
    double Ls(double t_);
    double Lm(double t_);
    double L(double t_);
    // Radius evolution of SNR and reverse shock
    double R_SNR(double t_);
    double R_RS(double t_);

    // Radius evolution of PWN
    double R_FREE(double t_);
    double P_FREE(double t_);
    // double R_REV(double t_);
    void R_REV(double* Rrev_, int num_);
    double P_REV(double t_, double* Rrev_);
    double R_ST(double t_, double R_t_s, double P_in_t_s);
    // double R_PWN(double t_);
    double P_PWN(double t_, double* Rrev_);

    // Radius evolution of TS
    double PIfun(double Rts_, double Rpwn_);
    double Hypgeom_2F1(int i, double t_);
    double Qfun(double t_);
    double fun(double t_, double Rts_, double Rpwn_);
    double R_TS(double t_, double Rpwn_);
    
    ofstream fout;

    cout << setprecision(2);
    cout << setiosflags(ios::fixed) << "tauc=" << tauc << "\tTage=" << Tage << "\ttau0=" << tau0 << resetiosflags(ios::fixed);
    cout << setiosflags(ios::scientific) << "\tLs0=" << Ls0 << "\tLm0=" << Lm0 << endl << resetiosflags(ios::scientific);

    double t[20000];
    for(int i=0; i<20000; i++)
        t[i] = i+1.0;


    // Evolution of SNR and reverse shock
    double Rsnr[20000];
    fout.open("Rsnr.csv");
    for(int i=0; i<20000; i++){
        Rsnr[i] = R_SNR(t[i]);
        fout << t[i] << "," << Rsnr[i] << endl;
    }
    fout.close();

    double Rrs[20000];
    fout.open("Rrs.csv");
    for(int i=0; i<20000; i++){
        Rrs[i] = R_RS(t[i]);
        fout << t[i] << "," << Rrs[i] << endl;
    }
    fout.close();

    // Evolution of PWN
    // 1. Free expansion phase
    double Rpwn[20000];
    for(int i=0; i<20000; i++){
        Rpwn[i] = R_FREE(t[i]);
    }

    for(int i=0; i<20000; i++){
        if( Rpwn[i] > Rrs[i] ){
            t_r = t[i];
            cout << setiosflags(ios::fixed) << "t_r=" << t_r << endl << resetiosflags(ios::fixed);
            break;
        }
    }
    
    // 2. Reverberation phase
    int num = round(20000.0-t_r)+2;
    double Rrev[num];
    R_REV(Rrev, num);

    int index1 = round(t_r);
    for(int i=index1; i<20000; i++){
        Rpwn[i] = Rrev[i-index1+2];
    }

    for(int i=0; i<20000; i++){
        double v_fs = (2.0/5.0)*sqrt(2.026*E_SN/rho_ISM)*pow(pow(R_SNR(t_ST), 5.0/2.0)+sqrt(2.026*E_SN/rho_ISM)*(t[i]-t_ST)*yr, -3.0/5.0);
        double P_Sed = 0.75*rho_ISM*pow(v_fs, 2.0);  // pressure at the SNR forward shock
        if( P_Sed < P_PWN(t[i], Rpwn) ){
            t_s = t[i];
            cout << setiosflags(ios::fixed) << "t_s=" << t_s << endl << resetiosflags(ios::fixed);
            break;
        }
    }

    // 3. Sedov-Taylor phase
    int index2 = round(t_s);
    for(int i=index2; i<20000; i++){
        Rpwn[i] = R_ST(t[i], Rpwn[index2-1], P_PWN(t_s, Rpwn));
    }

    fout.open("Rpwn.csv");
    for(int i=0; i<20000; i++)
        fout << t[i] << "," << Rpwn[i] << endl;
    fout.close();


    // Evolution of termination shock

    // double Rts[20000];
    // fout.open("Rts.csv");
    // for(int i=0; i<20000; i=i+2){
    //     Rts[i] = R_TS(t[i], Rpwn[i]);
    //     fout << t[i] << "," << Rts[i] << endl;
    //     cout << i << "\t" << "TS" << "\t" << Rts[i] <<endl;
    // }
    // fout.close();

    double Rts[20000];
    fout.open("Rts.csv");
    for(int i=0; i<20000; i++){
        Rts[i] = sqrt( L(t[i])/(4.0*pi*c*P_PWN(t[i], Rpwn)) );
        fout << t[i] << "," << Rts[i] << endl;
        // cout << i << "\t" << "TS" << "\t" << Rts[i] <<endl;
    }
    fout.close();


    // Evolution of pressure
    double Ppwn[20000];
    fout.open("Ppwn.csv");
    for(int i=0; i<20000; i++){
        Ppwn[i] = P_PWN(t[i], Rpwn);
        fout << t[i] << "," << Ppwn[i] << endl;
    }
    fout.close();


    // Evolution of the magnetic field
    double Einj[20000] = {0.0}, WB[20000] = {0.0}, B0[20000] = {0.0}, beta3 = -0.629, dt = yr;
    for(int i=1; i<20000; i++){
        Einj[i] = Einj[i-1] + L(t[i])*dt;
        WB[i] = eta_B*L(t[i])*dt - (Rpwn[i]/Rpwn[i-1]-2.0)*WB[i-1];
        // B0[i] = sqrt( 6.0*WB[i]/(pow(Rpwn[i], 3.0)-pow(Rts[i], 3.0)) );  // B=B0*pow(r/rts,beta), beta3=0
        B0[i] = sqrt( 2.0*(3.0+2.0*beta3)*pow(Rts[i], 2.0*beta3)*WB[i]/(pow(Rpwn[i], 3.0+2.0*beta3)-pow(Rts[i], 3.0+2.0*beta3)) );  // B=B0*pow(r/rts,beta3), beta3
    }
    cout << setiosflags(ios::scientific) << "Einj=" << Einj[int(Tage)] << "\tWB=" << WB[int(Tage)] << resetiosflags(ios::scientific);
    cout << setiosflags(ios::fixed) << "\t B_Tage=" << B0[int(Tage)]*1e6 << endl << resetiosflags(ios::fixed);
    fout.open("B0.csv");
    for(int i=1; i<20000; i++)
        fout << t[i] << "," << B0[i] << endl;
    fout.close();


    cout << "Done!" <<endl;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////FUNCTION////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
double Ls(double t_){  // spin-down luminosity
    return( Ls0*pow(1.0+t_/tau0, -n) );
}

double Lm(double t_){  // magnetic luminosity
    // return( Lm0*pow(t_/t0, -alpha) );
    return( Lm0*pow(1.0+t_/t0, -alpha) );
}

double L(double t_){
    // return( Ls(t_)+Lm(t_) );
    return( Ls(t_) );
}

double R_SNR(double t_){
    if(t_ <= t_ST)
        return( 1.12*R_ch*pow(t_/t_ch, 2.0/3.0) );
    else
        return( pow(pow(R_SNR(t_ST), 5.0/2.0)+sqrt(2.026*E_SN/rho_ISM)*(t_-t_ST)*yr, 2.0/5.0) );
}

double R_RS(double t_){
    if(t_ <= t_core)
        return(R_SNR(t_)/1.19);
    else
        return( (1.49-0.16*(t_-t_core)/t_ch-0.46*log(t_/t_core))*(t_/t_ch)*R_ch );
}

double R_FREE(double t_){
    double R0, s, R_FREE_temp;
    R0 = pow(pow(E_SN, 3.0)*pow(Ls0, 2.0)/pow(M_ej, 5.0), 1.0/10.0)*pow(t_*yr, 6.0/5.0);
    // R0 = pow(pow(E_SN, 3.0)*pow(Lm0+Ls0, 2.0)/pow(M_ej, 5.0), 1.0/10.0)*pow(t_*yr, 6.0/5.0);
    s = (t_/tau0)/(1.0+t_/tau0);
    R_FREE_temp = R0*pow(tau0/(tau0+t_), 6.0/5.0)/(1.0-s);

    double temp = 0.0;
    for(int j=0; j<5; j++)
        temp += coeff[j]*pow(s, j);
    return(R_FREE_temp*temp);
}

double P_FREE(double t_){
    double P1 = 0.0;
    double dt = yr;
    for(int i=1; i<=round(t_); i++)
        P1 += L(i)*R_FREE(i)*dt;
    P1 = P1/(4.0*pi*pow(R_FREE(t_), 4.0));
    return(P1);
}

// double R_REV(double t_){
//     int num1 = round(t_r), num2 = round(t_), num = num2-num1+2;
//     double dt = yr;
//     double R[num];
//     R[0] = R_FREE(t_r-1.0);
//     R[1] = R_FREE(t_r);
//     double rho_ej = 3.0*M_ej/(4.0*pi*pow(V0*t_r*yr, 3.0));  // density of the ejecta
//     double M_sw = (4.0*pi/3.0)*rho_ej*pow(R[1], 3.0);  // ejecta mass swept up by the pulsar wind nebula
//     double P_in, P_out;
//     for(int i=1; i<num-1; i++){
//         double temp = 0.0;
//         for(int j=0; j<i; j++){
//             temp += L(num1+j)*R[j+1]*dt;
//         }
//         P_in = P_FREE(t_r) + temp/(4.0*pi*pow(R[i], 4.0));
//         P_out = 0.037*E_SN/pow(R_SNR(num1+i), 3.0);  // 50% the Sedov solution
//         R[i+1] = 4.0*pi*pow(R[i], 2.0)*pow(dt, 2.0)/M_sw * (P_in-P_out) + 2.0*R[i] - R[i-1];
//     }
//     return( R[num-1] );
// }

void R_REV(double* Rrev_, int num_){
    double dt = yr;
    Rrev_[0] = R_FREE(t_r-1.0);
    Rrev_[1] = R_FREE(t_r);
    double rho_ej = 3.0*M_ej/(4.0*pi*pow(V0*t_r*yr, 3.0));  // density of the ejecta
    double M_sw = (4.0*pi/3.0)*rho_ej*pow(Rrev_[1], 3.0);  // ejecta mass swept up by the pulsar wind nebula
    double P_in, P_out;
    for(int i=1; i<num_-1; i++){
        double temp = 0.0;
        for(int j=0; j<i; j++){
            temp += L(t_r+j)*Rrev_[j+1]*dt;
        }
        P_in = (P_FREE(t_r)*(4.0*pi*pow(R_FREE(t_r), 4.0)) + temp)/(4.0*pi*pow(Rrev_[i], 4.0));
        // P_out = 0.037*E_SN/pow(R_SNR(t_r+i), 3.0);  // 50% the Sedov solution
        double v_fs = (2.0/5.0)*sqrt(2.026*E_SN/rho_ISM)*pow(pow(R_SNR(t_ST), 5.0/2.0)+sqrt(2.026*E_SN/rho_ISM)*(t_r+i-t_ST)*yr, -3.0/5.0);
        P_out = 0.25*(0.75*rho_ISM*pow(v_fs, 2.0));  // 25% the Sedov solution
        Rrev_[i+1] = 4.0*pi*pow(Rrev_[i], 2.0)*pow(dt, 2.0)/M_sw * (P_in-P_out) + 2.0*Rrev_[i] - Rrev_[i-1];
    }
}

double P_REV(double t_, double* Rrev_){
    double P1 = 0.0, P2 = 0.0;
    double dt = yr;
    P1 = P_FREE(t_r)*(4.0*pi*pow(R_FREE(t_r), 4.0));

    int index1 = round(t_r), index2 = round(t_);
    for(int i=index1; i<index2; i++)
        P2 += L(i)*Rrev_[i-index1+1]*dt;

    return( (P1+P2)/(4.0*pi*pow(Rrev_[index2-index1+1], 4.0)) );
}

double R_ST(double t_, double R_t_s, double P_in_t_s){
    double v_fs = (2.0/5.0)*sqrt(2.026*E_SN/rho_ISM)*pow(pow(R_SNR(t_ST), 5.0/2.0)+sqrt(2.026*E_SN/rho_ISM)*(t_-t_ST)*yr, -3.0/5.0);  // speed of the forward shock
    double P_out = 0.75*rho_ISM*pow(v_fs, 2.0);
    return( R_t_s*pow(P_in_t_s/P_out, 1.0/4.0) );
}

double P_PWN(double t_, double* Rpwn_){
    double Ppwn = 0.0;
    double dt = yr;
    int index = round(t_);
    for(int i=0; i<index; i++)
        Ppwn += L(i)*Rpwn_[i]*dt;
    Ppwn = Ppwn/(4.0*pi*pow(Rpwn_[index-1], 4.0));
    return(Ppwn);
}



// Numerical calculation of TS radius by N. Bucciantini, et al. (2004)
double PIfun(double Rts_, double Rpwn_){
    double y = sqrt(81.*sigma/2.0)*Rpwn_/Rts_;
    double tempt1 = pow( 1.0+pow(y, 2.0)+sqrt(pow(1+y*y, 2.0)-1.0), -1.0/3.0 );
    double tempt2 = pow( 1.0+pow(y, 2.0)+sqrt(pow(1+y*y, 2.0)-1.0), 1.0/3.0 );
    double G = 1.0+tempt1+tempt2;
    return( 27.0/pow(G, 4.0)*(2.0+3.0*pow(y, 2.0)/pow(G, 2.0)) );
}

double Hypgeom_2F1(int i, double t_){
    arb_t a1,b1,c1,z1,res1;
    arb_init(a1);arb_init(b1);arb_init(c1);arb_init(z1);arb_init(res1);
    arb_set_d(a1,n-1.+alphao+i);arb_set_d(b1,1.+alphao+i);arb_set_d(c1,2.+alphao+i);arb_set_d(z1,-t_/tau0);
    arb_hypgeom_2f1(res1,a1,b1,c1,z1,0,50);
    // arb_printd(res1,10);
    double res2 = arf_get_d(arb_midref(res1),ARF_RND_NEAR);
    arb_clear(a1);arb_clear(b1);arb_clear(c1);arb_clear(z1);arb_clear(res1);
    return(res2);
}

double Qfun(double t_){
    double tempt1 = 0.0, tempt2 = 0.0;
    for(int i=0; i<5; i++){
        tempt1 += coeff[i]*pow(t_/tau0,i)/(1.+alphao+i)*Hypgeom_2F1(i,t_);
        tempt2 += coeff[i]*pow(t_/tau0,i)/pow(1.+t_/tau0,i);
    }
    return(tempt1*(t_/tau0)*pow(1+t_/tau0,n-2.*(1.-alphao))*pow(tempt2,-2.0));
}

double fun(double t_, double Rts_, double Rpwn_){
    double PI_ycd,PI0,Qt,R0,factor;
    PI_ycd = PIfun(Rts_, Rpwn_);
    PI0 = PIfun(Rts_, Rts_);
    Qt = Qfun(t_);
    R0 = pow(pow(E_SN, 3.0)*pow(Ls0, 2.0)/pow(M_ej, 5.0), 1.0/10.0)*pow(t_*yr, 6.0/5.0);
    factor = 81.0*sigma/2.0 * tau0*yr*c/R0;
    return(81.0*sigma/2.0*pow(Rpwn_/Rts_, 2.0)*PI_ycd - factor*PI0*Qt);
}

double R_TS(double t_, double Rpwn_){  // radius of termination shock
    double Rts_;
    double Rmin = 3e-12, Rmax = Rpwn_;
    if(fun(t_,Rmin,Rpwn_)*fun(t_,Rmax,Rpwn_)>0.0){
        return(0.0);
    }
    do{
        Rts_ = (Rmin+Rmax)/2.0;
        if(fun(t_,Rmin,Rpwn_)*fun(t_,Rts_,Rpwn_)>0.0)
            Rmin=Rts_;
        else
            Rmax=Rts_;
    }while(fabs(fun(t_,Rts_,Rpwn_))>pow(10.0, -6.0));
    
    return(Rts_);
}
