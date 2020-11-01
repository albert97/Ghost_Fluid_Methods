#include <vector>
#include <random>
#include <iostream>
#include <string>
#include <fstream>
#include <regex>
#include <zconf.h>
#include <iostream>
#include <vector>
#include <random>
#include <cmath>

using namespace std;




/*--------------------------------------- set up the parameters --------------------------------------- */
int N = 400; // number of cells
double L = 1.0;   // domain length
double CFL = 0.9;
double tStop =0.0012;; //0.000312629;
double x1 = 0.05;
double x0 = 0.5;
double gam1 = 1.4;
double gam2 = 1.67;
double DensityL = 1.3333;
double DensityM = 1.0;
double DensityR = 0.1379;
double VelocityL =  0.3535 * sqrt(100000);
double VelocityM = 0.0;
double VelocityR= 0.0;
double PressureL= 1.5 * pow(10.0, 5);
double PressureM = 1.0 * pow(10.0, 5);
double PressureR= 1.0 * pow(10.0, 5);
double w= 0;
double Tolerance = 10.0e-6;



/*-----------------------------------------Level-set Equation-------------------------------------------*/

/*--------------------------Display fucntion level-set function ------------------------*/
void Display_levelset ( double const *U_data, int N){

    cout<<"\n"<<endl;

    for (int i=0;i<N;i++){

        cout<<U_data[i]<<"\t";
    }

    cout<<"\n"<<endl;
}
/*----------------------------------------------Sign function-----------------------------------------------*/

double sng(double Psi){

    if (Psi>=0){return 1;}
    if (Psi<0){ return -1;}

}


/*-------------------------------------Initial Condition level-set -----------------------------------------*/

void initial_condition_level_set (double dx, double *phi, int N){

    for (int i =0 ; i< N; i++){

        double x= i *dx;
        //cout<<"x\t"<<x<<endl;
        phi[i] = x - x0;
        //    cout<< "u\t"<<u[i]<<endl;

    }

}



/*--------------------------------------Boundary Condition level-set ---------------------------------------*/

double *boundary_condition_level_set( double *Property, int N ){

    double *Property_extend = new double [N+2];


    for (int i = 0; i<N; i++){

        Property_extend[i+1] = Property[i];

    }

    Property_extend[0] =  Property[0];
    Property_extend[N+1] = Property[N-1];


    return Property_extend;

}


/*--------------------------Level-set Equation----------------------------*/

double *Hamilton_Jacobi(int N, double *phi_extend, double dx, double dt, double *v){

    double *phi_new = new double [N];

    for (int i =0 ; i<N;i++){

        if (v[i] >=0){

            phi_new[i] = phi_extend[i+1] - dt/dx * v[i] * (phi_extend[i+1] - phi_extend[i]);

            //    cout<<u_new[i]<<"\t=\t"<<u_extend[i]<< "\t-\t"<< dt/dx<< "\t*\t"<< v<< "\t*\t"<<"\t ( \t" <<u_extend[i+1] << "\t- \t"<<u_extend[i] <<"\t)\t"<<endl;


        }

        if(v[0]<=0){

            phi_new[i] = phi_extend[i+1] - dt/dx * v[i] * (phi_extend[i+2] - phi_extend[i+1]);
            //    cout<<u_new[i]<<"\t=\t"<<u_extend[i]<< "\t-\t"<< dt/dx<< "\t*\t"<< v<< "\t*\t"<<"\t ( \t" <<u_extend[i+1] << "\t- \t"<<u_extend[i] <<"\t)\t"<<endl;

        }

    }

    return phi_new;

}

/* -----------------------------------------------Fast Sweeping----------------------------------------------*/

double *FastSweeping(int N, double *psi_extend, double dx){


    double *psi_new = new double [N];
    double psi_bar ;

    for (int i =0 ; i<N;i++){

        if (sng(psi_extend[i+1]) > 0 ){

            double psi_x ;
            //   cout<< "\tpsi_extend[i+2]\t"<<psi_extend [i+2]<<"\tpsi_extend[i]\t"<<psi_extend[i];
            if ( psi_extend [i+2]>= psi_extend[i]){
                psi_x =  psi_extend[i];
            }
            else{
                psi_x =  psi_extend[i+2];
            }

            psi_bar = dx + psi_x;

            //       cout << "\tpsi_x\t" <<psi_x<<"\tpsi_bar\t"<<psi_bar;

            if ( psi_bar > psi_extend[i+1]){

                psi_extend[i+1] = psi_bar;
            }
            //      cout << " \t psi_extend[i+1] \t" <<psi_extend[i+1];
        }


        if (sng(psi_extend[i+1]) < 0.0){
            //cout<< "\tpsi_extend[i+2]\t"<<psi_extend [i+2]<<"\tpsi_extend[i]\t"<<psi_extend[i];

            double psi_x ;
            if ( psi_extend [i+2]>= psi_extend[i]){
                psi_x =  psi_extend[i+2];
            }
            else{
                psi_x =  psi_extend[i];
            }

            psi_bar = -dx +psi_x;
            //      cout << "\tpsi_x\t" <<psi_x<<"\tpsi_bar\t"<<psi_bar;

            if ( psi_bar < psi_extend[i+1]){

                psi_extend[i+1] = psi_bar;
            }
            //       cout << " \t psi_extend[i+1] \t" <<psi_extend[i+1];
        }
    }



    for (int i=0; i< N; i++){

        psi_new[i]= psi_extend[i+1];
    }

    return psi_new;

}



/*---------------------------------------------Finding the interface--------------------------------------*/
int interface_finding(double *phi, int N) {

    int interface;


    for (int i = 0; i < N; i++) {
        //cout<<"phi[i]\t"<<phi[i]<<"sng(phi[i])\t"<<sng(phi[i])<<endl;
        if (sng(phi[i]) != sng(phi[i]+2) ) {

            interface = i+1 ;
        }
    }
    return interface;
}



/* ------------- inital condition and boundary condiitons  ----------- */

void initialConditions_ERS(int N, double dx, double *speed, double *density, double *pressure, double *energy){

    double EL = (PressureL/(DensityL * ( gam1 - 1.0)))+(VelocityL * VelocityL)/2.0 ;
    double EM = (PressureM/(DensityM * ( gam1 - 1.0)))+(VelocityM * VelocityM)/2.0 ;
    double ER = (PressureR/(DensityR * ( gam2 - 1.0)))+(VelocityR * VelocityR)/2.0 ;

    for (int i = 0; i < N; i++){

        double x=(i)*dx;

        if (x <= x0){

            if (x <= x1){
                speed[i] = VelocityL;
                density[i] = DensityL;
                pressure[i] = PressureL;
                energy[i] = EL;
            }

            else{

                speed[i] = VelocityM;
                density[i] = DensityM;
                pressure[i] = PressureM;
                energy[i] = EM;

            }
        }


        else{
            speed[i] = VelocityR;
            density[i] = DensityR;
            pressure[i] = PressureR;
            energy[i] = ER;
        }

    }


}

/* ----------This follwowing function allows you to print results on screen -------- */
void Display_ERS ( double *data,double *data2,double *data3, double *data4,int N){


    cout<<"\n"<<endl;

    for (int i=0;i<N;i++){

        cout<<data[i]<<"\t";
    }

    cout<<"\n"<<endl;

    for (int i=0;i<N;i++){

        cout<<data2[i]<<"\t";
    }

    cout<<"\n"<<endl;

    for (int i=0;i<N;i++){

        cout<<data3[i]<<"\t";
    }

    cout<<"\n"<<endl;

    for (int i=0;i<N;i++){

        cout<<data4[i]<<"\t";
    }

    cout<<"\n"<<endl;
}



double *boundaryConditions2( double *Property){

    double *Property_extend = new double [(N+2)];

    for (int i = 0; i<N; i++){

        Property_extend[i+1] = Property[i];
    }

    Property_extend[0] = Property[1];
    Property_extend[N+1] = Property[N-1];


    return Property_extend;
}



/* ----------This follwowing function allows you to print results on screen -------- */
void Display ( double **U_data, double **F_data, double *data,double *data2,double *data3, double *data4,int N){

    for (int i=0;i<3;i++){
        for (int j=0;j<N;j++){
            cout<<U_data[i][j]<<"\t";
        }
        cout<<endl;
    }

    cout<<"\n"<<endl;

    for (int i=0;i<3;i++){
        for (int j=0;j<N;j++){
            cout<<F_data[i][j]<<"\t";
        }
        cout<<endl;
    }

    cout<<"\n"<<endl;

    for (int i=0;i<N;i++){

        cout<<data[i]<<"\t";
    }

    cout<<"\n"<<endl;

    for (int i=0;i<N;i++){

        cout<<data2[i]<<"\t";
    }

    cout<<"\n"<<endl;

    for (int i=0;i<N;i++){

        cout<<data3[i]<<"\t";
    }

    cout<<"\n"<<endl;

    for (int i=0;i<N;i++){

        cout<<data4[i]<<"\t";
    }

    cout<<"\n"<<endl;
}

/*-------------Compute Ghost Fluid domain Boundary for second order Approximate solver ------------------*/

double **GF_Domain_boundaryConditions (int N, double**u) {

    double *u_extend = new double[(N + 4) * 3];
    double **U_extend = new double *[3];

    U_extend[0] = u_extend;

    U_extend[1] = u_extend + (N + 4);

    U_extend[2] = u_extend + (N + 4) * 2;

    for (int i = 0; i < N; i++) {

        U_extend[0][i + 2] = u[0][i];

        U_extend[1][i + 2] = u[1][i];

        U_extend[2][i + 2] = u[2][i];

    }


    U_extend[0][1] = u[0][1];
    U_extend[0][0] = u[0][1];
    U_extend[0][N + 2] = u[0][N - 1];
    U_extend[0][N + 3] = u[0][N - 1];

    U_extend[1][1] = u[1][1];
    U_extend[1][0] = u[1][1];
    U_extend[1][N + 2] = u[1][N - 1];
    U_extend[1][N + 3] = u[1][N - 1];

    U_extend[2][1] = u[2][1];
    U_extend[2][0] = u[2][1];
    U_extend[2][N + 2] = u[2][N - 1];
    U_extend[2][N + 3] = u[2][N - 1];

    return U_extend;
}

/* ------------------------------------------------ solvers ---------------------------------------------- */

/* --------------- Compute Gamma related constants --------------- */
/* --------------- Compute Gamma related constants --------------- */

    double G1(double gam){

        return (gam - 1.0)/(2.0 * gam);

    }
    double G2(double gam){

        return (gam + 1.0)/(2.0 * gam);

    }

    double G3(double gam){

        return 2.0 * gam/(gam - 1.0);

    }

    double G4(double gam){

        return 2.0/(gam - 1.0);

    }

    double G5(double gam){

        return 2.0/(gam + 1.0);

    }

    double G6(double gam){

        return (gam - 1.0)/(gam + 1.0);

    }

    double G7(double gam){

        return (gam - 1.0)/2.0;

    }

    double G8(double gam){

        return gam - 1.0;

    }

/* --------------- Pressure initial guess --------------- */

    double P_M(double DR, double DL, double UR, double UL, double PL, double PR, double CL, double CR, double GammaL, double GammaR){

        double PM;
        double QUSER = 2.0;
        double CUP  = 0.25 * ( DL + DR ) * ( CL + CR );
        double PPV  = 0.5 * ( PL + PR ) + 0.5 * ( UL - UR ) * CUP;
        PPV = max(0.0, PPV);
        double PMIN = min ( PL,  PR );
        double PMAX = max( PL,  PR );
        double QMAX = PMAX/PMIN;

        if (QMAX <= QUSER && PMIN <= PPV && PPV <= PMAX){

            //Select PVRS Riemann Solver
            PM = PPV;

        }
        else{

            if (PPV < PMIN){

                //Select Two_Shock Riemann Solver

                double PQ  = pow((PL/PR), G1(GammaL));
                double UM  = (PQ*UL/CL + UR/CR + G4(GammaL)*(PQ - 1.0))/(PQ/CL + 1.0/CR);
                double PTL = 1.0 + G7(GammaL)*(UL - UM)/CL;
                double PTR = 1.0 + G7(GammaL)*(UM - UR)/CR;
                PM  = 0.5*(pow((PL*PTL), G3(GammaL) ) + pow((PR*PTR), G3(GammaR)));
            }
            else{

                //Select Two_Shock Riemann Solver with PVRS as estimate
                double GEL = sqrt((G5(GammaL)/DL)/(G6(GammaL)*PL + PPV));
                double GER = sqrt((G5(GammaR)/DR)/(G6(GammaR)*PR + PPV));
                PM  = (GEL*PL + GER*PR - (UR - UL))/(GEL + GER);

            }
        }

        return PM;

    }



/* --------------- Exact Solution for pressure and velocity in star region is found --------------- */

// Solving f(P, WL, WR ) = fL ( P, WL ) + fR ( P, WR ) + UR - UL =0

    double ConstA_Func(double gam, double D){

        return 2.0 / ((gam+1.0)*D);

    }

    double ConstB_Func(double gam, double P){

        return (gam - 1.0) * P / (gam+1.0);

    }

    void Pressure_star(double &F, double &FD, double DK, double PK, double CK, double P_Guess , double gamk) {

        //Shock Wave

        if (P_Guess > PK) {

            double QRT = sqrt(ConstA_Func(gamk, DK) / (ConstB_Func(gamk, PK) + P_Guess));

            F = (P_Guess - PK) * QRT;

            //FL = ( P_Guess - PL) * pow(( ConstA_Func( gam, DL ))/(P_Guess + ConstB_Func(gam, PL)), 0.5);

            FD = (1.0 - 0.5 * (P_Guess - PK) / (ConstB_Func(gamk, PK) + P_Guess)) * QRT;
        }
            //Rarefraction Wave

        else {

            double PRAT = P_Guess / PK;

            F = G4(gamk) * CK * (pow(PRAT, G1(gamk)) - 1.0);
            // FL = ((2 * CL)/(gam - 1)) * (pow( P_Guess/PL, (gam-1)/(2 * gam)) - 1);

            FD = (1.0 / (DK * CK)) * pow(PRAT, -G2(gamk));

        }

    }

/*--------------------------Calculate the solution for pressure and velocity in the star Region ------------------------------*/

    double *Newton_Raphson( double DR, double DL,double UR, double UL, double PL, double PR, double CL, double CR, double gammaL, double gammaR, double Tolerance ){

        double *UPM = new double [2];
        double P_star;
        double dU;
        double FL, FR;
        double FDL, FDR;
        double U_star;
        int NI = 0;
        double CHANGE;
        double P_Guess;

        P_Guess = P_M( DR, DL, UR, UL, PL, PR, CL, CR, gammaL, gammaR);

        // cout<<"inital Pressure guesss is \t"<<P_Guess<<endl;

        dU = UR-UL;
        //cout<<"Pressure left\t "<<PL<<"Pressure right\t"<<PR<<endl;

        //Compute Pressure

        do{

            Pressure_star( FL, FDL, DL, PL, CL, P_Guess, gammaL );
            Pressure_star( FR, FDR, DR, PR, CR, P_Guess, gammaR );

            // cout<< "(\t"<<"\t FL \t"<<FL <<"\t+\t"<<"FR\t"<< FR <<"\t+\t"<<"dU\t"<< dU<<"\t)\t"<<"\t / \t"<<"\t(\t"<<"FDL \t"<<FDL<<"\t + \t"<<"\tFDR\t" << FDR<<"\t)\t"<<endl;
            P_star = P_Guess - ( FL + FR + dU)/(FDL + FDR);
            CHANGE = 2.0 * abs((P_star - P_Guess)/(P_star + P_Guess));
            NI ++;
            // cout<<"Numer of iteration is \t"<<NI<<"\t"<<"Change in Pressure is \t"<< CHANGE<< "\t"<<endl;
            if ( P_star < 0.0){

                P_star = Tolerance;
            }
            P_Guess = P_star;

        } while ( CHANGE > Tolerance && NI <20 );

        //Compute Velocity

        U_star = 0.5 * ( UL + UR + FR - FL);

        //  cout<<"the new velocity is \t "<<U_star<<endl;

        UPM[0] = U_star;
        UPM[1] = P_star;

        return UPM;
    }


    double *U_star_L  ( double DL, double PL, double *UPM, double gamL) {

        double *DUP = new double [3];
        double UM= UPM[0];
        double PM= UPM[1];
        double PML;

        //Sampling point lies to the left of the contact discontinuity
        if (PM <= PL){
            //left Rarefraction
            // Sampled point is Star Left state

            DUP[0] = DL*pow((PM/PL), (1.0/gamL));
            DUP[1] = UM;
            DUP[2] = PM;

        }
        else{

            //left shock
            PML = PM / PL;

            // Sampled point is Star Left state
            DUP[0] = DL*(PML + G6(gamL))/(PML * G6(gamL) + 1.0);
            DUP[1] = UM;
            DUP[2] = PM;


        }


        return DUP;


    }


    double *U_star_R(double DR,  double PR, double *UPM, double gamR) {

        double *DUP = new double [3];
        double UM= UPM[0];
        double PM= UPM[1];

        double PMR;

        if (PM > PR){

            //Right shock

            PMR = PM/PR;
            //Sampled point is Star Right state

            DUP[0] = DR * (PMR + G6(gamR))/(PMR * G6(gamR) + 1.0);
            DUP[1] = UM;
            DUP[2] = PM;

        }

        else{

            //Right Rarefraction

            //Sampled point is Star Right state

            DUP[0] = DR * pow((PM/PR), (1.0/gamR));
            DUP[1] = UM;
            DUP[2] = PM;

        }
        return DUP;

    }

/*----------------------------------------Ghost Fluid Boundary-----------------------------------------*/

void **Ghost_Fluid_Boundary(double DUP_StarL, double DUP_StarR, double *P1, double *P2, double const *P, double const *psi, int N, int N1) {

    //initialised both domain
    double *P1_temp = new double [N];
    double *P2_temp = new double [N];

    //cout<<P_extend[N1+1]<<endl;

    for ( int i =0; i< N; i++){

        P1_temp [i] = P[i];
        P2_temp [i] = P[i];

    }

    //Implement the ghost fluid values into ghost fluid doamin
    //  cout<< P1_extend[N1]<<endl;
    P1_temp[N1+1] = DUP_StarL;
    P2_temp[N1] = DUP_StarR;



    for (int i = 0; i<N; i++) {

        //Populate ghost cell for U_star_left
        if (sng(psi[i]) > 0.0) {

            P1_temp[i] = P1_temp[N1+1];

        }


        //Populate ghost cell for U_star_right

        if (sng(psi[i]) < 0.0) {

            P2_temp[i] = P2_temp[N1];

        }
    }


    for (int i = 0; i < N; i++) {

        P1[i] = P1_temp[i];
        P2[i] = P2_temp[i];
    }
}

void Riemann_based_Ghost_Fluid_Boundary(int N1, int N, double const *pressure, double const *density, double const *uspeed, double *density1, double *density2, double *uspeed1, double *uspeed2, double *pressure1, double *pressure2, double *phi, double **U1, double **U2 ) {

    double *UPM_temp ;
    double *DUP_R ;
    double *DUP_L ;


    double CL = sqrt(gam1 * (pressure[N1] / density[N1]));
    double CR = sqrt(gam2 * (pressure[N1+ 1] / density[N1+ 1]));

    UPM_temp = Newton_Raphson(density[N1 + 1], density[N1], uspeed[N1 + 1], uspeed[N1], pressure[N1], pressure[N1 + 1], CL, CR, gam1, gam2, Tolerance);

    DUP_L = U_star_L(density[N1], pressure[N1], UPM_temp, gam1);
    DUP_R = U_star_R(density[N1 + 1], pressure[N1 + 1], UPM_temp, gam2);


    //checker loop

    cout<< " intermediate velcoity   \t"<<  UPM_temp[0]<<endl;
    cout<< " intermediate pressure\t "<< UPM_temp[1]<<endl;
    cout<< " left star states density  \t"<<  DUP_L[0]<<endl;
    cout<< "right star state density \t "<< DUP_R[0]<<endl;


    //Populate the Intermediate state into all ghost fluid regions by Fast sweeping method

    Ghost_Fluid_Boundary (DUP_L[0], DUP_R[0], density1, density2, density, phi, N, N1);

    Ghost_Fluid_Boundary (DUP_L[1], DUP_R[1], uspeed1, uspeed2, uspeed, phi, N, N1);

    Ghost_Fluid_Boundary (DUP_L[2], DUP_R[2], pressure1, pressure2, pressure, phi, N, N1);

    for (int i =0 ; i<N; i++){

        U1[0][i] = density1[i];
        U1[1][i] = density1[i] * uspeed1[i];
        U1[2][i] = (pressure1[i] / ( gam1 - 1.0))+(uspeed1[i] * uspeed1[i]* density1[i])/2.0;
    }

    for (int i =0 ; i<N; i++){

        U2[0][i] = density2[i];
        U2[1][i] = density2[i] * uspeed2[i];
        U2[2][i] = (pressure2[i] / ( gam2 - 1.0))+(uspeed2[i] * uspeed2[i]* density2[i])/2.0;
    }

}



double *boundaryConditionsHLLC( double *Property){

    double *Property_extend = new double [(N+4)];

    for (int i = 0; i<N; i++){

        Property_extend[i+2] = Property[i];
    }

    Property_extend[1] = Property[1];
    Property_extend[0] = Property[1];
    Property_extend[N+2] = Property[N-1];
    Property_extend[N+3] = Property[N-1];

    return Property_extend;
}



double sound_speed (double pressure, double rho, double gam) {


    return sqrt(gam * (pressure / rho));


}


double dt_func( double dx,  double *speed_extend, double *pressure_extend, double *density_extend, double gam){

    double a;
    double dt;

    double cmax = 0, ctemp;

    //Finding max wave speed

    /* ??????? i am confused if it should be N+3 rather than N+2 ????????*/
    for (int i=1; i<N+3;i++){

        a = sound_speed(pressure_extend[i], density_extend[i], gam);

        ctemp = speed_extend[i] + a ;

        if (cmax<ctemp){

            cmax =ctemp;

        }

    }

    //Determine the value for dt

    dt = CFL * (dx/cmax);

    return dt;

}





/*----------------------------------------------------------------------------------------------------------------*/
/*                                                                                                                */
/*                                 Approximate Solver HLLC with MUSCL-HANCOCK                                     */
/*                                                                                                                */
/*                                                                                                                */
/*                                                                                                                */
/*----------------------------------------------------------------------------------------------------------------*/

/* --------------------------------- solvers ------------------------------------ */


/*--------------------Conservtive to primitive ---------------------*/


double *CtoP(double *US, double gam){

    double *w = new double[3];

    w[0] = US[0];
    w[1] = US[1] / US[0];
    w[2] = (US[2] - 0.5 * US[0] * pow(( US[1] / US[0]),2)) * (gam -1);

    return w;


}

//work out the slope_vector

double *U_to_F(double *U_value, double gam) {

    double * flux= new double[3];

    flux[0] = U_value[1];
    flux[1] = U_value[1]*(U_value[1]/U_value[0]) + (gam-1)*(U_value[2] - 0.5* U_value[0] * pow((U_value[1]/ U_value[0]),2.0));
    flux[2] = (U_value[1]/U_value[0]) * (U_value[2] + (gam-1)*( U_value[2] - 0.5 * U_value[0] * pow(( U_value[1]/U_value[0]),2.0)));

    return flux;

}
//Can do switch cases with Slop limitor

double min_bee( double r, double SlopeLR){

    double min_bee =0.0;

    if( r <= 0){ min_bee = 0.0; }
    if( r>=0 && r<=1 ){ min_bee = r; }
    if (r >=1 && SlopeLR >=1){ min_bee = 1.0; }
    if (r >=1 && SlopeLR <= 1){ min_bee = SlopeLR; }

    return min_bee;

}
double super_bee( double r, double SlopeLR){

    double superbee =0.0;

    if( r <= 0){ superbee = 0.0; }
    if( r>=0 && r<= 0.5 ){ superbee = 2.0*r; }
    if( r>=0.5 && r<=1 ){ superbee = 1.0; }
    if (r >=1 && r<=2 &&r<= SlopeLR ){ superbee = r; }
    if (r >=1 && SlopeLR <= 2 && SlopeLR <= r ){ superbee = SlopeLR; }
    if (r >=1 && SlopeLR >=2  &&  r>= 2 ){ superbee =2.0; }

    return superbee;

}





double **Muscl_Hancock(double w, double dx, double dt, double **U_list, double gam){

    // U_list is a three by three matrix that has the conservetive values for three ajucent points
    // calculate the slope_vector

    double di, di1, di2, Sdi;


    //store extrapolated values

    double **U_update = new double *[3];
    for (int i =0; i<3; i++){
        U_update[i] = new double[2];
    }

    // Store Ubar values

    double **U_MH = new double *[3];
    for (int i =0; i<3; i++){
        U_MH[i] = new double[2];
    }

    double r;
    double SlopeLR, SlopeLimiter;
    double beta =1;
    for (int i =0 ; i<3 ; i++){

        di1 = U_list[i][1] - U_list[i][0];

        di2 = U_list[i][2] - U_list[i][1];

        r = di1 / di2;

        SlopeLR = ( 2* beta )/ ((1-w) + (1+w)*r);

        SlopeLimiter = super_bee(r, SlopeLR);

        di = 0.5 * ( 1 + w ) * di1 + 0.5 * ( 1 - w ) * di2;

        Sdi = di * SlopeLimiter;
        //Boundary extrapolated left values

        U_update[i][0] = U_list[i][1] - 0.5 * Sdi;

        //     cout<<U_update[i][0]<<" \t = \t"<<U_list[i][1]<<"\t - \t" << 0.5 <<"\t*\t"<<di<<"\n";

        //Boundary extrapolated right values

        U_update[i][1] = U_list[i][1] + 0.5 * Sdi;
        //     cout<< U_update[i][1] <<" \t = \t"<<U_list[i][1]<< "\t+\t"<< 0.5<<"\t*\t"<< di<<"\n";

    }


    double *U_L= new double[3];

    U_L[0] = U_update[0][0] ;
    U_L[1] = U_update[1][0] ;
    U_L[2] = U_update[2][0] ;

    double *U_R= new double[3];

    U_R[0] = U_update[0][1] ;
    U_R[1] = U_update[1][1] ;
    U_R[2] = U_update[2][1] ;

    double *flux_L;

    double *flux_R;

    flux_L = U_to_F( U_L, gam );
    flux_R = U_to_F( U_R, gam );


    //NOTE assumptions Ur state vector to flux conversion is accurate

    for (int i = 0 ; i< 3; i++){

        //Ubar left
        U_MH[i][0] =  U_L[i] + 0.5 * ( dt/dx ) *( flux_L[i]- flux_R[i]);

        //  cout<<U_MH[i][0]<<"="<< U_L[i]<< "\t+\t" <<0.5 <<"\t*\t" <<"\t(\t"<< dt/dx<<" )" <<"*"<<"(" <<flux_L[i]<<"- "<<flux_R[i]<<")"<<"\n";
        //Ubar right
        U_MH[i][1] =  U_R[i] + 0.5 * ( dt/dx ) *( flux_L[i]- flux_R[i]);


    }



    return U_MH;
}

/* ---------------------------pressure based speed estimate --------------------------*/

double *pressure_based_wave_speed(double p_L, double p_R, double u_L, double u_R, double rho_L,double rho_R, double gam){

    double *wavespeed = new double[2];

    double a_L, a_R, a_bar;
    double rho_bar;
    a_L = sound_speed(p_L, rho_L, gam);
    a_R = sound_speed(p_R, rho_R, gam);
    a_bar = 0.5 * (a_L+a_R);
    rho_bar = 0.5 * (rho_L + rho_R);

    double Ppvrs, p_star;
    Ppvrs = 0.5 * (p_L + p_R) - 0.5 * (u_R - u_L) * rho_bar * a_bar;

    if (Ppvrs > 0){ p_star = Ppvrs;}
    else {p_star = 0; }

    double q_L, q_R;
    if (p_star <= p_L ){ q_L = 1; }
    else { q_L = sqrt(1.0 + ((gam +1 )/( 2 * gam)) * ( (p_star / p_L) -1)); }


    if (p_star <= p_R ){ q_R = 1; }
    else { q_R = sqrt(1.0 + ((gam +1 )/( 2 * gam)) * ( (p_star / p_R) -1)); }


    //S_L

    wavespeed[0] = u_L - a_L * q_L;

    //S_R

    wavespeed[1] = u_R + a_R * q_R;


    return wavespeed;


}


double *U_stark(double s_k, double s_star, double u_k,  double p_k, double rho_k, double e_k){

    double *U_star = new double[3];

    double ft= rho_k * ( (s_k - u_k)/(s_k - s_star));

    U_star[0]= ft *1.0;

    U_star[1]= ft * s_star;

    U_star[2]= ft * ((e_k /rho_k) + (s_star - u_k)*(s_star + (p_k/(rho_k * (s_k - u_k)))));

    return U_star;
}

double **U_update(double dx, double dt, double **U_extend, double gam){

    double S_L, S_R;
    double a_L, a_R; // If you want to use Roe average Wavespeed then a_L and a_R is there for you
    double S_star;

    double *wavespeed;

    double *U_hllc = new double [(N)*3];
    double **U_new = new double *[3];


    U_new[0] = U_hllc;
    U_new[1] = U_hllc + N;
    U_new[2] = U_hllc + N*2;

    double *f_force = new double [(N+1)*3];
    double **f_new = new double *[3];

    f_new[0] = f_force;
    f_new[1] = f_force + N+1;
    f_new[2] = f_force + (N+1)*2;

    double *Ustar_L ;
    double *Ustar_R ;


    double **U_temp = new double *[3];
    for (int i = 0; i<3 ; i++){
        U_temp[i] = new double [3];
    }


    double **U_temp2 = new double *[3];
    for (int i = 0 ; i<3; i++){

        U_temp2 [i] = new double [3];

    }

    double **U_MH1, **U_MH2;
    double *U_MH_L = new double [3];
    double *U_MH_R = new double [3];


    for (int i=1;i<N+2;i++){

        for (int j = 0 ; j<3; j++){

            U_temp[j][0] = U_extend[j][i-1];

            U_temp[j][1] = U_extend[j][i];

            U_temp[j][2] = U_extend[j][i+1];

        }

        U_MH1 = Muscl_Hancock( w, dx, dt, U_temp, gam);


        for (int j = 0 ; j<3; j++){
            U_temp2 [j][0] = U_extend [j][i];
            U_temp2 [j][1] = U_extend [j][i+1];
            U_temp2 [j][2] = U_extend [j][i+2];

        }

        U_MH2 = Muscl_Hancock(w, dx,dt, U_temp2, gam);


        for (int k =0; k< 3; k++){

            U_MH_L[k] = U_MH1[k][1];
            U_MH_R[k] = U_MH2[k][0];

        }

        double *WL, *WR ;
        WL = CtoP(U_MH_L, gam);
        WR = CtoP(U_MH_R, gam);


        double pressureL, pressureR, densityL, densityR,  speedL, speedR, energyL, energyR;

        densityL = WL[0];
        densityR = WR[0];
        speedL = WL[1];
        speedR = WR[1];
        pressureL = WL[2];
        pressureR = WR[2];
        energyL = (WL[2]) / (gam - 1) + 0.5 * WL[1] * WL[1] * WL[0];
        energyR = (WR[2]) / (gam - 1) + 0.5 * WR[1] * WR[1] * WR[0];

        wavespeed = pressure_based_wave_speed(pressureL, pressureR,  speedL, speedR, densityL, densityR, gam);

        S_L = wavespeed[0];
        S_R = wavespeed[1];

        double num =pressureR-pressureL + densityL * speedL * (S_L - speedL) - densityR * speedR * (S_R - speedR);

        double dem = densityL * (S_L - speedL) - densityR * (S_R - speedR);


        S_star= num/dem;


        if (S_L>=0){

            f_new[0][i-1] = U_MH_L[1];
            f_new[1][i-1] = U_MH_L[1] * (U_MH_L[1] / U_MH_L [0]) + (gam-1)*(U_MH_L[2] - 0.5 * U_MH_L[0] * pow((U_MH_L[1]/ U_MH_L[0]),2.0));
            f_new[2][i-1] = (U_MH_L[1]/U_MH_L[0]) * (U_MH_L[2] + (gam-1)*( U_MH_L[2] - 0.5 * U_MH_L[0] * pow(( U_MH_L[1] / U_MH_L[0]),2.0)));


        }

        if (S_L<=0 && S_star>=0){
            //    cout<< "second Condition"<<"\t";

            Ustar_L = U_stark(S_L, S_star, speedL, pressureL, densityL, energyL);

            //    cout<< "Ustar"<<Ustar_L[2]<<endl;

            f_new[0][i-1] = U_MH_L[1] + S_L * ( Ustar_L[0] - U_MH_L[0]);
            f_new[1][i-1] = (U_MH_L[1] * (U_MH_L[1] / U_MH_L [0]) + (gam-1)*(U_MH_L[2] - 0.5 * U_MH_L[0] * pow((U_MH_L[1]/ U_MH_L[0]),2.0))) + S_L * ( Ustar_L[1] - U_MH_L[1]);
            f_new[2][i-1] = ((U_MH_L[1]/U_MH_L[0]) * (U_MH_L[2] + (gam-1)*( U_MH_L[2] - 0.5 * U_MH_L[0] * pow(( U_MH_L[1] / U_MH_L[0]),2.0)))) + S_L * ( Ustar_L[2] -  U_MH_L[2]);
            //   cout<< "Ustar"<<f_new[1][i]<<endl;

        }

        if (S_R >= 0 && S_star <= 0){

            //        cout<< "third Condition"<<"\t";
            Ustar_R = U_stark(S_R, S_star, speedR, pressureR, densityR, energyR);

            f_new[0][i-1] = U_MH_R[1] + S_R * (Ustar_R[0] - U_MH_R[0]);
            f_new[1][i-1] = (U_MH_R[1]*( U_MH_R[1] / U_MH_R[0]) + (gam-1)*(U_MH_R[2] - 0.5* U_MH_R[0] * pow((U_MH_R[1] / U_MH_R[0]),2.0)))+ S_R * (Ustar_R[1] - U_MH_R[1]);
            f_new[2][i-1] = ((U_MH_R[1] / U_MH_R[0]) * (U_MH_R[2] + (gam-1)*( U_MH_R[2]- 0.5 * U_MH_R[0] * pow(( U_MH_R[1] / U_MH_R[0]),2.0)))) + S_R * (Ustar_R[2] - U_MH_R[2]);

        }

        if (S_R<=0){
            //         cout<< "fourth Condition"<<"\t";
            f_new[0][i-1] = U_MH_R[1];
            f_new[1][i-1] = U_MH_R[1]*( U_MH_R[1] / U_MH_R[0]) + (gam-1)*(U_MH_R[2] - 0.5* U_MH_R[0] * pow((U_MH_R[1] / U_MH_R[0]),2.0));
            f_new[2][i-1] = (U_MH_R[1] / U_MH_R[0]) * (U_MH_R[2] + (gam-1)*( U_MH_R[2]- 0.5 * U_MH_R[0] * pow(( U_MH_R[1] / U_MH_R[0]),2.0)));

        }
    }

    //  cout<<U_new[0][0] <<"="<<U_extend[0][2] <<"-"<< dt/dx<<"x"<<f_new[0][1]<<"-"<<f_new[0][0]<<"\n";
    for(int i = 0; i<N; i++){

        for (int j=0; j<3; j++){

            U_new[j][i] = U_extend[j][i+2] - (dt/dx) * ( f_new[j][i+1]-f_new[j][i]);


        }
    }

    return U_new ;

}



//Claim plotfile from fstream and get ready for output
fstream plotfile;

int main(void){


    // state vector

    double *u = new double[N*3];
    double **U = new double *[3];

    //DENSITY rho

    U[0]= u;

    //RHO*u

    U[1] = u + N;

    //TOTAL SPECIFIC ENERGY

    U[2] = u  + N * 2;

    double *f= new double[N*3];
    double **flux = new double *[3];

    //rho*u

    flux[0] = f;

    //rho*u**2 + p

    flux[1] = f + N;

    //u*(E+P)

    flux[2] = f + N * 2;


    double **U1 = new double*[3] ;
    double **U2= new double*[3];
    for (int i=0; i< 3; i++){

        U1[i] = new double [N];
        U2[i] = new double [N];

    }

    double dx = L/N;

    double *uspeed = new double[N];
    double *density = new double[N];
    double *pressure = new double[N];
    double *energy = new double[N];
    double *internal_energy = new double[N];


    // Level-set function

    double *phi=new double[N];
    double *phi_extend;

    initial_condition_level_set(dx, phi, N );

    Display_levelset(phi,N);

    //Initial conditions for Exact Solver

    initialConditions_ERS(N, dx, uspeed, density , pressure, energy);


    Display_ERS( uspeed, density, pressure, energy,N);


    //Split the primitive variables in parts depending on the interfae

    double *density1= new double[N];
    // double *density1_extend;
    double *density2 = new double[N];
    //double *density2_extend;
    double *uspeed1 = new double[N];
    //double *uspeed1_extend;
    double *uspeed2 = new double[N];
    //double *uspeed2_extend;
    double *pressure1 = new double[N];
    //double *pressure1_extend;
    double *pressure2 = new double[N];
    //double *pressure2_extend;
    double *energy1 = new double[N];
    double *energy2 = new double[N];


    double *uspeed_extend ;
    double *density_extend ;
    double *pressure_extend ;
    double *energy_extend ;



    double **flux1 = new double *[3];
    double **flux2 = new double *[3];

    for (int i=0; i< 3; i++){

        flux1[i] = new double [N];
        flux2[i] = new double [N];

    }


    double **U1_extend ;
    double **U2_extend;

    //  double **flux_extend ;

    double *uspeed1_extend2 ;
    double *density1_extend2 ;
    double *pressure1_extend2 ;
    double *uspeed2_extend2 ;
    double *density2_extend2 ;
    double *pressure2_extend2 ;


    double t=0.00;
    double dt , dt1, dt2 ;

    int N1;
    N1 = interface_finding(phi, N);
    cout<< N1<<endl;
    do{

        //Implement the Ghost fluid boundaryes
        Riemann_based_Ghost_Fluid_Boundary( N1, N,  pressure, density, uspeed,  density1,  density2,  uspeed1,  uspeed2,  pressure1,  pressure2, phi, U1,  U2);

        cout<<"Display the first ghost fluid domains after implementing the first ghost fluid boundary conditions "<<endl;
        Display(U1, flux1, uspeed1, density1, pressure1, energy1,N);

        cout<<"Display the second ghost fluid domains after implementing the first ghost fluid boundary conditions "<<endl;
        Display(U2, flux2, uspeed2, density2, pressure2, energy2,N);



        //Implement the domain boundarys
        U1_extend = GF_Domain_boundaryConditions(N, U1);

        U2_extend = GF_Domain_boundaryConditions(N, U2);

        uspeed1_extend2 = boundaryConditionsHLLC(uspeed1);

        density1_extend2 = boundaryConditionsHLLC(density1);

        pressure1_extend2 = boundaryConditionsHLLC(pressure1);

        uspeed2_extend2 = boundaryConditionsHLLC(uspeed2);

        density2_extend2 = boundaryConditionsHLLC(density2);

        pressure2_extend2 = boundaryConditionsHLLC(pressure2);


        //Fimd the smallest time step
        dt1 =dt_func(dx, uspeed1_extend2, pressure1_extend2, density1_extend2, gam1);
        dt2 =dt_func(dx, uspeed2_extend2, pressure2_extend2, density2_extend2, gam2);

        if (dt1 >dt2){ dt = dt2;}
        else{ dt = dt1;}


        //Update Level set function using level set equation  and reinilisasition

        phi_extend = boundary_condition_level_set(phi, N);

        phi = Hamilton_Jacobi(N, phi_extend, dx,dt, uspeed);

        Display_levelset(phi,N);

        phi_extend = boundary_condition_level_set(phi, N);

        phi = FastSweeping(N, phi_extend, dx);

        Display_levelset(phi,N);

        //Find the location of the interface for the new time step
        N1 = interface_finding(phi, N);
        cout<<"interface location \t"<<N1<<endl;


        //Call the approximate solver

        U1= U_update(dx, dt, U1_extend, gam1);
        U2= U_update(dx, dt, U2_extend, gam2);

        cout<<"Preint the state vector of two materials \t"<<endl;

        Display(U1,U2, uspeed, density, pressure, energy,N);


        for (int i =0; i<N; i++){

            density1[i] = U1[0][i] ;

            uspeed1[i] = U1[1][i] / U1[0][i];

            pressure1[i] = (U1[2][i] - 0.5 * U1[0][i] * pow(( U1[1][i] / U1[0][i]),2)) * (gam1 -1.0);

            energy1[i] = U1[2][i] ;

            flux1[0][i] =  U1[0][i] * U1[1][i] / U1[0][i];

            flux1[1][i] =  U1[0][i] * (U1[1][i] / U1[0][i]) * (U1[1][i] / U1[0][i]) + ((U1[2][i] - 0.5 * U1[0][i] * pow(( U1[1][i] / U1[0][i]),2)) * (gam1 -1.0));

            flux1[2][i] = ( U1[1][i] / U1[0][i]) * (U1[2][i] + ((U1[2][i] - 0.5 * U1[0][i] * pow(( U1[1][i] / U1[0][i]),2)) * (gam1 -1.0)));

            density2 [i] = U2[0][i] ;

            uspeed2 [i] = U2[1][i] / U2[0][i];

            pressure2[i] = (U2[2][i] - 0.5 * U2[0][i] * pow((U2[1][i] / U2[0][i]),2)) * (gam2 -1.0);

            energy2[i] = U2[2][i] ;

            flux2[0][i] =  U2[0][i] * U2[1][i] / U2[0][i];

            flux2[1][i] =  U2[0][i] * (U2[1][i] / U2[0][i]) * (U2[1][i] / U2[0][i]) + ((U2[2][i] - 0.5 * U2[0][i] * pow(( U2[1][i] / U2[0][i]),2)) * (gam2 -1.0));

            flux2[2][i] = ( U2[1][i] / U2[0][i]) * (U2[2][i] + ((U2[2][i] - 0.5 * U2[0][i] * pow(( U2[1][i] / U2[0][i]),2)) * (gam2 -1.0)));

        }
        for (int i =0; i<N1+1; i++){

            pressure [i] = pressure1[i];
            uspeed [i] = uspeed1[i];
            density [i] = density1[i];

        }

        for (int i =N1+1; i<N; i++){

            pressure [i] = pressure2[i];
            uspeed [i] = uspeed2[i];
            density [i] = density2[i];

        }

        for (int i = 0; i< N1+1;i++){

            for (int j=0; j<3;j++){

                U[j][i] = U1[j][i];
                flux[j][i] = flux1[j][i];

            }
        }

        for (int i = N1+1; i< N;i++){

            for (int j=0; j<3;j++){

                U[j][i] = U2[j][i];
                flux[j][i] = flux2[j][i];

            }
        }



        cout<<"2. Display the first ghost fluid domains after implementing the first ghost fluid boundary conditions "<<endl;
        Display(U1, flux1, uspeed1, density1, pressure1, energy1,N);

        cout<<"2. Display the second ghost fluid domains after implementing the first ghost fluid boundary conditions "<<endl;
        Display(U2, flux2, uspeed2, density2, pressure2, energy2,N);

        cout<<"Display the Two material cimbined vectors  "<<endl;
        Display(U, flux, uspeed, density, pressure, energy,N);


        t += dt;

    }while(t<tStop);
    for (int i =0 ; i<N1+1;i++){

        internal_energy [i] = pressure[i]/(density[i]*(gam1 -1));
    }
    for (int i =N1+1 ; i<N;i++){

        internal_energy [i] = pressure[i]/(density[i]*(gam2 -1));
    }
    plotfile.open("plotSOD_MM_400.dat", fstream::out);

    if(plotfile.fail()){

        cout<<"Error opening file"<<endl;
        exit(0);
    }

    for (int i=0;i<N;i++){

        double x=(i)*dx;

        plotfile<<x<<"\t"<<uspeed[i]<<"\t"<<density[i]<<"\t"<<pressure[i]<<"\t"<<internal_energy[i]<<"\t"<<phi[i]<<"\t"<<endl;
    }


    return 0;

}