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
/*
int N = 400; // number of cells
double L = 1.0;   // domain length
double CFL = 0.9;
double tStop =0.2;
double x0 = 0.5;
double gam=1.4;           // advenction velocity
double DL= 1.0;
double UL= 0.0;
double PL= 1.0;
double DR= 0.125;
double UR= 0.0;
double PR= 0.1;
double w= 0;
double Tolerance = 10.0e-6;

*/
//test 2
/*
int N = 400; // number of cells
double L = 1.0;   // domain length
double CFL = 0.9;
double tStop = 0.15;
double x0 = 0.5;
double gam=1.4;           // advenction velocity
double DL= 1.0;
double UL= -2.0;
double PL= 0.4;
double DR= 1.0;
double UR= 2.0;
double PR= 0.4;
double w= 0.0;
double Tolerance = 10.0e-6;
*/
//test 3

int N = 400; // number of cells
double L = 1.0;   // domain length
double CFL = 0.9;
double tStop = 0.012;
double x0 = 0.5;
double gam=1.4;           // advenction velocity
double DL= 1.0;
double UL= 0;
double PL= 1000.0;
double DR= 1.0;
double UR= 0.0;
double PR= 0.01;
double w= 0.0;
double Tolerance = 10.0e-6;


// Test 4
/*
int N = 100; // number of cells
double L = 1.0;   // domain length
double CFL = 0.9;
double tStop = 0.035;
double x0 = 0.4;
double gam=1.4;           // advenction velocity
double DL= 5.99924;
double UL= 19.5975;
double PL= 460.894;
double DR= 5.99242;
double UR= -6.19633;
double PR= 46.0950;
double w= 0.0;
double Tolerance = 10.0e-6;
*/

/*----------------------------------Function for sound speed calculation---------------------------------*/

double sound_speed (double pressure, double rho) {

    return sqrt(gam * (pressure / rho));

}

/*----------------------------------Function for time-step calculation---------------------------------*/


double dt_func( double dx,  double *speed_extend, double *pressure_extend, double *density_extend){

    double a;
    double dt;

    double cmax = 0, ctemp;

    //Finding max wave speed

    /* ??????? i am confused if it should be N+3 rather than N+2 ????????*/
    for (int i=1; i<N+3;i++){

        a = sound_speed(pressure_extend[i], density_extend[i]);

        ctemp = speed_extend[i] + a ;

        if (cmax<ctemp){

            cmax =ctemp;

        }

    }

    //Determine the value for dt

    dt = CFL * (dx/cmax);

    return dt;

}


/*-----------------------------------------Level-set Equation-------------------------------------------*/

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

/*----------------------------------------------Sign function-----------------------------------------------*/

double sng(double Psi){

    if (Psi>0){return 1;}
    if (Psi==0){ return 0;}
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


/*------------------------Reimplementation of level_set Equation for Riemann exact solver for simnple advection ------------------*/

double *reimplement_levelset(double F_time, double advect_v, double dx, int N){

    double x0_new = x0 + advect_v * F_time;

    double *phi_new = new double [N];

    for (int i =0 ; i< N; i++){

        double x= i *dx;
        //cout<<"x\t"<<x<<endl;
        phi_new[i] = x - x0_new;
        //    cout<< "u\t"<<u[i]<<endl;

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
            cout<< "\tpsi_extend[i+2]\t"<<psi_extend [i+2]<<"\tpsi_extend[i]\t"<<psi_extend[i];
            if ( psi_extend [i+2]>= psi_extend[i]){
                psi_x =  psi_extend[i];
            }
            else{
                psi_x =  psi_extend[i+2];
            }

            psi_bar = dx + psi_x;

            cout << "\tpsi_x\t" <<psi_x<<"\tpsi_bar\t"<<psi_bar;

            if ( psi_bar > psi_extend[i+1]){

                psi_extend[i+1] = psi_bar;
            }
            cout << " \t psi_extend[i+1] \t" <<psi_extend[i+1];
        }


        if (sng(psi_extend[i+1]) < 0.0){

            cout<< "\tpsi_extend[i+2]\t"<<psi_extend [i+2]<<"\tpsi_extend[i]\t"<<psi_extend[i];

            double psi_x ;
            if ( psi_extend [i+2]>= psi_extend[i]){
                psi_x =  psi_extend[i+2];
            }
            else{
                psi_x =  psi_extend[i];
            }

            psi_bar = -dx +psi_x;
            cout << "\tpsi_x\t" <<psi_x<<"\tpsi_bar\t"<<psi_bar;

            if ( psi_bar < psi_extend[i+1]){

                psi_extend[i+1] = psi_bar;
            }
            cout << " \t psi_extend[i+1] \t" <<psi_extend[i+1];
        }
    }



    for (int i=0; i< N; i++){

        psi_new[i]= psi_extend[i+1];
    }

    return psi_new;

}




/*---------------------------------------------Finding the interface--------------------------------------*/

int interface_finding(double *phi, int N){

    int interface;

    for (int i = 0; i<(N-1); i++){

        if ( sng(phi[i])!=sng(phi[i+1])){

            interface = i;

        }
    }

    return interface ;
}




/*----------------------------------------Ghost Fluid Boundary-----------------------------------------*/

void **Ghost_Fluid_Boundary(double DUP_StarL, double DUP_StarR, double *P1, double *P2, double const *P_extend, double const *psi_extend, int N, int N1) {

    //initialised both domain
    double *P1_extend = new double [N+2];
    double *P2_extend = new double [N+2];

    //cout<<P_extend[N1+1]<<endl;

    for ( int i =0; i< N+2; i++){

        P1_extend [i] = P_extend[i];
        P2_extend [i] = P_extend[i];

    }

    //Implement the ghost fluid values into ghost fluid doamin
  //  cout<< P1_extend[N1]<<endl;
    P1_extend[N1+1] = DUP_StarL;
    P2_extend[N1] = DUP_StarR;


    //Populate ghost cell for U_star_left
    for (int i = N1+1; i<N;i++) {

        if (sng(psi_extend[i + 1]) > 0.0) {

            P1_extend[i+1] = P1_extend[i];

        }


        if (sng(psi_extend[i + 1]) < 0.0) {

            P1_extend[i+1] = P1_extend[i];

        }
    }

    //Populate ghost cell for U_star_right
    for (int i = N1; i >= 0; i--) {

        if (sng(psi_extend[i - 1]) > 0) {

            P2_extend[i-1] = P2_extend[i];

        }


        if (sng(psi_extend[i - 1]) < 0.0) {

            P2_extend[i-1] = P2_extend[i];

        }

    }
    
    for (int i = 0; i < N+2; i++) {

        P1[i] = P1_extend[i];
        P2[i] = P2_extend[i];
    }
}


/* ------------------------------------------------ solvers ---------------------------------------------- */

/* --------------- Compute Gamma related constants --------------- */

double G1 = (gam - 1.0)/(2.0 * gam);
double G2 = (gam + 1.0)/(2.0 * gam);
double G3 = 2.0 * gam/(gam - 1.0);
double G4 = 2.0/(gam - 1.0);
double G5 = 2.0/(gam + 1.0);
double G6 = (gam - 1.0)/(gam + 1.0);
double G7 = (gam - 1.0)/2.0;
double G8 = gam - 1.0;



/* --------------- Pressure initial guess --------------- */

double P_M(double DR, double DL, double UR, double UL, double PL, double PR, double CL, double CR){

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

            double PQ  = pow((PL/PR), G1);
            double UM  = (PQ*UL/CL + UR/CR + G4*(PQ - 1.0))/(PQ/CL + 1.0/CR);
            double PTL = 1.0 + G7*(UL - UM)/CL;
            double PTR = 1.0 + G7*(UM - UR)/CR;
            PM  = 0.5*(pow((PL*PTL), G3 ) + pow((PR*PTR), G3));
        }
        else{

            //Select Two_Shock Riemann Solver with PVRS as estimate
            double GEL = sqrt((G5/DL)/(G6*PL + PPV));
            double GER = sqrt((G5/DR)/(G6*PR + PPV));
            PM  = (GEL*PL + GER*PR - (UR - UL))/(GEL + GER);

        }
    }

    return PM;

}



/* --------------- Exact Solution for pressure and velocity in star region is found --------------- */

// Solving f(P, WL, WR ) = fL ( P, WL ) + fR ( P, WR ) + UR - UL =0

double ConstA_Func(double gam, double D){

    return 2.0 / ((gam+1)*D);

}

double ConstB_Func(double gam, double P){

    return (gam - 1) * P / (gam+1);

}

void Pressure_star(double &F, double &FD, double DK, double PK, double CK, double P_Guess ) {

    //Shock Wave

    if (P_Guess > PK) {

        double QRT = sqrt(ConstA_Func(gam, DK) / (ConstB_Func(gam, PK) + P_Guess));

        F = (P_Guess - PK) * QRT;

        //FL = ( P_Guess - PL) * pow(( ConstA_Func( gam, DL ))/(P_Guess + ConstB_Func(gam, PL)), 0.5);

        FD = (1.0 - 0.5 * (P_Guess - PK) / (ConstB_Func(gam, PK) + P_Guess)) * QRT;
    }
        //Rarefraction Wave

    else {

        double PRAT = P_Guess / PK;

        F = G4 * CK * (pow(PRAT, G1) - 1.0);
        // FL = ((2 * CL)/(gam - 1)) * (pow( P_Guess/PL, (gam-1)/(2 * gam)) - 1);

        FD = (1.0 / (DK * CK)) * pow(PRAT, -G2);

    }

}


/*--------------------------Calculate the solution for pressure and velocity in the star Region ------------------------------*/

double *Newton_Raphson( double DR, double DL,double UR, double UL, double PL, double PR, double CL, double CR, double Tolerance ){

    double *UPM = new double [2];
    double P_star;
    double dU;
    double FL, FR;
    double FDL, FDR;
    double U_star;
    int NI = 0;
    double CHANGE;
    double P_Guess;

    P_Guess = P_M( DR, DL, UR, UL, PL, PR, CL, CR);

    // cout<<"inital Pressure guesss is \t"<<P_Guess<<endl;

    dU = UR-UL;
    //cout<<"Pressure left\t "<<PL<<"Pressure right\t"<<PR<<endl;

    //Compute Pressure

    do{

        Pressure_star( FL, FDL, DL, PL, CL, P_Guess );
        Pressure_star( FR, FDR, DR, PR, CR, P_Guess );

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


/*----------------------------------------------------------------------------------------------------------------*/
/*                                                                                                                */
/*       To Sample the solution throughout the wave pattern. Pressure PM and velocity UM in the                   */
/*              Star Region are known. Sampling is performed in terms of the 'speed' S = X/T.                     */
/*                                          Sampled values are D, U, P                                            */
/*                                                                                                                */
/*----------------------------------------------------------------------------------------------------------------*/


double *U_star_L  ( double DL, double PL, double *UPM) {

    double *DUP = new double [3];
    double UM= UPM[0];
    double PM= UPM[1];
    double PML;

    //Sampling point lies to the left of the contact discontinuity
    if (PM <= PL){
        //left Rarefraction
        // Sampled point is Star Left state

        DUP[0] = DL*pow((PM/PL), (1.0/gam));
        DUP[1] = UM;
        DUP[2] = PM;

        }
    else{

        //left shock
        PML = PM / PL;

        // Sampled point is Star Left state
        DUP[0] = DL*(PML + G6)/(PML * G6 + 1.0);
        DUP[1] = UM;
        DUP[2] = PM;


        }


    return DUP;


}


double *U_star_R(double DR,  double PR, double *UPM) {

    double *DUP = new double [3];
    double UM= UPM[0];
    double PM= UPM[1];

    double PMR;

    if (PM > PR){

        //Right shock

        PMR = PM/PR;
        //Sampled point is Star Right state

        DUP[0] = DR * (PMR + G6)/(PMR * G6 + 1.0);
        DUP[1] = UM;
        DUP[2] = PM;

    }

    else{

        //Right Rarefraction

        //Sampled point is Star Right state

        DUP[0] = DR * pow((PM/PR), (1.0/gam));
        DUP[1] = UM;
        DUP[2] = PM;

    }
    return DUP;

}



double  *Primitive_varaible( double DR, double DL,double UR, double UL, double PL, double PR, double CL, double CR, double S, double const *UPM) {

    double *DUP = new double [3];
    double UM= UPM[0];
    double PM= UPM[1];
    double SHL, SHR, SL, SR, STL, STR;
    double CML, CMR, PML, PMR;
    double C;

    if(S <= UM){

        //Sampling point lies to the left of the contact discontinuity
        if (PM <= PL){

            //left Rarefraction
            SHL = UL - CL;

            if (S <= SHL ){

                //Sampled point is left data state
                //    cerr<<"I left data state"<<endl;
                DUP[0] = DL;
                DUP[1] = UL;
                DUP[2] = PL;
            }
            else{

                CML = CL* pow((PM/PL), G1);
                STL = UM - CML;

                if (S > STL){

                    //Sampled point is Star Left state
                    //  cout<<"I am at rarefraction left star state "<<endl;

                  //  cout<< DL<< "*"<<"pow"<<"("<<"("<<PM<<"/"<<PL<<")"<<"."<<"("<<1.0<<"/"<<gam<<")"<<")";
                    DUP[0] = DL*pow((PM/PL), (1.0/gam));
                    DUP[1] = UM;
                    DUP[2] = PM;
                    //  cout<<"PM"<<PM<<endl;
                    // cout<<"PL"<<PL<<endl;
                  //  cout<<"density"<<DUP[0]<<endl;
                }
                else{

                    //     cerr<<"I left fan"<<endl;
                    //Sampled point is inside left fan

                    C = G5*(CL + G7*(UL - S));
                    DUP[0] = DL* pow((C/CL), G4);
                    DUP[1] = G5*(CL + G7*UL + S);
                    DUP[2] = PL*pow((C/CL), G3);

                    //    cerr<<"pressure being updated "<<DUP[2] <<"\t=\t"<< "\tPL\t"<<PL<<"\t*\t"<<"pow\t"<<"(\t"<<"\t(\t"<<"\tC\t"<<C<<"\t/\t"<<"CL"<<CL<<"\t)\t"<<"\t,\t"<<"\tG3\t"<< G3<<"\t)\t"<<endl;

                }
            }

        }

            //left shock

        else{

            PML = PM / PL;
            SL  = UL - CL * sqrt(G2 * PML + G1);

            if( S <= SL ){

                //Sampled point is left data state

                DUP[0] = DL;
                DUP[1] = UL;
                DUP[2] = PL;

            }
            else{

                // Sampled point is Star Left state
                //   cout<<"I am at shock left star state "<<endl;
                DUP[0] = DL*(PML + G6)/(PML * G6 + 1.0);
                DUP[1] = UM;
                DUP[2] = PM;


                //   cout<<"PM"<<PM<<endl;
                // cout<<"PL"<<PL<<endl;
                //   cout<<"density"<<DUP[0]<<endl;

            }
        }

    }

        //Sampling point lies to the right of the contact discontinuity
    else{

        if (PM > PR){

            //Right shock

            PMR = PM/PR;
            SR  = UR + CR * sqrt(G2 * PMR + G1);

            if ( S >= SR ){

                //Sample point is right data state

                DUP[0] = DR;
                DUP[1] = UR;
                DUP[2] = PR;

            }
            else{

                //Sampled point is Star Right state
                //    cout<<"I am at right right shock star state"<<endl;
                DUP[0] = DR * (PMR + G6)/(PMR * G6 + 1.0);
                DUP[1] = UM;
                DUP[2] = PM;

            }

        }

        else{

            //Right Rarefraction
            SHR = UR + CR;

            if (S >= SHR ){

                //Sampled point is right data state

                DUP[0] = DR;
                DUP[1] = UR;
                DUP[2] = PR;

            }
            else{

                CMR = CR * pow((PM/PR), G1);
                STR = UM + CMR;

                if (S <= STR ){

                    //Sampled point is Star Right state

                    //    cout<<"I am at right reareraction star state"<<endl;
                    DUP[0] = DR * pow((PM/PR), (1.0/gam));
                    DUP[1] = UM;
                    DUP[2] = PM;



                }
                else{

                    //Sampled point is inside left fan

                    C = G5*(CR - G7*(UR - S));
                    DUP[0] = DR * pow((C/CR),G4);
                    DUP[1] = G5 * (-CR + G7*UR + S);
                    DUP[2] = PR * pow((C/CR), G3);
                }

            }

        }
    }

    return DUP;
}

/* ------------- inital condition and boundary condiitons  ----------- */

void initialConditions_ERS(int N, double dx, double *speed, double *density, double *pressure, double *energy){

    double EL = (PL/(DL * ( gam - 1.0)))+(UL * UL)/2.0 ;
    double ER = (PR/(DR * ( gam - 1.0)))+(UR * UR)/2.0 ;

    for (int i = 0; i < N; i++){

        double x=(i)*dx+dx/2.0;

        if (x<x0){
            speed[i] = UL;
            density[i] = DL;
            pressure[i] = PL;
            energy[i]=EL;
        }
        else{
            speed[i] = UR;
            density[i] = DR;
            pressure[i] = PR;
            energy[i] = ER;
        }

    }

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


//<double>( *FiniteDifferenceSolver)(int, vector<double>, double, double, double);
//switch()

//Claim plotfile from fstream and get ready for output
fstream plotfile;




int main(void){

    double dx = L/N;

    double *uspeed = new double[N];
    double *density = new double[N];
    double *pressure = new double[N];
    double *energy = new double[N];
    double *internal_energy = new double[N];


    // Level-set function

    double *phi=new double[N];
    double *phi_extend ;

    initial_condition_level_set(dx, phi, N );


    //Initial conditions for Exact Solver

    initialConditions_ERS(N, dx, uspeed,density , pressure, energy);


    //Split the primitive variables in parts depending on the interfae

    double *density1= new double[N];
    double *density1_extend = new double[N+2];
    double *density2 = new double[N];
    double *density2_extend = new double[N+2];
    double *uspeed1 = new double[N];
    double *uspeed1_extend = new double[N+2];
    double *uspeed2 = new double[N];
    double *uspeed2_extend = new double[N+2];
    double *pressure1 = new double[N];
    double *pressure1_extend = new double[N+2];
    double *pressure2 = new double[N];
    double *pressure2_extend = new double[N+2];


//    Display(U, flux, uspeed,density, pressure, energy,N);

//    Display_Ghost_fluid( U1, temp1);
//    Display_Ghost_fluid( U2, temp2);


    double *uspeed_extend ;
    double *density_extend ;
    double *pressure_extend ;
    double *energy_extend ;

    int N1;

    double t=0;
    double dt ;

    N1 = interface_finding(phi, N);
    cout<< N1<<endl;

    double XPOS ;
    double S ;

    double *UPM_temp ;

    Display_ERS( uspeed, density, pressure, energy,N);

    uspeed_extend = boundaryConditions2(uspeed);

    uspeed1_extend = boundaryConditions2(uspeed1);

    uspeed1_extend = boundaryConditions2(uspeed2);

    density_extend = boundaryConditions2(density);

    density1_extend = boundaryConditions2(density1);

    density2_extend = boundaryConditions2(density2);

    pressure_extend = boundaryConditions2(pressure);

    pressure1_extend = boundaryConditions2(pressure1);

    pressure2_extend = boundaryConditions2(pressure2);

    energy_extend = boundaryConditions2(energy);

 //   dt =dt_func(dx, uspeed_extend, pressure_extend, density_extend);


    /* ----------------level_set function implementation--------------*/
    phi_extend = boundary_condition_level_set(phi, N);

    //
    //Stationary Boundary therefore no need to consider evolution equation
    //

 //   phi = Hamilton_Jacobi(N, phi_extend, dx, dt, uspeed);
  //  phi_extend = boundary_condition_level_set(phi, N);


    //Obtain the intermediate states

    double *DUP_R = new double [3];
    double *DUP_L = new double [3];

    for (int i = 0; i<N ; i++ ){

        if (pressure_extend[i] != pressure_extend[i+1] || density_extend[i] != density_extend[i+1] || uspeed_extend[i] != uspeed_extend[i+1]){

            double CL = sqrt(gam * (pressure_extend[i]/density_extend[i]));
            double CR = sqrt(gam * (pressure_extend[i+1]/density_extend[i+1]));

            UPM_temp = Newton_Raphson( density_extend[i+1], density_extend[i], uspeed_extend[i+1], uspeed_extend[i], pressure_extend[i], pressure_extend[i+1], CL, CR, Tolerance);

            DUP_L = U_star_L(density_extend[i],  pressure_extend[i], UPM_temp);
            DUP_R = U_star_R(density_extend[i+1], pressure_extend[i+1], UPM_temp);

            }

    }


    //checker loop

    cout<< " left star states density  \t"<<  DUP_L[0]<<endl;
    cout<< "right star state density \t "<< DUP_R[0]<<endl;



    //Populate the Intermediate state into all ghost fluid regions by Fast sweeping method

    Ghost_Fluid_Boundary (DUP_L[0], DUP_R[0], density1_extend, density2_extend, density_extend, phi_extend, N, N1);

    Ghost_Fluid_Boundary (DUP_L[1], DUP_R[1], uspeed1_extend, uspeed2_extend, uspeed_extend, phi_extend, N, N1);

    Ghost_Fluid_Boundary (DUP_L[2], DUP_R[2], pressure1_extend, pressure2_extend , pressure_extend, phi_extend, N, N1);


    phi=reimplement_levelset(tStop, UPM_temp[0], dx,  N);
    N1 = interface_finding(phi, N);
    cout<< N1<<endl;

    for ( int i =0; i<N1; i++){

        double *DUP_temp1;
        double CL1 = sqrt(gam * (pressure1_extend[i]/density1_extend[i]));
        double CR1 = sqrt(gam * (pressure1_extend[i+1]/density1_extend[i+1]));

        //  cout<<"\tCL\t"<<CL<<endl;
        XPOS = (i-0.5) *dx;
        S = (XPOS -x0)/tStop;

        //cout<<" wave speed \t"<<S<<endl;
        DUP_temp1 = Primitive_varaible(density1_extend[i+1], density1_extend[i],uspeed1_extend[i+1], uspeed1_extend[i], pressure1_extend[i], pressure1_extend[i+1], CL1, CR1, S, UPM_temp);

        density_extend[i+1]= DUP_temp1 [0];
        uspeed_extend[i+1] = DUP_temp1 [1];
        pressure_extend[i+1] = DUP_temp1 [2];

    }

    for (int i = N1; i<N; i++){

        double *DUP_temp2;
        double CL2 = sqrt(gam * (pressure2_extend[i]/density2_extend[i]));
        double CR2 = sqrt(gam * (pressure2_extend[i+1]/density2_extend[i+1]));

        //  cout<<"\tCL\t"<<CL<<endl;
        XPOS = (i-0.5) *dx;
        S = (XPOS -x0)/tStop;


        DUP_temp2 = Primitive_varaible(density2_extend[i+1], density2_extend[i],uspeed2_extend[i+1], uspeed2_extend[i], pressure2_extend[i], pressure2_extend[i+1], CL2, CR2, S, UPM_temp);

        density_extend[i+1]= DUP_temp2 [0];
        uspeed_extend[i+1] = DUP_temp2 [1];
        pressure_extend[i+1] = DUP_temp2 [2];

    }



        // double distance = (UPM_temp[0] * tStop) / dx;
   // cout << "distance travelled "<<distance;
   // N1 = N1 +  18;


    for (int i =0; i<N;i++){

        density[i]= density_extend[i+1];
        uspeed[i]= uspeed_extend[i+1];
        pressure[i]= pressure_extend[i+1];
    }

    for (int i =0 ; i<N;i++){

        internal_energy [i] = pressure[i]/(density[i]*(gam -1));
    }


    Display_ERS( uspeed, density, pressure, internal_energy,N);

    plotfile.open("plotmixed_material_Riemann_400.dat", fstream::out);

    if(plotfile.fail()){

        cout<<"Error opening file"<<endl;
        exit(0);
    }

    for (int i=0;i<N;i++){

        double x=(i-1)*dx;

        plotfile<<x<<"\t"<<uspeed[i]<<"\t"<<density[i]<<"\t"<<pressure[i]<<"\t"<<internal_energy[i]<<"\t" <<phi[i]<<endl;
    }



    delete[] uspeed;
    delete[] density;
    delete[] pressure;
    delete[] energy;
    delete[] uspeed_extend;
    delete[] density_extend;
    delete[] pressure_extend;
    delete[] energy_extend;



    return 0;


}
