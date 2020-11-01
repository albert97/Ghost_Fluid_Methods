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
#include <float.h>
#include <math.h>

/*--------------------------------------- set up the parameters --------------------------------------- */
using namespace std;
int N = 100; // number of cells
double L = 1.0;   // domain length
double W = 1.0;   // domain width
double CFL = 0.9;
double tStop = 0.0014; //
double s1 = 0.05; //Intial Shock Location
double IFrad = (35.0/180.0) * M_PI;
double x2 = 0.4;
double x3 = 0.6;
double gam1 = 1.4;
double gam2 = 1.67;
double DensityL = 1.3333;
double DensityM1 = 1.0;
double DensityM2 = 0.1379;
double DensityR = 1;
double VelocityNL = 0.3535 * sqrt(100000);
double VelocityNM1 = 0.0;
double VelocityNM2 = 0.0;
double VelocityNR= 0.0;
double VelocityTL = 0.0;
double VelocityTM1 = 0.0;
double VelocityTM2 = 0.0;
double VelocityTR= 0.0;
double PressureL= 1.5 * pow(10.0, 5);
double PressureM1 = 1.0 * pow(10.0, 5);
double PressureM2 = 1.0 * pow(10.0, 5);
double PressureR= 1.0 * pow(10.0, 5);
double w= 0;
double Tolerance = 10.0e-6;


/*-------------------------------------Initial Condition level-set -----------------------------------------*/

void initial_condition_level_set (double dx, vector<vector<double>> & phi, int N){

    double phi1_temp, phi2_temp;

    for (int i =0 ; i< N; i++){

        double tempangle = (double (i) * (W/N)/tan(IFrad));

        // cout<< " tempangle \t"<< tempangle <<endl;

        for (int j =0 ; j<N; j++){

            double x = j *dx + 0.5*dx ;

            //  cout<< "X\t"<<x<<endl;

            double x2_new = x2 - tempangle;
            double x3_new = x3 - tempangle;

            phi1_temp = x2_new - x ;

            phi2_temp = x - x3_new ;

            phi[i][j] = max( phi1_temp, phi2_temp);


        }


    }


}

/*----------------------------------------------Sign function-----------------------------------------------*/

double sng(double Psi){

    if (Psi>=0){return 1;}
    if (Psi<0){ return -1;}

}

/*----------------------------------------------Sign function 2 for fast sweeping -----------------------------------------------*/

double sngFS(double Psi){

    if (Psi>=0){return 1;}
    if (Psi==0){return 0;}
    if (Psi<0){ return -1;}

}


/*-----------------------------------------Level-set Equation-------------------------------------------*/



/*--------------------------Level-set Equation----------------------------*/

vector<vector<double>>Hamilton_Jacobi(int N, vector<vector<double>> const &phi_extend, double dx, double dt, vector<vector<double>>const & v, vector<vector<double>>const &u){

    vector<vector<double>> phi_new(N, vector<double> (N, 0.0));

    double phidx;
    double phidy;

    for (int i = 0 ; i<N; i++){

        for (int j=0; j<N; j++){

            if (u[i][j] >= 0 && v[i][j] >= 0){

                phidx = (phi_extend[i+1][j+1] - phi_extend[i+1][j]);
                phidy = (phi_extend[i+1][j+1] - phi_extend[i][j+1]);

                // cout<<phi_new[i][j]<<"\t=\t"<<phi_extend[i+1][j+1]<< "\t-\t"<< dt/dx<< "\t*\t"<< v[i][j]<< "\t*\t"<<"\t ( \t" <<phi_extend[i+1][j+1] << "\t- \t"<<phi_extend[i+1][j] <<"\t)\t"<<endl;

            }
            else{
                if(u[i][j] < 0 && v[i][j] < 0){

                    phidx =(phi_extend[i+1][j+2] - phi_extend[i+1][j+1]);
                    phidy =(phi_extend[i+2][j+1] - phi_extend[i+1][j+1]);

                }
                else{
                    if(u[i][j] >= 0 && v[i][j] <0){

                        phidx = (phi_extend[i+1][j+1] - phi_extend[i+1][j]);
                        phidy = (phi_extend[i+2][j+1] - phi_extend[i+1][j+1]);
                    }
                    else{
                        phidx =(phi_extend[i+1][j+2] - phi_extend[i+1][j+1]);
                        phidy = (phi_extend[i+1][j+1] - phi_extend[i][j+1]);

                    }
                }

                //    cout<<phi_new[i][j]<<"\t=\t"<<phi_extend[i+1][j+1]<< "\t-\t"<< dt/dx<< "\t*\t"<< v[i][j]<< "\t*\t"<<"\t ( \t" <<phi_extend[i+1][j+2] << "\t- \t"<<phi_extend[i+1][j+1] <<"\t)\t"<<endl;

            }

            phi_new[i][j] = phi_extend[i+1][j+1] - dt/dx * (u[i][j] * phidx + v[i][j]* phidy);

        }
    }

    return phi_new;

}

/*--------------------------Display fucntion level-set function ------------------------*/
void Display_levelset (vector<vector<double>>const &U_data, int N){


    cout<<"\n"<<endl;

    for (int i=N-1;i >=0;i--){

        for(int j=0; j<N;j++){

            cout<<U_data[i][j]<<"\t";
        }
        cout<<"\n"<<endl;

    }


}




void initialConditions_ERS(int N, double dx, vector<vector<double>>& uspeed, vector<vector<double>>& vspeed, vector<vector<double>>& density, vector<vector<double>>& pressure, vector<vector<double>>&energy){

    //
    //Split the Normal and tangential component into velocity in x and y direction
    //

    double VelocityUL,VelocityUM1,VelocityUM2,VelocityUR;
    double VelocityVL,VelocityVM1,VelocityVM2,VelocityVR;

    VelocityUL = cos(0.5 * M_PI - IFrad) * VelocityNL + cos(IFrad) * VelocityTL;
    VelocityUM1 = cos(0.5 * M_PI - IFrad) * VelocityNM1 + cos(IFrad) * VelocityTM1;
    VelocityUM2 = cos(0.5 * M_PI - IFrad) * VelocityNM2 + cos(IFrad) * VelocityTM2;
    VelocityUR = cos(0.5 * M_PI - IFrad) *  VelocityNR + cos(IFrad) *  VelocityTR;
    VelocityVL = sin(0.5 * M_PI - IFrad) * VelocityNL + sin(IFrad) * VelocityTL;
    VelocityVM1 = sin(0.5 * M_PI - IFrad) * VelocityNM1 + sin(IFrad) * VelocityTM1;
    VelocityVM2 = sin(0.5 * M_PI - IFrad) * VelocityNM2 + sin(IFrad) * VelocityTM2;
    VelocityVR = sin(0.5 * M_PI - IFrad) *  VelocityNR + sin(IFrad) *  VelocityTR;

    double EL = (PressureL/(DensityL * ( gam1 - 1.0)))+((VelocityUL * VelocityUL)+ (VelocityVL * VelocityVL))/2.0 ;
    double EM1 = (PressureM1/(DensityM1 * ( gam1 - 1.0)))+((VelocityUM1 * VelocityUM1) + (VelocityVM1 * VelocityVM1))/2.0 ;
    double EM2 = (PressureM2/(DensityM2 * ( gam2 - 1.0)))+((VelocityUM2 * VelocityUM2)+ (VelocityVM2 * VelocityVM2))/2.0 ;
    double ER = (PressureR/(DensityR * ( gam1 - 1.0)))+ ((VelocityUR * VelocityUR) + (VelocityVR * VelocityVR) )/2.0 ;

    for (int i = 0; i < N; i++) {

        double tempangle = (double(i) * (W / N) / tan(IFrad));

        for (int j = 0; j < N; j++) {

            double x = j * dx + 0.5*dx;
            double x2_new = x2 - tempangle;
            double x3_new = x3 - tempangle;
            double s1_new = s1 - tempangle;

            if (x2_new >= 0.0 && x3_new >= 0.0) {

                if (x < x3_new) {

                    if (x < x2_new) {

                        if ( s1_new >= 0 && x <= s1_new) {

                            uspeed[i][j] = VelocityUL;
                            vspeed[i][j] = VelocityVL;
                            density[i][j] = DensityL;
                            pressure[i][j] = PressureL;
                            energy[i][j] = EL * DensityL;

                        } else {

                            uspeed[i][j] = VelocityUM1;
                            vspeed[i][j] = VelocityVM1;
                            density[i][j] = DensityM1;
                            pressure[i][j] = PressureM1;
                            energy[i][j] = EM1 * DensityM1 ;


                        }
                    } else {

                        uspeed[i][j] = VelocityUM2;
                        vspeed[i][j] = VelocityVM2;
                        density[i][j] = DensityM2;
                        pressure[i][j] = PressureM2;
                        energy[i][j] = EM2 * DensityM2;

                    }
                }
                else {

                    uspeed[i][j] = VelocityUR;
                    vspeed[i][j] = VelocityVR;
                    density[i][j] = DensityR;
                    pressure[i][j] = PressureR;
                    energy[i][j] = ER * DensityR;
                }

            } else {

                if (x2_new < 0.0 && x3_new >= 0.0) {

                    if (x <= x3_new) {

                        uspeed[i][j] = VelocityUM2;
                        vspeed[i][j] = VelocityVM2;
                        density[i][j] = DensityM2;
                        pressure[i][j] = PressureM2;
                        energy[i][j] = EM2 * DensityM2;

                    } else {

                        uspeed[i][j] = VelocityUR;
                        vspeed[i][j] = VelocityVR;
                        density[i][j] = DensityR;
                        pressure[i][j] = PressureR;
                        energy[i][j] = ER * DensityR;

                    }
                }

                else {

                    uspeed[i][j] = VelocityUR;
                    vspeed[i][j] = VelocityVR;
                    density[i][j] = DensityR;
                    pressure[i][j] = PressureR;
                    energy[i][j] = ER * DensityR;

                }

            }


        }

    }
}



/*--------------------------------------Boundary Condition level-set ---------------------------------------*/

vector<vector<double>>boundary_condition_level_set( vector<vector<double>>const &Property, int N ){


    vector<vector<double>>  Property_extend(N+2, vector<double> (N+2, 0.0));


    for (int i = 0; i<N; i++){

        for(int j =0;j<N;j++){

            Property_extend[i+1][j+1] = Property[i][j];

        }

    }


    //fill the first and last row
    for (int j=0; j<N;j++){

        Property_extend[0][j+1] =  Property[0][j];
        Property_extend[N+1][j+1] =  Property[N-1][j];
    }

    //fill the first and last Column
    for (int i=0; i<N;i++){

        Property_extend[i+1][0] =  Property[i][0];
        Property_extend[i+1][N+1] =  Property[i][N-1];
    }


    Property_extend[0][0] = Property_extend[0][1];
    Property_extend[0][N+1] = Property_extend[0][N];
    Property_extend[N+1][0] = Property_extend[N+1][1];
    Property_extend[N+1][N+1] = Property_extend[N+1][N];

    return Property_extend;

}



/* ----------This follwowing function allows you to print results on screen -------- */


void Display_ERS ( vector<vector<double>> const &data, vector<vector<double>> const &data2, vector<vector<double>> const &data3, vector<vector<double>>const &data4, vector<vector<double>> const &data5, int N){


    cout<<"\n"<<endl;



    for (int i=N-1;i >=0;i--){

        for(int j=0; j<N;j++){

            cout<<data[i][j]<<"\t";
        }
        cout<<"\n"<<endl;

    }


    cout<<"\n"<<endl;


    for (int i=N-1;i >=0;i--){

        for(int j=0; j<N;j++){

            cout<<data2[i][j]<<"\t";
        }
        cout<<"\n"<<endl;

    }



    cout<<"\n"<<endl;



    for (int i=N-1;i >=0;i--){

        for(int j=0; j<N;j++){

            cout<<data3[i][j]<<"\t";
        }
        cout<<"\n"<<endl;

    }



    cout<<"\n"<<endl;


    for (int i=N-1;i >=0;i--){

        for(int j=0; j<N;j++){

            cout<<data4[i][j]<<"\t";
        }
        cout<<"\n"<<endl;

    }



    cout<<"\n"<<endl;


    for (int i=N-1;i >=0;i--){

        for(int j=0; j<N;j++){

            cout<<data5[i][j]<<"\t";
        }
        cout<<"\n"<<endl;

    }
}
void Display_2DarrayN4(vector<vector<double>>data){

    for (int i=N+3;i >=0;i--){

        for(int j=0; j<N+4;j++){

            cout<<data[i][j]<<"\t";

        }
        cout<<"\n"<<endl;

    }
}

void Display_2DarrayN(vector<vector<double>> const &data){

    for (int i=N-1;i >=0;i--){

        for(int j=0; j<N;j++){

            cout<<data[i][j]<<"\t";

        }
        cout<<"\n"<<endl;

    }
}

vector<vector<double>>boundaryConditions2( vector<vector<double>> const &Property){


    vector<vector<double>> Property_extend(N+2, vector<double> (N+2, 0.0));


    for (int i = 0; i<N; i++){

        for(int j = 0; j< N; j++){

            Property_extend[i+1][j+1] = Property[i][j];

        }

    }



    //Fill the first and last row
    for(int j=0; j<N; j++){

        Property_extend[0][j+1] = Property[0][j];
        Property_extend[N+1][j+1] = Property[N-1][j];
    }


    //Fill the first and last column
    for(int i=0; i<N; i++){

        Property_extend[i+1][0] = Property[i][0];
        Property_extend[i+1][N+1] = Property[i][N-1];
    }

    //Asiign values to the four corners of the extended 2D array

    Property_extend[0][0] = Property_extend[0][1];
    Property_extend[0][N+1] = Property_extend[0][N];
    Property_extend[N+1][0] = Property_extend[N+1][1];
    Property_extend[N+1][N+1] = Property_extend[N+1][N];

    return Property_extend;
}


/*---------------------------------------------Finding the interface in X direction --------------------------------------*/


vector<vector<int>> interface_findingX(vector<vector<double>>const &phi, int N){


    vector<vector<int>>interface(2, vector<int> (N, 0));


    for (int i = 0; i<N; i++){

        int k =0 ;

        for(int j =0; j<N-2; j++){

            if (phi[i][0] <=0){

                interface[0][i] = 0;
                if (sng(phi[i][j]) != sng(phi[i][j+1]) ) {

                    interface[1][i] = j+1 ;

                }
            }
            else{

                //
                //If the loop has reached the end then there is no interface
                //
                if( j == N-3) {

                    interface[0][i] = 0;
                    interface[1][i] = 0;

                }
                else{

                    if (sng(phi[i][j]) != sng(phi[i][j+1])) {

                        interface[k][i] = j+1 ;

                        k++;
                        if (k == 2){

                            break;
                        }
                    }

                }

            }

        }

    }


    return interface;

}


/*---------------------------------------------Finding the interface in Y direction --------------------------------------*/


vector<vector<int>>interface_findingY(vector<vector<double>>const &phi, int N){


    vector<vector<int>>interface(2, vector<int> (N, 0));


    for (int i = 0; i<N; i++){

        int k =0 ;

        //  cout<<"I am wroking here1"<<endl;

        for(int j =0; j<N-2; j++){


     //       cout<< "\tj\t"<<j<<"\t"<<phi[j][i]<<"\t"<<phi[j+1][i]<<"\t"<<phi[j+2][i]<<"\t"<<sng(phi[j][i])<<"\t"<<sng(phi[j+1][i])<<"\t"<<sng(phi[j+2][i])<<"\t";

            if(phi[0][i] >=0){

                if (sng(phi[j][i]) != sng(phi[j+1][i]) ) {

                    interface[k][i] = j+1 ;

                    k++;
                }

            }
            else{

                interface[0][i] = 0;
                if (sng(phi[j][i]) != sng(phi[j+1][j]) ) {

                    interface[1][i] = j+1 ;
                }

            }

        }
        //    cout<<"i\t"<<i<<endl;

        //    cout<<interface[0][i]<<"\t"<<interface[1][i]<<endl;

    }


    return interface;

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

        //    cout<<  F <<"= "<<G4(gamk) <<"*"<< CK<< ""<< "("<<pow(PRAT, G1(gamk)) <<"-"<< 1.0<<")"<<endl;

        FD = (1.0 / (DK * CK)) * pow(PRAT, -G2(gamk));

    }

}


/*--------------------------Calculate the solution for pressure and velocity in the star Region ------------------------------*/

vector<double>Newton_Raphson( double DR, double DL,double UR, double UL, double PL, double PR, double CL, double CR, double gammaL, double gammaR, double Tolerance ){

    vector<double>UPM(2);
    double P_star;
    double dU;
    double FL=0.0, FR=0.0;
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

        //  cout<< "FL\t"<<FL<<endl;

        Pressure_star( FR, FDR, DR, PR, CR, P_Guess, gammaR );
        // cout<<"FR\t"<<FR << "FL\t"<<FL<<endl;

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
    if (abs( FL) < 0.000001){
        FL= 0.0;
    }

    U_star = 0.5 * ( UL + UR + FR - FL);

    if (abs(U_star) < 0.00000000001 ){
        U_star =0.0;
    }

    //   cout<<"FR\t"<<FR << "FL\t"<<FL<<endl;

    //  cout<< U_star <<"=" <<0.5<< "*"<< "(" <<UL<< "+"<< UR<< "+"<< FR<< "-"<< FL<<")"<<endl;

    // cout<<"the new velocity is \t "<<U_star<< "\tUR\t"<<UR<<"\tUL\t"<<UL<<endl;
    // cout<<"the new pressire is \t "<<P_star<<"\tPl\t"<<PL<<"\tPR\t"<<PR<<endl;

    UPM[0] = U_star;
    UPM[1] = P_star;

    return UPM;
}


vector<double> U_star_L ( double DL, double PL, vector<double>const &UPM, double gamL) {

    vector<double>DUP(3);
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


vector<double>U_star_R(double DR,  double PR, vector<double>const &UPM, double gamR) {

    vector<double>DUP(3);
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




double bilinear_interpolation( int X1, int X2, int Y1, int Y2, vector<vector<double>>const &property){

    //Linear interpolation

    //  cout<<"\t Bottom left\t"<< property[X1][Y1]<<"\t Bottom right \t"<<property[X2][Y1]<<"\tTop left\t"<<property[X1][Y2]<<"\t Top right \t"<<property[X2][Y2]<<endl;
    double fxy1=  0.5 * property[X1][Y1] + 0.5 * property[X2][Y1];
    double fxy2= 0.5 * property[X1][Y2] + 0.5 * property[X2][Y2];

    //Linear in Y direction and Bilinear as a whole

    double fxy;

    fxy = 0.5 * fxy1 + 0.5* fxy2;

    return fxy;
}




void Reimann_boundary_extrapolation(int X, int Y,int N, double &PL, double &PR, const double dx,  vector<vector<double>>const &property, vector<vector<double>> const &property_extend,  vector<vector<double>>const &phi, const double *Normal){

    //Calcualte the positon of the interface with the cell

    vector<double>XP(2);

    //     cout<<"\tX\t"<<X<<"\tY\t"<<Y<<endl;

    //x is the row number while y is your column number

    double x  = double(X) / double(N) + 0.005;
    double y = double(Y) /double(N) + 0.005;

    //     cout<<"\tx\t"<<x<<"\ty\t"<<y<<endl;

    XP[0] = y - Normal[0]*phi[X][Y];
    XP[1] = x - Normal[1]*phi[X][Y];

//    cout<<"location of the interface"<<"\t x \t"<< XP[0]<<"\t y\t"<< XP[1]<<endl;
    //Find the location of your extrapolated cell


    vector<double>XP1(2);
    vector<double>XP2(2);

    for (int i =0; i<2;i++){


        //Right materal


        // XP1[0] is the column extraction (i.e in x direction) while XP[1] is the row extractipn in y directiom
        XP1[i] = XP[i] + 1.5 * dx * Normal[i];

        //Left materal

        XP2[i] = XP[i] - 1.5 * dx * Normal[i];



    }

    double XP1_temp1 =  XP1[0] * N; // The y value
    double XP1_temp2 =  XP1[1] * N; // The x value

    double XP2_temp1 =  XP2[0] * N; // The y value
    double XP2_temp2 =  XP2[1] * N; // The x value

//    cout<< XP1_temp1<<"\t"<<XP1_temp2<<"\t"<< XP2_temp1<<"\t"<<XP2_temp2<<endl;

    int XL1,XL2,YL1,YL2;
    int XR1,XR2,YR1,YR2;

    if ( Normal[0] >= 0 && Normal[1]>= 0){

        XR1 = int(XP1_temp1);
        YR1 = int(XP1_temp2) ;

        XR2 = int(XP1_temp1)+1;
        YR2 = int(XP1_temp2)+1;

        //     cout<<"right interpolated points \t" <<"XR1\t"<<XR1<<"YR1\t"<<YR1<<"XR2\t"<<XR2<<"YR2\t"<<YR2<<"\n";

/*
    cout<<"Location of the right cell"<<endl;
    cout<<"location of the bottom right cell"<<"\t XR1\t"<<XR1<<"\t YR1 \t"<<YR1<<endl;
    cout<<"location of the bottom right cell"<<"\t XR2\t"<<XR2<<"\t YR2 \t"<<YR2<<endl;
*/

        if(XR1 > 0 && XR2 > 0 && YR1 > 0 && YR2 > 0 && XR1 <N && XR2 <N && YR1 < N && YR2 < N){


            PR=  bilinear_interpolation( YR1, YR2, XR1, XR2, property);

            //               cout<<"PR BI\t"<<PR<<endl;


        }
        else {

            PR = bilinear_interpolation( YR1+2, YR2+2, XR1+2, XR2+2, property_extend);

            //     PR = property[X][Y + 1];
            //       cout<<"PR NotBI\t"<<PR<<endl;


        }


        XL1 = int(XP2_temp1);
        YL1 = int(XP2_temp2);

        XL2 = int(XP2_temp1) - 1;
        YL2 = int(XP2_temp2) - 1;

/*
    cout<<"Location of the left cell"<<endl;
    cout<<"location of the bottom left cell"<<"\t XL1\t"<<XL1<<"\t YL1 \t"<<YL1<<endl;
    cout<<"location of the bottom left cell"<<"\t XL2\t"<<XL2<<"\t YL2 \t"<<YL2<<endl;
*/

        //    cout<<"left interpolated points \t" <<"XL1\t"<<XL1<<"YL1\t"<<YL1<<"XL2\t"<<XL2<<"YL2\t"<<YL2<<"\n";

        if(XL1 > 0 && XL2 > 0 && YL1 >0 && YL2 > 0&& XL1 <N && XL2 <N && YL1 < N && YL2 < N){


            PL =  bilinear_interpolation( YL1, YL2, XL1, XL2, property);

        }
        else{

            PL= bilinear_interpolation(YL1+2, YL2+2, XL1+2, XL2+2, property_extend);

            //    PL = property[X][Y-1];
        }
        //              cout<<"PL\t"<<PL<<endl;



    }
    else{

        XL1 = int(XP1_temp1);
        YL1 = int(XP1_temp2) ;

        XL2 = int(XP1_temp1)-1;
        YL2 = int(XP1_temp2)-1;

/*
        cout<<"Location of the left cell"<<endl;
        cout<<"location of the bottom right cell"<<"\t XL1 \t"<<XL1<<"\t YL1 \t"<<YL1<<endl;
        cout<<"location of the bottom right cell"<<"\t XL2 \t"<<XL2<<"\t YL2 \t"<<YL2<<endl;
*/

        // cout<<"left interpolated points \t" <<"XL1\t"<<XL1<<"YL1\t"<<YL1<<"XL2\t"<<XL2<<"YL2\t"<<YL2<<"\n";
        // cout<< "i am working here3"<<endl;

        if(XL1 > 0 && XL2 > 0 && YL1 > 0 && YL2 > 0 && XL1 <N && XL2 <N && YL1 < N && YL2 < N){

            PL =  bilinear_interpolation( YL1, YL2, XL1, XL2, property);

        }
        else {

            PL =  bilinear_interpolation( YL1+2, YL2+2, XL1+2, XL2+2, property_extend);

            // cout<< "i am working here4"<<endl;
            //cout<<"X-1 \t"<< X<<"\t Y \t"<<Y-1<<endl;
            //       PL= property[X][Y -1];

        }

        //     cout<<"PL\t"<<PL<<endl;

        XR1 = int(XP2_temp1);
        YR1 = int(XP2_temp2);

        XR2 = int(XP2_temp1) + 1;
        YR2 = int(XP2_temp2) + 1;

/*
        cout<<"Location of the right cell"<<endl;
        cout<<"location of the bottom left cell"<<"\t XR1\t"<<XR1<<"\t YR1 \t"<<YR1<<endl;
        cout<<"location of the bottom left cell"<<"\t XR2\t"<<XR2<<"\t YR2 \t"<<YR2<<endl;
*/
        //       cout<<"right interpolated points \t" <<"XR1\t"<<XR1<<"YR1\t"<<YR1<<"XR2\t"<<XR2<<"YR2\t"<<YR2<<"\n";


        if(XR1 > 0 && XR2 > 0 && YR1 >0 && YR2 > 0 && XR1 <N && XR2 <N && YR1 < N && YR2 < N){

            PR=  bilinear_interpolation( YR1, YR2, XR1, XR2, property);

        }
        else{


            PR=  bilinear_interpolation( YR1+2, YR2+2, XR1+2, XR2+2, property_extend);

            //  cout<<"X \t"<< X<<"\t Y+1 \t"<<Y+1<<endl;
            //  PR = property[X][Y+1];

        }

        //   cout<<"PR\t"<<PR<<endl;



    }


}



vector<vector<double>>transmisive_boundary_condition(vector<vector<double>>Property){

    vector<vector<double>> Property_extend(N+4, vector<double> (N+4, 0.0));


    for (int i =0 ; i< N; i++){
        for (int j=0; j<N; j++){

            Property_extend[i+2][j+2] = 0.0;

        }
    }

    for (int i =0 ; i< N; i++){
        for (int j=0; j<N; j++){

            Property_extend[i+2][j+2] = Property[i][j];

        }
    }


    //fill the first two and last two rolls
    for (int j = 0 ; j<N; j++){

        Property_extend[0][j+2] = Property[0][j];
        Property_extend[1][j+2] = Property[0][j];
        Property_extend[N+2][j+2] = Property[N-1][j];
        Property_extend[N+3][j+2] = Property[N-1][j];

    }

    //fill the first two and last two columns
    for (int i = 0; i<N; i++){

        Property_extend[i+2][0] = Property[i][0];
        Property_extend[i+2][1] = Property[i][0];
        Property_extend[i+2][N+2] = Property[i][N-1];
        Property_extend[i+2][N+3] = Property[i][N-1];

    }


    //The Copying is done in a way thay the four cell values at each conor are not alighed to any values, but zeros


    //Asiign values to the four corners of the extended 2D array

    //Bottom left
    Property_extend[0][0] =  Property_extend[2][2];
    Property_extend[0][1] =  Property_extend[2][2];
    Property_extend[1][1] =  Property_extend[2][2];
    Property_extend[1][0] =  Property_extend[2][2];
    //Top right
    Property_extend[N+3][N+3] =  Property_extend[N+1][N+1];
    Property_extend[N+3][N+2] =  Property_extend[N+1][N+1];
    Property_extend[N+2][N+2] =  Property_extend[N+1][N+1];
    Property_extend[N+2][N+3] =  Property_extend[N+1][N+1];
    //Bottom Right
    Property_extend[0][N+3] =  Property_extend[3][N+1];
    Property_extend[0][N+2] =  Property_extend[3][N+1];
    Property_extend[1][N+2] =  Property_extend[3][N+1];
    Property_extend[1][N+3] =  Property_extend[3][N+1];
    //TopLeft
    Property_extend[N+3][0] =  Property_extend[N+1][2];
    Property_extend[N+3][1] =  Property_extend[N+1][2];
    Property_extend[N+2][0] =  Property_extend[N+1][2];
    Property_extend[N+2][1] = Property_extend[N+1][2];


    int X1 = 2;
    int X2 = 3;
    int Y1 = 2;
    int Y2 = 2;

    //
    //Bottom two rows
    //


    for (int j =0; j<N; j++){

        Property_extend[1][j+X1] = Property[0][j];
        Property_extend[0][j+X2] = Property[0][j];

    }

    //
    //left two columns
    //

    for (int i =0; i<N; i++){

        Property_extend[i+Y1][1] = Property[i][1];
        Property_extend[i+Y2][0] = Property[i][0];

    }
    return Property_extend;


}




void Riemann_based_Ghost_Fluid_Boundary(vector<int>IFX, vector<vector<int>>IFY, int const N, double const dx, vector<vector<double>>pressure,  vector<vector<double>>density, vector<vector<double>>uspeed, vector<vector<double>>vspeed, vector<vector<double>>normalvelocity, vector<vector<double>>pressure_extend, vector<vector<double>>density_extend, vector<vector<double>>normalvelocity_extend, vector<vector<double>>phi, vector<vector<double>>& density1, vector<vector<double>>& density2, vector<vector<double>>& uspeed1, vector<vector<double>>& uspeed2, vector<vector<double>>& vspeed1, vector<vector<double>>& vspeed2, vector<vector<double>>& pressure1, vector<vector<double>>& pressure2, const double *Normal1, const double *Normal2) {

    vector<double>UPM_temp1;
    vector<double>DUP_R1;
    vector<double>DUP_L1 ;
    vector<double>UPM_temp2;
    vector<double>DUP_R2;
    vector<double>DUP_L2 ;
    double pressureL;
    double pressureR;
    double densityL;
    double densityR;
    double nspeedL;
    double nspeedR;


    for (int i =0; i< N; i++){

        for (int j =0; j<N; j++){

            density1[i][j] = density[i][j];
            density2[i][j] = density[i][j];
            uspeed1[i][j] = uspeed[i][j];
            uspeed2[i][j] = uspeed[i][j];
            vspeed1[i][j] = vspeed[i][j];
            vspeed2[i][j] = vspeed[i][j];
            pressure1[i][j] = pressure[i][j];
            pressure2[i][j] = pressure[i][j];

        }
    }



    for (int i =0; i< N;i++){



        if (IFY[0][i]==0 && IFY[1][i] > 0  && ((double(IFY[1][i]) - double(IFY[0][i])) != (double(IFY[1][0]) - double(IFY[0][0])))){


            //        cout << "single interface" << endl;

            Reimann_boundary_extrapolation(IFX[i], IFY[1][i],  N, pressureL, pressureR, dx, pressure, pressure_extend, phi, Normal2);

            //    cout<<"\tPressure\t"<<pressureL<<"\tPressurel\t"<<pressureR<<endl;

            Reimann_boundary_extrapolation(IFX[i], IFY[1][i],  N, densityL, densityR, dx, density, density_extend, phi, Normal2);

            //       cout<<"\tdensityL\t"<<densityL<<"\tdensityR\t"<<densityR<<endl;

            Reimann_boundary_extrapolation(IFX[i], IFY[1][i],  N, nspeedL, nspeedR, dx, normalvelocity, normalvelocity_extend, phi, Normal2);

            //    cout<<"\tuspeedL\t"<<uspeedL<<"\tuspeedR\t"<<nspeedR<<endl;

            //Calculate the intermediate state at the first boundary
            double CL1 = sqrt(gam1 * (pressureL / densityL));
            double CR1 = sqrt(gam2 * (pressureR / densityR));

            UPM_temp1 = Newton_Raphson(densityR, densityL, nspeedR, nspeedL, pressureL, pressureR, CL1, CR1, gam2, gam1, Tolerance);


            DUP_L1 = U_star_L(densityL, pressureL, UPM_temp1, gam2);
            DUP_R1 = U_star_R(densityR, pressureR, UPM_temp1, gam1);

            //         cout << "\t DensityL1 \t" << DUP_L1[0] << "\t" << "DensityR1 \t" << DUP_R1[0] << endl;



       //     cout.precision(17);
       //     cout<<"DUP_R1[0]\t"<<DUP_R1[0]<<"\tDUP_L1[0]\t"<<DUP_L1[0]<< " \tDUP_R1[1]\t"<<DUP_R1[1]<<"\tDUP_L1[1]\t"<< DUP_L1[1]<<"\tDUP_R1[2]\t"<<DUP_R1[2]<<"\tDUP_L1[2]\t"<<DUP_L1[2]<<"\n";

            density1[i][IFY[1][i] -1 ] = DUP_R1[0];
            density2[i][IFY[1][i]] = DUP_L1[0];

            uspeed1[i][IFY[1][i] -1] = DUP_R1[1];
            uspeed2[i][IFY[1][i]] = DUP_L1[1];

            pressure1[i][IFY[1][i] - 1] = DUP_R1[2];
            pressure2[i][IFY[1][i]] = DUP_L1[2];

        }

        else{

            if (IFY[0][i] >= 0 && IFY[1][i] > 0){


                //         cout<<"Two interface \t"<<endl;

                if (IFY[0][i] == 0){

                    Reimann_boundary_extrapolation(IFX[i-1], IFY[0][i-1], N, pressureL, pressureR, dx, pressure, pressure_extend, phi, Normal1);

                    //           cout << "\tPressureL\t" << pressureL << "\tPressureR\t" << pressureR << endl;

                    Reimann_boundary_extrapolation(IFX[i-1], IFY[0][i-1], N, densityL, densityR, dx, density, density_extend,  phi, Normal1);

                    //           cout << "\tdensityL\t" << densityL << "\tdensityR\t" << densityR << endl;

                    Reimann_boundary_extrapolation(IFX[i-1], IFY[0][i-1], N, nspeedL, nspeedR, dx, normalvelocity, normalvelocity_extend, phi, Normal1);

                    //         cout << "\tuspeedL\t" << uspeedL<< "\tnspeedR\t" << nspeedR << endl;

                    double CL1 = sqrt(gam1 * (pressureL / densityL));
                    double CR1 = sqrt(gam2 * (pressureR / densityR));

                    UPM_temp1 = Newton_Raphson(densityR, densityL, nspeedR, nspeedL, pressureL, pressureR, CL1, CR1, gam1, gam2, Tolerance);


                    DUP_L1 = U_star_L(densityL, pressureL, UPM_temp1, gam1);

                    DUP_R1 = U_star_R(densityR, pressureR, UPM_temp1, gam2);


                }

                else{

                    cout<<"\n";

                    Reimann_boundary_extrapolation(IFX[i], IFY[0][i], N, pressureL, pressureR, dx, pressure, pressure_extend, phi, Normal1);

                    //        cout << "\tPressureL\t" << pressureL << "\tPressureR\t" << pressureR << endl;

                    Reimann_boundary_extrapolation(IFX[i], IFY[0][i], N, densityL, densityR, dx, density, density_extend, phi, Normal1);

//                        cout << "\tdensityL\t" << densityL << "\tdensityR\t" << densityR << endl;

                    Reimann_boundary_extrapolation(IFX[i], IFY[0][i], N, nspeedL, nspeedR, dx, normalvelocity, normalvelocity_extend, phi, Normal1);

                    //              cout << "\tnspeedL\t" << nspeedL << "\tnspeedR\t" << nspeedR << endl;

                    double CL1 = sqrt(gam1 * (pressureL / densityL));
                    double CR1 = sqrt(gam2 * (pressureR / densityR));

                    UPM_temp1 = Newton_Raphson(densityR, densityL, nspeedR, nspeedL, pressureL, pressureR, CL1, CR1, gam1, gam2, Tolerance);


                    DUP_L1 = U_star_L(densityL, pressureL, UPM_temp1, gam1);

                    DUP_R1 = U_star_R(densityR, pressureR, UPM_temp1, gam2);

                }

                //Calculate the intermediate state at the first boundary



                //checker loop


                //Calculate the intermediate state at the second boundary

                //Due to the normal to the seond interface is the opposite to that of the first one, Your left and right values are the wrong way round

//                cout << "Second Interface \t" << endl;

                Reimann_boundary_extrapolation(IFX[i], IFY[1][i], N, pressureL, pressureR, dx, pressure, pressure_extend, phi, Normal2);

//                cout << "\tPressurer\t" << pressureL << "\tPressurel\t" << pressureR << endl;

                Reimann_boundary_extrapolation(IFX[i], IFY[1][i], N, densityL, densityR, dx, density, density_extend, phi, Normal2);

//                cout << "\tdensityR\t" << densityL << "\tdensityL\t" << densityR << endl;

                Reimann_boundary_extrapolation(IFX[i], IFY[1][i], N, nspeedL, nspeedR, dx, normalvelocity, normalvelocity_extend,  phi, Normal2);

//                 cout << "\tnspeedR\t" << nspeedL << "\tnspeedL\t" << nspeedR << endl;


                double CL2 = sqrt(gam1 * (pressureL / densityL));
                double CR2 = sqrt(gam2 * (pressureR / densityR));

                UPM_temp2 = Newton_Raphson(densityR, densityL, nspeedR, nspeedL, pressureL, pressureR, CL2, CR2, gam2, gam1, Tolerance);

                DUP_L2 = U_star_L(densityL, pressureL, UPM_temp2, gam2);

                DUP_R2 = U_star_R(densityR, pressureR, UPM_temp2, gam1);


                if (IFY[0][i] == 0){

               //     cout.precision(17);
               //     cout<<"DUP_R1[0]\t"<<DUP_R1[0]<<"\tDUP_L1[0]\t"<<DUP_L1[0]<< " \tDUP_R1[1]\t"<<DUP_R1[1]<<"\tDUP_L1[1]\t"<< DUP_L1[1]<<"\tDUP_R1[2]\t"<<DUP_R1[2]<<"\tDUP_L1[2]\t"<<DUP_L1[2]<<"\n";


                    density1[i][IFY[0][i]] = DUP_L1[0];
                    //                cout <<"new density\t"<<IFY[0][i]<<'\t'<< DUP_L1[0]<<endl;
                    density1[i][IFY[1][i]-1] = DUP_R2[0];
                    density2[i][IFY[1][i]] = DUP_L2[0];


                    uspeed1[i][IFY[0][i]] = DUP_L1[1] ;
                    uspeed1[i][IFY[1][i] -1] = DUP_R2[1];
                    uspeed2[i][IFY[1][i]] =  DUP_L2[1];
                    //    cout<<" uspeed2[i][IFY[1][i]] 1 \t"<< uspeed2[i][IFY[1][i]]<<endl;

                    pressure1[i][IFY[0][i]] = DUP_L1[2];
                    pressure1[i][IFY[1][i]-1] = DUP_R2[2];
                    pressure2[i][IFY[1][i]] = DUP_L2[2];

                }

                else {
                  //  cout.precision(17);

                  //  cout<<"DUP_R1[0]\t"<<DUP_R1[0]<<"\tDUP_L1[0]\t"<<DUP_L1[0]<< " \tDUP_R1[1]\t"<<DUP_R1[1]<<"\tDUP_L1[1]\t"<< DUP_L1[1]<<"\tDUP_R1[2]\t"<<DUP_R1[2]<<"\tDUP_L1[2]\t"<<DUP_L1[2]<<"\n";

                    density1[i][IFY[0][i]] = DUP_L1[0];
                    //         cout << "new density\t" << IFY[0][i] << '\t' << DUP_L1[0] << endl;
                    density1[i][IFY[1][i] - 1] = DUP_R2[0];
                    density2[i][IFY[0][i] - 1] = DUP_R1[0];
                    density2[i][IFY[1][i]] = DUP_L2[0];


                    uspeed1[i][IFY[0][i]] = DUP_L1[1];
                    uspeed1[i][IFY[1][i] - 1] = DUP_R2[1];
                    uspeed2[i][IFY[0][i] - 1] = DUP_R1[1];
                    uspeed2[i][IFY[1][i]] = DUP_L2[1];


                    pressure1[i][IFY[0][i]] = DUP_L1[2];
                    pressure1[i][IFY[1][i] - 1] = DUP_R2[2];
                    pressure2[i][IFY[0][i] - 1] = DUP_R1[2];
                    pressure2[i][IFY[1][i]] = DUP_L2[2];

                }

            }
            else{
                break;
            }

        }

    }

}



void Combinevelocityvectors( vector<vector<int>>const &IFY, vector<vector<double>> const &tangentialvelocity , vector<vector<double>>& uspeed1, vector<vector<double>>& uspeed2, vector<vector<double>>& vspeed1, vector<vector<double>>& vspeed2){


    for (int i =0; i< N;i++){


        //Single interface

        if (IFY[0][i]==0 && IFY[1][i] > 0  && ((double(IFY[1][i]) - double(IFY[0][i])) != (double(IFY[1][0]) - double(IFY[0][0])))) {


            double NormU1IF2 = uspeed1[i][IFY[1][i] -1];
            double NormU2IF2 = uspeed2[i][IFY[1][i]];

            //Second interface in both material

            uspeed1[i][IFY[1][i] -1] = NormU1IF2 * cos(0.5 * M_PI - IFrad) + tangentialvelocity[i][IFY[1][i] - 1] * cos(IFrad) ;
            uspeed2[i][IFY[1][i]] = NormU2IF2 * cos(0.5 * M_PI - IFrad) + tangentialvelocity[i][IFY[1][i]] * cos(IFrad);

            vspeed1[i][IFY[1][i] -1]  =  NormU1IF2 * sin(0.5 * M_PI - IFrad) + tangentialvelocity[i][IFY[1][i] -1] * sin(IFrad) ;
            vspeed2[i][IFY[1][i]]  =  NormU2IF2* sin(0.5 * M_PI - IFrad) + tangentialvelocity[i][IFY[1][i]]  * sin(IFrad);
          //  cout.precision(17);

        }
        else{

            if (IFY[0][i]>=0 && IFY[1][i] > 0 ){

                double NormU1IF1 = uspeed1[i][IFY[0][i]];
                double NormU1IF2 = uspeed1[i][IFY[1][i] -1];

                double NormU2IF1 = uspeed2[i][IFY[0][i] -1];
                double NormU2IF2 = uspeed2[i][IFY[1][i]];

                //first interface in both material
                uspeed1[i][IFY[0][i]] = NormU1IF1 * cos(0.5 * M_PI - IFrad) + tangentialvelocity[i][IFY[0][i]] * cos(IFrad) ;
                uspeed2[i][IFY[0][i] -1] = NormU2IF1 * cos(0.5 * M_PI - IFrad) + tangentialvelocity[i][IFY[0][i] -1] * cos(IFrad);

                vspeed1[i][IFY[0][i]] = NormU1IF1  * sin(0.5 * M_PI - IFrad) + tangentialvelocity[i][IFY[0][i]] * sin(IFrad) ;
                vspeed2[i][IFY[0][i] -1] =  NormU2IF1 * sin(0.5 * M_PI - IFrad) + tangentialvelocity[i][IFY[0][i] -1] * sin(IFrad);

                //Second interface in both material

                uspeed1[i][IFY[1][i] -1] = NormU1IF2 * cos(0.5 * M_PI - IFrad) + tangentialvelocity[i][IFY[1][i] - 1] * cos(IFrad) ;
                uspeed2[i][IFY[1][i]] = NormU2IF2 * cos(0.5 * M_PI - IFrad) + tangentialvelocity[i][IFY[1][i]] * cos(IFrad);

                vspeed1[i][IFY[1][i] -1]  =  NormU1IF2 * sin(0.5 * M_PI - IFrad) + tangentialvelocity[i][IFY[1][i] -1] * sin(IFrad) ;
                vspeed2[i][IFY[1][i]]  =  NormU2IF2* sin(0.5 * M_PI - IFrad) + tangentialvelocity[i][IFY[1][i]]  * sin(IFrad);

              //  cout.precision(17);
              //  cout<<"uspeed1[i][IFY[1][i] -1]\t"<<uspeed1[i][IFY[1][i] -1]<<" uspeed2[i][IFY[1][i]]\t "<< uspeed2[i][IFY[1][i]]<<"\n";

            }
            else{
                break;
            }

        }

    }



}




/*----------------------------------------Ghost Fluid Boundary-----------------------------------------*/

//Ghost FLuid Forward to pouplate both forward regions in two materials
void Ghost_Fluid_Boundary_Forward(double dx, vector<double> &P1, vector<double>&P2, vector<double> const &P1_extend, vector<double> const & P2_extend, vector<vector<double>>phi_extend, int N,int row) {

    //initialised both domain

    vector<double>P1_temp_extend(N+2);
    vector<double>P2_temp_extend(N+2);

    //cout<<P_extend[N1+1]<<endl;


    for (int i = 0; i < N + 2; i++) {

        P1_temp_extend[i] = P1_extend[i];
        P2_temp_extend[i] = P2_extend[i];

    }
    //  cout<<"I am working here 2"<<endl;

    for (int j = 0; j < N; j++) {

        //   cout<<"I am working here 3"<<endl;

        if (sngFS(phi_extend[row][j + 1]) > 0.0) {

            if (sngFS(phi_extend[row][j + 2]) != sngFS(phi_extend[row][j])) {

                P2_temp_extend[j + 1] = P2_temp_extend[j + 1];
                //        cout << "P2_temp_extend[j + 1] \t"<< P2_temp_extend[j + 1] << endl;

            } else {

                double P_bar_temp2;

                //Second material Forward

                double P_x2;

                //        cout<<"I am working here 4"<<endl;

                if (((phi_extend[row][j + 2] - phi_extend[row][j + 1]) / dx) >= 0) {

                    //    cout<<"I am working here 1"<<endl;
                    //  if the grad of level set is greater than zero, we  use forward sweep

                    //       cout << "two X values \t" << P2_temp_extend[j + 2] << "\t" << P2_temp_extend[j] << endl;

                    if (P2_temp_extend[j + 2] >= P2_temp_extend[j]) {

                        P_x2 = P2_temp_extend[j];

                    } else {

                        P_x2 = P2_temp_extend[j + 2];

                    }
                    P_bar_temp2 =  P_x2 ;

                    if (P_bar_temp2 < P2_temp_extend[j + 1] && ((-P_bar_temp2  + P2_temp_extend[j + 1]) / (P2_temp_extend[j + 1])) > 0.0001 ) {

                        P2_temp_extend[j + 1] = P_bar_temp2;

                    }
                    //      cout<< "new Value Material 2 \t"<< P2_temp_extend[j + 1]<<endl;

                }
            }
        }
        else {

            //First Material forward

            //  cout<<"phi_extend[row][j + 1]\t"<<phi_extend[row][j + 1]<<endl;
            //   cout<<" P1_temp_extend[j + 1]\t"<< P1_temp_extend[j + 1]<<endl;
            //  cout<<"I am working here 4"<<endl;
            //   cout<< "j+1\t"<<j+1<<endl;
            if (sngFS(phi_extend[row][j + 1]) < 0.0) {

                if (sngFS(phi_extend[row][j + 2]) != sngFS(phi_extend[row][j])) {

                    P1_temp_extend[j + 1] = P1_temp_extend[j + 1];

                    //         cout << "\t P1_temp_extend[j + 1] \t" << P1_temp_extend[j + 1] << "\t" << endl;

                } else {

                    double P_x1;
                    double P_bar_temp1;

                    if (((phi_extend[row][j + 2] - phi_extend[row][j + 1]) / dx) < 0) {

                        //  if the grad of level set is less than zero, we  use forward sweep
                        if (P1_temp_extend[j + 2] >= P1_temp_extend[j]) {

                            P_x1 = P1_temp_extend[j+2];

                        } else {

                            P_x1 = P1_temp_extend[j];

                        }

                        //        cout<<"P_x1\t"<<P_x1<<endl;
                        P_bar_temp1 = P_x1 ;

                        if (P_bar_temp1 > P1_temp_extend[j + 1] &&  ( (P_bar_temp1  - P1_temp_extend[j + 1]) / (P1_temp_extend[j + 1])) > 0.0001) {

                            P1_temp_extend[j + 1] = P_bar_temp1;

                        }
                        //          cout<< "new Value Material 1 \t"<< P1_temp_extend[j + 1]<<endl;

                    }
                }
            }
        }
    }


    for (int i = 0; i < N; i++) {

        P1[i] = P1_temp_extend[i+1];
        P2[i] = P2_temp_extend[i+1];

    }
    //   delete[] P1_temp_extend;
    //  delete[] P2_temp_extend;

}

// Ghost fluid backward to populate two materials's backward regions

void Ghost_Fluid_Boundary_Backward(double dx, vector<double> &P1, vector<double> &P2, vector<double> const &P1_extend, vector<double> const &P2_extend, vector<vector<double>>phi_extend, int N, int row) {


    //initialised both domain
    vector<double>P1_temp_extend(N+2);
    vector<double>P2_temp_extend(N+2);

    //cout<<P_extend[N1+1]<<endl;

    for (int i = 0; i < N + 2; i++) {

        P1_temp_extend[i] = P1_extend[i];
        P2_temp_extend[i] = P2_extend[i];

    }


    for (int j = N + 1; j >1; j--) {

        //First Material

        //    cout<<phi_extend[row][j - 1]<<endl;
        //      cout<<"First Material \t"<<P1_temp_extend[j - 1]<<endl;
        if (phi_extend[row][j - 1] < 0.0) {

            if (sngFS(phi_extend[row][j - 2]) != sngFS(phi_extend[row][j])) {

                P1_temp_extend[j - 1] = P1_temp_extend[j - 1];

                //     cout<<"P1_temp_extend[j - 1]  \t"<<P1_temp_extend[j - 1]<<endl;
                //     cout<<"P2_temp_extend[j - 1]  \t"<<P2_temp_extend[j - 1]<<endl;


            } else {

                double P_x1;
                double P_bar_temp_BW1;

                if (((phi_extend[row][j - 1] - phi_extend[row][j - 2]) / dx) >= 0) {

                    //First Material backward Sweep

                    //        cout<<" P1_temp_extend[j - 2] \t"<<P1_temp_extend[j - 2]<<"\tP1_temp_extend[j]\t"<<P1_temp_extend[j]<<endl;

                    if (P1_temp_extend[j - 2] >= P1_temp_extend[j]) {

                        P_x1 = P1_temp_extend[j-2];

                    } else {

                        P_x1 = P1_temp_extend[j];
                    }


                    P_bar_temp_BW1 =  P_x1;

                    if (P_bar_temp_BW1 > P1_temp_extend[j - 1] && ( (P_bar_temp_BW1  - P1_temp_extend[j - 1]) / (P1_temp_extend[j - 1])) > 0.0001) {

                        P1_temp_extend[j - 1] = P_bar_temp_BW1;

                    }

                    //    cout<<" new Value Material 1 \t "<<P1_temp_extend[j - 1]<<endl;

                }

            }

        }
            //  cout<<" new Value Material 1 \t "<<P1_temp_extend[j - 1]<<endl;
            //Second material
        else{

            if (phi_extend[row][j - 1] > 0.0) {

                //   cout<< j-1<<"\t"<<P2_temp_extend[j - 1]<<"\t"<<endl;

                if (sngFS(phi_extend[row][j - 2]) != sngFS(phi_extend[row][j])) {

                    P2_temp_extend[j - 1] = P2_temp_extend[j - 1];

                } else {

                    double P_x2;
                    double P_bar_temp_BW2;

                    //      cout<< phi_extend[j - 1] <<"\t-\t"<<phi_extend[j - 2]<<endl;

                    if (((phi_extend[row][j - 1] - phi_extend[row][j - 2]) / dx) <= 0) {

                        //    cout<<"P1_temp_extend[j - 2] \t"<<P1_temp_extend[j - 2]<<"\t P1_temp_extend[j]\t"<<P1_temp_extend[j]<<endl;

                        if (P2_temp_extend[j - 2] >= P2_temp_extend[j]) {

                            P_x2 = P2_temp_extend[j];

                        } else {

                            P_x2 = P2_temp_extend[j-2];

                        }


                        P_bar_temp_BW2 = P_x2 ;


                        if (P_bar_temp_BW2 < P2_temp_extend[j - 1] &&  ((-P_bar_temp_BW2  + P2_temp_extend[j - 1]) / (P2_temp_extend[j - 1])) > 0.0001) {

                            P2_temp_extend[j - 1] = P_bar_temp_BW2;

                        }
                        //      cout<<" new Value Material 2"<< P2_temp_extend[j - 1]<<endl;
                    }
                }
            }

        }

        //. cout<<" new Value Material 2 \t"<< P2_temp_extend[j - 1]<<"'\tnew Value Material 1 \t"<< P1_temp_extend[j - 1]<<endl;

    }

    for (int i = 0; i < N; i++) {

        P1[i] = P1_temp_extend[i+1];
        P2[i] = P2_temp_extend[i+1];

//        cout<<"P1 \t"<< P1[i]<<"\tP2\t"<< P2[i] <<endl;

    }

    //  delete[] P1_temp_extend;
    //delete[] P2_temp_extend;

}

double *solution_qudartic(double a, double b, double c){

    double * solution = new double[2];
    double temp = pow(b, 2) - 4 * a * c;
    if (temp >= 0){

        solution[0] = (-b + sqrt(temp))/ (2*a);
        solution[1] = (-b - sqrt(temp))/ (2*a);

    }
    else{

        solution[0] = -10e20;
        solution[1] = -10e20;
        cerr<<"Quadratic equation solving error sqrt of zero"<<endl;

    }
    return solution;

}


void Ghost_Fluid_Boundary_YMaterial1(double dx, vector<vector<double>> &P1, vector<vector<double>> &P1_extend, vector<vector<double>> const &phi_extend, int N) {

    //First material Forward

    for (int j = 0; j < N; j++) {
//        cout<<"I am working here 0"<<endl;

        for (int i = 0; i<N; i++) {

            //       cout<<"I am working here 1"<<endl;
            //   cout<<phi_extend[i+1][j + 1]<<"\t"<<sngFS(phi_extend[i+1][j + 1])<<endl;

            if (sngFS(phi_extend[i + 1][j + 1]) < 0.0) {

                if (sngFS(phi_extend[i + 2][j + 1]) != sngFS(phi_extend[i][j + 1])) {

                    //If the points are immediately next to the interface, the values of the points are kept constant

                    P1_extend[i + 1][j + 1] = P1_extend[i + 1][j + 1];

                    //            cout << "P1_extend[i + 1][j + 1] \t"<< P1_extend[i + 1][j + 1] << endl;

                } else {
                    //      cout<<"I am working here 2"<<endl;
                    double *P_bar_temp1;


                    double P_x1, P_y1;

                    double P_bar;




                    if ((phi_extend[i+2][j+1] - phi_extend[i+1][j+1])/dx < 0){
                        //   cout<<"I am working here 3"<<endl;
                        //  cout<< "i\t" <<i<<endl;



                        /*
                        P_x1 = max(P1_extend[i+2][j+1], P1_extend[i][j+1]);
                        P_y1 = max( P1_extend[i+1][j+2], P1_extend[i+1][j]);

                        // cout<<"P_y1 \t"<<P_y1 <<"\tP_x1\t"<<P_x1;
                        double temp = pow((-2.0*(P_x1 + P_y1)), 2) - 4.0 * 2.0 * (P_y1 * P_y1 + P_x1 * P_x1);

                        if(temp < 0){

                            //   cout<<"illdefined"<<endl;
                            if (P_y1> P_x1){
                                P_bar =P_y1;
                            }
                            else{
                                P_bar =P_x1;
                            }

                            //   cout<<"Old velue\t"<<P1_extend[i + 1][j + 1];
                            if (P_bar > P1_extend[i + 1][j + 1] ){

                                P1_extend[i + 1][j + 1] = P_bar;

                            }
                            //   cout<<"New velue\t"<<P1_extend[i + 1][j + 1];
                        }else{

                            P_bar_temp1 = solution_qudartic(2.0, (-2.0*(P_x1 + P_y1)), (P_y1*P_y1 + P_x1 * P_x1));

                            if (P_bar_temp1[0] ==-10e20){
                                cout<<"i\t"<<i<<"\t j \t"<<j<<endl;

                            }
                            //    cout<<P_bar_temp1[0]<< "\t"<< P_bar_temp1[1]<<"\n";

                            if (P_bar_temp1[0] > P_bar_temp1[1]){
                                P_bar =P_bar_temp1[0];
                            }
                            else{
                                P_bar =  P_bar_temp1[1];

                            }

*/
                        P_y1 = max( P1_extend[i+1][j+2], P1_extend[i+1][j]);

                        P_bar =   P_y1;

                        if (P_bar > P1_extend[i + 1][j + 1]){

                            P1_extend[i + 1][j + 1] = P_bar;

                        }


                        //          cout<<"P1_extend[i + 1][j + 1]\t"<< P1_extend[i + 1][j + 1];

                    }

                }

            }

        }

    }


    for (int j = 0; j < N; j++) {

        for (int i = N + 1; i > 1; i--) {


            if (sngFS(phi_extend[i - 1][j + 1]) < 0) {

                if (sngFS(phi_extend[i - 2][j + 1]) != sngFS(phi_extend[i][j + 1])) {

                    P1_extend[i - 1][j + 1] = P1_extend[i - 1][j + 1];

                } else {

                    double P_x1, P_y1;
                    double *P_bar_temp1;
                    double P_bar;

                    if ((phi_extend[i - 2][j + 1] - phi_extend[i - 1][j + 1]) / dx <= 0) {

                        /*
                        P_y1 = max( P1_extend[i-2][j+1], P1_extend[i][j+1]);

                        P_x1 = max(P1_extend[i-1][j+2], P1_extend[i-1][j]);

                        double temp = pow((-2.0*(P_x1 + P_y1)), 2) - 4.0*2.0*(P_y1 * P_y1 + P_x1 * P_x1);

                        if (temp < 0 ){

                            P_bar = max(P_y1, P_x1);

                            if (P_bar > P1_extend[i - 1][j + 1]){

                                P1_extend[i - 1][j + 1] = P_bar;

                            }

                        }
                        else{

                            P_bar_temp1 = solution_qudartic(2.0,  -2.0*(P_x1 + P_y1), (P_y1*P_y1 + P_x1 * P_x1));

                            P_bar = max(P_bar_temp1[0], P_bar_temp1[1]);

                            if (P_bar > P1_extend[i - 1][j + 1]){

                                P1_extend[i - 1][j + 1] = P_bar;

                            }

                        }
                        */

                        P_y1 = max(P1_extend[i - 2][j + 1], P1_extend[i][j + 1]);
                        P_bar = P_y1;
                        if (P_bar > P1_extend[i - 1][j + 1]) {

                            P1_extend[i - 1][j + 1] = P_bar;


                        }

                    }

                }

            }

        }
    }



    for (int i = 0; i<N; i++){

        for (int j =0; j<N; j++){

            P1[i][j] = P1_extend[i+1][j+1];

        }
    }

}


void Ghost_Fluid_Boundary_YMaterial2(double dx, vector<vector<double>>& P2, vector<vector<double>> &P2_extend, vector<vector<double>> const &phi_extend, int N) {

    //second material Forward

    for (int j = 0; j < N; j++) {

        for (int i = 0; i<N; i++) {


            //   cout<<phi_extend[i+1][j + 1]<<"\t"<<sngFS(phi_extend[i+1][j + 1])<<endl;

            if (sngFS(phi_extend[i + 1][j + 1]) > 0.0) {

                if (sngFS(phi_extend[i + 2][j + 1]) != sngFS(phi_extend[i][j + 1])) {

                    //If the points are immediately next to the interface, the values of the points are kept constant

                    P2_extend[i + 1][j + 1] = P2_extend[i + 1][j + 1];
                    //     cout << "P1_extend[i + 1][j + 1] \t"<< P2_extend[i + 1][j + 1] << endl;

                } else {

                    //    cout<<"I am working here 3"<<endl;
                    double *P_bar_temp2;

                    double P_bar;

                    double P_x2, P_y2;

                    if ((phi_extend[i+2][j+1] - phi_extend[i+1][j+1])/dx >=0){
                        //       cout<<"I am working here 4"<<endl;

                        /*
                        P_y2 = min(P2_extend[i+2][j+1], P2_extend[i][j+1]);

                        P_x2 = min(P2_extend[i+1][j+2], P2_extend[i+1][j]);


                        double temp = pow((-2.0*(P_x2 + P_y2)), 2) - 4.0*2.0*(P_y2 * P_y2 + P_x2 * P_x2);

                        if ( temp < 0){

                            P_bar = min(P_x2, P_y2);

                            if (P_bar < P2_extend[i + 1][j + 1] ){

                                P2_extend[i + 1][j + 1] = P_bar;
                            }


                        }
                        else{

                            P_bar_temp2 = solution_qudartic(2.0,  -2.0*(P_x2 + P_y2), (P_y2 * P_y2 + P_x2 * P_x2));

                            P_bar = min(P_bar_temp2[0], P_bar_temp2[1]);

                            if (P_bar < P2_extend[i + 1][j + 1]){

                                P2_extend[i + 1][j + 1] = P_bar;
                            }

                        }

                        */
                        P_y2 = min(P2_extend[i+2][j+1], P2_extend[i][j+1]);

                        P_bar = P_y2;

                        if (P_bar < P2_extend[i + 1][j + 1]){

                            P2_extend[i + 1][j + 1] = P_bar;
                        }




                    }
                    //    cout<<"I am working here 4"<<endl;

                }

            }

        }

    }


    for (int j = 0; j < N; j++) {

        for (int i = N+1 ; i >1 ; i--) {


            if (sngFS(phi_extend[i-1][j+1])>0){

                //  cout <<"I am working here 2"<<endl;

                if (sngFS(phi_extend[i-2][j+1])!= sngFS(phi_extend[i][j+1])){

                    P2_extend[i-1][j+1] = P2_extend[i-1][j+1];

                }

                else{


                    double P_x2, P_y2;

                    double *P_bar_temp2;

                    double P_bar;

                    if ((phi_extend[i-2][j+1] - phi_extend[i-1][j+1])/dx >= 0){

                        /*
                        P_y2 = min( P2_extend[i-2][j+1], P2_extend[i][j+1]);

                        P_x2 = min(P2_extend[i-1][j+2], P2_extend[i-1][j]);
                        //      cout<<"P_y2 \t"<<P_y2 <<"\tP_x1\t"<<P_x2;

                        //       cout<<"I AM WPKTIN HERE 1"<<endl;
                        double temp = pow(( -2.0*(P_x2 + P_y2)), 2) - 4.0 * 2.0 * (P_y2 * P_y2 + P_x2 * P_x2);

                        if (temp< 0){

                            P_bar = min(P_x2, P_y2);

                            if (P_bar < P2_extend[i - 1][j + 1] ){

                                P2_extend[i - 1][j + 1] = P_bar;

                            }


                        }
                        else{

                            P_bar_temp2 = solution_qudartic(2.0,  -2.0*(P_x2 + P_y2), (P_y2 * P_y2 + P_x2 * P_x2));

                            P_bar = min(P_bar_temp2[0], P_bar_temp2[1]);

                            if (P_bar < P2_extend[i - 1][j + 1]){

                                P2_extend[i - 1][j + 1] = P_bar;

                            }


                        }
                        */
                        //    cout<<"I AM WPKTIN HERE 2"<<endl;

                        P_y2 = min( P2_extend[i-2][j+1], P2_extend[i][j+1]);

                        P_bar = P_y2;

                        if (P_bar < P2_extend[i - 1][j + 1]){

                            P2_extend[i - 1][j + 1] = P_bar;

                        }
                    }

                }

            }

        }
    }

    //  cout<<"I AM WPKTIN HERE 3"<<endl;

//    cout <<"I am workinh here 5"<<endl;
    for (int i = 0; i<N; i++){
        for (int j =0; j<N; j++){

            P2[i][j] = P2_extend[i+1][j+1];

            //   cout<<   P2[i][j]<< "="<< P2_extend[i+1][j+1]<<endl;

        }
    }


}

/*
void Ghost_Fluid_Boundary_YMaterial1(double dx, vector<vector<double>> &P1, vector<vector<double>> &P1_extend, vector<vector<double>> const &phi_extend, int N) {

    //First material Forward

    for (int j = 0; j < N; j++) {
//        cout<<"I am working here 0"<<endl;

        for (int i = 0; i < N; i++) {

            //       cout<<"I am working here 1"<<endl;
            //   cout<<phi_extend[i+1][j + 1]<<"\t"<<sngFS(phi_extend[i+1][j + 1])<<endl;

            if (sngFS(phi_extend[i + 1][j + 1]) < 0.0) {

                if (sngFS(phi_extend[i + 2][j + 1]) != sngFS(phi_extend[i][j + 1])) {

                    //If the points are immediately next to the interface, the values of the points are kept constant

                    P1_extend[i + 1][j + 1] = P1_extend[i + 1][j + 1];

                    //            cout << "P1_extend[i + 1][j + 1] \t"<< P1_extend[i + 1][j + 1] << endl;

                } else {
                    //      cout<<"I am working here 2"<<endl;
                    double *P_bar_temp1;


                    double P_x1, P_y1;

                    double P_bar;


                    if ((phi_extend[i + 2][j + 1] - phi_extend[i + 1][j + 1]) / dx < 0) {
                        //   cout<<"I am working here 3"<<endl;
                        //  cout<< "i\t" <<i<<endl;


                        P_x1 = max(P1_extend[i + 2][j + 1], P1_extend[i][j + 1]);
                        P_y1 = max(P1_extend[i + 1][j + 2], P1_extend[i + 1][j]);

                        // cout<<"P_y1 \t"<<P_y1 <<"\tP_x1\t"<<P_x1;
                        double temp = pow((-2.0 * (P_x1 + P_y1)), 2) - 4.0 * 2.0 * (P_y1 * P_y1 + P_x1 * P_x1);

                        if (temp < 0) {

                            //   cout<<"illdefined"<<endl;

                            P_bar = max(P_y1 , P_x1);

                            //   cout<<"Old velue\t"<<P1_extend[i + 1][j + 1];
                            if (P_bar > P1_extend[i + 1][j + 1]) {

                                P1_extend[i + 1][j + 1] = P_bar;

                            }
                            //   cout<<"New velue\t"<<P1_extend[i + 1][j + 1];
                        } else {

                            P_bar_temp1 = solution_qudartic(2.0, (-2.0 * (P_x1 + P_y1)), (P_y1 * P_y1 + P_x1 * P_x1));

                            P_bar = max(P_bar_temp1[0], P_bar_temp1[1]);

                            if (P_bar > P1_extend[i + 1][j + 1]) {

                                P1_extend[i + 1][j + 1] = P_bar;

                            }

                            //          cout<<"P1_extend[i + 1][j + 1]\t"<< P1_extend[i + 1][j + 1];

                        }

                    }

                }

            }

        }
    }



    for (int j = 0; j < N; j++) {

        for (int i = N + 1; i > 1; i--) {


            if (sngFS(phi_extend[i - 1][j + 1]) < 0) {

                if (sngFS(phi_extend[i - 2][j + 1]) != sngFS(phi_extend[i][j + 1])) {

                    P1_extend[i - 1][j + 1] = P1_extend[i - 1][j + 1];

                } else {

                    double P_x1, P_y1;
                    double *P_bar_temp1;
                    double P_bar;

                    if ((phi_extend[i - 2][j + 1] - phi_extend[i - 1][j + 1]) / dx <= 0) {


                        P_y1 = max( P1_extend[i-2][j+1], P1_extend[i][j+1]);

                        P_x1 = max(P1_extend[i-1][j+2], P1_extend[i-1][j]);

                        double temp = pow((-2.0*(P_x1 + P_y1)), 2) - 4.0*2.0*(P_y1 * P_y1 + P_x1 * P_x1);

                        if (temp < 0 ){

                            P_bar = max(P_y1, P_x1);

                            if (P_bar > P1_extend[i - 1][j + 1]){

                                P1_extend[i - 1][j + 1] = P_bar;

                            }

                        }
                        else{

                            P_bar_temp1 = solution_qudartic(2.0,  -2.0*(P_x1 + P_y1), (P_y1*P_y1 + P_x1 * P_x1));

                            P_bar = max(P_bar_temp1[0], P_bar_temp1[1]);

                            if (P_bar > P1_extend[i - 1][j + 1]){

                                P1_extend[i - 1][j + 1] = P_bar;

                            }

                        }

                    }

                }

            }

        }

    }



    for (int i = 0; i<N; i++){

        for (int j =0; j<N; j++){

            P1[i][j] = P1_extend[i+1][j+1];

        }
    }

}


void Ghost_Fluid_Boundary_YMaterial2(double dx, vector<vector<double>>& P2, vector<vector<double>> &P2_extend, vector<vector<double>> const &phi_extend, int N) {

    //second material Forward

    for (int j = 0; j < N; j++) {

        for (int i = 0; i < N; i++) {


            //   cout<<phi_extend[i+1][j + 1]<<"\t"<<sngFS(phi_extend[i+1][j + 1])<<endl;

            if (sngFS(phi_extend[i + 1][j + 1]) > 0.0) {

                if (sngFS(phi_extend[i + 2][j + 1]) != sngFS(phi_extend[i][j + 1])) {

                    //If the points are immediately next to the interface, the values of the points are kept constant

                    P2_extend[i + 1][j + 1] = P2_extend[i + 1][j + 1];
                    //     cout << "P1_extend[i + 1][j + 1] \t"<< P2_extend[i + 1][j + 1] << endl;

                } else {

                    //    cout<<"I am working here 3"<<endl;
                    double *P_bar_temp2;

                    double P_bar;

                    double P_x2, P_y2;

                    if ((phi_extend[i + 2][j + 1] - phi_extend[i + 1][j + 1]) / dx >= 0) {
                        //       cout<<"I am working here 4"<<endl;


                        P_y2 = min(P2_extend[i + 2][j + 1], P2_extend[i][j + 1]);

                        P_x2 = min(P2_extend[i + 1][j + 2], P2_extend[i + 1][j]);


                        double temp = pow((-2.0 * (P_x2 + P_y2)), 2) - 4.0 * 2.0 * (P_y2 * P_y2 + P_x2 * P_x2);

                        if (temp < 0) {

                            P_bar = min(P_x2, P_y2);

                            if (P_bar < P2_extend[i + 1][j + 1]) {

                                P2_extend[i + 1][j + 1] = P_bar;
                            }


                        } else {

                            P_bar_temp2 = solution_qudartic(2.0, -2.0 * (P_x2 + P_y2), (P_y2 * P_y2 + P_x2 * P_x2));

                            P_bar = min(P_bar_temp2[0], P_bar_temp2[1]);

                            if (P_bar < P2_extend[i + 1][j + 1]) {

                                P2_extend[i + 1][j + 1] = P_bar;
                            }

                        }


                    }
                    //    cout<<"I am working here 4"<<endl;

                }

            }

        }

    }


    for (int j = 0; j < N; j++) {

        for (int i = N + 1; i > 1; i--) {


            if (sngFS(phi_extend[i - 1][j + 1]) > 0) {

                //  cout <<"I am working here 2"<<endl;

                if (sngFS(phi_extend[i - 2][j + 1]) != sngFS(phi_extend[i][j + 1])) {

                    P2_extend[i - 1][j + 1] = P2_extend[i - 1][j + 1];

                } else {


                    double P_x2, P_y2;

                    double *P_bar_temp2;

                    double P_bar;

                    if ((phi_extend[i - 2][j + 1] - phi_extend[i - 1][j + 1]) / dx >= 0) {


                        P_y2 = min(P2_extend[i - 2][j + 1], P2_extend[i][j + 1]);

                        P_x2 = min(P2_extend[i - 1][j + 2], P2_extend[i - 1][j]);
                        //      cout<<"P_y2 \t"<<P_y2 <<"\tP_x1\t"<<P_x2;

                        //       cout<<"I AM WPKTIN HERE 1"<<endl;
                        double temp = pow((-2.0 * (P_x2 + P_y2)), 2) - 4.0 * 2.0 * (P_y2 * P_y2 + P_x2 * P_x2);

                        if (temp < 0) {

                            P_bar = min(P_x2, P_y2);

                            if (P_bar < P2_extend[i - 1][j + 1]) {

                                P2_extend[i - 1][j + 1] = P_bar;

                            }


                        } else {

                            P_bar_temp2 = solution_qudartic(2.0, -2.0 * (P_x2 + P_y2), (P_y2 * P_y2 + P_x2 * P_x2));

                            P_bar = min(P_bar_temp2[0], P_bar_temp2[1]);

                            if (P_bar < P2_extend[i - 1][j + 1]) {

                                P2_extend[i - 1][j + 1] = P_bar;

                            }
                        }
                    }
                }
            }
        }
    }



    for (int i = 0; i<N; i++){
        for (int j =0; j<N; j++){

            P2[i][j] = P2_extend[i+1][j+1];

            //   cout<<   P2[i][j]<< "="<< P2_extend[i+1][j+1]<<endl;

        }
    }


}




void Ghost_Fluid_Boundary_XMaterial1(double dx, vector<vector<double>>&P1, vector<vector<double>>P1_extend, vector<vector<double>>phi_extend, int N) {

    //First material Forward

    for (int j = 0; j < N; j++) {
//        cout<<"I am working here 0"<<endl;

        for (int i = 0; i<N; i++) {

            //       cout<<"I am working here 1"<<endl;
            //   cout<<phi_extend[i+1][j + 1]<<"\t"<<sngFS(phi_extend[i+1][j + 1])<<endl;

            if (sngFS(phi_extend[i + 1][j + 1]) < 0.0) {

                if (sngFS(phi_extend[i + 1][j + 2]) != sngFS(phi_extend[i+1][j ])) {

                    //If the points are immediately next to the interface, the values of the points are kept constant

                    P1_extend[i + 1][j + 1] = P1_extend[i + 1][j + 1];

        //            cout << "P1_extend[i + 1][j + 1] \t"<< P1_extend[i + 1][j + 1] << endl;

                } else {
                    //      cout<<"I am working here 2"<<endl;
                    double *P_bar_temp1;

                    double P_x1, P_y1;

                    double P_bar;

                    if ((phi_extend[i+1][j+2] - phi_extend[i+1][j+1])/dx < 0){
                        //   cout<<"I am working here 3"<<endl;
                      //  cout<< "i\t" <<i<<endl;
                        P_y1 = max(P1_extend[i+2][j+1], P1_extend[i][j+1]);

                        P_x1 = max(P1_extend[i+1][j+2], P1_extend[i+1][j]);

                     // cout<<"P_y1 \t"<<P_y1 <<"\tP_x1\t"<<P_x1;
                        double temp = pow((-2.0*(P_x1 + P_y1)), 2) - 4.0 * 2.0 * (P_y1 * P_y1 + P_x1 * P_x1);

                        if(temp < 0){

                            //   cout<<"illdefined"<<endl;
                            P_bar = max( P_y1,  P_x1);

                         //   cout<<"Old velue\t"<<P1_extend[i + 1][j + 1];
                            if (P_bar > P1_extend[i + 1][j + 1] ){

                                P1_extend[i + 1][j + 1] = P_bar;

                            }
                        //   cout<<"New velue\t"<<P1_extend[i + 1][j + 1];
                        }else{

                            P_bar_temp1 = solution_qudartic(2.0, (-2.0*(P_x1 + P_y1)), (P_y1*P_y1 + P_x1 * P_x1));

                            if (P_bar_temp1[0] ==-10e20){
                                cout<<"i\t"<<i<<"\t j \t"<<j<<endl;

                            }
                        //    cout<<P_bar_temp1[0]<< "\t"<< P_bar_temp1[1]<<"\n";

                            P_bar = max(P_bar_temp1[0], P_bar_temp1[1]);

                            if (P_bar > P1_extend[i + 1][j + 1]){

                                P1_extend[i + 1][j + 1] = P_bar;

                            }


                        }
              //          cout<<"P1_extend[i + 1][j + 1]\t"<< P1_extend[i + 1][j + 1];

                    }

                }

            }

        }

    }


    for (int i = 0; i < N; i++) {

        for (int j = N+1 ; j > 1; j--) {


            if (sngFS(phi_extend[i+1][j-1])<0){

                if (sngFS(phi_extend[i+1][j-2])!= sngFS(phi_extend[i+1][j])){

                    P1_extend[i+1][j-1] = P1_extend[i+1][j-1];

                }
                else{

                    double P_x1, P_y1;
                    double *P_bar_temp1;
                    double P_bar;

                    if ((phi_extend[i+1][j-2] - phi_extend[i+1][j-1])/dx < 0){

                        P_y1 = max( P1_extend[i+2][j-1], P1_extend[i][j-1]);

                        P_x1 = max(P1_extend[i+1][j-2], P1_extend[i+1][j]);

                        double temp = pow((-2.0*(P_x1 + P_y1)), 2) - 4.0*2.0*(P_y1 * P_y1 + P_x1 * P_x1);

                        if (temp < 0 ){

                            P_bar = max(P_y1, P_x1);

                            if (P_bar > P1_extend[i + 1][j - 1]){

                                P1_extend[i + 1][j - 1] = P_bar;

                            }

                        }
                        else{

                            P_bar_temp1 = solution_qudartic(2.0,  -2.0*(P_x1 + P_y1), (P_y1*P_y1 + P_x1 * P_x1));

                            P_bar = max(P_bar_temp1[0], P_bar_temp1[1]);

                            if (P_bar > P1_extend[i + 1][j - 1]){

                                P1_extend[i + 1][j - 1] = P_bar;

                            }

                        }


                    }

                }

            }

        }

    }


    for (int i = 0; i<N; i++){

        for (int j =0; j<N; j++){

            P1[i][j] = P1_extend[i+1][j+1];

        }
    }

}


void Ghost_Fluid_Boundary_XMaterial2(double dx, vector<vector<double>>& P2, vector<vector<double>> P2_extend, vector<vector<double>> const phi_extend, int N) {

    //second material Forward

    for (int j = 0; j < N; j++) {

        for (int i = 0; i<N; i++) {


            //   cout<<phi_extend[i+1][j + 1]<<"\t"<<sngFS(phi_extend[i+1][j + 1])<<endl;

            if (sngFS(phi_extend[i + 1][j + 1]) > 0.0) {

                if (sngFS(phi_extend[i + 1][j + 2]) != sngFS(phi_extend[i + 1][j ])) {

                    //If the points are immediately next to the interface, the values of the points are kept constant

                    P2_extend[i + 1][j + 1] = P2_extend[i + 1][j + 1];
               //     cout << "P1_extend[i + 1][j + 1] \t"<< P2_extend[i + 1][j + 1] << endl;

                } else {

                //    cout<<"I am working here 3"<<endl;
                    double *P_bar_temp2;

                    double P_bar;

                    double P_x2, P_y2;

                    if ((phi_extend[i+1][j+2] - phi_extend[i+1][j+1])/dx >=0){
                 //       cout<<"I am working here 4"<<endl;

                        P_y2 = min(P2_extend[i+2][j+1], P2_extend[i][j+1]);

                        P_x2 = min(P2_extend[i+1][j+2], P2_extend[i+1][j]);


                        double temp = pow((-2.0*(P_x2 + P_y2)), 2) - 4.0*2.0*(P_y2 * P_y2 + P_x2 * P_x2);

                        if ( temp < 0){

                            P_bar = min(P_x2, P_y2);

                            if (P_bar < P2_extend[i + 1][j + 1] ){

                                P2_extend[i + 1][j + 1] = P_bar;
                            }


                        }
                        else{

                            P_bar_temp2 = solution_qudartic(2.0,  -2.0*(P_x2 + P_y2), (P_y2 * P_y2 + P_x2 * P_x2));

                            P_bar = min(P_bar_temp2[0], P_bar_temp2[1]);

                            if (P_bar < P2_extend[i + 1][j + 1]){

                                P2_extend[i + 1][j + 1] = P_bar;
                            }

                        }


                    }
                //    cout<<"I am working here 4"<<endl;

                }

            }

        }

    }


    for (int i = 0; i < N; i++) {

        for (int j = N+1 ; j >1 ; j--) {


            if (sngFS(phi_extend[i+1][j-1])>0){

                 //  cout <<"I am working here 2"<<endl;

                if (sngFS(phi_extend[i+1][j-2])!= sngFS(phi_extend[i+1][j])){

                    P2_extend[i+1][j-1] = P2_extend[i+1][j-1];

                }

                else{


                    double P_x2, P_y2;

                    double *P_bar_temp2;

                    double P_bar;

                    if ((phi_extend[i+1][j-2] - phi_extend[i+1][j-1])/dx > 0){

                        P_y2 = min( P2_extend[i+2][j-1], P2_extend[i][j-1]);

                        P_x2 = min(P2_extend[i+1][j-2], P2_extend[i+1][j]);
                  //      cout<<"P_y2 \t"<<P_y2 <<"\tP_x1\t"<<P_x2;

                //       cout<<"I AM WPKTIN HERE 1"<<endl;
                        double temp = pow(( -2.0*(P_x2 + P_y2)), 2) - 4.0 * 2.0 * (P_y2 * P_y2 + P_x2 * P_x2);

                        if (temp< 0){

                            P_bar = min(P_x2, P_y2);

                            if (P_bar < P2_extend[i + 1][j - 1] ){

                                P2_extend[i + 1][j - 1] = P_bar;

                            }


                        }
                        else{

                            P_bar_temp2 = solution_qudartic(2.0,  -2.0*(P_x2 + P_y2), (P_y2 * P_y2 + P_x2 * P_x2));

                            P_bar = min(P_bar_temp2[0], P_bar_temp2[1]);

                            if (P_bar < P2_extend[i + 1][j - 1] && (( -P_bar + P2_extend[i + 1][j - 1])/ P2_extend[i + 1][j - 1] > 0.00001)){

                                P2_extend[i + 1][j - 1] = P_bar;

                            }


                        }
                    //    cout<<"I AM WPKTIN HERE 2"<<endl;


                    }

                }

            }

        }
    }

  //  cout<<"I AM WPKTIN HERE 3"<<endl;

//    cout <<"I am workinh here 5"<<endl;
    for (int i = 0; i<N; i++){
        for (int j =0; j<N; j++){

            P2[i][j] = P2_extend[i+1][j+1];

             //   cout<<   P2[i][j]<< "="<< P2_extend[i+1][j+1]<<endl;

        }
    }


}




*/

void Primitive_to_Conservative (int N, double gamA, double gamB, vector<vector<double>>const density1, vector<vector<double>> const density2, vector<vector<double>> const uspeed1, vector<vector<double>> const uspeed2, vector<vector<double>> const vspeed1, vector<vector<double>> const vspeed2, vector<vector<double>> const pressure1, vector<vector<double>> const pressure2, vector<vector<double>> &U1D,  vector<vector<double>> &U1MU, vector<vector<double>> &U1MV, vector<vector<double>> &U1E, vector<vector<double>> &U2D, vector<vector<double>> &U2MU, vector<vector<double>> &U2MV, vector<vector<double>> &U2E){


    for (int i =0 ; i<N; i++){

        for (int j =0; j<N; j++){

            U1D[i][j] = density1[i][j];
            U1MU[i][j] = density1[i][j] * uspeed1[i][j];
            U1MV[i][j] = density1[i][j] * vspeed1[i][j];
            U1E[i][j] = (pressure1[i][j] / ( gamA - 1.0))+(((uspeed1[i][j] * uspeed1[i][j]) + (vspeed1[i][j] * vspeed1[i][j]))* density1[i][j])/2.0;
            U2D[i][j] = density2[i][j];
            U2MU[i][j] = density2[i][j] * uspeed2[i][j];
            U2MV[i][j] = density2[i][j] * vspeed2[i][j];
            U2E[i][j] = (pressure2[i][j] / ( gamB - 1.0))+(((uspeed2[i][j] * uspeed2[i][j]) + (vspeed2[i][j] * vspeed2[i][j])) * density2[i][j])/2.0;

            // cout<< U2E[i][j] <<"\t=\t"<<"(\t"<<pressure2[i][j]<< "\t/\t" <<"\t(\t"<< gamB <<"\t-\t"<< "1.0"<<"\t)\t"<<"\t)\t"<<"\t+\t"<<"\t(\t"<<uspeed2[i][j]<< "\t*\t"<< uspeed2[i][j] <<"\t*\t"<< density2[i][j]<<"\t)\t"<<"\t/\t"<<"2.0"<<endl;

        }
    }


}


/*-------------Compute Ghost Fluid domain Boundary for second order Approximate solver ------------------*/

vector<vector<double>> GF_Domain_boundaryConditions (int N, vector<vector<double>>const &u) {


    vector<vector<double>>U_extend(N+4, vector<double> (N+4, 0));


    for (int i = 0; i < N; i++) {

        for (int j =0 ; j < N; j++){

            U_extend[i + 2][j + 2] = u[i][j];

        }

    }


    //fill the first two and last two rolls
    for (int j = 0 ; j<N; j++){

        U_extend[0][j+2] = u[0][j];
        U_extend[1][j+2] = u[0][j];
        U_extend[N+2][j+2] = u[N-1][j];
        U_extend[N+3][j+2] = u[N-1][j];

    }

    //fill the first two and last two columns
    for (int i = 0; i<N; i++){

        U_extend[i+2][0] = u[i][0];
        U_extend[i+2][1] = u[i][0];
        U_extend[i+2][N+2] = u[i][N-1];
        U_extend[i+2][N+3] = u[i][N-1];

    }

    //The Copying is done in a way thay the four cell values at each conor are not alighed to any values, but zeros


    //Asiign values to the four corners of the extended 2D array

    //Bottom left
    U_extend[0][0] = U_extend[2][2];
    U_extend[0][1] = U_extend[2][2];
    U_extend[1][1] = U_extend[2][2];
    U_extend[1][0] = U_extend[2][2];
    //Top right
    U_extend[N+3][N+3] = U_extend[N+1][N+1];
    U_extend[N+3][N+2] = U_extend[N+1][N+1];
    U_extend[N+2][N+2] = U_extend[N+1][N+1];
    U_extend[N+2][N+3] = U_extend[N+1][N+1];
    //Bottom Right
    U_extend[0][N+3] = U_extend[3][N+1];
    U_extend[0][N+2] = U_extend[3][N+1];
    U_extend[1][N+2] = U_extend[3][N+1];
    U_extend[1][N+3] = U_extend[3][N+1];
    //TopLeft
    U_extend[N+3][0] = U_extend[N+1][2];
    U_extend[N+3][1] = U_extend[N+1][2];
    U_extend[N+2][0] = U_extend[N+1][2];
    U_extend[N+2][1] = U_extend[N+1][2];


    return U_extend;

}


double sound_speed (double pressure, double rho, double gam) {


    return sqrt(gam * (pressure / rho));


}


double dt_func( double dx, vector<double> const &speed_extend, vector<double> const &pressure_extend, vector<double> const &density_extend, double const gam){

    double a;
    double dt;

    double cmax = 0, ctemp;

    //Finding max wave speed

    /* ??????? i am confused if it should be N+3 rather than N+2 ????????*/
    for (int i=0; i<N+4;i++){

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

double dt_f(double dx, double dy, vector<vector<double>>const &uspeed1_extend2, vector<vector<double>>const &vspeed1_extend2,  vector<vector<double>>const &pressure1_extend2, vector<vector<double>>const &density1_extend2, double gam){

    double dtNx = 10.0;
    double dtNy = 10.0;

    //
    //FInd the min dt in x direction
    //

    for(int i =0 ; i<N+4; i++){

        double dt1_tempx;
        dt1_tempx =dt_func(dx, uspeed1_extend2[i], pressure1_extend2[i], density1_extend2[i], gam);
        if (dt1_tempx < dtNx){

            dtNx = dt1_tempx;

        }

    }

    //
    //Find the min dt in y direction
    //
    for(int i = 0 ; i<N; i++) {

        double dt1_tempy;

        vector<double>vspeed1y_extend(N+4);
        vector<double>pressure1y_extend (N+4);
        vector<double>density1y_extend (N+4);

        for (int j = 0; j < N + 4; j++) {

            vspeed1y_extend[j] = vspeed1_extend2[j][i + 2];
            pressure1y_extend[j] = pressure1_extend2[j][i + 2];
            density1y_extend[j] = density1_extend2[j][i + 2];

        }

        dt1_tempy =dt_func(dy, vspeed1y_extend, pressure1y_extend, density1y_extend, gam);
        if (dt1_tempy < dtNy){

            dtNy = dt1_tempy;

        }

    }

    cout<<"dtNx\t"<<dtNx<<"dtNy\t"<<dtNy<<endl;
    //
    //Compare the two dt and find the minimium
    //

    if ( dtNx < dtNy){
        return dtNx;
    }
    else{
        return  dtNy;
    }



}





/*----------------------------------------------------------------------------------------------------------------*/
/*                                                                                                                */
/*                                 Approximate Solver HLLC with MUSCL-HANCOCK                                     */
/*                                                                                                                */
/*                                                                                                                */
/*                                                                                                                */
/*----------------------------------------------------------------------------------------------------------------*/

/* --------------------------------- solvers ------------------------------------ */



vector<double>CtoP(vector<double> US, double gam){


    vector<double>w(4);

    w[0] = US[0];       // rho
    w[1] = US[1] / US[0];           //velocity in x direction
    w[2] = US[2] / US[0];           //velocity in y direction
    w[3] = (US[3] - 0.5 * US[0] * (pow(( US[1] / US[0]),2) + pow(( US[2] / US[0]),2))) * (gam -1);          //pressure


    return w;


}

//work out the slope_vector

vector<double>U_to_F(vector<double> const& U_value, double gam, bool xsweep) {

    vector<double>flux(4);
    if (xsweep){

        flux[0] = U_value[1];
        flux[1] = U_value[1]*(U_value[1]/U_value[0]) + (gam-1)*(U_value[3] - 0.5* U_value[0] * (pow((U_value[1]/ U_value[0]),2.0) + pow((U_value[2]/ U_value[0]),2.0)));
        flux[2] = U_value[1]* (U_value[2] / U_value[0]);
        flux[3] = (U_value[1]/U_value[0]) * (U_value[3] + (gam-1)*( U_value[3] - 0.5 * U_value[0] * (pow(( U_value[1]/U_value[0]),2.0) + pow((U_value[2]/ U_value[0]),2.0))));

    }
    else{

        flux[0] = U_value[2];
        flux[1] = U_value[2]* (U_value[1] / U_value[0]);
        flux[2] = U_value[2]*(U_value[2]/U_value[0]) + (gam-1)*(U_value[3] - 0.5* U_value[0] * (pow((U_value[1]/ U_value[0]),2.0) + pow((U_value[2]/ U_value[0]),2.0)));
        flux[3] = (U_value[2]/U_value[0]) * (U_value[3] + (gam-1)*( U_value[3] - 0.5 * U_value[0] * (pow(( U_value[2]/U_value[0]),2.0) + pow(( U_value[1]/U_value[0]),2.0))));


    }


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





vector<vector<double>>Muscl_Hancock( double dx, double dt, vector<vector<double>> const &U_list, double gam, bool xsweep){

    // U_list is a three by three matrix that has the conservetive values for three ajucent points
    // calculate the slope_vector

    double di, di1, di2, Sdi;


    //store extrapolated values

    vector<vector<double>>U_update(4, vector<double> (2, 0));


    // Store Ubar values



    vector<vector<double>>U_MH (4, vector<double> (2, 0));


    double r;
    double SlopeLR, SlopeLimiter;
    double beta =1;
    double TOL = 1E-10;
    for (int i =0 ; i<4 ; i++){

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

    vector<double>U_L(4);

    U_L[0] = U_update[0][0] ;
    U_L[1] = U_update[1][0] ;
    U_L[2] = U_update[2][0] ;
    U_L[3] = U_update[3][0] ;


    vector<double>U_R(4);

    U_R[0] = U_update[0][1] ;
    U_R[1] = U_update[1][1] ;
    U_R[2] = U_update[2][1] ;
    U_R[3] = U_update[3][1] ;


    vector<double>flux_L;
    vector<double>flux_R;


    flux_L = U_to_F( U_L, gam, xsweep);
    flux_R = U_to_F( U_R, gam, xsweep);

    //   cout<< flux_R[2]<<endl;
    //NOTE assumptions Ur state vector to flux conversion is accurate

    for (int i = 0 ; i< 4; i++){

        //Ubar left
        U_MH[i][0] =  U_L[i] + 0.5 * ( dt/dx ) *( flux_L[i]- flux_R[i]);

        //    cout<<  U_L[2]<<endl;
        //      cout<<U_MH[i][0]<<"="<< U_L[i]<< "\t+\t" <<0.5 <<"\t*\t" <<"\t(\t"<< dt/dx<<" )" <<"*"<<"(" <<flux_L[i]<<"- "<<flux_R[i]<<")"<<"\n";
        //Ubar right
        U_MH[i][1] =  U_R[i] + 0.5 * ( dt/dx ) *( flux_L[i]- flux_R[i]);


    }
    return U_MH;
}

/* ---------------------------pressure based speed estimate --------------------------*/

vector<double>pressure_based_wave_speed(double p_L, double p_R, double u_L, double u_R, double rho_L,double rho_R, double gam){


    vector<double>wavespeed(2);

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


vector<double>U_stark(double s_k, double s_star, double u_k,  double p_k, double rho_k, double e_k , double v_k, bool xsweep){


    vector<double>U_star(4);
    if(xsweep){

        double ft= rho_k * ( (s_k - u_k)/(s_k - s_star));

        U_star[0]= ft *1.0;

        U_star[1]= ft * s_star;

        U_star[2]= ft * v_k;

        U_star[3]= ft * ((e_k /rho_k) + (s_star - u_k)*(s_star + (p_k/(rho_k * (s_k - u_k)))));

    }
    else{

        double ft= rho_k * ( (s_k - v_k)/(s_k - s_star));

        U_star[0]= ft *1.0;

        U_star[1]= ft * u_k;

        U_star[2]= ft * s_star;

        U_star[3]= ft * ((e_k /rho_k) + (s_star - v_k)*(s_star + (p_k/(rho_k * (s_k - v_k)))));

    }



    return U_star;
}


bool AreSame(double a, double b)
{
    return fabs(a - b) < numeric_limits<double>::epsilon();
}

bool doubletest(double a, double b){

    if (a == b){

        return true;
    }
    else{
        //   cout<<"a\t"<<a<<"\tb\t"<<b <<"\t";
        return false;
    }
}

vector<vector<double>>U_update(double dx, double dt, vector<double> const & UND_extend, vector<double> const & UNMU_extend, vector<double> const & UNMV_extend, vector<double> const & UNE_extend, double gam, bool xsweep) {

    double S_L, S_R;
    double a_L, a_R;        // If you want to use Roe average Wavespeed then a_L and a_R is there for you
    double S_star;


    vector<double>wavespeed;
    vector<vector<double>>U_new(4, vector<double> (N, 0));

    vector<vector<double>>f_new(4, vector<double> (N+1, 0));

    vector<double> Ustar_L;
    vector<double> Ustar_R;

    vector<vector<double>>U_temp(4, vector<double> (3, 0));

    vector<vector<double>>U_temp2(4, vector<double> (3, 0));

    vector<vector<double>> U_MH1;
    vector<vector<double>> U_MH2;
    vector<double >U_MH_L(4);
    vector<double >U_MH_R(4);

    vector<vector<double>>U_extend(4, vector<double> (N+4, 0));



    for (int i = 0; i < N + 4; i++) {

        U_extend[0][i] = UND_extend[i];
        U_extend[1][i] = UNMU_extend[i];
        U_extend[2][i] = UNMV_extend[i];
        U_extend[3][i] = UNE_extend[i];

    }


    for (int i = 1; i < N + 2; i++) {

        for (int j = 0; j < 4; j++) {


            U_temp[j][0] = U_extend[j][i - 1];

            U_temp[j][1] = U_extend[j][i];

            U_temp[j][2] = U_extend[j][i + 1];


        }


        U_MH1 = Muscl_Hancock( dx, dt, U_temp, gam, xsweep);


        for (int j = 0; j < 4; j++) {

            U_temp2[j][0] = U_extend[j][i];

            U_temp2[j][1] = U_extend[j][i + 1];

            U_temp2[j][2] = U_extend[j][i + 2];

        }


        U_MH2 = Muscl_Hancock( dx, dt, U_temp2, gam, xsweep);


        for (int k = 0; k < 4; k++) {

            U_MH_L[k] = U_MH1[k][1];
            //      cout<<U_MH1[2][1]<<endl;
            U_MH_R[k] = U_MH2[k][0];


        }

        vector<double>WL;
        vector<double>WR;

        WL = CtoP(U_MH_L, gam);
        WR = CtoP(U_MH_R, gam);


        double pressureL, pressureR, densityL, densityR, speedL, speedR, vspeedL, vspeedR, energyL, energyR;


        if (xsweep) {


            densityL = WL[0];
            densityR = WR[0];
            speedL = WL[1];
            speedR = WR[1];
            vspeedL = WL[2];
            vspeedR = WR[2];
            pressureL = WL[3];
            pressureR = WR[3];
            energyL = (WL[3]) / (gam - 1) + 0.5 * (WL[1] * WL[1] +  WL[2] * WL[2])* WL[0];
            energyR = (WR[3]) / (gam - 1) + 0.5 * (WR[1] * WR[1] +  WR[2] * WR[2])* WR[0];

            wavespeed = pressure_based_wave_speed(pressureL, pressureR, speedL, speedR, densityL, densityR, gam);

            S_L = wavespeed[0];
            S_R = wavespeed[1];

            double num = pressureR - pressureL + densityL * speedL * (S_L - speedL) - densityR * speedR * (S_R - speedR);

            double dem = densityL * (S_L - speedL) - densityR * (S_R - speedR);

            S_star = num / dem;

            //      cout<<" \tS_R\t"<<S_R<<" \tS_L\t"<<S_L<<"\tpressureL\t"<<pressureL << " \t pressureR \t"<<pressureR << "\t vspeedL \t"<< speedL<< " \tvspeedR\t"<<speedR<<"\tdensityL \t" <<densityL<<"\tDensityR\t"<<densityR<<endl;


            //   cout<< " U_MH_L[0] \t"<< U_MH_L[0]<< " U_MH_L[1] \t"<< U_MH_L[1]<<"\t  U_MH_L[2] \t"<< U_MH_L[2]<<"\t  U_MH_L[3] \t"<< U_MH_L[3]<<endl;

            if (S_L >= 0) {

                f_new[0][i - 1] = U_MH_L[1];
                f_new[1][i - 1] = U_MH_L[1] * (U_MH_L[1] / U_MH_L[0]) +
                                  (gam - 1) * (U_MH_L[3] - 0.5 * U_MH_L[0] * (pow((U_MH_L[1] / U_MH_L[0]), 2.0) +  pow((U_MH_L[2] / U_MH_L[0]), 2.0)));
                f_new[2][i - 1] = U_MH_L[1] * U_MH_L[2] / U_MH_L[0];
                //       cout<< "  f_new[1][i - 1]\t"<< f_new[1][i - 1]<<"\t"<< f_new[2][i - 1]<<endl;
                f_new[3][i - 1] = (U_MH_L[1] / U_MH_L[0]) * (U_MH_L[3] + (gam - 1) * (U_MH_L[3] - 0.5 * U_MH_L[0] * (pow((U_MH_L[1] /U_MH_L[0]),2.0) + pow((U_MH_L[2] / U_MH_L[0]), 2.0))));

            }

            if (S_L < 0 && S_star >= 0) {

                //    cout<< "second Condition"<<"\t";
                Ustar_L = U_stark(S_L, S_star, speedL, pressureL, densityL, energyL, vspeedL, xsweep);

                //    cout<< "Ustar"<<Ustar_L[2]<<endl;

                f_new[0][i - 1] = U_MH_L[1] + S_L * (Ustar_L[0] - U_MH_L[0]);
                f_new[1][i - 1] = (U_MH_L[1] * (U_MH_L[1] / U_MH_L[0]) + (gam - 1) * (U_MH_L[3] - 0.5 * U_MH_L[0] * (pow((U_MH_L[1] / U_MH_L[0]), 2.0) +  pow((U_MH_L[2] / U_MH_L[0]), 2.0)))) + S_L * (Ustar_L[1] - U_MH_L[1]);

                //  cout<< f_new[1][i - 1]<< "\t=\t" <<"(\t"<<U_MH_L[0]<<"\t"<<U_MH_L[1]<<"\t"<<U_MH_L[2]<<"\t"<<U_MH_L[3]<<"\t"<<S_L<< "\t*\t"<<"\t(\t"<< Ustar_L[1]<<"\t-\t"<< U_MH_L[1];

                f_new[2][i - 1] = U_MH_L[1] * U_MH_L[2] / U_MH_L[0] + S_L * (Ustar_L[2] - U_MH_L[2]);
                //          cout<< "  f_new[1][i - 1]\t"<< f_new[1][i - 1]<<"\t"<< f_new[2][i - 1]<<endl;
                f_new[3][i - 1] = ((U_MH_L[1] / U_MH_L[0]) * (U_MH_L[3] + (gam - 1) * (U_MH_L[3] - 0.5 * U_MH_L[0] *  (pow((U_MH_L[1] / U_MH_L[0]), 2.0) +  pow((U_MH_L[2] / U_MH_L[0]), 2.0))))) + S_L * (Ustar_L[3] - U_MH_L[3]);
                //   cout<< "Ustar"<<f_new[1][i]<<endl;
                //   delete[] Ustar_L;
            }

            if (S_R >= 0 && S_star < 0) {

                //        cout<< "third Condition"<<"\t";
                Ustar_R = U_stark(S_R, S_star, speedR, pressureR, densityR, energyR, vspeedR, xsweep);

                f_new[0][i - 1] = U_MH_R[1] + S_R * (Ustar_R[0] - U_MH_R[0]);
                f_new[1][i - 1] = (U_MH_R[1] * (U_MH_R[1] / U_MH_R[0]) + (gam - 1) * (U_MH_R[3] - 0.5 * U_MH_R[0] *  (pow((U_MH_R[1] / U_MH_R[0]), 2.0) +  pow((U_MH_R[2] / U_MH_R[0]), 2.0)))) + S_R * (Ustar_R[1] - U_MH_R[1]);
                f_new[2][i - 1] = U_MH_R[1] * U_MH_R[2] / U_MH_R[0] + S_R * (Ustar_R[2] - U_MH_R[2]);
                //      cout<< "  f_new[1][i - 1]\t"<< f_new[1][i - 1]<<"\t"<< f_new[2][i - 1]<<endl;
                f_new[3][i - 1] = ((U_MH_R[1] / U_MH_R[0]) * (U_MH_R[3] + (gam - 1) * (U_MH_R[3] - 0.5 * U_MH_R[0] *  (pow((U_MH_R[1] / U_MH_R[0]), 2.0) +  pow((U_MH_R[2] / U_MH_R[0]), 2.0))))) +S_R * (Ustar_R[3] - U_MH_R[3]);
                //delete[] Ustar_R;
            }

            if (S_R <  0) {

                //         cout<< "fourth Condition"<<"\t";
                f_new[0][i - 1] = U_MH_R[1];
                f_new[1][i - 1] = U_MH_R[1] * (U_MH_R[1] / U_MH_R[0]) +
                                  (gam - 1) * (U_MH_R[3] - 0.5 * U_MH_R[0] *  (pow((U_MH_R[1] / U_MH_R[0]), 2.0) +  pow((U_MH_R[2] / U_MH_R[0]), 2.0)));
                f_new[2][i - 1] = U_MH_R[1] * U_MH_R[2] / U_MH_R[0];
                //         cout<< "  f_new[1][i - 1]\t"<< f_new[1][i - 1]<<"\t"<< f_new[2][i - 1]<<endl;
                f_new[3][i - 1] = (U_MH_R[1] / U_MH_R[0]) * (U_MH_R[3] + (gam - 1) * (U_MH_R[3] - 0.5 * U_MH_R[0] *  (pow((U_MH_R[1] / U_MH_R[0]), 2.0) +  pow((U_MH_R[2] / U_MH_R[0]), 2.0))));

            }
        } else {

            //
            //If the sweep is in y direction
            //

            densityL = WL[0];
            densityR = WR[0];
            speedL = WL[1];
            speedR = WR[1];
            vspeedL = WL[2];
            vspeedR = WR[2];
            pressureL = WL[3];
            pressureR = WR[3];
            energyL = (WL[3]) / (gam - 1) + 0.5 * (WL[2] * WL[2] + WL[1] * WL[1])* WL[0];
            energyR = (WR[3]) / (gam - 1) + 0.5 * (WR[2] * WR[2] + WR[1] * WR[1]) * WR[0];

            wavespeed = pressure_based_wave_speed(pressureL, pressureR, vspeedL, vspeedR, densityL, densityR, gam);

            S_L = wavespeed[0];

            S_R = wavespeed[1];

            //         cout<<" \tS_R\t"<<S_R<<" \tS_L\t"<<S_L<<"\tpressureL\t"<<pressureL << " \t pressureR \t"<<pressureR << "\t vspeedL \t"<< vspeedL<< " \tvspeedR\t"<<vspeedR<<"\tdensityL \t" <<densityL<<"\tDensityR\t"<<densityR<<endl;


            double num =
                    pressureR - pressureL + densityL * vspeedL * (S_L - vspeedL) - densityR * vspeedR * (S_R - vspeedR);

            double dem = densityL * (S_L - vspeedL) - densityR * (S_R - vspeedR);


            S_star = num / dem;


            if (S_L >= 0) {

                f_new[0][i - 1] = U_MH_L[2];
                f_new[1][i - 1] = U_MH_L[2] * U_MH_L[1] / U_MH_L[0];
                f_new[2][i - 1] = U_MH_L[2] * (U_MH_L[2] / U_MH_L[0]) +
                                  (gam - 1) * (U_MH_L[3] - 0.5 * U_MH_L[0] * (pow((U_MH_L[2] / U_MH_L[0]), 2.0) + pow((U_MH_L[1] / U_MH_L[0]), 2.0)));
                f_new[3][i - 1] = (U_MH_L[2] / U_MH_L[0]) * (U_MH_L[3] + (gam - 1) * (U_MH_L[3] - 0.5 * U_MH_L[0] * (pow((U_MH_L[2] / U_MH_L[0]), 2.0) + pow((U_MH_L[1] / U_MH_L[0]), 2.0))));

            }

            if (S_L <= 0 && S_star >= 0) {
                //    cout<< "second Condition"<<"\t";
                Ustar_L = U_stark(S_L, S_star, speedL, pressureL, densityL, energyL, vspeedL, xsweep);

                //     cout<< "Ustar_L"<<Ustar_L[2]<<endl;

                f_new[0][i - 1] = U_MH_L[2] + S_L * (Ustar_L[0] - U_MH_L[0]);
                f_new[1][i - 1] = U_MH_L[2] * U_MH_L[1] / U_MH_L[0] + S_L * (Ustar_L[1] - U_MH_L[1]);
                //  cout<<f_new[1][i -1]<< "\t=\t" << U_MH_L[2] <<"\t*\t"<< U_MH_L[1]<< "\t/\t"<< U_MH_L[0] <<"\t+\t"<< S_L<< "\t*\t"<< "\t(\t"<<Ustar_L[1] <<"\t-\t"<< U_MH_L[1]<<"\t)\t"<<endl;
                f_new[2][i - 1] = (U_MH_L[2] * (U_MH_L[2] / U_MH_L[0]) + (gam - 1) * (U_MH_L[3] - 0.5 * U_MH_L[0] * (pow((U_MH_L[2] / U_MH_L[0]), 2.0) + pow((U_MH_L[1] / U_MH_L[0]), 2.0)))) + S_L * (Ustar_L[2] - U_MH_L[2]);
                f_new[3][i - 1] = ((U_MH_L[2] / U_MH_L[0]) * (U_MH_L[3] + (gam - 1) * (U_MH_L[3] - 0.5 * U_MH_L[0] * (pow((U_MH_L[2] / U_MH_L[0]), 2.0) + pow((U_MH_L[1] / U_MH_L[0]), 2.0))))) + S_L * (Ustar_L[3] - U_MH_L[3]);
                //   cout<< "Ustar"<<f_new[1][i]<<endl;

                //    delete[] Ustar_L;
            }

            if (S_R >= 0 && S_star <= 0) {

                //        cout<< "third Condition"<<"\t";
                Ustar_R = U_stark(S_R, S_star, speedR, pressureR, densityR, energyR, vspeedR, xsweep);

                f_new[0][i - 1] = U_MH_R[2] + S_R * (Ustar_R[0] - U_MH_R[0]);
                f_new[1][i - 1] = U_MH_R[2] * U_MH_R[1] / U_MH_R[0] + S_R * (Ustar_R[1] - U_MH_R[1]);

                f_new[2][i - 1] = (U_MH_R[2] * (U_MH_R[2] / U_MH_R[0]) + (gam - 1) * (U_MH_R[3] - 0.5 * U_MH_R[0] * (pow((U_MH_R[2] / U_MH_R[0]), 2.0) + pow((U_MH_R[1] / U_MH_R[0]), 2.0)))) + S_R * (Ustar_R[2] - U_MH_R[2]);
                f_new[3][i - 1] = ((U_MH_R[2] / U_MH_R[0]) * (U_MH_R[3] + (gam - 1) * (U_MH_R[3] - 0.5 * U_MH_R[0] *(pow((U_MH_R[2] / U_MH_R[0]), 2.0) + pow((U_MH_R[1] / U_MH_R[0]), 2.0))))) + S_R * (Ustar_R[3] - U_MH_R[3]);

                //    delete[] Ustar_R;
            }

            if (S_R <= 0) {
                //         cout<< "fourth Condition"<<"\t";
                f_new[0][i - 1] = U_MH_R[2];
                f_new[1][i - 1] = U_MH_R[2] * U_MH_R[1] / U_MH_R[0];
                //        cout<< "  f_new[1][i - 1]"<< f_new[1][i - 1]<<endl;
                f_new[2][i - 1] = U_MH_R[2] * (U_MH_R[2] / U_MH_R[0]) + (gam - 1) * (U_MH_R[3] - 0.5 * U_MH_R[0] * (pow((U_MH_R[2] / U_MH_R[0]), 2.0) + pow((U_MH_R[1] / U_MH_R[0]), 2.0)));
                f_new[3][i - 1] = (U_MH_R[2] / U_MH_R[0]) * (U_MH_R[3] + (gam - 1) * (U_MH_R[3] - 0.5 * U_MH_R[0] * (pow((U_MH_R[2] / U_MH_R[0]), 2.0) + pow((U_MH_R[1] / U_MH_R[0]), 2.0))));

            }

        }

        //  delete[] WL;
        //  delete[] WR;


    }


    //freeMemory(U_temp, 4);
    // freeMemory(U_temp2,4);
    // delete[] U_MH_L;
    //delete[] U_MH_R;

    //cout<<U_new[1][0] <<"="<<U_extend[1][2] <<"-"<< dt/dx<<"x"<<f_new[1][1]<<"-"<<f_new[1][0]<<"\n";

    //cout<<U_new[1][4] <<"="<<U_extend[1][6] <<"-"<< dt/dx<<"x"<<f_new[1][5]<<"-"<<f_new[1][4]<<"\n";


    for (int i = 0; i < N; i++) {

        for (int j = 0; j < 4; j++) {

            U_new[j][i] = U_extend[j][i + 2] - (dt / dx) * (f_new[j][i + 1] - f_new[j][i]);

            // cout<<U_new[j][i] <<"="<<U_extend[j][i + 2] <<"-"<< dt/dx<<"x"<<f_new[j][i + 1]<<"-"<<f_new[j][i]<<"\n";

            if (U_new[j][i] < 0.00000000000001) {
                U_new[j][i] = 0.0;
            }

        }

    }

    return U_new;

}




void Precision_bounding(vector<vector<double>> &Property,int N_cell){

    for(int i =0; i<N_cell; i++){
        for(int j =0; j<N_cell-1;j++){

            if ((abs(Property[i][j+1] - Property[i][j])/Property[i][j]) < 0.000000000000001){

                Property[i][j+1] = 100000.0;

            }

        }
    }

}




void Precision_checking(vector<vector<double>> &Property,int N_cell){

    for(int i =0; i<N_cell; i++){
        for(int j =0; j<N_cell;j++){

            if (abs(Property[i][j])  < 0.00000000001){

                Property[i][j] = 0.0;

            }

        }
    }


}



fstream plotfile;

int main(void) {

    // state vector

    vector<vector<double>> UD(N, vector<double> (N, 0.0));
    vector<vector<double>> UMU(N, vector<double> (N, 0.0));
    vector<vector<double>> UMV(N, vector<double> (N, 0.0));
    vector<vector<double>> UE(N, vector<double> (N, 0.0));
    vector<vector<double>> FD(N, vector<double> (N, 0.0));
    vector<vector<double>> FMV(N, vector<double> (N, 0.0));
    vector<vector<double>> FE(N, vector<double> (N, 0.0));

    vector<vector<double>> U1D(N, vector<double> (N, 0.0));
    vector<vector<double>> U1MU(N, vector<double> (N, 0.0));
    vector<vector<double>> U1MV(N, vector<double> (N, 0.0));
    vector<vector<double>> U1E(N, vector<double> (N, 0.0));
    vector<vector<double>> U2D(N, vector<double> (N, 0.0));
    vector<vector<double>> U2MU(N, vector<double> (N, 0.0));
    vector<vector<double>> U2MV(N, vector<double> (N, 0.0));
    vector<vector<double>> U2E(N, vector<double> (N, 0.0));

    vector<vector<double>> U1(N, vector<double> (N, 0.0));
    vector<vector<double>> U2(N, vector<double> (N, 0.0));


    double dx = L / N;
    double dy = W / N;

    vector<vector<double>> uspeed(N, vector<double> (N, 0.0));
    vector<vector<double>> vspeed(N, vector<double> (N, 0.0));
    vector<vector<double>> density(N, vector<double> (N, 0.0));
    vector<vector<double>> pressure(N, vector<double> (N, 0.0));
    vector<vector<double>> energy(N, vector<double> (N, 0.0));
    vector<vector<double>> internal_energy (N, vector<double> (N, 0.0));

    vector<vector<double>> phi(N, vector<double> (N, 0.0));

    // Level-set function


    initial_condition_level_set(dx, phi, N);

    Display_levelset(phi, N);


    initialConditions_ERS(N, dx, uspeed, vspeed, density, pressure, energy);


    Display_ERS(uspeed, vspeed, density, pressure, energy, N);


    vector<vector<double>> phi_extend(N, vector<double> (N, 0.0));
    phi_extend = boundary_condition_level_set(phi, N);
    Display_levelset(phi_extend, N + 2);

    //
    //Split the primitive variables in parts depending on the interface
    //

    vector<vector<double>> density1(N, vector<double> (N, 0.0));
    vector<vector<double>> density2(N, vector<double> (N, 0.0));
    vector<vector<double>> uspeed1(N, vector<double> (N, 0.0));
    vector<vector<double>> uspeed2(N, vector<double> (N, 0.0));
    vector<vector<double>> vspeed1(N, vector<double> (N, 0.0));
    vector<vector<double>> vspeed2(N, vector<double> (N, 0.0));
    vector<vector<double>> pressure1(N, vector<double> (N, 0.0));
    vector<vector<double>> pressure2(N, vector<double> (N, 0.0));
    vector<vector<double>> energy1(N, vector<double> (N, 0.0));
    vector<vector<double>> energy2(N, vector<double> (N, 0.0));


    vector<vector<double>> density1_extend(N+2, vector<double> (N+2, 0.0));
    vector<vector<double>> density2_extend(N+2, vector<double> (N+2, 0.0));
    vector<vector<double>> uspeed1_extend(N+2, vector<double> (N+2, 0.0));
    vector<vector<double>> uspeed2_extend(N+2, vector<double> (N+2, 0.0));
    vector<vector<double>> vspeed1_extend(N+2, vector<double> (N+2, 0.0));
    vector<vector<double>> vspeed2_extend(N+2, vector<double> (N+2, 0.0));
    vector<vector<double>> pressure1_extend(N+2, vector<double> (N+2, 0.0));
    vector<vector<double>> pressure2_extend(N+2, vector<double> (N+2, 0.0));
    vector<vector<double>> energy1_extend(N+2, vector<double> (N+2, 0.0));
    vector<vector<double>> energy2_extend(N+2, vector<double> (N+2, 0.0));


    vector<vector<double>> uspeed_extend(N+4, vector<double> (N+4, 0.0));
    vector<vector<double>> vspeed_extend(N+4, vector<double> (N+4, 0.0));
    vector<vector<double>> density_extend(N+4, vector<double> (N+4, 0.0));
    vector<vector<double>> pressure_extend(N+4, vector<double> (N+4, 0.0));
    vector<vector<double>> energy_extend(N+4, vector<double> (N+4, 0.0));
    vector<vector<double>> normalvelocity_extend(N+4, vector<double> (N+4, 0.0));


    vector<vector<double>> U1D_extend;
    vector<vector<double>> U1MU_extend;
    vector<vector<double>> U1MV_extend;
    vector<vector<double>> U1E_extend;
    vector<vector<double>> U2D_extend;
    vector<vector<double>> U2MU_extend;
    vector<vector<double>> U2MV_extend;
    vector<vector<double>> U2E_extend;

    vector<vector<double>> uspeed1_extend2;
    vector<vector<double>> vspeed1_extend2 ;
    vector<vector<double>> density1_extend2;
    vector<vector<double>> pressure1_extend2;
    vector<vector<double>> uspeed2_extend2;
    vector<vector<double>> vspeed2_extend2;
    vector<vector<double>> density2_extend2;
    vector<vector<double>> pressure2_extend2;

   double t = 0;
   double dt, dt1, dt2;

    //Find the location of the interface


    vector<vector<int>>interfacey;
    vector<vector<int>>interfacexin;

    vector<int>interfacex(N);
    for (int i = 0; i < N; i++) {

        interfacex[i] = i;

    }




    interfacey = interface_findingY(phi, N);
    interfacexin = interface_findingX(phi, N);


    for (int i = 0; i < N; i++) {

        cout << interfacey[0][i] << "\t, \t" << interfacey[1][i] << "\t";

    }


    cout<<"\n";
    for (int i = 0; i < N; i++) {

        cout << interfacexin[0][i] << "\t, \t" << interfacexin[1][i] << "\t";

    }

    //Define the Normal to the Level set fucntion
    //First Nromal
    double *Normal1 = new double[2];
    Normal1[0] = -cos(0.5 * M_PI - IFrad);
    Normal1[1] = -sin(0.5 * M_PI - IFrad);

    double *Normal2 = new double[2];
    Normal2[0] = cos(0.5 * M_PI - IFrad);
    Normal2[1] = sin(0.5 * M_PI - IFrad);



    do{

        //  Display_ERS (uspeed2, vspeed2, density2, pressure2, energy2, N);
        //
        //split the velocity vector into normal and tangential components relative to the interface.
        //

        vector<vector<double>> normalvelocity(N, vector<double> (N, 0.0));
        vector<vector<double>> tangentialvelocity(N, vector<double> (N, 0.0));


        for (int i =0 ; i<N; i++){
            for (int j =0; j<N; j++){

                normalvelocity[i][j] = uspeed[i][j] * cos(0.5 * M_PI -  IFrad) + vspeed[i][j] * cos(IFrad);

                tangentialvelocity[i][j] = uspeed[i][j] * sin(0.5 * M_PI -  IFrad) + vspeed[i][j] * sin(IFrad);


            }

        }


        pressure_extend= transmisive_boundary_condition(pressure);
        density_extend = transmisive_boundary_condition(density);
        normalvelocity_extend = transmisive_boundary_condition(normalvelocity);


        Riemann_based_Ghost_Fluid_Boundary( interfacex, interfacexin, N, dx, pressure, density, uspeed, vspeed, normalvelocity,  pressure_extend,  density_extend, normalvelocity_extend,  phi, density1, density2, uspeed1, uspeed2, vspeed1, vspeed2, pressure1, pressure2, Normal1, Normal2);

        Combinevelocityvectors(interfacexin, tangentialvelocity, uspeed1, uspeed2, vspeed1, vspeed2);




        cout<<"Display materal1\t"<<endl;


//        Display_ERS(density1, uspeed1, vspeed1, pressure1, energy1, N );

        cout<<"Display materal2\t"<<endl;

  //      Display_ERS(density2, uspeed2, vspeed2, pressure2, energy2, N );

        cout<<"\n";

        density1_extend = boundaryConditions2(density1);
        density2_extend = boundaryConditions2(density2);
        uspeed1_extend = boundaryConditions2(uspeed1);
        uspeed2_extend = boundaryConditions2(uspeed2);
        vspeed1_extend = boundaryConditions2(vspeed1);
        vspeed2_extend = boundaryConditions2(vspeed2);
        pressure1_extend = boundaryConditions2(pressure1);
        pressure2_extend = boundaryConditions2(pressure2);
        for (int i = 1 ;i< N+2; i++){
            density1_extend[i][0] = density1_extend[i-1][1];
            density2_extend[i][0] = density2_extend[i-1][1];
            uspeed1_extend[i][0] = uspeed1_extend[i-1][1];
            uspeed2_extend[i][0] = uspeed2_extend[i-1][1];
            vspeed1_extend[i][0] = vspeed1_extend[i-1][1];
            vspeed2_extend[i][0] = vspeed2_extend[i-1][1];
            pressure1_extend[i][0] = pressure1_extend[i-1][1];
            pressure2_extend[i][0] = pressure2_extend[i-1][1];
        }


        for (int i =0; i<N; i++){

            Ghost_Fluid_Boundary_Forward(dx, density1[i], density2[i], density1_extend[i+1], density2_extend[i+1], phi_extend, N, i+1);
            Ghost_Fluid_Boundary_Forward(dx, uspeed1[i], uspeed2[i], uspeed1_extend[i+1], uspeed2_extend[i+1], phi_extend, N, i+1);
            Ghost_Fluid_Boundary_Forward(dx, vspeed1[i], vspeed2[i], vspeed1_extend[i+1], vspeed2_extend[i+1], phi_extend, N, i+1);
            Ghost_Fluid_Boundary_Forward(dx, pressure1[i], pressure2[i], pressure1_extend[i+1], pressure2_extend[i+1], phi_extend, N, i+1);

        }



        cout<<"Display material 1 after extrapolation 1\t"<<endl;

 //       Display_ERS(density1, uspeed1, vspeed1, pressure1, energy1, N );

        density1_extend = boundaryConditions2(density1);
        density2_extend = boundaryConditions2(density2);
        uspeed1_extend = boundaryConditions2(uspeed1);
        uspeed2_extend = boundaryConditions2(uspeed2);
        vspeed1_extend = boundaryConditions2(vspeed1);
        vspeed2_extend = boundaryConditions2(vspeed2);
        pressure1_extend = boundaryConditions2(pressure1);
        pressure2_extend = boundaryConditions2(pressure2);
        cout<<"\n"<<endl;
        //    cout<<"I am working here2"<<endl;


        for (int i =0; i<N; i++){

            Ghost_Fluid_Boundary_Backward(dx, density1[i], density2[i], density1_extend[i+1], density2_extend[i+1], phi_extend, N, i+1);
            Ghost_Fluid_Boundary_Backward(dx, uspeed1[i], uspeed2[i], uspeed1_extend[i+1], uspeed2_extend[i+1], phi_extend, N, i+1);
            Ghost_Fluid_Boundary_Backward(dx, vspeed1[i], vspeed2[i], vspeed1_extend[i+1], vspeed2_extend[i+1], phi_extend, N, i+1);
            Ghost_Fluid_Boundary_Backward(dx, pressure1[i], pressure2[i], pressure1_extend[i+1], pressure2_extend[i+1], phi_extend, N, i+1);

        }



        cout<<"Display material2 extrapolation 1\t"<<endl;

     //   Display_ERS(density2, uspeed2, vspeed2, pressure2, energy2, N );

        cout<<"\n";


        density1_extend = boundaryConditions2( density1);
        density2_extend = boundaryConditions2( density2);
        uspeed1_extend = boundaryConditions2(uspeed1);
        uspeed2_extend = boundaryConditions2(uspeed2);
        vspeed1_extend = boundaryConditions2(vspeed1);
        vspeed2_extend = boundaryConditions2(vspeed2);
        pressure1_extend = boundaryConditions2(pressure1);
        pressure2_extend = boundaryConditions2(pressure2);


        Ghost_Fluid_Boundary_YMaterial1( dx, density1, density1_extend , phi_extend, N);
        Ghost_Fluid_Boundary_YMaterial1( dx, uspeed1, uspeed1_extend , phi_extend, N);
        Ghost_Fluid_Boundary_YMaterial1( dx, vspeed1, vspeed1_extend , phi_extend,  N);
        Ghost_Fluid_Boundary_YMaterial1( dx, pressure1, pressure1_extend , phi_extend, N);



        Ghost_Fluid_Boundary_YMaterial2( dx, density2, density2_extend , phi_extend, N);
        Ghost_Fluid_Boundary_YMaterial2( dx, uspeed2, uspeed2_extend , phi_extend, N);
        Ghost_Fluid_Boundary_YMaterial2( dx, vspeed2, vspeed2_extend , phi_extend, N);
        Ghost_Fluid_Boundary_YMaterial2( dx, pressure2, pressure2_extend , phi_extend, N);

/*

        density1_extend = boundaryConditions2( density1);
        density2_extend = boundaryConditions2( density2);
        uspeed1_extend = boundaryConditions2(uspeed1);
        uspeed2_extend = boundaryConditions2(uspeed2);
        vspeed1_extend = boundaryConditions2(vspeed1);
        vspeed2_extend = boundaryConditions2(vspeed2);
        pressure1_extend = boundaryConditions2(pressure1);
        pressure2_extend = boundaryConditions2(pressure2);


        Ghost_Fluid_Boundary_XMaterial1( dx, density1, density1_extend , phi_extend, N);
        Ghost_Fluid_Boundary_XMaterial1( dx, uspeed1, uspeed1_extend , phi_extend, N);
        Ghost_Fluid_Boundary_XMaterial1( dx, vspeed1, vspeed1_extend , phi_extend,  N);
        Ghost_Fluid_Boundary_XMaterial1( dx, pressure1, pressure1_extend , phi_extend, N);




        Ghost_Fluid_Boundary_XMaterial2( dx, density2, density2_extend , phi_extend, N);
        Ghost_Fluid_Boundary_XMaterial2( dx, uspeed2, uspeed2_extend , phi_extend, N);
        Ghost_Fluid_Boundary_XMaterial2( dx, vspeed2, vspeed2_extend , phi_extend, N);
        Ghost_Fluid_Boundary_XMaterial2( dx, pressure2, pressure2_extend , phi_extend, N);


        density1_extend = boundaryConditions2( density1);
        density2_extend = boundaryConditions2( density2);
        uspeed1_extend = boundaryConditions2(uspeed1);
        uspeed2_extend = boundaryConditions2(uspeed2);
        vspeed1_extend = boundaryConditions2(vspeed1);
        vspeed2_extend = boundaryConditions2(vspeed2);
        pressure1_extend = boundaryConditions2(pressure1);
        pressure2_extend = boundaryConditions2(pressure2);


        Ghost_Fluid_Boundary_YMaterial1( dx, density1, density1_extend , phi_extend, N);
        Ghost_Fluid_Boundary_YMaterial1( dx, uspeed1, uspeed1_extend , phi_extend, N);
        Ghost_Fluid_Boundary_YMaterial1( dx, vspeed1, vspeed1_extend , phi_extend,  N);
        Ghost_Fluid_Boundary_YMaterial1( dx, pressure1, pressure1_extend , phi_extend, N);


        cout<<"Display materal1\t"<<endl;

        Display_ERS(density1, uspeed1, vspeed1, pressure1, energy1, N);


        Ghost_Fluid_Boundary_YMaterial2( dx, density2, density2_extend , phi_extend, N);
        Ghost_Fluid_Boundary_YMaterial2( dx, uspeed2, uspeed2_extend , phi_extend, N);
        Ghost_Fluid_Boundary_YMaterial2( dx, vspeed2, vspeed2_extend , phi_extend, N);
        Ghost_Fluid_Boundary_YMaterial2( dx, pressure2, pressure2_extend , phi_extend, N);

        cout<<"Display materal2\t"<<endl;

        Display_ERS( energy2,density2, uspeed2, vspeed2, pressure2, N);

*/

        Primitive_to_Conservative(N,gam1, gam2, density1, density2, uspeed1, uspeed2, vspeed1, vspeed2, pressure1, pressure2, U1D, U1MU, U1MV, U1E, U2D, U2MU, U2MV, U2E) ;


        //Implement the domain boundarys

        uspeed1_extend2 = GF_Domain_boundaryConditions(N,uspeed1);
        vspeed1_extend2 = GF_Domain_boundaryConditions(N,vspeed1);
        density1_extend2 = GF_Domain_boundaryConditions(N, density1);
        pressure1_extend2 = GF_Domain_boundaryConditions(N,pressure1);
        uspeed2_extend2 = GF_Domain_boundaryConditions(N,uspeed2);
        vspeed2_extend2 = GF_Domain_boundaryConditions(N,vspeed2);
        density2_extend2 = GF_Domain_boundaryConditions(N,density2);
        pressure2_extend2 = GF_Domain_boundaryConditions(N,pressure2);
        dt1 = dt_f(dx, dy, uspeed1_extend2, vspeed1_extend2,  pressure1_extend2, density1_extend2, gam1);

        dt2 =dt_f(dx, dy, uspeed2_extend2, vspeed2_extend2, pressure2_extend2, density2_extend2, gam2);

        if (dt1 >dt2){ dt = dt2;}

        else{ dt = dt1;}

        cout<<"\n";
        cout<<"dt1\t"<<dt1<<"dt2\t"<<dt2<<"dt\t" <<dt<<endl;
        //Update Level set function using level set equation  and reinilisasition


        phi_extend = boundary_condition_level_set(phi, N);



        phi = Hamilton_Jacobi(N, phi_extend, dx,dt, vspeed, uspeed);

       Display_2DarrayN( phi);


        interfacey = interface_findingY(phi, N);
        interfacexin = interface_findingX(phi, N);

        cout<<"Interfaces are"<<endl;
        for (int i =0 ;i<N; i++){

            cout<<interfacexin[0][i]<<"\t, \t"<<interfacexin[1][i]<<"\t";


        }
        cout<<"Interfaces are"<<endl;
        for (int i =0 ;i<N; i++){

            cout<<interfacey[0][i]<<"\t, \t"<<interfacey[1][i]<<"\t";


        }






        U1D_extend = GF_Domain_boundaryConditions(N, U1D);

        U1MU_extend = GF_Domain_boundaryConditions(N, U1MU);

        U1MV_extend = GF_Domain_boundaryConditions(N, U1MV);

        U1E_extend = GF_Domain_boundaryConditions(N, U1E);

        U2D_extend = GF_Domain_boundaryConditions(N, U2D);

        U2MU_extend = GF_Domain_boundaryConditions(N, U2MU);

        U2MV_extend = GF_Domain_boundaryConditions(N, U2MV);

        U2E_extend = GF_Domain_boundaryConditions(N, U2E);


        //--------------------------------------------------------//
        //                                                        //
        //Perform dimensional Spliting                            //
        //Let's start from x direction                            //
        //                                                        //
        //--------------------------------------------------------//

        /*--------------------------X Sweep using inital conditions -----------------------*/


        for(int i = 0 ; i<N; i++){

            //Call the approximate solver


/*
            if(i == 3){

                for(int i =0; i<N; i++){
                    if(U2MU_extend[1][i+2]!=0){

                        cout<<"notzero"<<endl;

                    }
                    if( U2D_extend [1][i+2]!=0){

                        cout<<"notzero"<<endl;

                    }
                }

            }
*/
            U1= U_update(dx, dt, U1D_extend[i+2], U1MU_extend[i+2],  U1MV_extend[i+2], U1E_extend[i+2],  gam1, true);

            U2= U_update(dx, dt, U2D_extend[i+2], U2MU_extend[i+2],  U2MV_extend[i+2], U2E_extend[i+2], gam2, true);
       //     cout<<"\n";
            if(i == 3){

                for(int i =0; i<N; i++){

                    cout<<U2[1][i]<<"\t";
                }

            }


            for (int k =0; k<N; k++){

                U1D[i][k] = U1[0][k];
                U1MU[i][k] = U1[1][k];
                U1MV[i][k] = U1[2][k];
                U1E[i][k] = U1[3][k];
                U2D[i][k] = U2[0][k];
                U2MU[i][k] = U2[1][k];
                U2MV[i][k] = U2[2][k];
                // cout<<U2MV[i][k] <<"\t=\t"<< U2[2][k]<<"\t i \t "<<i<<"\t k \t "<<k<<endl;
                U2E[i][k] = U2[3][k];

            }

        }




    //    cout<<" After xsweep update Material 1"<<endl;

   //     Display_ERS ( U1D, U1MU, U1MV, U1E, U2D, N);

   //     cout<<" After x sweep update Material 2"<<endl;

    //    Display_ERS ( U2D, U2MU, U2MV, U2E, U1D, N);




        U1D_extend = GF_Domain_boundaryConditions(N, U1D);

        U1MU_extend = GF_Domain_boundaryConditions(N, U1MU);

        U1MV_extend = GF_Domain_boundaryConditions(N, U1MV);

        U1E_extend = GF_Domain_boundaryConditions(N, U1E);

        U2D_extend = GF_Domain_boundaryConditions(N, U2D);

        U2MU_extend = GF_Domain_boundaryConditions(N, U2MU);

        U2MV_extend = GF_Domain_boundaryConditions(N, U2MV);

        U2E_extend = GF_Domain_boundaryConditions(N, U2E);

        /*--------------------------Y Sweep using inital conditions -----------------------*/


        for(int i = 0 ; i<N; i++){

            vector<double> U1Dy_extend(N+4);
            vector<double> U1MUy_extend(N+4);
            vector<double> U1MVy_extend(N+4);
            vector<double> U1Ey_extend(N+4);
            vector<double> U2Dy_extend(N+4);
            vector<double>U2MUy_extend (N+4);
            vector<double>U2MVy_extend (N+4);
            vector<double>U2Ey_extend (N+4);

            for (int  j = 0; j<N+4; j++){

                U1Dy_extend[j] = U1D_extend[j][i+2];
                U1MUy_extend[j] = U1MU_extend[j][i+2];
                U1MVy_extend[j] = U1MV_extend[j][i+2];
                U1Ey_extend[j] = U1E_extend[j][i+2];
                U2Dy_extend[j] = U2D_extend[j][i+2];
                U2MUy_extend[j] = U2MU_extend[j][i+2];
                U2MVy_extend[j] = U2MV_extend[j][i+2];
                U2Ey_extend[j] = U2E_extend[j][i+2];


            }


            //Call the approximate solver
            U1= U_update(dx, dt, U1Dy_extend, U1MUy_extend, U1MVy_extend, U1Ey_extend, gam1, false);

            U2= U_update(dx, dt, U2Dy_extend, U2MUy_extend, U2MVy_extend, U2Ey_extend, gam2, false);


            for (int k =0; k<N; k++){

                U1D[k][i] = U1[0][k];
                U1MU[k][i] = U1[1][k];
                U1MV[k][i]= U1[2][k];
                U1E[k][i] = U1[3][k];

                U2D[k][i]= U2[0][k];
                U2MU[k][i] = U2[1][k];
                U2MV[k][i] = U2[2][k];
                U2E[k][i] = U2[3][k];
            }


        }
        for (int i =0; i<N; i++){

            for (int j=0; j<N;j++){

                density1[i][j] = U1D[i][j] ;

                uspeed1[i][j] = U1MU[i][j] / U1D[i][j];

                vspeed1[i][j] = U1MV[i][j] / U1D[i][j];

                pressure1[i][j] = (U1E[i][j] - 0.5 * U1D[i][j] * (pow(( U1MU[i][j] / U1D[i][j]),2) + pow(( U1MV[i][j] / U1D[i][j]),2))) * (gam1 -1.0);

                energy1[i][j] = U1E[i][j] ;

                density2[i][j] = U2D[i][j] ;

                uspeed2[i][j] = U2MU[i][j] / U2D[i][j];

                vspeed2[i][j] = U2MV[i][j] / U2D[i][j];

                pressure2[i][j] = (U2E[i][j] - 0.5 * U2D[i][j] * (pow(( U2MU[i][j] / U2D[i][j]),2) + pow(( U2MV[i][j] / U2D[i][j]),2))) * (gam2 -1.0);

                energy2[i][j] = U2E[i][j] ;


            }

            //    cout<<  vspeed2[i][0] << "\t=\t"<< U2MV[i][0] <<"\t/\t"<< U2D[i][0]<<endl;

        }


        for (int i = 0 ; i< N;i++){

            // cout<<"Interface location \t "<<interfacey[0][i]+1<<endl;

            for (int j =0; j < interfacexin[0][i] ; j++){

                pressure[i][j]= pressure1[i][j];
                uspeed[i][j] = uspeed1[i][j];
                vspeed[i][j] = vspeed1[i][j];
                density[i][j] = density1[i][j];
                energy[i][j] = energy1[i][j];

            }

            //cout<<"Interface location \t "<<interfacey[1][i]+1<<endl;

            for (int j =  interfacexin[0][i]; j<interfacexin[1][i]; j++){

                pressure[i][j]= pressure2[i][j];
                uspeed[i][j] = uspeed2[i][j];
                vspeed[i][j] = vspeed2[i][j];
                density[i][j] = density2[i][j];
                energy[i][j] = energy2[i][j];

            }

            for (int j = interfacexin[1][i]; j<N; j++){

                pressure[i][j]= pressure1[i][j];
                uspeed[i][j] = uspeed1[i][j];
                vspeed[i][j] = vspeed1[i][j];
                density[i][j] = density1[i][j];
                energy[i][j] = energy1[i][j];
            }



        }


        Precision_bounding(pressure,N);
   //     Precision_checking(uspeed,N);
   //     Precision_checking(vspeed,N);


/*

        cout<<"\n"<<endl;

        cout<<"Display materal1\t"<<endl;

        Display_ERS(density1, uspeed1, vspeed1, pressure1, energy1, N );

        cout<<"Display materal2\t"<<endl;

        Display_ERS(density2, uspeed2, vspeed2, pressure2, energy2, N );

        cout<<"\n";
*/




        cout<<"Final velocity Profile"<<endl;

        Display_2DarrayN(uspeed);

        cout<<"\n";

        cout<<"Final velocity2 Profile"<<endl;

        Display_2DarrayN(vspeed);

        cout<<"\n";
        cout<<"Final density Profile"<<endl;

        Display_2DarrayN( density);

        cout<<"\n";
        cout<<"Final pressure Profile"<<endl;

        Display_2DarrayN(pressure);

        cout<<"\n";

       t+=dt;

        cout<<"dt\t"<<dt<<"\tt\t"<<t<<endl;


     //   t=tStop;

        //
        //Free all memory
        //

    }while( t<tStop);

    plotfile.open("2D_Angled_interface_density.dat", fstream::out);

    if(plotfile.fail()){

        cout<<"Error opening file"<<endl;
        exit(0);
    }

    for (int i=0;i<N;i++){

        double x= i*dx;

        for (int j =0; j<N;j++){

            double y= j*dx;

            plotfile<<x<<"\t"<<y<<"\t"<<density[i][j]<<"\t"<<endl;
        }


    }

    return 0;


}



