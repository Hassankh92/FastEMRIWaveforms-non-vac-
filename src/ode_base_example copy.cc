#include <math.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_ellint.h>
#include <algorithm>

#include <hdf5.h>
#include <hdf5_hl.h>

#include <complex>
#include <cmath>

#include "Interpolant.h"
#include "global.h"
#include "Utility.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>
#include <iomanip> // std::setprecision
#include <cstring>

#include "dIdt8H_5PNe10.h"
#include "ode.hh"

#define pn5_Y
#define pn5_citation1 Pn5_citation
__deriv__ void pn5(double *pdot, double *edot, double *Ydot,
                   double *Omega_phi, double *Omega_theta, double *Omega_r,
                   double epsilon, double a, double p, double e, double Y, double *additional_args)
{
    // evaluate ODEs

    // the frequency variables are pointers!
    double x = Y_to_xI(a, p, e, Y);
    KerrGeoCoordinateFrequencies(Omega_phi, Omega_theta, Omega_r, a, p, e, x);

    int Nv = 10;
    int ne = 10;
    *pdot = epsilon * dpdt8H_5PNe10(a, p, e, Y, Nv, ne);

    // needs adjustment for validity
    Nv = 10;
    ne = 8;
    *edot = epsilon * dedt8H_5PNe10(a, p, e, Y, Nv, ne);

    Nv = 7;
    ne = 10;
    *Ydot = epsilon * dYdt8H_5PNe10(a, p, e, Y, Nv, ne);
}

// Initialize flux data for inspiral calculations
void load_and_interpolate_flux_data(struct interp_params *interps, const std::string &few_dir)
{

    // Load and interpolate the flux data
    std::string fp = "few/files/FluxNewMinusPNScaled_fixed_y_order.dat";
    fp = few_dir + fp;
    ifstream Flux_file(fp);

    if (Flux_file.fail())
    {
        throw std::runtime_error("The file FluxNewMinusPNScaled_fixed_y_order.dat did not open sucessfully. Make sure it is located in the proper directory (Path/to/Installation/few/files/).");
    }

    // Load the flux data into arrays
    string Flux_string;
    vector<double> ys, es, Edots, Ldots;
    double y, e, Edot, Ldot;
    while (getline(Flux_file, Flux_string))
    {

        stringstream Flux_ss(Flux_string);

        Flux_ss >> y >> e >> Edot >> Ldot;

        ys.push_back(y);
        es.push_back(e);
        Edots.push_back(Edot);
        Ldots.push_back(Ldot);
    }

    // Remove duplicate elements (only works if ys are perfectly repeating with no round off errors)
    sort(ys.begin(), ys.end());
    ys.erase(unique(ys.begin(), ys.end()), ys.end());

    sort(es.begin(), es.end());
    es.erase(unique(es.begin(), es.end()), es.end());

    Interpolant *Edot_interp = new Interpolant(ys, es, Edots);
    Interpolant *Ldot_interp = new Interpolant(ys, es, Ldots);

    interps->Edot = Edot_interp;
    interps->Ldot = Ldot_interp;
}

// Class to carry gsl interpolants for the inspiral data
// also executes inspiral calculations
SchwarzEccFlux::SchwarzEccFlux(std::string few_dir)
{
    interps = new interp_params;

    // prepare the data
    // python will download the data if
    // the user does not have it in the correct place
    load_and_interpolate_flux_data(interps, few_dir);
    // load_and_interpolate_amp_vec_norm_data(&amp_vec_norm_interp, few_dir);
}

#define SchwarzEccFlux_num_add_args 0
#define SchwarzEccFlux_spinless
#define SchwarzEccFlux_equatorial
#define SchwarzEccFlux_file1 FluxNewMinusPNScaled_fixed_y_order.dat
__deriv__ void SchwarzEccFlux::deriv_func(double *pdot, double *edot, double *xdot,
                                          double *Omega_phi, double *Omega_theta, double *Omega_r,
                                          double epsilon, double a, double p, double e, double x, double *additional_args)
{
    if ((6.0 + 2. * e) > p)
    {
        *pdot = 0.0;
        *edot = 0.0;
        *xdot = 0.0;
        return;
    }

    SchwarzschildGeoCoordinateFrequencies(Omega_phi, Omega_r, p, e);
    *Omega_theta = *Omega_phi;

    double y1 = log((p - 2. * e - 2.1));

    // evaluate ODEs, starting with PN contribution, then interpolating over remaining flux contribution

    double yPN = pow((*Omega_phi), 2. / 3.);

    double EdotPN = (96 + 292 * Power(e, 2) + 37 * Power(e, 4)) / (15. * Power(1 - Power(e, 2), 3.5)) * pow(yPN, 5);
    double LdotPN = (4 * (8 + 7 * Power(e, 2))) / (5. * Power(-1 + Power(e, 2), 2)) * pow(yPN, 7. / 2.);

    double Edot = -epsilon * (interps->Edot->eval(y1, e) * pow(yPN, 6.) + EdotPN);

    double Ldot = -epsilon * (interps->Ldot->eval(y1, e) * pow(yPN, 9. / 2.) + LdotPN);

    *pdot = (-2 * (Edot * Sqrt((4 * Power(e, 2) - Power(-2 + p, 2)) / (3 + Power(e, 2) - p)) * (3 + Power(e, 2) - p) * Power(p, 1.5) + Ldot * Power(-4 + p, 2) * Sqrt(-3 - Power(e, 2) + p))) / (4 * Power(e, 2) - Power(-6 + p, 2));

    // handle e = 0.0
    if (e > 0.)
    {
        *edot = -((Edot * Sqrt((4 * Power(e, 2) - Power(-2 + p, 2)) / (3 + Power(e, 2) - p)) * Power(p, 1.5) *
                       (18 + 2 * Power(e, 4) - 3 * Power(e, 2) * (-4 + p) - 9 * p + Power(p, 2)) +
                   (-1 + Power(e, 2)) * Ldot * Sqrt(-3 - Power(e, 2) + p) * (12 + 4 * Power(e, 2) - 8 * p + Power(p, 2))) /
                  (e * (4 * Power(e, 2) - Power(-6 + p, 2)) * p));
    }
    else
    {
        *edot = 0.0;
    }

    *xdot = 0.0;
}

// destructor
SchwarzEccFlux::~SchwarzEccFlux()
{

    delete interps->Edot;
    delete interps->Ldot;
    delete interps;
}




// void load_and_interpolate_flux_data_kerr_circ(struct interp_params *interps, const std::string &few_dir)
// {

//     // Load and interpolate the flux data
//     std::string fp = "/few/files/Kerr_flux_minus_PN1_kerr_circ.dat";
//     fp = few_dir + fp;
//     ifstream Flux_file(fp);

//     if (Flux_file.fail())
//     {
//         throw std::runtime_error("The file Kerr_flux_minus_PN1_kerr_circ did not open sucessfully. Make sure it is located in the proper directory (Path/to/Installation/few/files/).");
//     }

//     // Load the flux data into arrays
//     string Flux_string;
//     vector<double> ys, as, EdotsKerr;
//     double y, a, EdotKerr;
//     while (getline(Flux_file, Flux_string))
//     {

//         stringstream Flux_ss(Flux_string);

//         Flux_ss >> y >> a >> EdotKerr;

//         ys.push_back(y);
//         as.push_back(a);
//         EdotsKerr.push_back(EdotKerr);
//     }
//     // printf("in the kerr load %1.6e, %1.6e\n", a, y);

//     // Remove duplicate elements (only works if ys are perfectly repeating with no round off errors)
//     sort(ys.begin(), ys.end());
//     ys.erase(unique(ys.begin(), ys.end()), ys.end());

//     sort(as.begin(), as.end());
//     as.erase(unique(as.begin(), as.end()), as.end());

//     // notice that if you resort ys and a you have to change also Edots

//     Interpolant *EdotKerr_interp = new Interpolant(ys, as, EdotsKerr);

//     interps->EdotKerrCirc = EdotKerr_interp;
// }

// // KerrCircFlux_Has

// KerrCircFlux::KerrCircFlux(std::string few_dir)
// {
//     interps = new interp_params;
//     load_and_interpolate_flux_data_kerr_circ(interps, few_dir);
// }

// double KerrCircFlux::EdotPN(double r, double a)
// {
//     double y = pow(1. / (sqrt(r * r * r) + a), 2. / 3.);
//     double res = 6.4 * pow(y, 5);
//     return res;
// }

// #define KerrCircFlux_num_add_args 0
// #define KerrCircFlux_equatorial
// #define KerrCircFlux_circular
// #define KerrCircFlux_file1 Kerr_flux_minus_PN1_kerr_circ.dat
// __deriv__
// void KerrCircFlux::deriv_func(double *pdot, double *edot, double *xdot,
//                                                        double *Omega_phi, double *Omega_theta, double *Omega_r,
//                                                        double epsilon, double a, double p, double e, double x, double *additional_args)
// {
//     double p_sep = get_separatrix(a, e, x);
//     if (p_sep > p)
//     {
//         *pdot = 0.0;
//         *edot = 0.0;
//         *xdot = 0.0;
//         return;
//     }

//     double u = log((p - p_sep + 3.9));

//     KerrGeoCoordinateFrequencies(Omega_phi, Omega_theta, Omega_r, a, p, e, x);
//     *Omega_theta = *Omega_phi;
//     *Omega_r = *Omega_phi;
//     // evaluate ODEs, starting with PN contribution, then interpolating over remaining flux contribution

//     double yPNKerr = pow(1. / (sqrt(p * p * p) + a), 2. / 3.);

//     double EdotKerr = epsilon * (interps->EdotKerrCirc->eval(u, a) * pow(yPNKerr, 6.) + EdotPN(p, a));
//     double LdotKerr = EdotKerr / *Omega_phi;

//     double dL_dp = (-3 * Power(a, 3) + Power(a, 2) * (8 - 3 * p) * Sqrt(p) + (-6 + p) * Power(p, 2.5) + 3 * a * p * (-2 + 3 * p)) / (2. * Power(2 * a + (-3 + p) * Sqrt(p), 1.5) * Power(p, 1.75));
//     *pdot = -LdotKerr / dL_dp;

//     *edot = 0.0;
//     *xdot = 0.0;
// }

// // destructor
// KerrCircFlux::~KerrCircFlux()
// {

//     delete interps->EdotKerrCirc;
//     delete interps;
// }



void load_and_interpolate_flux_data_kerr_circ(struct interp_params *interps, const std::string &few_dir)
{

    // Load and interpolate the flux data
    std::string fp = "/few/files/Kerr_flux_minus_PN1_kerr_circ.dat";
    fp = few_dir + fp;
    ifstream Flux_file(fp);

    if (Flux_file.fail())
    {
        throw std::runtime_error("The file Kerr_flux_minus_PN1_kerr_circ did not open sucessfully. Make sure it is located in the proper directory (Path/to/Installation/few/files/).");
    }
    // printf("seccessfully opened the file\n");

    // Load the flux data into arrays
    string Flux_string;
    vector<double> ys, as, Edots_Kerr;
    double y, a, Edot_Kerr;
    while (getline(Flux_file, Flux_string))
    {

        stringstream Flux_ss(Flux_string);

        Flux_ss >> y >> a >> Edot_Kerr;

        ys.push_back(y);
        as.push_back(a);
        Edots_Kerr.push_back(Edot_Kerr);
        // printf("Edots_Kerr %1.6e\n", Edot_Kerr);
    }
    // printf("in the kerr load %1.6e, %1.6e\n", a, y);

    // Remove duplicate elements (only works if ys are perfectly repeating with no round off errors)
    sort(ys.begin(), ys.end());
    ys.erase(unique(ys.begin(), ys.end()), ys.end());

    sort(as.begin(), as.end());
    as.erase(unique(as.begin(), as.end()), as.end());

    // notice that if you resort ys and a you have to change also Edots

    Interpolant *Edot_Kerr_interp = new Interpolant(ys, as, Edots_Kerr);

    interps->Edot_Kerr = Edot_Kerr_interp;
}

// KerrCircFlux_Has

KerrCircFlux::KerrCircFlux(std::string few_dir)
{
    interps = new interp_params;
    load_and_interpolate_flux_data_kerr_circ(interps, few_dir);
}

double KerrCircFlux::EdotPN(double r, double a)
{
    double y = pow(1. / (sqrt(r * r * r) + a), 2. / 3.);
    double res = 6.4 * pow(y, 5);
    return res;
}

#define KerrCircFlux_num_add_args 0
#define KerrCircFlux_equatorial
#define KerrCircFlux_circular
#define KerrCircFlux_file1 Kerr_flux_minus_PN1_kerr_circ.dat
__deriv__
void KerrCircFlux::deriv_func(double *pdot, double *edot, double *xdot,
                                                       double *Omega_phi, double *Omega_theta, double *Omega_r,
                                                       double epsilon, double a, double p, double e, double x, double *additional_args)
{
    double p_sep = get_separatrix(a, e, x);
    if (p_sep > p)
    {
        *pdot = 0.0;
        *edot = 0.0;
        *xdot = 0.0;
        return;
    }

    double u = log((p - p_sep + 3.9));

    KerrGeoCoordinateFrequencies(Omega_phi, Omega_theta, Omega_r, a, p, e, x);
    *Omega_theta = *Omega_phi;
    *Omega_r = *Omega_phi;
    // evaluate ODEs, starting with PN contribution, then interpolating over remaining flux contribution

    double yPNKerr = pow(1. / (sqrt(p * p * p) + a), 2. / 3.);

    double Edot_Kerr = epsilon * (interps->Edot_Kerr->eval(u, a) * pow(yPNKerr, 6.) + EdotPN(p, a));
    double Ldot_Kerr = Edot_Kerr / *Omega_phi;

    double dL_dp = (-3 * Power(a, 3) + Power(a, 2) * (8 - 3 * p) * Sqrt(p) + (-6 + p) * Power(p, 2.5) + 3 * a * p * (-2 + 3 * p)) / (2. * Power(2 * a + (-3 + p) * Sqrt(p), 1.5) * Power(p, 1.75));
    *pdot = -Ldot_Kerr / dL_dp;

    *edot = 0.0;
    *xdot = 0.0;
}

// destructor
KerrCircFlux::~KerrCircFlux()
{

    delete interps->Edot_Kerr;
    delete interps;
}







// Disk effect with the new relativistic Kerr
MigTorqKerrCircFlux::MigTorqKerrCircFlux(std::string few_dir)
{
    interps = new interp_params;
    load_and_interpolate_flux_data_kerr_circ(interps, few_dir);
}


double MigTorqKerrCircFlux::EdotPN(double r, double a)
{
    double y = pow(1./(sqrt(r*r*r) + a), 2./3.) ;
    double res = 6.4*pow(y,5); 
    return res;
}


#define MigTorqKerrCircFlux_num_add_args 2
#define MigTorqKerrCircFlux_equatorial
#define MigTorqKerrCircFlux_circular
#define MigTorqKerrCircFlux_file1 Kerr_flux_minus_PN1_kerr_circ.dat
__deriv__
void MigTorqKerrCircFlux::deriv_func(double* pdot, double* edot, double* xdot,
                      double* Omega_phi, double* Omega_theta, double* Omega_r,
                      double epsilon, double a, double p, double e, double x, double* additional_args)
{
    double p_sep = get_separatrix(a, e, x);
    if (p_sep > p)
    {
        *pdot = 0.0;
        *edot = 0.0;
        *xdot = 0.0;
        return;
    }

    double u = log(p - p_sep + 3.9);
  
  //  /*Disk effect */
    double A = additional_args[0];
    double nr = additional_args[1];

    double LdotPN_circ = +32/5.* pow(p, -7./2.);
    double Ldot_disk = epsilon * A * Power((p/10.),nr) * LdotPN_circ;
   
    double yPNKerr = pow(1./(sqrt(p*p*p) + a), 2./3.) ;

    KerrGeoCoordinateFrequencies(Omega_phi, Omega_theta, Omega_r, a, p, e, x);
     *Omega_theta = *Omega_phi;
        *Omega_r = *Omega_phi;
  

    double EdotKerr = epsilon*(interps->EdotKerrCirc->eval(u, a)*pow(yPNKerr,6.) + EdotPN(p,a)) ;
    double LdotKerr = EdotKerr/ (*Omega_phi) ;
    
    double dL_dp = (-3*Power(a,3) + Power(a,2)*(8 - 3*p)*Sqrt(p) + (-6 + p)*Power(p,2.5) + 3*a*p*(-2 + 3*p))/(2.*Power(2*a + (-3 + p)*Sqrt(p),1.5)*Power(p,1.75));
    *pdot =  -(LdotKerr + Ldot_disk  )/dL_dp;

    *edot = 0.0;

    *xdot = 0.0;
}

// destructor
MigTorqKerrCircFlux::~MigTorqKerrCircFlux()
{

    delete interps->EdotKerrCirc;
    delete interps;

}






void load_and_interpolate_flux_data_cloud(struct interp_params *interps, const std::string& few_dir){

	// Load and interpolate the flux data
    std::string fp = "/few/files/Relativistic_cloud_data.dat";
    fp = few_dir + fp;
	ifstream Flux_file(fp);

    if (Flux_file.fail())
    {
        throw std::runtime_error("The file Relativistic_cloud_data.dat did not open sucessfully. Make sure it is located in the proper directory (Path/to/Installation/few/files/).");
    }

	// Load the flux data into arrays
	string Flux_string;
	vector<double> ps, alphas, Edots_cloud_Hor00, Edots_cloud_Inf00, Edots_cloud_Hor11, Edots_cloud_Inf11, Edots_cloud_Hor22, Edots_cloud_Inf22, Edots_cloud_Hor33, Edots_cloud_Inf33, Edots_cloud_Hor44, Edots_cloud_Inf44, Edots_cloud_Hor_all, Edots_cloud_Inf_all;
	double p, alpha, Edot_cloud_hor00, Edot_cloud_inf00, Edot_cloud_hor11, Edot_cloud_inf11, Edot_cloud_hor22, Edot_cloud_inf22, Edot_cloud_hor33, Edot_cloud_inf33, Edot_cloud_hor44, Edot_cloud_inf44, Edot_cloud_hor_all, Edot_cloud_inf_all;;
	while(getline(Flux_file, Flux_string)){

		stringstream Flux_ss(Flux_string);

        Flux_ss  >> p >> alpha >> Edot_cloud_inf00 >> Edot_cloud_hor00 >> Edot_cloud_inf11 >> Edot_cloud_hor11 >> Edot_cloud_inf22 >> Edot_cloud_hor22 >> Edot_cloud_inf33 >> Edot_cloud_hor33 >> Edot_cloud_inf44 >> Edot_cloud_hor44 >> Edot_cloud_inf_all >> Edot_cloud_hor_all;


		ps.push_back(p);
		alphas.push_back(alpha);

        Edots_cloud_Hor_all.push_back(Edot_cloud_hor_all);
        Edots_cloud_Inf_all.push_back(Edot_cloud_inf_all);



	}
	// Remove duplicate elements (only works if ys are perfectly repeating with no round off errors)
	sort( ps.begin(), ps.end() );
	ps.erase( unique( ps.begin(), ps.end() ), ps.end() );

	sort( alphas.begin(), alphas.end() );
	alphas.erase( unique( alphas.begin(), alphas.end() ), alphas.end() );

    // notice that if you resort ys and a you have to change also Edots


    Interpolant *Edot_cloud_hor_interp_all = new Interpolant(ps, alphas, Edots_cloud_Hor_all);
    Interpolant *Edot_cloud_inf_interp_all = new Interpolant(ps, alphas, Edots_cloud_Inf_all);


    interps->Edot_cloud_hor_all = Edot_cloud_hor_interp_all;
    interps->Edot_cloud_inf_all = Edot_cloud_inf_interp_all;


}



CloudKerrCircFlux::CloudKerrCircFlux(std::string few_dir)
{
    interps = new interp_params;
    load_and_interpolate_flux_data_kerr_circ(interps, few_dir);
    load_and_interpolate_flux_data_cloud(interps, few_dir);

}


double CloudKerrCircFlux::EdotPN(double r, double a)
{
    double y = pow(1./(sqrt(r*r*r) + a), 2./3.) ;
    double res = 6.4*pow(y,5); 
    return res;
}


#define CloudKerrCircFlux_num_add_args 2
#define CloudKerrCircFlux_circular
#define CloudKerrCircFlux_equatorial
#define CloudKerrCircFlux_file1 Kerr_flux_minus_PN1_kerr_circ.dat
#define CloudKerrCircFlux_file2 Relativistic_cloud_data.dat

__deriv__
void CloudKerrCircFlux::deriv_func(double* pdot, double* edot, double* xdot,
                  double* Omega_phi, double* Omega_theta, double* Omega_r,
                  double epsilon, double a, double p, double e, double x, double* additional_args)
{

    double p_sep = get_separatrix(a, e, x); 
    

    if (p_sep > p)
    {
        *pdot = 0.0;
        *edot = 0.0;
        *xdot = 0.0;
        return;
    }
    
    double u = log((p - p_sep + 3.9));
    double yPN = pow(1./(sqrt(p*p*p) + a), 2./3.) ;

    
    // additional_args[3] = -1;
    // the cloud effect calculations:  
    double alpha = additional_args[0];  
    double Mb_M = additional_args[1];
    double cloud_horizon_flux_switch = additional_args[2];
    cloud_horizon_flux_switch = 1.0;

    double cloud_Edot_Hor;// = 0.0;
    double cloud_Edot_Inf;// = 0.0;


 if(alpha==0){
         cloud_Edot_Hor = 0.0;
         cloud_Edot_Inf = 0.0;}
 else{
    if (p < 6.1){
        cloud_Edot_Hor = 0.0;
        cloud_Edot_Inf = 0.0;
        }
    else{     
        cloud_Edot_Hor = Mb_M * interps->Edot_cloud_hor_all->eval(p, alpha);
        cloud_Edot_Inf = Mb_M * interps->Edot_cloud_inf_all->eval(p, alpha);
        }
    }



    KerrGeoCoordinateFrequencies(Omega_phi, Omega_theta, Omega_r, a, p, e, x);
     *Omega_theta = *Omega_phi;
        *Omega_r = *Omega_phi;
    // evaluate ODEs, starting with PN contribution, then interpolating over remaining flux contribution

    double EdotKerr = (interps->EdotKerrCirc->eval(u, a)*pow(yPN,6.) + EdotPN(p,a));
    double Edot = epsilon*(EdotKerr + cloud_horizon_flux_switch * cloud_Edot_Hor + cloud_Edot_Inf) ;
    double Ldot = Edot/ *Omega_phi;
    
    double dL_dp = (-3*Power(a,3) + Power(a,2)*(8 - 3*p)*Sqrt(p) + (-6 + p)*Power(p,2.5) + 3*a*p*(-2 + 3*p))/(2.*Power(2*a + (-3 + p)*Sqrt(p),1.5)*Power(p,1.75));
    *pdot =  -Ldot/dL_dp;

    *edot = 0.0;
    *xdot = 0.0;



}

// destructor
CloudKerrCircFlux::~CloudKerrCircFlux()
{

    delete interps->EdotKerrCirc;
    delete interps->Edot_cloud_hor_all; 
    delete interps->Edot_cloud_inf_all;
    delete interps;


}





