/**
 * First attempt at building solar system model with Mercury perturbations at regular intervals
 *
 *  TO DO:
 *    x Change planet locations to astronomical units
 *    x Implement normal distribution generator
 *    x Add Mercury perturbations to heartbeat function
 *    o Implement multiple simulations
 *    o Arrays for ecc max and min
 *    o Implement sorted stratified resampling method
 *    o Splitting and killing
 *
 *
 *    UNITS:
 *        distance: Astronomical Unit
 *        time    : Earth year/2pi
 *        mass    : Solar Mass
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

#define SOLAR_MASS 1.988544e30     // Solar Mass in kg
#define AU 149597870700         // Astronomical Unit in m

double ss_pos[10][3] =
{
    {-0.0071371792,  -0.0027959975,   0.0002062985},   // Sun
    {-0.1372307845,  -0.4500836156,  -0.0243920085},   // Mercury
    {-0.7254394755,  -0.0354503057,	  0.0412204805},   // Venus
    {-0.1842959633,   0.9644233550,	  0.0002051592},   // Earth
    { 1.3835787426,  -0.0162123156,	 -0.0342613643},   // Mars
    { 3.9940399820,   2.9357801095,	 -0.1015789781},   // Jupiter
    { 6.3992716804,   6.5671936194,	 -0.3688701277},   // Saturn
    { 14.4247194352, -13.7371174663, -0.2379353893},   // Uranus
    { 16.8049097889, -24.9945588868,  0.1274291784},   // Neptune
    {-9.8824894091,	 -27.9615926237,  5.8506522150},   // Pluto
};

double ss_vel[10][3] =
{
    { 0.0003126630,	-0.0004305821, -0.0000054844},   // Sun
    { 1.2423933957,	-0.3752679646, -0.1446310572},   // Mercury
    { 0.0467091823,	-1.1802411832, -0.0188085105},   // Venus
    {-0.9997460569,	-0.1843565013, -0.0000040752},   // Earth
    { 0.0393485692,	 0.8824412195,  0.0175302689},   // Mars
    {-0.2652545529,	 0.3741287103,  0.0043881376},   // Jupiter
    {-0.2492122285,	 0.2257229652,  0.0059791254},   // Saturn
    { 0.1559974780,	 0.1549397540, -0.0014455449},   // Uranus
    { 0.1502521714,	 0.1028649997, -0.0055809444},   // Neptune
    { 0.1763813425,	-0.0898258456, -0.0414075629},   // Pluto
};

double ss_mass[10] =            // Masses relative to Solar Mass
{
    1.,                         // Sun
    3.302e23/SOLAR_MASS,        // Mercury
    48.685e23/SOLAR_MASS,       // Venus
    6.0477246e24/SOLAR_MASS,    // Earth
    6.4185e23/SOLAR_MASS,       // Mars
    1898.13e24/SOLAR_MASS,      // Jupiter
    5.68319e26/SOLAR_MASS,      // Saturn
    86.8103e24/SOLAR_MASS,      // Uranus
    102.41e24/SOLAR_MASS,       // Neptune
    1.4639248e22/SOLAR_MASS,    // Pluto
};


void heartbeat(struct reb_simulation* r, double merc_ecc_max, double merc_ecc_min);
double tmax;

double gaussian(){
    /*
    Compute a Guassian random variable using the Marsaglia (Box-Muller) method
    */
    double u, v, s, z;
    s = 1.1;
    while (s < 1.e-5 || s >= 1.){
      u = 2*((double) rand()/(double) RAND_MAX) - 1.;
      v = 2*((double) rand()/(double) RAND_MAX) - 1.;

      s = pow(u, 2.) + pow(v, 2.);
    }

    z = u * pow((-2.*log(s))/s, .5);

    return z;
}

struct reb_simulation* init_sim(){
    struct reb_simulation* r = reb_create_simulation();
    // Setup constants
    // **** CHECK UNITS
    r->dt             = pow(65., .5)*2*M_PI/365.25; // Corresponds to ~8.062 days
    tmax              = 5e6*2*M_PI;            // 5 Myr
    r->G              = 1.;               // in AU^3 / SM / (year/2pi)^2
    r->ri_whfast.safe_mode     = 0;        // Turn off safe mode. Need to call reb_integrator_synchronize() before outputs.
    r->ri_whfast.corrector     = 11;        // 11th order symplectic corrector
    r->integrator        = REB_INTEGRATOR_WHFAST;
    r->heartbeat        = heartbeat;
    r->exact_finish_time = 1; // Finish exactly at tmax in reb_integrate(). Default is already 1.
    //r->integrator        = REB_INTEGRATOR_IAS15;        // Alternative non-symplectic integrator

    // Initial conditions
    for (int i=0;i<10;i++){
        struct reb_particle p = {0};
        p.x  = ss_pos[i][0];         p.y  = ss_pos[i][1];         p.z  = ss_pos[i][2];
        p.vx = ss_vel[i][0];         p.vy = ss_vel[i][1];         p.vz = ss_vel[i][2];
        p.m  = ss_mass[i];
        reb_add(r, p);
    }

    // Initial Gaussian perturbation to Mercury x-coord
    struct reb_particle merc = r->particles[1];
    merc.x += (0.38/AU)*gaussian();

    reb_move_to_com(r);

    return r;
}

void heartbeat(struct reb_simulation* r, double merc_ecc_max, double merc_ecc_min){
    if (reb_output_check(r, 10000.)){           // Display (default heartbeat function)
        reb_output_timing(r, tmax);
        reb_integrator_synchronize(r);
        /*
        FILE* f = fopen("energy.txt","a");
        double e = reb_tools_energy(r);
        fprintf(f,"%e %e\n",r->t, fabs((e-e_init)/e_init));
        fclose(f);
        */
    }
    if (reb_output_check(r, 5e5*2*M_PI)){         // Perturb Mercury x-coord every 0.5 Myr
      struct reb_particle merc = r->particles[1];
      double pert = 0.38*gaussian();
      merc.x += pert/AU;
      printf("\nPerturbed Mercury's x-coordinate by %f m\n", pert);
    }
    if (reb_output_check(r, 1e6*2*M_PI)){         // Splitting and killing every 1 Myr
      double merc_ecc_range = merc_ecc_max - merc_ecc_min;

    }
}

int main(int argc, char* argv[]){

    int N = 1;
    if (argc == 1){
      printf("Initialising a single simulation\n\n");
    } else if (argc == 2){
      N = atoi(argv[1]);
      printf("Initialising %d simulations\n\n", N);
    } else {
      printf("Incorrect input arguments: aborting\n");
      return 1;
    }

    struct reb_simulation** sims = malloc(N*sizeof(struct reb_simulation*));
    struct rebx_extras** rebx = malloc(N*sizeof(struct rebx_extras*));
    for (int i = 0; i < N; i++){
      printf("Initialising simulation %d\n", i+1);
      sims[i] = init_sim();

      rebx[i] = rebx_attach(sims[i]);
      // Could also add "gr" or "gr_full" here.  See documentation for details.
      struct rebx_force* gr = rebx_load_force(rebx[i], "gr");
      rebx_add_force(rebx[i], gr);
      // Have to set speed of light in right units (set by G & initial conditions).  Here we use default units of AU/(yr/2pi)
      rebx_set_param_double(rebx[i], &gr->ap, "c", 10065.32);
    }

    // Get initial values of Mercury's eccentricity
    /*
    struct reb_particle merc = sim->particles[1];
    struct reb_orbit merc_orb = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
    double merc_ecc_max = merc_orb.e;
    double merc_ecc_min = merc_orb.e;

    //double tmax = 5.e-1;
    reb_integrate(sim, tmax);
    //rebx_free(rebx);    // this explicitly frees all the memory allocated by REBOUNDx
    reb_free_simulation(sim);
    */
    for (int i = 0; i < N; i++){
      printf("Freeing simulation %d\n", i+1);

      rebx_free(rebx[i]);
      reb_free_simulation(sims[i]);
    }
    free(sims);

    printf("\n ====END==== \n");
}
