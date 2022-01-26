/**
 * Brute force method.
 * Multiple simulations with Mercury perturbation, but no DMC

 ******** REQUIRES MODIFIED rebound.h FILE WITH NEW SIMULATION PARAMETER sim_id ********

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
#include <omp.h>
#include "rebound.h"
#include "reboundx.h"

#define SOLAR_MASS 1.988544e30    // Solar Mass in kg
#define AU 149597870700           // Astronomical Unit in m

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


void heartbeat(struct reb_simulation* r);

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

double r_hill(struct reb_simulation* r, int planet_id){
    /*
    Compute the Hill radius of a specified planet in simulation r
    */
    struct reb_orbit o = reb_tools_particle_to_orbit(r->G, r->particles[planet_id], r->particles[0]);
    double m_star   = r->particles[0].m;
    double m_planet = r->particles[planet_id].m;
    return o.a*(1.-o.e)*pow(m_planet/(3.*m_star), 1./3.);
}

struct reb_simulation* init_sim(int sim_id){
    /*
    Initialise simulation with ID sim_id
    */
    struct reb_simulation* r = reb_create_simulation();

    // Setup constants
    r->sim_id               = sim_id;
    r->collision            = REB_COLLISION_DIRECT;
    r->collision_resolve    = reb_collision_resolve_halt;
    r->dt                   = pow(65., .5)*2*M_PI/365.25; // Corresponds to ~8.062 days
    r->G                    = 1.;         // in AU^3 / SM / (year/2pi)^2
    r->ri_whfast.safe_mode  = 0;          // Turn off safe mode. Need to call reb_integrator_synchronize() before outputs.
    r->ri_whfast.corrector  = 11;         // 11th order symplectic corrector
    r->integrator           = REB_INTEGRATOR_WHFAST;
    r->heartbeat            = heartbeat;
    r->exact_finish_time    = 1;          // Finish exactly at tmax in reb_integrate(). Default is already 1.

    // Initial conditions
    for (int i=0;i<10;i++){
        struct reb_particle p = {0};
        p.x  = ss_pos[i][0];         p.y  = ss_pos[i][1];         p.z  = ss_pos[i][2];
        p.vx = ss_vel[i][0];         p.vy = ss_vel[i][1];         p.vz = ss_vel[i][2];
        p.m  = ss_mass[i];

        if (idx >= 1){p.r = r_hill(r, idx);}

        reb_add(r, p);
    }
    reb_move_to_com(r);
    return r;
}

void heartbeat(struct reb_simulation* r){
    if (reb_output_check(r, 5e5*2*M_PI)){         // Perturb Mercury x-coord every 0.5 Myr
      for (int idx = 1; idx < 10; idx++){         // Update Hill radius of all planets
        r->particles[idx].r = r_hill(r, idx);
      }

      double pert = 100*gaussian();
      r->particles[1].x += pert/AU;
      printf("\nPerturbed Mercury's x-coordinate by %f m\n", pert);
    }
    if (reb_output_check(r, 10000.)){           // Synchronise integrator
        reb_integrator_synchronize(r);
    }
    if (reb_output_check(r, 50000.)){
        struct reb_orbit merc_orb = reb_tools_particle_to_orbit(r->G, r->particles[1], r->particles[0]);

        char id_str[4];
        sprintf(id_str, "%d", r->sim_id);
        char filename[64];
        sprintf(filename, "sim_%s_ecc_bf.csv", id_str);
        FILE* fpt;

        fpt = fopen(filename, "a");
        fprintf(fpt, "%f, ", merc_orb.e);
        fclose(fpt);
    }
}

int main(int argc, char* argv[]){

    // Get inputs ==============================================================
    int N;
    if (argc == 1){
      N = 1;
    } else if (argc == 2){
      N = atoi(argv[1]);
    } else {
      printf("Incorrect input arguments: aborting\n");
      return 1;
    }

    double tmax = 2e8*2*M_PI;     // Final simulation time

    // Initialise simulations ==================================================
    struct reb_simulation** sims = malloc(N*sizeof(struct reb_simulation*));
    struct rebx_extras** rebx = malloc(N*sizeof(struct rebx_extras*));

    for (int i = 0; i < N; i++){
      printf("Initialising simulation %d\n", i+1);
      sims[i] = init_sim(i);

      rebx[i] = rebx_attach(sims[i]);
      struct rebx_force* gr = rebx_load_force(rebx[i], "gr");
      rebx_add_force(rebx[i], gr);
      rebx_set_param_double(rebx[i], &gr->ap, "c", 10065.32);   // Set speed of light in units AU/(yr/2pi)
    }

    // Integrate simulations ===================================================
    printf("============ Starting simulations ============");
#pragma omp parallel num_threads(8)
{
#pragma omp for
    for (int idx = 0; idx < N; idx++){
      printf("\n\nIntegrating simulation %d\n", idx+1);
      reb_integrate(sims[idx], tmax);
    }
}

    // Free simulations ========================================================
    for (int idx = 0; idx < N; idx++){
      printf("\nFreeing simulation %d", idx+1);

      rebx_free(rebx[idx]);
      reb_free_simulation(sims[idx]);
    }
    free(sims);

    printf("\n============ END ============\n");
    return 0;
}
