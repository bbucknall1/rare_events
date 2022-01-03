/**
 * Parallel implementation using OpenMP

 ******** REQUIRES MODIFIED rebound.h FILE WITH NEW SIMULATION PARAMETERS ********

 *
 *  TO DO:
 *    x Change planet locations to astronomical units
 *    x Implement normal distribution generator
 *    x Add Mercury perturbations to heartbeat function
 *    x Implement multiple simulations
 *    x **NO LONGER NEEDED** Arrays for ecc max and min
 *    x Update eccentricities at each heartbeat call?
 *    x **NO LONGER NEEDED** Array for particle weights
 *    x REWEIGHTING
 *    x Implement sorted stratified resampling method
 *    x Splitting and killing - use identity splitting function V(x) = theta(x) = eccentricity range or constant multiple?
 *        can't copy simulations in place!!!!
 *    x Copy rebx to copied simulations
 *    x Manually copy new simulation parameters (weights, etc...)
 *    x Hill radii
 *    o Custom collision detection function that outputs to file?
 *    x Output eccentricities to file?
 *    o Find a way to make weights diverge more
 *
 *    o Tidy up and optimise DMC part of main()
 *    o Make sim indexing consistent (start at 0 or 1)
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


void heartbeat(struct reb_simulation* r);
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

double r_hill(struct reb_simulation* r, int planet_id){
    struct reb_orbit o = reb_tools_particle_to_orbit(r->G, r->particles[planet_id], r->particles[0]);
    double m_star   = r->particles[0].m;
    double m_planet = r->particles[planet_id].m;
    return o.a*(1.-o.e)*pow(m_planet/(3.*m_star), 1./3.);
}

int collision_resolve_halt_print(struct reb_simulation* const r, struct reb_collision c){
    printf("Collision occurred in simulation %d!\tHalting!\n", r->sim_id);

    char id_str[4];
    sprintf(id_str, "%d", r->sim_id);

    char filename[64];
    sprintf(filename, "sim_%s_ecc_dmc.csv", id_str);

    FILE* fpt;
    fpt = fopen(filename, "a");

    fprintf(fpt, "Collision occurred between planets %d and %d. Halting\n", c.p1, c.p2);

    fclose(fpt);

    r->sim_weight = 0.;
    r->status = REB_EXIT_COLLISION;
    return 0;
}

struct reb_simulation* init_sim(int sim_id){
    struct reb_simulation* r = reb_create_simulation();
    // Setup constants
    r->sim_id = sim_id;
    r->sim_weight = 1.;     // 1 = exp(0) .... eq (4)
    r->prev_V = 0.;

    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = collision_resolve_halt_print;

    r->dt             = pow(65., .5)*2*M_PI/365.25; // Corresponds to ~8.062 days
    //tmax              = 5e6*2*M_PI;            // 5 Myr
    r->G              = 1.;               // in AU^3 / SM / (year/2pi)^2
    r->ri_whfast.safe_mode     = 0;        // Turn off safe mode. Need to call reb_integrator_synchronize() before outputs.
    r->ri_whfast.corrector     = 11;        // 11th order symplectic corrector
    r->integrator        = REB_INTEGRATOR_WHFAST;
    r->heartbeat        = heartbeat;
    r->exact_finish_time = 1; // Finish exactly at tmax in reb_integrate(). Default is already 1.
    //r->integrator        = REB_INTEGRATOR_IAS15;        // Alternative non-symplectic integrator

    // Initial conditions
      for (int idx = 0; idx < 10; idx++){
          struct reb_particle p = {0};
          p.x  = ss_pos[idx][0];         p.y  = ss_pos[idx][1];         p.z  = ss_pos[idx][2];
          p.vx = ss_vel[idx][0];         p.vy = ss_vel[idx][1];         p.vz = ss_vel[idx][2];
          p.m  = ss_mass[idx];

          if (idx >= 1){p.r = r_hill(r, idx);}

          reb_add(r, p);
      }

    reb_move_to_com(r);

    return r;
}

void heartbeat(struct reb_simulation* r){
    if (reb_output_check(r, 1e5*2*M_PI)){        // TEMP // Update Hill radii and perturb Mercury x-coord every 10 Myr
      for (int idx = 1; idx < 10; idx++){
        r->particles[idx].r = r_hill(r, idx);
      }

      //double pert = 0.38*gaussian();
      double pert = 100*gaussian();
      r->particles[1].x += pert/AU;
      printf("\nPerturbed Mercury's x-coordinate by %f m\n", pert);
    }

    if (reb_output_check(r, 10000.)){           // Display (default heartbeat function)
        //reb_output_timing(r, tmax);
        reb_integrator_synchronize(r);

        // Update max and min eccentricities
        struct reb_orbit merc_orb = reb_tools_particle_to_orbit(r->G, r->particles[1], r->particles[0]);
        if (merc_orb.e > r->merc_ecc_max){
          r->merc_ecc_max = merc_orb.e;
        }
        if (merc_orb.e < r->merc_ecc_min){
          r->merc_ecc_min = merc_orb.e;
        }
    }

    if (reb_output_check(r, 50000.)){
        struct reb_orbit merc_orb = reb_tools_particle_to_orbit(r->G, r->particles[1], r->particles[0]);
        char id_str[4];
        sprintf(id_str, "%d", r->sim_id);

        char filename[64];
        sprintf(filename, "sim_%s_ecc_dmc.csv", id_str);

        FILE* fpt;
        fpt = fopen(filename, "a");

        fprintf(fpt, "%f, ", merc_orb.e);

        fclose(fpt);
    }
}

void sort_sims(double* thetas, struct reb_simulation** sims, int N){
    // Use insertion sort to order simulations in increasing values of theta (eccentricity range)
    int j;
    double temp_theta;
    struct reb_simulation* temp_sim;

    for (int i = 1; i < N; i++){
      temp_theta = thetas[i];
      temp_sim   = sims[i];
      j = i - 1;
      while (j >= 0 && thetas[j] > temp_theta){
        thetas[j+1] = thetas[j];
        sims[j+1]   = sims[j];
        j--;
      }
      thetas[j+1] = temp_theta;
      sims[j+1]   = temp_sim;
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

    // Initialise simulations ==================================================
    struct reb_simulation** sims = malloc(N*sizeof(struct reb_simulation*));
    struct rebx_extras** rebx = malloc(N*sizeof(struct rebx_extras*));

    double avg_weight = 1.;

    struct reb_orbit merc_orb;

    for (int i = 0; i < N; i++){
      printf("Initialising simulation %d\n", i+1);
      sims[i] = init_sim(i);

      // Get initial values of Mercury's eccentricity in each simulation
      merc_orb = reb_tools_particle_to_orbit(sims[i]->G, sims[i]->particles[1], sims[i]->particles[0]);
      sims[i]->merc_ecc_max = merc_orb.e;
      sims[i]->merc_ecc_min = merc_orb.e;

      rebx[i] = rebx_attach(sims[i]);
      // Could also add "gr" or "gr_full" here.  See documentation for details.
      struct rebx_force* gr = rebx_load_force(rebx[i], "gr");
      rebx_add_force(rebx[i], gr);
      // Have to set speed of light in right units (set by G & initial conditions).  Here we use default units of AU/(yr/2pi)
      rebx_set_param_double(rebx[i], &gr->ap, "c", 10065.32);
    }
    printf("============ Starting simulations ============");

    // Integrate simulations ===================================================
    double times[5] = {2e5*2*M_PI, 4e8*2*M_PI, 6e8*2*M_PI, 8e8*2*M_PI, 10e8*2*M_PI};   // TEMP  // Max time 1Gyr

    double resampling_bnds[N];
    double total_sum_weights;
    double partial_sum_weights;

    for (int i = 0; i < 5; i++){                  // i is resampling iteration
      // ======================== Integrate simulations ========================
#pragma omp parallel num_threads(8)
{
#pragma omp for
      for (int idx = 0; idx < N; idx++){            // loop over simulations
        if (sims[idx]->status != REB_EXIT_COLLISION){
          int thread_id = omp_get_thread_num();
          printf("\n\nIntegrating simulation %d on thread %d until resampling time %d\n", sims[idx]->sim_id + 1, thread_id, i+1);
          reb_integrate(sims[idx], times[i]);
        }
      }
}
      printf("\nAll simulations are now at time %f\n", times[i]);

      // =========================== 2a: Reweighting ===========================
      // Create array of 'thetas' (Eccentricity range)
      double theta;
      double thetas[N];
      double new_V;
      double new_weights[N];

      // Sum total weights
      total_sum_weights = 0.;
      for (int idx = 0; idx < N; idx++){
        total_sum_weights += sims[idx]->sim_weight;
      }
      avg_weight = total_sum_weights/N;

      for (int idx = 0; idx < N; idx++){
        theta = sims[idx]->merc_ecc_max - sims[idx]->merc_ecc_min;
        new_V = 10*theta;
        new_weights[idx] = avg_weight*exp(new_V - sims[idx]->prev_V);     // eq (5)
        sims[idx]->prev_V = new_V;
        thetas[idx] = theta;
        printf("Simulation %d had max and min eccentricity of (%f, %f), so theta = %f\n", sims[idx]->sim_id+1, sims[idx]->merc_ecc_max, sims[idx]->merc_ecc_min, theta);
      }

      for (int idx = 0; idx < N; idx++){
        sims[idx]->sim_weight = new_weights[idx];
      }

      // =========================== 2b: Resampling ============================

      // Sort sims into increasing thetas
      sort_sims(thetas, sims, N);

      // Relabel simulations based on their new ordering
      for (int idx = 0; idx < N; idx++){
        //sims[idx]->sim_id = idx;
        printf("(After sorting) Sim %d now has theta %f and weight %f\n", sims[idx]->sim_id, thetas[idx], sims[idx]->sim_weight);
      }

      // Sum total weights
      total_sum_weights = 0.;
      for (int idx = 0; idx < N; idx++){
        total_sum_weights += sims[idx]->sim_weight;
      }
      avg_weight = total_sum_weights/N;

      // Find bounds for sorted stratified resampling
      partial_sum_weights = 0.;
      for (int idx = 0; idx < N-1; idx++){
        partial_sum_weights += sims[idx]->sim_weight;
        resampling_bnds[idx] = partial_sum_weights/total_sum_weights;
        printf("Resampling bound %d is %f\n", idx+1, resampling_bnds[idx]);
      }
      resampling_bnds[N-1] = 1.;
      printf("Resampling bound %d is %f\n", N, resampling_bnds[N-1]);

      // ============================ Split & Kill =============================
      struct reb_simulation** sims_temp = malloc(N*sizeof(struct reb_simulation*));
      for (int j = 0; j < N; j++){
        double Qarg = ((double) j + ((double) rand()/(double) RAND_MAX))/((double) N);
        printf("j = %d\tQarg = %f\n", j, Qarg);

        for (int idx = 0; idx < N; idx++){
          if (Qarg < 1.){
            sims_temp[j] = reb_copy_simulation(sims[idx]);
            // Reset function pointers
            sims_temp[j]->heartbeat = heartbeat;
            sims_temp[j]->sim_id = j;
            sims_temp[j]->sim_weight = sims[idx]->sim_weight;
            sims_temp[j]->prev_V = sims[idx]->prev_V;
            printf("Simulation %d is now a copy of simulation %d\n", sims[j]->sim_id, idx);
            break;
          }
        }
      }

      // Then reset max and min eccentricites to current
      for (int idx = 0; idx < N; idx++){
        sims[idx] = reb_copy_simulation(sims_temp[idx]);
        // Reset function pointers and custom parameters
        sims[idx]->heartbeat = heartbeat;
        sims[idx]->sim_id = sims_temp[idx]->sim_id;                // sim_id keeps track of starting id.
        sims[idx]->sim_weight = sims_temp[idx]->sim_weight;
        sims[idx]->prev_V = sims_temp[idx]->prev_V;

        // Reattach reboundx
        rebx[idx] = rebx_attach(sims[idx]);
        // Could also add "gr" or "gr_full" here.  See documentation for details.
        struct rebx_force* gr = rebx_load_force(rebx[idx], "gr");
        rebx_add_force(rebx[idx], gr);
        // Have to set speed of light in right units (set by G & initial conditions).  Here we use default units of AU/(yr/2pi)
        rebx_set_param_double(rebx[idx], &gr->ap, "c", 10065.32);

        // Free temporary simulation
        reb_free_simulation(sims_temp[idx]);

        // Reset max and min eccentricies
        merc_orb = reb_tools_particle_to_orbit(sims[idx]->G, sims[idx]->particles[1], sims[idx]->particles[0]);
        sims[idx]->merc_ecc_max = merc_orb.e;
        sims[idx]->merc_ecc_min = merc_orb.e;
      }

      free(sims_temp);
    }

    // Free simulations ========================================================
    for (int i = 0; i < N; i++){
      printf("\nFreeing simulation %d", i+1);

      rebx_free(rebx[i]);
      reb_free_simulation(sims[i]);
    }
    free(sims);

    printf("\n ====END==== \n");
    return 0;
}
