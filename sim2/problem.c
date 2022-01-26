/**
 * sim2: Same algorithm as sim1, but using a toy solar system model

 ******** REQUIRES MODIFIED rebound.h FILE WITH NEW SIMULATION PARAMETERS ********

 *
 *  TO DO:
 *
 *    o Tidy up and optimise DMC part of main()
 *    x Edit copy_sim function to include custom parameters - or inline my own?
 *    x Check paper for details on how to retreive unbiased probability - Am I using the right weight?
 *    x In resampling: Selection of x-coordinate needs to exclude halted simulations
 *    x What about using average eccentricity as theta in this case?
 *        The other planets are also unstable, so it doesn't seem appropriate anymore to
*         focus solely on the innermost. Maybe we should also perturb the others?
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

void heartbeat(struct reb_simulation* r);
double tmax;

inline void my_copy_sim(struct reb_simulation* source, struct reb_simulation* dest){
    /*
    Copy simulation stored in source to destination.
    Required to copy function pointers and custom variables.
    */
    // Use REBOUNDs default copy function
    dest = reb_copy_simulation(source);
    // Set function pointers
    dest->heartbeat = heartbeat;
    dest->collision = REB_COLLISION_DIRECT;
    dest->collision_resolve = collision_resolve_halt_print;
    // Copy custom variables
    dest->sim_id = source->sim_id;
    dest->sim_weight = source->sim_weight;
    dest->prev_V = source->prev_V;
}

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

int collision_resolve_halt_print(struct reb_simulation* const r, struct reb_collision c){
    /*
    Custom collision resolutiomn function that halts the simulation, prints information about
    the collision to file, and sets simulation weight to 0 (ensuring that it does not get copied)
    */
    printf("Collision occurred in simulation %d!\tHalting!\n", r->sim_id);
    printf("Planets had Hill radii %f and %f\n", r->particles[c.p1].r, r->particles[c.p2].r);

    char id_str[4];
    sprintf(id_str, "%d", r->sim_id);

    char filename[64];
    sprintf(filename, "sim_%s_ecc_dmc.csv", id_str);

    FILE* fpt;
    fpt = fopen(filename, "a");

    fprintf(fpt, "Collision occurred between planets %d and %d. Halting.\nSimulation had weight %f\n", c.p1, c.p2, r->sim_weight);

    fclose(fpt);

    r->sim_weight = 0.;
    r->merc_ecc_max = 0.;
    r->merc_ecc_min = 0.;
    r->status = REB_EXIT_COLLISION;
    return 0;
}

struct reb_simulation* init_sim(int sim_id){
    struct reb_simulation* r = reb_create_simulation();
    // Setup constants
    r->sim_id               = sim_id;
    r->sim_weight           = 1.;       // 1 = exp(0) .... eq (4)
    r->prev_V               = 0.;
    r->collision            = REB_COLLISION_DIRECT;
    r->collision_resolve    = collision_resolve_halt_print;
    r->dt                   = pow(65., .5)*2*M_PI/365.25;   // Corresponds to ~8.062 days
    r->G                    = 1.;       // in AU^3 / SM / (year/2pi)^2
    r->ri_whfast.safe_mode  = 0;        // Turn off safe mode. Need to call reb_integrator_synchronize() before outputs.
    r->ri_whfast.corrector  = 11;       // 11th order symplectic corrector
    r->integrator           = REB_INTEGRATOR_WHFAST;
    r->heartbeat            = heartbeat;
    r->exact_finish_time    = 1;        // Finish exactly at tmax in reb_integrate(). Default is already 1.

    // Initial conditions - Model unstable system
    struct reb_particle star = {0};
    star.m = 1;
    star.r = 0;                            // Star is pointmass
    reb_add(r, star);

    for (int idx=0; idx<9; idx++){
        double a = 1e4*(1.+(double)idx/(double)(8));        // semi major axis
        double v = sqrt(1./a);                     // velocity (circular orbit)
        struct reb_particle planet = {0};
        planet.m = 1e-4;
        planet.r = 0.5*r_hill(r, idx+1);                    // Set planet radius to hill radius
                                    // A collision is recorded when planets get within their hill radius
                                    // The hill radius of the particles might change, so it should be recalculated after a while
        planet.lastcollision = 0;
        planet.x = a;
        planet.vy = v;
        reb_add(r, planet);
    }
    reb_move_to_com(r);
    return r;
}

void heartbeat(struct reb_simulation* r){
    if (reb_output_check(r, 1e6*2*M_PI)){         // Update Hill radii and perturb Mercury x-coord every 10 Myr
      for (int idx = 1; idx < 10; idx++){
        r->particles[idx].r = 0.5*r_hill(r, idx);
      }

      double pert = 100*gaussian();
      r->particles[1].x += pert;
      printf("\nPerturbed Mercury's x-coordinate by %f m\n", pert);
    }

    if (reb_output_check(r, 10000.)){           // Display (default heartbeat function)
        reb_output_timing(r, tmax);
        reb_integrator_synchronize(r);

        // Update max and min eccentricities
        struct reb_orbit orb;
        double ecc = 0.;
        for (int i = 1; i < 10; i++){
          orb = reb_tools_particle_to_orbit(r->G, r->particles[i], r->particles[0]);
          ecc += orb.e;
        }
        ecc = ecc/9.;
        if (ecc > r->merc_ecc_max){
          r->merc_ecc_max = ecc;
        }
        if (ecc < r->merc_ecc_min){
          r->merc_ecc_min = ecc;
        }
    }

    if (reb_output_check(r, 50000.)){
        struct reb_orbit orb;
        double ecc = 0.;
        for (int i = 1; i < 10; i++){
            orb = reb_tools_particle_to_orbit(r->G, r->particles[i], r->particles[0]);
            ecc += orb.e;
        }
    ecc = ecc/9.;

    char id_str[4];
    sprintf(id_str, "%d", r->sim_id);
    char filename[64];
    sprintf(filename, "sim_%s_ecc_dmc.csv", id_str);
    FILE* fpt;

    fpt = fopen(filename, "a");
    fprintf(fpt, "%f, ", ecc);
    fclose(fpt);
    }
}

void sort_sims(double* thetas, struct reb_simulation** sims, int N){
    /*
    Use insertion sort to order simulations in increasing values of theta (eccentricity range)
    */
    int j;
    double temp_theta;
    struct reb_simulation* temp_sim;

    for (int i = 1; i < N; i++){
      temp_theta = thetas[i];
      my_copy_sim(sims[i], temp_sim);           // copy sims[i] to temp_sim
      /*
      temp_sim = reb_copy_simulation(sims[i]);
      // Reset function pointers and custom variables
      temp_sim->heartbeat = heartbeat;
      temp_sim->sim_id = sims[i]->sim_id;
      temp_sim->sim_weight = sims[i]->sim_weight;
      temp_sim->prev_V = sims[i]->prev_V;
      temp_sim->collision = REB_COLLISION_DIRECT;
      temp_sim->collision_resolve = collision_resolve_halt_print;
      */

      j = i - 1;
      while (j >= 0 && thetas[j] > temp_theta){
        thetas[j+1] = thetas[j];
        my_copy_sim(sims[j], sims[j+1]);
        /*
        sims[j+1] = reb_copy_simulation(sims[j]);
        // Reset function pointers and custom variables
        sims[j+1]->heartbeat = heartbeat;
        sims[j+1]->sim_id = sims[j]->sim_id;
        sims[j+1]->sim_weight = sims[j]->sim_weight;
        sims[j+1]->prev_V = sims[j]->prev_V;
        sims[j+1]->collision = REB_COLLISION_DIRECT;
        sims[j+1]->collision_resolve = collision_resolve_halt_print;
        */

        j--;
      }
      thetas[j+1] = temp_theta;
      my_copy_sim(temp_sim, sims[j+1]);
      /*
      sims[j+1] = reb_copy_simulation(temp_sim);
      // Reset function pointers and custom variables
      sims[j+1]->heartbeat = heartbeat;
      sims[j+1]->sim_id = temp_sim->sim_id;
      sims[j+1]->sim_weight = temp_sim->sim_weight;
      sims[j+1]->prev_V = temp_sim->prev_V;
      sims[j+1]->collision = REB_COLLISION_DIRECT;
      sims[j+1]->collision_resolve = collision_resolve_halt_print;
      */
    }
}

int main(int argc, char* argv[]){

    // Get inputs ==============================================================
    int N;              // Number of simulations
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

    double avg_weight = 1.;

    for (int idx = 0; idx < N; idx++){
      printf("Initialising simulation %d\n", idx);
      sims[idx] = init_sim(idx);

      struct reb_orbit orb;
      double ecc = 0.;
      for (int i = 1; i < 10; i++){
        orb = reb_tools_particle_to_orbit(sims[idx]->G, sims[idx]->particles[i], sims[idx]->particles[0]);
        ecc += orb.e;
      }
      ecc = ecc/9.;
      sims[idx]->merc_ecc_max = ecc;
      sims[idx]->merc_ecc_min = ecc;
    }
    printf("============ Starting simulations ============");

    // Integrate simulations ===================================================
    double times[5] = {2e6*2*M_PI, 4e6*2*M_PI, 6e6*2*M_PI, 8e6*2*M_PI, 10e6*2*M_PI};     // Max time 0.1Gyr

    double resampling_bnds[N];
    double total_sum_weights;
    double partial_sum_weights;

    int num_halted;

    for (int i = 0; i < 5; i++){                  // i is resampling iteration
      // ======================== Integrate simulations ========================
      for (int idx = 0; idx < N; idx++){
        if (sims[idx]->status != REB_EXIT_COLLISION){
          printf("\n\nIntegrating simulation %d until resampling time %f\n", sims[idx]->sim_id, times[i]);
	        reb_integrate(sims[idx], times[i]);
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
        new_V = theta;
        new_weights[idx] = avg_weight*exp(new_V - sims[idx]->prev_V);     // eq (5)
        sims[idx]->prev_V = new_V;
        thetas[idx] = theta;
        printf("Simulation %d had max and min eccentricity of (%f, %f), so theta = %f\n", sims[idx]->sim_id, sims[idx]->merc_ecc_max, sims[idx]->merc_ecc_min, theta);
      }

      num_halted = 0;
      for (int idx = 0; idx < N; idx++){
        if (sims[idx]->status != REB_EXIT_COLLISION){
          sims[idx]->sim_weight = new_weights[idx];
        } else {
          num_halted++;
        }
      }

      // =========================== 2b: Resampling ============================

      // Sort sims into increasing thetas
      sort_sims(thetas, sims, N);
      for (int idx = 0; idx < N; idx++){
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
      for (int j = num_halted; j < N; j++){   // Exclude halted sims (since they have theta = 0 they are at beginning of array after sorting)
        double Qarg = ((double) (j-num_halted) + ((double) rand()/(double) RAND_MAX))/((double) N-num_halted);
        printf("j = %d\tQarg = %f\n", j, Qarg);

        for (int idx = 0; idx < N; idx++){
          if (Qarg < resampling_bnds[idx]){
              my_copy_sim(sims[idx], sims_temp[j]);
              sims_temp[j]->sim_id = sims[j]->sim_id;       // Retains starting id
              /*
              sims_temp[j] = reb_copy_simulation(sims[idx]);
              // Reset function pointers and custom variables, but keep sim_id
              sims_temp[j]->heartbeat = heartbeat;
              sims_temp[j]->sim_id = sims[j]->sim_id;         // Retains starting id
              sims_temp[j]->sim_weight = sims[idx]->sim_weight;
              sims_temp[j]->prev_V = sims[idx]->prev_V;
              sims_temp[j]->collision = REB_COLLISION_DIRECT;
              sims_temp[j]->collision_resolve = collision_resolve_halt_print;
              */
              printf("Simulation %d is now a copy of simulation %d\n", sims[j]->sim_id, sims[idx]->sim_id);
              break;
          }
        }
      }

      // Then reset max and min eccentricites to current
      for (int idx = num_halted; idx < N; idx++){
        if (sims[idx]->status != REB_EXIT_COLLISION){     // Only copy if not halted
            my_copy_sim(sims_temp[idx], sims[idx]);
            /*
            sims[idx] = reb_copy_simulation(sims_temp[idx]);
            // Reset function pointers and custom parameters
            sims[idx]->heartbeat = heartbeat;
            sims[idx]->sim_id = sims_temp[idx]->sim_id;                // sim_id keeps track of starting id.
            sims[idx]->sim_weight = sims_temp[idx]->sim_weight;
            sims[idx]->prev_V = sims_temp[idx]->prev_V;
            sims[idx]->collision = REB_COLLISION_DIRECT;
            sims[idx]->collision_resolve = collision_resolve_halt_print;
            */
        }
        // Free temporary simulation
        reb_free_simulation(sims_temp[idx]);

        struct reb_orbit orb;
        double ecc = 0.;
        for (int i = 1; i < 10; i++){
          orb = reb_tools_particle_to_orbit(sims[idx]->G, sims[idx]->particles[i], sims[idx]->particles[0]);
          ecc += orb.e;
        }
        ecc = ecc/9.;
        sims[idx]->merc_ecc_max = ecc;
        sims[idx]->merc_ecc_min = ecc;
      }

      free(sims_temp);
    }

    // Free simulations ========================================================
    for (int i = 0; i < N; i++){
      printf("\nFreeing simulation %d", i);

      reb_free_simulation(sims[i]);
    }
    free(sims);

    printf("\n ====END==== \n");
    return 0;
}
