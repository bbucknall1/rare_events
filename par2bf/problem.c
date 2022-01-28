/**
 * par2bf: Parallel implementation of model in sim2, but parallel and no DMC

 ******** REQUIRES MODIFIED rebound.h FILE WITH NEW SIMULATION PARAMETERS ********

 *
 *  TO DO:
 *
 *    o Tidy up and optimise DMC part of main()
 *    o Edit copy_sim function to include custom parameters - or inline my own?
 *    x Check paper for details on how to retreive unbiased probability - Am I using the right weight?
 *    x In resampling: Selection of x-coordinate needs to exclude halted simulations
 *    o What about using average eccentricity as theta in this case?
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

#define SOLAR_MASS 1.988544e30     // Solar Mass in kg
#define AU 149597870700         // Astronomical Unit in m

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
    r->sim_id = sim_id;

    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = collision_resolve_halt_print;

    r->dt             = pow(65., .5)*2*M_PI/365.25; // Corresponds to ~8.062 days
    r->G              = 1.;               // in AU^3 / SM / (year/2pi)^2
    r->ri_whfast.safe_mode     = 0;        // Turn off safe mode. Need to call reb_integrator_synchronize() before outputs.
    r->ri_whfast.corrector     = 11;        // 11th order symplectic corrector
    r->integrator        = REB_INTEGRATOR_WHFAST;
    r->heartbeat        = heartbeat;
    r->exact_finish_time = 1; // Finish exactly at tmax in reb_integrate(). Default is already 1.

    // Initial conditions - Model unstable system
    struct reb_particle star = {0};
    star.m = 1;
    star.r = 0;                            // Star is pointmass
    reb_add(r, star);

    for (int idx=0; idx<9; idx++){
        double a = 6e4*(1.+(double)idx/(double)(8));        // semi major axis
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

      double pert = 1000*gaussian();
      r->particles[1].x += pert;
      printf("\nPerturbed Mercury's x-coordinate by %f m\n", pert);
    }

    if (reb_output_check(r, 10000.)){           // Display (default heartbeat function)
        //reb_output_timing(r, tmax);
        reb_integrator_synchronize(r);

        if (reb_output_check(r, 50000.)){
          struct reb_orbit orb;
          double ecc = 0.;
          for (int i = 1; i < 10; i++){
            orb = reb_tools_particle_to_orbit(r->G, r->particles[i], r->particles[0]);
            ecc += orb.e;
          }
          ecc = ecc/9.;
          //struct reb_orbit merc_orb = reb_tools_particle_to_orbit(r->G, r->particles[3], r->particles[0]);
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

    srand((unsigned int)time(NULL));

    // Initialise simulations ==================================================
    struct reb_simulation** sims = malloc(N*sizeof(struct reb_simulation*));

    for (int idx = 0; idx < N; idx++){
      printf("Initialising simulation %d\n", idx);
      sims[idx] = init_sim(idx);
    }
    printf("============ Starting simulations ============");

    // Integrate simulations ===================================================
    double tmax = 1e8*2*M_PI;     // Max time 0.1Gyr

      // ======================== Integrate simulations ========================
    for (int idx = 0; idx < N; idx++){
      if (sims[idx]->status != REB_EXIT_COLLISION){
        printf("\n\nIntegrating simulation %d\n", sims[idx]->sim_id);
	      reb_integrate(sims[idx], tmax);
        printf("\nSimulation %d complete!\n", sims[idx]->sim_id);
      }
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
