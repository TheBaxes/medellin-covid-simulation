#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "omp.h"
#include "grid.h"

//
//  benchmarking program
//
int main( int argc, char **argv )
{   
    
	
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h             to see this help\n" );
        printf( "-n <int>       to set the number of particles\n" );
        printf( "--hours <int>  to set the number of hours that will be simulated\n" );
        printf( "--stats        to specify if stats will be collected during run\n" );
        printf( "-o <filename>  to specify the output file name\n" );
        return 0;
    }

     int n = read_int( argc, argv, "-n", 100 );
    int PARAM_HOURS_TO_SIMULATE = read_int( argc, argv, "--hours", 1 );
    int NSTEPS = (int) floor((VEL_HUMAN/MOVEMENT_DELTA)*PARAM_HOURS_TO_SIMULATE*3600);
    
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    floydSolution floyd;
    std::vector<int> interesNodes;

    double simulation_time = read_timer( );
    init_particles( n, particles, floyd, interesNodes);
    simulation_time = read_timer( ) - simulation_time;
    printf("n = %d, init particles time = %g seconds\n", n, simulation_time);
    
    int nPOI = interesNodes.size();
    //int nNodes = sizeof(floyd.cost[0]) / sizeof(floyd.cost[0][0]); 
     
    char *savename = read_string( argc, argv, "-o", NULL );
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    if(fsave){
        // save params of the run     
        fprintf(fsave, "n_particles=%d,n_poi=%d,nsteps=%d,init_infect=%f,going_out=%f,prob_infection=%f\n", n, nPOI, NSTEPS, PARAM_RATE_INIT_INFECT, PARAM_GOING_OUT ,PARAM_PROB_INFECTION);
        fprintf(fsave,"ITERATIONS,N_INFECTED\n");
    }
    //Create grid
    int sizeGrid = (n/cutoff) +1; 
    grid_t grid;
    grid_init(grid, sizeGrid);
    for(int i=0; i < n; ++i){
        grid_add(grid, &particles[i]);
    }

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    #pragma omp parallel private(dmin) 
    {
        numthreads = omp_get_num_threads();
        for( int step = 0; step < 1000; step++ ){
        navg = 0;
        davg = 0.0;
	    dmin = 1.0;
        //
        //  compute all forces
        //
        #pragma omp for reduction (+:navg) reduction(+:davg)
        for( int i = 0; i < n; i++ ) 
        {
            particles[i].ax = particles[i].ay = 0;

            // 
            int gx = grid_coord(particles[i].x);
            int gy = grid_coord(particles[i].y);

            for(int x = max(gx - 1, 0); x <= min(gx + 1, sizeGrid-1); x++){
                for(int y = max(gy - 1, 0); y <= min(gy + 1, sizeGrid-1); y++){
                    linkedlist_t * curr = grid.grid[x * grid.size + y];
                    while(curr != 0)
                    {
                        apply_force(particles[i], *(curr->value),&dmin,&davg,&navg);
                        curr = curr->next;
                    }
                }
            }
        }
        
        //  move particles
        //

        #pragma omp for
        for( int i = 0; i < n; i++ ){
            int gc = grid_coord_flat(grid.size, particles[i].x, particles[i].y);
            move( particles[i] );

            //check if the particle chaged its position on the grid
            if(gc != grid_coord_flat(grid.size, particles[i].x, particles[i].y)){
                // then it is removed and added again
                if (! grid_remove(grid, &particles[i], gc))
                {
                    fprintf(stdout, "Error: Failed to remove particle '%p'.\n", &particles[i]);
                    exit(3);
                }
                grid_add(grid, &particles[i]);
            }

        } 
  
        if( find_option( argc, argv, "-no" ) == -1 ) 
        {
          //
          //  compute statistical data
          //
          #pragma omp master
          if (navg) { 
            absavg += davg/navg;
            nabsavg++;
          }

          #pragma omp critical
	      if (dmin < absmin) absmin = dmin; 
		
          //
          //  save if necessary
          //
          #pragma omp master
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
      }
    }

    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d,threads = %d, simulation time = %g seconds", n,numthreads, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -the minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");
    
    //
    // Printing summary data
    //
    if( fsum)
        fprintf(fsum,"%d %d %g\n",n,numthreads,simulation_time);

    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );

    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
