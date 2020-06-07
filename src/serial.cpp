#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"

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

    /**
    int e = 11;
    printf("Particle %d is going %d -> %d\nPATH: ", e, particles[e].id_initial_node, particles[e].id_object_node );
    printPath(floyd.path, particles[e].id_initial_node, particles[e].id_object_node);
    printf("\nInformation:\n");
    printf("Infection: %d\nGoingOut: %d\nSocial Discipline is: %f\n", particles[e].infected, particles[e].goingout, particles[e].socialDiscipline);
    printf("\n\n");
    printf("Current Node: %d, next_node: %d, meter_to_move: %f\n", particles[e].id_present_node, floyd.path[particles[e].id_object_node][particles[e].id_present_node], particles[e].meters_to_move);
    **/
    //
    //  simulate a number of time steps
    //
    simulation_time = read_timer( );
	
    for( int step = 0; step < NSTEPS; step++ )
    {
        //
        // Compute infection
        for( int i = 0; i < n; ++i )
        {
            for (int j = 0; j < n; ++j )
				apply_interaction( particles[i], particles[j]);
        }
 
        //
        //  move particles
        //
        for( int i = 0; i < n; i++ ) 
            move( particles[i], floyd, interesNodes );		
        

        //printf("Current Node: %d, next_node: %d, meter_to_move: %f\n", particles[e].id_present_node, floyd.path[particles[e].id_object_node][particles[e].id_present_node], particles[e].meters_to_move);


        if( find_option( argc, argv, "--stats" ) != -1 )
        {
          //
          // Computing statistical data
          //
          if((step%SAVEFREQ) == 0){
            int count_infected = 0;
            for(int i = 0; i < n; ++i){
                if(particles[i].infected) count_infected++;
            }   
            if(fsave) fprintf(fsave,"%d,%d\n", step, count_infected);
            else{
                printf("ITER %d - Infected: %d\n", step, count_infected);
            }
          }
        }
    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf("n = %d, simulation time = %g seconds\n", n, simulation_time);
    //printf("Particle %d is from %d to %d\nInfection: %d\nGoingOut: %d\nSocial Discipline is: %f\n", e, particles[e].id_present_node, particles[e].id_object_node, particles[e].infected, particles[e].goingout, particles[e].socialDiscipline); 

    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
