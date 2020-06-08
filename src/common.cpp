#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"
#include <iostream>
#include <vector>
#include <tuple>

unsigned int nNodes;

//
//  timer
//
double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
// Helper to create Array in C
//
float** createArray(int m, int n){
    float* values = (float*) calloc(m*n, sizeof(float));
    float** rows = (float**) malloc(n*sizeof(float*));
    for (int i=0; i<n; ++i)
    {
        rows[i] = values + i*m;
    }
    return rows;
}

int** createArrayInt(int m, int n){
    int* values = (int*) calloc(m*n, sizeof(int));
    int** rows = (int**) malloc(n*sizeof(int*));
    for (int i=0; i<n; ++i)
    {
        rows[i] = values + i*m;
    }
    return rows;
}

//
// Initialize the graph matrix
//
float** init_graph(){
    int u,v;
    float length;
    FILE *fp = fopen("graph.csv", "r");
    if(!fp){
        printf("ERROR: graph.csv not found!!\n");
        exit(1);
    }
    // Read number of nodes
    fscanf(fp, "%d", &nNodes);
    float **graph = createArray(nNodes, nNodes);
    // Read each line of graph.csv
    while(fscanf(fp, "%d, %d, %f", &u, &v, &length) != EOF){
        graph[u][v] = length;
        graph[v][u] = length;
        //printf("%d -> %d = %f\n", u, v, length);
    }
    return graph;
}

void getPathVector(int** path, std::vector<int> vector, int v,int u){
    if (path[v][u] == v)
		return;
    getPathVector(path, vector, v, path[v][u]);
    vector.push_back(path[v][u]);
}

std::vector<int> getPath(int** path, int v, int u){
    std::vector<int> vector;
    getPathVector(path, vector,v, u);
    return vector;
}

// Recursive Function to print path of given
// vertex u from source vertex v
void printPath(int** path, int v, int u)
{
	if (path[v][u] == v)
		return;

	printPath(path, v, path[v][u]);
	std::cout << path[v][u] << " ";
}

// Function to print the shortest cost with path
// information between all pairs of vertices
void printSolution(float** cost, int** path)
{
	for (int v = 0; v < nNodes; v++)
	{
		for (int u = 0; u < nNodes; u++)
		{
			if (u != v && path[v][u] != -1)
			{
				std::cout << "Shortest Path from " << v << " -> " << u << " is (" << v << " ";
				printPath(path, v, u);
				std::cout << u << ")" << std::endl;
			}
		}
	}
}


//
// Floyd implementation to calculate distances for every node in the path
//
// Function to run Floyd-Warshell algorithm
floydSolution floydWarshell(float **adjMatrix)
{
	// cost[] and parent[] stores shortest-path
	// (shortest-cost/shortest route) information
    float **cost = createArray(nNodes, nNodes);
    int **path = createArrayInt(nNodes, nNodes);
	//int cost[nNodes][nNodes], path[nNodes][nNodes];
    
	// initialize cost[] and parent[]
	for (int v = 0; v < nNodes; v++)
	{
		for (int u = 0; u < nNodes; u++)
		{
			// initally cost would be same as weight
			// of the edge
			cost[v][u] = adjMatrix[v][u];

			if (v == u)
				path[v][u] = 0;
			else if (cost[v][u] != INFINITY)
				path[v][u] = v;
			else
				path[v][u] = -1;
		}
	}

	// run Floyd-Warshell
	for (int k = 0; k < nNodes; k++)
	{
		for (int v = 0; v < nNodes; v++)
		{
			for (int u = 0; u < nNodes; u++)
			{
				// If vertex k is on the shortest path from v to u,
				// then update the value of cost[v][u], path[v][u]

				if (cost[v][k] != INFINITY && cost[k][u] != INFINITY
					&& cost[v][k] + cost[k][u] < cost[v][u])
				{
					cost[v][u] = cost[v][k] + cost[k][u];
					path[v][u] = path[k][u];
				}
			}

			// if diagonal elements become negative, the
			// graph contains a negative weight cycle
			if (cost[v][v] < 0)
			{
				printf("ERROR: Negative distance found in during Floyd execution, exiting...\n");
				exit(1);
			}
		}
	}

	// Print the shortest path between all pairs of vertices
	//printSolution(cost, path);

    floydSolution solution;
    solution.cost = cost;
    solution.path = path;
    return solution;
}

std::tuple<std::vector<int>, std::vector<int> > init_nodes(){
    FILE *fp = fopen("nodes_info.csv", "r");
    if(!fp){
        printf("ERROR: nodes_info.csv not found!!\n");
        exit(1);
    }

    std::vector<int> startNodes;
    std::vector<int> interestNodes;
    
    unsigned int  id;
    unsigned int  isInteres;

    while(fscanf(fp, "%u, %u", &id, &isInteres) != EOF){
        if(isInteres){
            interestNodes.push_back(id);
        }
        else {
            startNodes.push_back(id);
        }
    } 
    //std::tuple<std::vector<int>, std::vector<int> > 
    return std::make_tuple(startNodes, interestNodes);
}

void particle_is_going_out(particle_t &p){
    if(((double) rand() / (RAND_MAX)) < PARAM_GOING_OUT*(p.socialDiscipline/100)){
            p.goingout = true;
    }
    else{
            p.goingout = false;
    }
}

//
//  Initialize the particle positions and velocities
//
void init_particles( int n, particle_t *p, floydSolution &solution,  std::vector<int> &interesNodes )
{    
    float **graph = init_graph();
   
    for(int i = 0; i<nNodes; ++i){
        for(int j=0; j<nNodes; ++j){
            if(graph[i][j] == 0 && i != j) graph[i][j] = INFINITY; 
        }
    }
    solution = floydWarshell(graph);
    
    std::tuple<std::vector<int>, std::vector<int> > nodes = init_nodes();
    std::vector<int> startNodes;
    std::tie(startNodes, interesNodes) = nodes;

    srand( (unsigned int)time( NULL ) );
        
    int nNodesInteres = interesNodes.size();
    int nStartNodes = nNodes - nNodesInteres;

     if(n>nStartNodes){
        printf("ERROR, number of persons must be less than the number of available nodes to start.\nBut %d persons > %d nodes\n", n, nStartNodes);
        exit(2);
    }


    int *startNodesShuffle = (int*)malloc( nStartNodes * sizeof(int) );
    for( int i = 0; i < nStartNodes; i++)
        startNodesShuffle[i] = i;

    int *interesNodesShuffle = (int*)malloc( nNodesInteres * sizeof(int) );
    for( int i = 0; i < nNodesInteres; i++)
        interesNodesShuffle[i] = i;
         

    for( int i = 0; i < n; ++i ) 
    {
        //
        //  Assing start nodes to person
        //  Ensure each person has a node but one node has at most one person
        //
        int j = rand()%(nStartNodes-i); 
        int k = startNodesShuffle[j];
        startNodesShuffle[j] = startNodesShuffle[nStartNodes-i-1];
        p[i].id_initial_node = startNodes[k];
        p[i].id_present_node = startNodes[k];
        //printf("Particle %d init in node %d\n", i,k);
        
        //
        //  Assing Interes nodes to person
        //  Ensure each person has a node of POI to go
        //
        //  First ensure each POI has a person assigned
        if(nNodesInteres-i > 0){
            j = rand()%(nNodesInteres-i); 
            k = interesNodesShuffle[j];
            interesNodesShuffle[j] = interesNodesShuffle[nNodesInteres-i-1];
            p[i].id_object_node = interesNodes[k];
        }
        // Once that is done, assing the rest of persons to POI uniformly ramdom
        else{
           k = rand()%nNodesInteres;
           p[i].id_object_node = interesNodes[k];
        }
          
        // Compute path
        // p[i].path = getPath(solution.path ,p[i].id_present_node, p[i].id_object_node);

        // Is person infected? 
        if(((double) rand() / (RAND_MAX)) < PARAM_RATE_INIT_INFECT){
            p[i].infected = true;
            //printf("Particle %d infected\n", i);
        }

        // Social discipline of the person, a number of 0 to 100 with uniform distribution
        float valDiscipline = rand()%100;
        if (valDiscipline>50.0){
           p[i].socialDiscipline = 95.0;
        }
        //else if(valDiscipline<30.0){
        //   p[i].socialDiscipline = 50.0;
        //}
        else{
          p[i].socialDiscipline = 5.0;
        }


        // Is the person going out of their node? 
        particle_is_going_out(p[i]);
        

    }
    //DEBUG
    free( startNodesShuffle );
    free( interesNodesShuffle );
     
}


//
//  interact two particles
//
void apply_interaction( particle_t &p_a, particle_t &p_b)
{
    //check if particles are in the same node, and thus could interact
    if( p_a.id_present_node != p_b.id_present_node)
        return;
    
    /// INFECTION
    // the only case we will check is if a is infected and b is not
    if( !p_a.infected || p_b.infected)
    return;

    double prob_infection_b = ((p_a.socialDiscipline*p_b.socialDiscipline)/10000)*PARAM_PROB_INFECTION;
    if(((double) rand() / (RAND_MAX)) < prob_infection_b){
            p_b.infected = true;
    }
    return;

}

void move( particle_t &p, floydSolution &floyd, std::vector<int> &interesNodes )
{   
    //Check if I wanna go out
    if(!p.goingout){
        particle_is_going_out(p);
        return;
    }
    // If the particle reach the object node
    if(p.id_present_node == p.id_object_node){
      if(p.id_present_node == p.id_initial_node){

         /// Tengo otro objetivo y se calcula de nuevo la probabilidad de salir     
         p.id_object_node = interesNodes[rand()%interesNodes.size()];
         // Is going out this time?
         particle_is_going_out(p);
         return;
       }
      
      // TODO: wait
      // New object is the initial
      p.id_object_node = p.id_initial_node;
    }

    // DO MOVE
    p.meters_to_move += MOVEMENT_DELTA;
    int next = floyd.path[p.id_object_node][p.id_present_node];
    int distance = floyd.cost[p.id_present_node][next];
    if (distance > p.meters_to_move){
       return;
    }
    p.meters_to_move -= distance;
    p.id_present_node = next;
    return; 
}

//
//  command line option processing
//
int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return atoi( argv[iplace+1] );
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}

