#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__
#include <vector>
#include <tuple>
inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

//
//  saving parameters
//
const float VEL_HUMAN = 1.39; // 1.38 m/s
const float MOVEMENT_DELTA = 80; //80 Meters -> the length of a block
const int SAVEFREQ = 5;
const double PARAM_RATE_INIT_INFECT = 0.05; // ratio between infected and not infected in the beginning of the simulation
const double PARAM_GOING_OUT = 0.4; // prob a given particle is going out
const double PARAM_PROB_INFECTION = 0.01; // prob of infection in a interaction

//
// Solution of floyd algorithm
//
typedef struct 
{
  float **cost;
  int **path;
} floydSolution;

//
// particle data structure: a person
//
typedef struct 
{ 
  // Moving related vars
  int id_present_node;
  int id_initial_node;
  int id_object_node;
  double meters_to_move;
  // Virus related vars
  bool infected;
  bool goingout;
  double socialDiscipline;
} particle_t;

//
//  timing routines
//
double read_timer( );
float** createArray(int m, int n);

//
//  simulation routines
//
void init_particles( int n, particle_t *p, floydSolution &solution,  std::vector<int> &interesNodes );
void apply_interaction( particle_t &p_a, particle_t &p_b);
void move( particle_t &p, floydSolution &floyd, std::vector<int> &interesNodes );

// graph stuff
void floydWarshell(int graph);
float** init_graph();
std::vector<int> getPath(float **path, int v, int u);
std::tuple<std::vector<int>, std::vector<int> > init_nodes();

void printPath(int** path, int v, int u);

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#endif
