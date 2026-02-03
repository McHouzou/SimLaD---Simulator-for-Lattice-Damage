/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************* COMPILE COMMAND ************ */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
//clang++ -O3 -march=native -fopenmp -o lf_exec lattice_fire.cpp



/*
TODO:
- Tidy up input scheme so that all geometric parameters are fed through input files.
*/


#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <unistd.h>
#include <time.h>
#include <random>
#include <unordered_map>
#include <cstdlib>
#include <cstdio>
#include <omp.h>


// ---- Compile-time constants and convenience macros ----
# define DIM 2               // Spatial dimension (2D)
# define PI M_PI             // Use standard library pi
# define _x 0                // Index for x component
# define _y 1                // Index for y component
//# define VRBS //Verbose run with messages for debugging

// ---- FIRE 2.0 parameters (relaxation algorithm) ----
const double alpha_init = 0.1; //0.1
const double dt_init = 0.0002; //0.002 = dt * 0.02
const double dt_min = dt_init*0.02;
const double dt_max = 15*dt_init; //10*dt_init
const double f_inc = 1.1;
const double f_dec = 0.1;
const double alpha_decay = 0.99;
const double dx_tol = 1e-2; //Minimum strain step for stretchToFracture function
//const double force_threshold = 0.0025;
const double force_threshold = 0.0025; //For micromechanical study need smallest error you can get.
const int max_iterations = 10000;
const int Npnegmax = 2000;

// ---- Geometrical parameters ----
//Triangular
const double aa = 3.8;    //Area of cross section (mm^2)
//const double ii = 0.01583; //Second moment of area (mm^4)
//const double yy = 0.5;     //Thickness of beams (mm)


using namespace std;

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ Program Classes ************* */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
// ---- Core data structures ----
// Node: position, velocity, force, boundary flags, and original index
struct node {
  public:
    double r[DIM];						// Current position
    double v[DIM];						// Velocity (FIRE integration variable)
    double f[DIM];						// Force accumulator (reset each step)
    int boundary[DIM];        //Boundary condition: 0 free, 1 fixed, 2 prescribed displacement
    int index;                //Index of node in sim.vnode vector
};
  // Bond: connectivity + axial parameters + bending stress bookkeeping
  struct bond {
  public:
    node *n1;                 //Pointer to node 1
    node *n2;                 //Pointer to node 2
    double k;                 //Stiffness
    double d0;                //Equilibrium length
    double fa;                //Axial stress
    double fb[2] = {0.0, 0.0};                //Bending stress contributions on each side
    double famax;             //Max axial stress (failure threshold)
    double fbmax;             //Max bending stress (failure threshold)
};
  // Bend: three-node angle constraint
  struct bend {
  public:
    node *n1;                 //Pointer to node 1
    node *n2;                 //POinter to node 2
    node *n3;                 //Pointer to node 3
    double k;                 //Stiffness
    double a0;                //Equilibrium angle (radians)
};
struct sim{
  /*Think about having a variable dx that gets doubled every time a double fracture is not
  reached, and is halved every time there's a double fracture. */
  public:
    double dx;                                //Strain-step (default small change in strain)
    double f_thresh, term_threshold;          //Threshold for FIRE termination, according to Leo, not FIRE paper.
    double kbe, ii, yy;								        //Bending Stiffness, second moment of area, and beam thickness - read through file.
    int n_nodes,n_extfs,n_bends,n_bonds;	    //Number of objects

    //Vectors of Nodes, Bonds and Bends.
    vector<node> vnode;
    vector<bond> vbond;
    vector<bend> vbend;

    //Variables for bending and axial energy calculation.
    double Eax, Ebe;
};
struct in{
  FILE *inp;
  string inputDir;
};
struct out{
  FILE *posF, *bondF, *systemE, *bforce;
  string outputDir;
  int outsteps;
};


/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************** SUBROUTINES *************** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

void stateHelp( ) {
/*
 * State the code usage
*/
	printf("Usage:\n");
	printf("\t-o\t[path to output data files]\tdefault='./'\n");
	printf("\t-i\t[path to input data files]\tdefault='./'\n");
	printf("\t-h\t[this help menu]\n\n");
}

double dotprod(const double *x, const double *y) {
/*
 * Returns the dot product of two vectors.
*/
	double r=0.0;
	for(int d=0; d<DIM; d++ ) r+=x[d]*y[d];
	return r;
}

double safe_acos(const double value) {
    if (value<=-1.0) return PI;
    else if (value>=1.0) return 0.0;
		else return acos(value);
}

bool all_same(vector<double> v, double tol){
  /*
  Return true if all of the elements are the same up to a certain tolerance. False otherwise
  */
  double first_el = v[0];
  for (int i = 0; i < v.size(); i++){
    if (abs(v[0] - v[i]) > tol)
      return 0;
  }
  return 1;
}

void find_max_indeces(vector<double> v, vector<int> &max_indeces){
  /*Returns the indeces to the 2 largest elements in the vector v
  This is intended to be used to find the largest two thetas.*/
  double max1 = -10, max2 = -10;
  int ind[2] = {-1,-1};

  for (int i = 0; i < v.size(); i++){
    if (abs(v[i]) > max2){
      if (abs(v[i]) > max1){
        max2 = max1;
        max1 = abs(v[i]);
        ind[1] = ind[0];
        ind[0] = i;
      }
      else{
        max2 = abs(v[i]);
        ind[1] = i;
      }
    }
  }
  max_indeces.push_back(ind[0]);
  max_indeces.push_back(ind[1]);
}

void copy_node_vectors(vector<node> &v1, vector<node> &v2){
  /*
  Take node vector 1 and copy it into vector node 2.
  Used when storing temporary sim before FIRE update in the stretchToFracture function.
  */
  v2.clear();
  for(int i = 0; i < v1.size(); i++){
    v2.push_back(v1[i]);
  }
}

void copy_ind_vectors(vector<int> &v1, vector<int> &v2){
  /*
  Take node vector 1 and copy it into vector 2.
  */
  v2.clear();
  for(int i = 0; i < v1.size(); i++){
    v2.push_back(v1[i]);
  }
}



/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ***************** INPUT ****************** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
void readSim(in *inf, out *outf, sim *sim){
  /*
  Reads simulation parameters from file sim.inp
  */
  string inFile = inf->inputDir + "sim.inp";
  char descriptive_str[1000];
  int arg_count;
  float f;

  inf->inp = fopen(inFile.c_str(), "r");
  if (!inf->inp){ //If file couldn't be opened
    printf("Error:\tFile %s could not be opened.\n", inFile.c_str());
    exit(EXIT_FAILURE);
  }
  //Read strain-step size
  arg_count = fscanf(inf->inp, "%f %s", &f, descriptive_str);
  sim->dx = f;
  //Read bending stiffness
  arg_count = fscanf(inf->inp, "%f %s", &f, descriptive_str);
  sim->kbe = f;
  //Read beam thickness
  arg_count = fscanf(inf->inp, "%f %s", &f, descriptive_str);
  sim->yy = f;
  //Read second moment of area
  arg_count = fscanf(inf->inp, "%f %s", &f, descriptive_str);
  sim->ii = f;
  /*
  //Read noise amplitude
  arg_count = fscanf(inf->inp, "%f %s", &f, descriptive_str);
  sim->noiseAmp = f;
  */
  fclose( inf->inp );
}

void readNodes(in *inf, sim *sim){
  /*
  Reads all the nodes from file nodes.inp
  */
  string inFile = inf->inputDir + "nodes.inp";
  char descriptive_str[1000];
  int arg_count, nnodes, boundary_boolx, boundary_booly;
  float x,y;

  inf->inp = fopen(inFile.c_str(), "r");
  if (!inf->inp){ //If file couldn't be opened
    printf("Error:\tFile %s could not be opened.\n", inFile.c_str());
    exit(EXIT_FAILURE);
  }

  //Read Number of Nodes
  arg_count = fscanf(inf->inp, "%d %s", &nnodes, descriptive_str);
  sim->n_nodes = nnodes;

  //Read Node Information and push into node vector
  for(int i = 0; i < nnodes; i++){
    arg_count = fscanf(inf->inp, "%f %f %d %d %s", &x, &y, &boundary_boolx, &boundary_booly, descriptive_str);
    node n; n.r[0] = x; n.r[1] = y; n.boundary[0] = boundary_boolx; n.boundary[1] = boundary_booly; n.index = i;
    sim->vnode.push_back(n);
    for(int d=0; d<DIM; d++) { //Initialize velocities and forces of nodes to zero.
			sim->vnode[i].v[d]=0.0;
			sim->vnode[i].f[d]=0.0;
		}
  }
  fclose(inf -> inp);
  sim->f_thresh = sim->n_nodes*pow(force_threshold,2);
  sim->term_threshold = 0.5; //Termination criterion for complete fracture (this is a stress value scaled by famax)
}

void readBonds(in *inf, sim *sim){
  /*
  Reads all the bonds from file bonds.inp
  */
  string inFile = inf->inputDir + "bonds.inp";
  char descriptive_str[1000];
  int arg_count, nbonds, in1, in2;
  float k, d0, famax, fbmax;

  inf->inp = fopen(inFile.c_str(), "r");
  if (!inf->inp){ //If file couldn't be opened
    printf("Error:\tFile %s could not be opened.\n", inFile.c_str());
    exit(EXIT_FAILURE);
  }

  //Read Number of Bonds
  arg_count = fscanf(inf->inp, "%d %s", &nbonds, descriptive_str);
  sim->n_bonds = nbonds;

  //Read Node Information and push into node vector
  for(int i = 0; i < nbonds; i++){
    arg_count = fscanf(inf->inp, "%d %d %f %f %f %f %s", &in1, &in2, &k, &d0, &famax, &fbmax, descriptive_str);
    bond b; b.n1 = &(sim->vnode[in1]); b.n2 = &(sim->vnode[in2]);
    b.k = k; b.d0 = d0; b.famax = famax; b.fbmax = fbmax;
    sim->vbond.push_back(b);
  }
  fclose(inf -> inp);
}

void readBends(in *inf, sim *sim){
  /*
  Reads all the bends from file bends.inp
  */
  string inFile = inf->inputDir + "bends.inp";
  char descriptive_str[1000];
  int arg_count, nbends, in1, in2, in3;
  float k, a0, fbmax;

  inf->inp = fopen(inFile.c_str(), "r");
  if (!inf->inp){ //If file couldn't be opened
    printf("Error:\tFile %s could not be opened.\n", inFile.c_str());
    exit(EXIT_FAILURE);
  }

  //Read Number of Bonds
  arg_count = fscanf(inf->inp, "%d %s", &nbends, descriptive_str);
  sim->n_bends = nbends;

  //Read Node Information and push into node vector
  for(int i = 0; i < nbends; i++){
    arg_count = fscanf(inf->inp, "%d %d %d %f %f %s", &in1, &in2, &in3, &k, &a0, descriptive_str);
    bend be; be.n1 = &(sim->vnode[in1]); be.n2 = &(sim->vnode[in2]); be.n3 = &(sim->vnode[in3]);
    be.k = k; be.a0 = a0*PI/180.0; //Bend equilibrium angle in radians
    sim->vbend.push_back(be);
  }
  fclose(inf -> inp);
}





/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ***************** OUTPUT ***************** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
void openOutput(FILE **outF, string outputDir, string fname){
  string strBuff = outputDir+fname;
  *outF = fopen(strBuff.c_str(), "w+");
  if(!*outF){ //file couldn't be opened
    printf("Error:\tFile '%s' could not be opened.\n", strBuff.c_str());
    exit(EXIT_FAILURE);
  }
}

void outputPos(FILE *outF, sim sim){
  /*
  Prints the instantaneous position for all nodes
  */
  fprintf(outF, "%d N_nodes\n", sim.n_nodes);
  //fprintf(outF, "%f\n", strain);
  for(int i = 0; i < sim.n_nodes; i++){
    fprintf(outF, "%f %f\n", sim.vnode[i].r[0], sim.vnode[i].r[1]); //Output x,y position of node
  }
}

void outputBonds(FILE *outF, sim &sim){
  /*
  Prints the instantaneous bond information where each row contains the indices of corresponding nodes,
  the axial and bending stresses in MPa.
  */
  //fprintf(outF, "%d\n", sim.n_bonds);
  int n1ind=0, n2ind=0;
  for(int i = 0; i < sim.n_bonds; i++){
    n1ind = sim.vbond[i].n1->index;
    n2ind = sim.vbond[i].n2->index;
    fprintf(outF, "%d %d %f %f %f\n", n1ind, n2ind, sim.vbond[i].fa/sim.vbond[i].famax, abs(sim.vbond[i].fb[0]/sim.vbond[i].fbmax), abs(sim.vbond[i].fb[1]/sim.vbond[i].fbmax)); //Node indeces for bond with index i
  }
  fprintf(outF, "\n");
}

void outputaxialStress(FILE *outF, double strain, sim sim){
  /*
  Prints the instantaneous axial and bending stresses in all bonds
  */
  double axial,bending, r1[DIM],r2[DIM];
  fprintf(outF, "%d\n", sim.n_bonds);
  fprintf(outF, "%f\n", strain);
  /*
  for(int i = 0; i<sim.n_bonds; i++){
    r1 = (sim.vbond[i].n1)->r; r2 = (sim.vbond[i].n2)->r;
    axial = sim.vbond[i].k*(sqrt((r1[0] - r2[0])^2 + (r1[1] - r2[1])^2) - sim.vbond[i].d0) //compression is taken to be negative
    fprintf(outF, "%f\n", axial);
  }
  */
}

void outputbendingStress(FILE *outF, double strain, sim sim){}

void outputEnergy(FILE *outF, sim sim){
  /*
  Output total energy in lattice, with axial first and bending after.
  */
  fprintf(outF, "%f %f Eaxial_Ebending\n", sim.Eax, sim.Ebe);
}

void outputForceDisp(FILE *outF, sim sim){
  /*
  Prints the instantaneous total force applied by the nodes of the lattice on the boundary (in the x and y directions separately)

  ASSUMES MODE I FRACTURE, SO DISPLACEMENTS ARE MEASURED IN Y
  */
  vector<int> index_vec;
  double ymax = 0.0, ymin = 1000000.0;    //Coordinates of lattice boundary - changes as it gets pulled. (mm)

  for(int i = 0; i < sim.n_nodes; i++){
    if (sim.vnode[i].boundary[1] != 0){ //If node is on y constraints (i.e. on the boundary)
      index_vec.push_back(i);
      if (sim.vnode[i].r[1] > ymax)  ymax = sim.vnode[i].r[1];
      if (sim.vnode[i].r[1] < ymin)  ymin = sim.vnode[i].r[1];
    }
  }
  fprintf(outF, "%f Lattice_Height\n", ymax - ymin);
  for(int j = 0; j < index_vec.size(); j++){
    fprintf(outF, "%d %f %f ind_Fx_Fy\n", sim.vnode[index_vec[j]].index, sim.vnode[index_vec[j]].f[0], sim.vnode[index_vec[j]].f[1]);
  }
}


/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ******* DYNAMICS/RELAXATION ************** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
void bondForce(bond b){
  /*
  Find/add force on each node connected to this bond.
  */
  double R12, r12[DIM], fcomponent;
  for(int d = 0; d<DIM; d++) r12[d] = b.n2->r[d] - b.n1->r[d];
  R12 = sqrt(dotprod(r12,r12));
  for(int d = 0; d<DIM; d++){
     fcomponent = b.k * (R12 - b.d0) * (r12[d]/R12);
     b.n1->f[d] += fcomponent;
     b.n2->f[d] -= fcomponent;
  }
}

void bendForce(bend be){
  /*
  Find/add force on each node of bend
  */
  double r12[DIM], r23[DIM], n12[DIM], n23[DIM], R12, R23, theta;
  for (int d = 0; d<DIM; d++){
    r12[d] = be.n2->r[d] - be.n1->r[d];
    r23[d] = be.n3->r[d] - be.n2->r[d];
  }
  R12 = sqrt(dotprod(r12, r12));
  R23 = sqrt(dotprod(r23, r23));
  for (int d = 0; d<DIM; d++){
    n12[d] = r12[d]/R12;
    n23[d] = r23[d]/R23;
  }

  //Keep and eye out for improved algorithms to find angle
  double c = dotprod(n12, n23);
  double n23ortho[DIM] = {n23[1], -1*n23[0]};
  double s = dotprod(n12, n23ortho);
  theta = atan2(s, c); //Sensitive to negative angles as well.

  //Calculate moment
  double moment = be.k * (be.a0 - theta); //M = -k * dtheta, where theta is measured between r12 and r23 a-clockwise

  //Update and add bending forces to nodes
  double n12ortho[DIM] = {-1*n12[1], n12[0]};
  n23ortho[0] = -1*n23[1]; n23ortho[1] = n23[0];
  double fcomponent1, fcomponent2;
  for(int d = 0; d<DIM; d++){
     fcomponent1 = moment * n12ortho[d] / R12;
     fcomponent2 = moment * n23ortho[d] / R23;
     be.n1->f[d] += fcomponent1;
     be.n3->f[d] += fcomponent2;
     be.n2->f[d] -= (fcomponent1 + fcomponent2);
  }
}

void computeForces(sim &sim){
    int n;
    //Zero all forces
    n = sim.n_nodes;
    for(int i = 0; i < n; i++){
      for(int d = 0; d < DIM; d++) sim.vnode[i].f[d] = 0.0; //Zero all forces
    }
    //Bonding Force Calculation
    n = sim.n_bonds;
    for(int i = 0; i < n; i++) bondForce(sim.vbond[i]);
    //Bending Force Calculation
    n = sim.n_bends;
    for(int i = 0; i < n; i++) bendForce(sim.vbend[i]);
}

void FIRE2(sim &sim){
  /* Improved FIRE algorithm, with FIRE 2.0 scheme and implicit euler implemented*/
  double fcomp, vcomp, Fmag, Vmag, Fvproj = 0.0;
  double alpha = alpha_init;
  double dt = dt_init;
  int i, Nppos = 0, Npneg = 0, n = sim.n_nodes;

  computeForces(sim); //Initial force calculation

  // Zero all velocities
  for (int j = 0; j < n; j++) {
      fill(begin(sim.vnode[j].v), end(sim.vnode[j].v), 0.0);
  }

  //MAIN LOOP
  for (i = 0; i < max_iterations; i++){
    if (Fvproj > 0){ //If there is improvement (power is positive, so sliding down hill)
      dt = min(dt * f_inc, dt_max);
      alpha *= alpha_decay;
    }
    else { //If moving uphill
      //printf("Going Backwards at iter = %d. Fmag = %f\n", i, Fmag);
      //dt = max(dt * f_dec, dt_min);
      dt = dt_init;
      alpha = alpha_init;

      //This corrects the position to 1/2 step, before zeroing v, if P < 0.
      for (int j = 0; j < n; j++) {
        auto &node = sim.vnode[j];
        for (int d = 0; d < DIM; d++) {
          if (node.boundary[d] == 0) { //Only move free nodes
            node.r[d] -= 0.5 * dt * node.v[d];
            node.v[d] = 0.0;
          }
        }
      }
    }

    //Calculate updated F, x, v, and projections and then check for convergence.
    //Compute dot(F,v), and use latter to determine if progress is being made
    Fmag = Vmag = 0.0;
    for (int j = 0; j < n; j++) {
      auto &node = sim.vnode[j];
      for (int d = 0; d < DIM; d++) {
        if (node.boundary[d] == 0) { //Only move free nodes
          double fcomp = node.f[d];
          node.v[d] += fcomp * dt; // Implicit Euler update
          double vcomp = node.v[d];
          Fmag += fcomp * fcomp;
          Vmag += vcomp * vcomp;
        }
      }
    }

    Fmag = sqrt(Fmag); Vmag = sqrt(Vmag);
    if (Fmag < sim.f_thresh) break; //Check for Convergence.
    //if (Fmag < force_threshold) break; //Check for Convergence.

    //CHANGE THIS TO INCLUDE OTHER NODE CONDITONALS - BOUNDARY MOTION, CONSTRAINT FOR EACH DOF.
    double one_minus_alpha = 1.0 - alpha;
    for (int j = 0; j < n; j++) {
      auto &node = sim.vnode[j];
      for (int d = 0; d < DIM; d++) {
        if (node.boundary[d] == 0) { //Only move free nodes
          node.v[d] = one_minus_alpha * node.v[d] + alpha * node.f[d] * (Vmag / Fmag);
          node.r[d] += dt * node.v[d];
        }
      }
    }

    computeForces(sim); //Compute forces at the end of each step, to be ready for next iteration.

    // Compute F·v projection
    Fvproj = 0.0;
    for (int j = 0; j < n; j++) {
      auto &node = sim.vnode[j];
      for (int d = 0; d < DIM; d++) {
        Fvproj += node.f[d] * node.v[d];
      }
    }
  }
  #ifdef VRBS
  printf("\tConverged after %d iterations. Fmag = %f\n", i, Fmag);
  #endif
}

void moveBoundary(sim &sim, int mode){
  /*
  Move the boundary of the lattice.
  mode 1: Affine deformation, where assumed Mode I loading - improved convergence.
  mode 2: Only nodes that are to be displaced (boundary_bool == 2) are moved. Can cater for either I or II (and beyond) 
          fracture, depending on input file.
  */
  int n = sim.n_nodes;
  double maxheight = 0.0, dx = sim.dx;
  switch (mode){

    //Affine deformation
    case 1:
      for(int i = 0; i<n; i++)
        if (sim.vnode[i].r[1] > maxheight)
          maxheight = sim.vnode[i].r[1];
      for(int i = 0; i<n; i++)
        sim.vnode[i].r[1] += dx * sim.vnode[i].r[1] / maxheight;
      break;

    //Only move dofs that are to be displaced (boundary flag = 2)
    case 2:
      //This boundary movement is correct, but not the affine one - slows down covergence by some sizeable fraction of 100%.
      for (int i=0; i<n; i++){
        for (int d=0; d<DIM; d++){
          if (sim.vnode[i].boundary[d] == 2){ //Only move nodes that are to be displaced
            sim.vnode[i].r[d] += dx;
          }
        }
      }
      break;

  }
}

void calcEnergy(sim &sim){
  /*
  Calculate the bending and axial energy of the lattice, which can then
  be written in a file using the outputEnergy routine above.

  save in sim.Eax and sim.Ebe variables.
  */

  //Zero energies
  sim.Eax = 0; sim.Ebe = 0;

  //These variables are used and refreshed in the following two for-loops
  double r12[DIM], r23[DIM], n12[DIM], n23[DIM], R12, R23, theta;

  //Axial Energy Calculation
  for (int i = 0; i < sim.n_bonds; i++){
    for(int d = 0; d<DIM; d++)
      r12[d] = sim.vbond[i].n2->r[d] - sim.vbond[i].n1->r[d];

    //Add the axial energy of each bond: 1/2 * k * dx**2
    sim.Eax += 0.5 * sim.vbond[i].k * pow((sqrt(dotprod(r12,r12)) - sim.vbond[i].d0), 2);
  }

  //Bending energy calculation
  for (int i = 0; i < sim.n_bends; i++){
    for (int d = 0; d<DIM; d++){
      r12[d] = sim.vbend[i].n2->r[d] - sim.vbend[i].n1->r[d];
      r23[d] = sim.vbend[i].n3->r[d] - sim.vbend[i].n2->r[d];
    }
    R12 = sqrt(dotprod(r12, r12));
    R23 = sqrt(dotprod(r23, r23));
    for (int d = 0; d<DIM; d++){
      n12[d] = r12[d]/R12;
      n23[d] = r23[d]/R23;
    }
    //Find bend angle
    double c = dotprod(n12, n23);
    double n23ortho[DIM] = {n23[1], -1*n23[0]};
    double s = dotprod(n12, n23ortho);
    theta = atan2(s, c); //Sensitive to negative angles as well.

    //Add the bending energy for each bend: 1/2 * k * dtheta**2
    sim.Ebe += 0.5 * sim.vbend[i].k * pow((theta - sim.vbend[i].a0), 2);
  }

}


/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ******* FRACTURE + TOPOLOGY CHANGES ****** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

double thetaOfBend(bend be){
  double r12[DIM], r23[DIM], n12[DIM], n23[DIM], R12, R23;

  for (int d = 0; d<DIM; d++){
    r12[d] = be.n2->r[d] - be.n1->r[d];
    r23[d] = be.n3->r[d] - be.n2->r[d];
  }
  R12 = sqrt(dotprod(r12, r12));
  R23 = sqrt(dotprod(r23, r23));
  for (int d = 0; d<DIM; d++){
    n12[d] = r12[d]/R12;
    n23[d] = r23[d]/R23;
  }

  //Keep an eye out for improved algorithms to find angle.
  double c = dotprod(n12, n23);
  double n23ortho[DIM] = {n23[1], -1*n23[0]};
  double s = dotprod(n12, n23ortho);
  return atan2(s, c); //Sensitive to negative angles as well.
}

void axialStress(bond &b){
  /*
  Add axial stresses to member corresponding to bond b
  */
  double R12, r12[DIM];
  for(int d = 0; d<DIM; d++) r12[d] = (b.n2)->r[d] - (b.n1)->r[d];
  R12 = sqrt(dotprod(r12,r12));
  b.fa = b.k * (R12 - b.d0) / aa; //Divide force by cross sectional area to get the stress.
}

void bendingStress(int i,  sim &sim){
  /*
  Add bending stresses to members corresponding to bend with index i
  */
  bend be = sim.vbend[i];
  //Calculate moment
  double moment = be.k * (thetaOfBend(be) - be.a0); //M = c * dtheta, where theta is measured between r12 and r23 anticlockwise.

  //Find which members to add this bending stress to. The ones that have end nodes (n1,n2) AND (n2,n3).
  int n = sim.n_bonds;
  for (int j = 0; j<n; j++){
    if (sim.vbond[j].n1 == be.n2){
      if (sim.vbond[j].n2 == be.n3)
        sim.vbond[j].fb[0] += 0.5 * moment * sim.yy / sim.ii; //Bending stress according to beam thickness and second moment of area
      else if (sim.vbond[j].n2 == be.n1)
        sim.vbond[j].fb[0] -= 0.5 * moment * sim.yy / sim.ii;
    }
    if (sim.vbond[j].n2 == be.n2){
      if (sim.vbond[j].n1 == be.n3)
        sim.vbond[j].fb[1] += 0.5 * moment * sim.yy / sim.ii; //Bending stress according to beam thickness and second moment of area
      else if (sim.vbond[j].n1 == be.n1)
        sim.vbond[j].fb[1] -= 0.5 * moment * sim.yy / sim.ii;
    }
  }
}

int checkCatastrophe(sim &sim){
  /*
  Check if the lattice has been broken fully (check if there is no resultant axial/bending stress in the members).
  */

  /*
  //Find the maximum stress in any bond, for debugging purposes.
  double max_stress = 0.0;
  for (int i = 0; i < sim.n_bonds; i++){
    double bond_stress = abs(sim.vbond[i].fa)/sim.vbond[i].famax + max(abs(sim.vbond[i].fb[0]), abs(sim.vbond[i].fb[1]))/sim.vbond[i].fbmax;
    if (bond_stress > max_stress){
      max_stress = bond_stress;
    }
  }
  printf("Maximum stress in lattice = %f\n", max_stress);
  */

  for (int i = 0; i < sim.n_bonds; i++){
    if (abs(sim.vbond[i].fa)/sim.vbond[i].famax + max(abs(sim.vbond[i].fb[0]), abs(sim.vbond[i].fb[1]))/sim.vbond[i].fbmax > sim.term_threshold){
      //printf("Stress exceeded threshold: σ = %f > %f\n", abs(sim.vbond[i].fa)/sim.vbond[i].famax + max(abs(sim.vbond[i].fb[0]), abs(sim.vbond[i].fb[1]))/sim.vbond[i].fbmax, sim.term_threshold);
      //printf("At bond %d\n", i);
      //printf("with end-nodes %d %d\n", sim.vbond[i].n1->index, sim.vbond[i].n2->index);
      return 0;
    }
  }

  /*
  //Or check that the boundary is sufficiently relaxed.
  for (int i = 0; i < sim.n_nodes; i++){
    if ((sim.vnode[i].boundary[1] != 0) && (abs(sim.vnode[i].f[1]) > sim.term_threshold)){ //If on the boundary the vertical force is larger than threshold:
      printf("Force on Boundary exceeded. F = %f\n", abs(sim.vnode[i].f[1]));
      return 0;
    }
  }
  */
  return 1;
}

int stressCheck(sim &sim, bool avalance) {
    /*
    Identifies node that has exceeded stress threshold if it's just one.
    Retunrs appropriate values to indicate no fracture or multiple ones, in which case
    re-relaxation is needed.

    Parameters:
    - sim: Reference to the simulation object.
    - avalance: If true, selects the node with the highest stress in case of multiple failures.

    Returns:
    - Index of the failed node.
    - -1 if no node has failed.
    - -2 if multiple nodes have failed and require further refinement.
    */

    vector<int> vbroken_ind;     // Indices of failed nodes
    vector<double> vbroken_stress; // Corresponding stress values

    // Reset all bond stresses before recalculating
    for (auto &bond : sim.vbond) {
        bond.fa = 0.0; // Axial stress reset
        bond.fb[0] = bond.fb[1] = 0.0; // Bending stress reset
    }

    // Compute axial and bending stress
    for (auto &bond : sim.vbond) axialStress(bond);
    for (int i = 0; i < sim.n_bends; i++) bendingStress(i, sim);

    // Identify failed nodes
    for (const auto &bond : sim.vbond) {
        for (int j = 0; j < 2; j++) { // Check both sides of the bond
            double stress_ratio = abs(bond.fa / bond.famax) + abs(bond.fb[j] / bond.fbmax);
            if (stress_ratio >= 1.0) {
                node* failed_node_ptr = (j == 0) ? bond.n1 : bond.n2;
                vbroken_ind.push_back(failed_node_ptr->index);
                vbroken_stress.push_back(stress_ratio);
            }
        }
    }

    if (vbroken_ind.empty()) return -1; // No nodes failed

    #ifdef VRBS
        printf("\t%ld Nodes have reached failure:\n", vbroken_ind.size());
        for (size_t i = 0; i < vbroken_ind.size(); i++)
            printf("\t-Node %d reached failure, stress ratio: %f\n", vbroken_ind[i], vbroken_stress[i]);
    #endif

    // Handle post-fracture scenario where we choose the most stressed node
    if (avalance) {
        int max_stress_index = max_element(vbroken_stress.begin(), vbroken_stress.end()) - vbroken_stress.begin();
        return vbroken_ind[max_stress_index];
    }

    // If only one node failed, determine next step
    if (vbroken_ind.size() == 1) return vbroken_ind[0];

    // If all stress values are nearly the same, pick a random node
    if (all_same(vbroken_stress, 1e-5)) {
        #ifdef VRBS
            printf("Found within tolerance, choosing one randomly.\n");
        #endif
        return vbroken_ind[rand() % vbroken_ind.size()];
    }

    // Requires further relaxation step if more than one reached fracture + avalance == False
    return -2;
}

int oldstretchToFracture(sim &sim){
  /*
  Identifies the first node that reaches the failure criterion under increasing strain.

  - Follow "de Waal et al. (2024)" for subsequent node splitting when failure occurs.
  - Saves a copy of the lattice if multiple failures require rollback and relaxation.

  Returns:
  - Index of the first node that breaks.
  */

  vector<node> temp_node_vec;

  int ind_broke = stressCheck(sim, false); // Initial stress check
  double original_dx = sim.dx; // Store original strain step size
  int iteration_count = 0;

  while (ind_broke < 0) {
      copy_node_vectors(sim.vnode, temp_node_vec); // Backup current state
      moveBoundary(sim, 1); // Apply strain increment
      FIRE2(sim); // Relaxation algorithm
      ind_broke = stressCheck(sim, false);

      if (ind_broke == -2) { // Multiple asymmetric nodes failed simultaneously
          sim.dx /= 2; // Reduce strain step size
          #ifdef VRBS
              printf("Multiple nodes reached failure at iteration %d. Reducing strain step to %.4f\n\n", iteration_count, sim.dx);
          #endif
          copy_node_vectors(temp_node_vec, sim.vnode); // Restore previous state
        }
        iteration_count++;
    }

    printf("\n-> Found breaking node (%d) at iteration %d\n", ind_broke, iteration_count);

    sim.dx = original_dx; // Restore original strain step size
    return ind_broke;
}

int stretchToFracture(sim &sim) {
    /*
    Finds the minimal displacement that breaks exactly one bond, approaching from above.
    Uses a binary-search style refinement on the applied strain.

    Returns:
    - Index of the first node that breaks.
    */

    vector<node> safe_state;       // last configuration before any node broke
    vector<node> failed_state;     // last known state with exactly one node broken

    int ind_broke = stressCheck(sim, false); // Initial stress check
    double original_dx = sim.dx; // Store original strain step size

    bool have_failed = false; // Flag to indicate if critical strain has been crossed
    int iteration_count = 0;

    copy_node_vectors(sim.vnode, safe_state); // initial safe state

    while (true) {
        // Start from the last safe state
        copy_node_vectors(safe_state, sim.vnode);
        moveBoundary(sim, 1); // Apply strain increment
        FIRE2(sim); // Relaxation algorithm
        ind_broke = stressCheck(sim, false);

        if (ind_broke == -2) {
            // Too far: multiple nodes broke → rollback (by not updating safe_state) and halve dx
            sim.dx *= 0.5;
        }
        else if (ind_broke < 0) {
            // No failure → update safe_state
            copy_node_vectors(sim.vnode, safe_state);
        }
        else if (ind_broke >= 0) {
            // Exactly one node broke
            copy_node_vectors(sim.vnode, failed_state);

            if (sim.dx <= dx_tol) {
                // Converged from above → done
                break;
            }
            // else: refine further in next iteration
            sim.dx *= 0.5;
        }
        iteration_count++;
    }

    printf("\n-> Found breaking node (%d), with dx = %f at iteration %d\n", ind_broke, sim.dx, iteration_count);

    sim.dx = original_dx; // Restore original strain step size
    copy_node_vectors(failed_state, sim.vnode); // Restore failed config.

    return ind_broke;
}

void broken_nodes(sim &sim, vector<int> &node_group1, vector<int> &node_group2, vector<int> &bend_ind, vector<int> max_indeces){
  /*
  Sort bend_ind, starting from max_bend. Clock or anti-clockwise does not matter.

  Fill up the node_group vector containing the group of nodes to be broken off.
  */
  vector<int> sorted_bend_ind; //Use temporarily and copy into bend_ind later.
  int current_bend_ind = bend_ind[max_indeces[0]];
  int current_node_ind = sim.vbend[current_bend_ind].n1->index;
  bool fill_group1 = true; //Determines whether group 1 or group 2 of nodes is being filled.

  node_group1.push_back(current_node_ind);
  sorted_bend_ind.push_back(current_bend_ind);

  //Sort bend_ind and fill up node_group vectors.
  int i = 0;
  while ((sorted_bend_ind.size() < bend_ind.size()) && (i < bend_ind.size())) {
      if (current_bend_ind != bend_ind[i]){
         if (sim.vbend[bend_ind[i]].n1->index == current_node_ind){
           current_bend_ind = bend_ind[i];
           sorted_bend_ind.push_back(current_bend_ind);
           current_node_ind = sim.vbend[bend_ind[i]].n3->index;
           if (fill_group1) node_group1.push_back(current_node_ind);
           else node_group2.push_back(current_node_ind);
           i = -1;
         }
         else if (sim.vbend[bend_ind[i]].n3->index == current_node_ind){
           current_bend_ind = bend_ind[i];
           sorted_bend_ind.push_back(current_bend_ind);
           current_node_ind = sim.vbend[bend_ind[i]].n1->index;
           if (fill_group1) node_group1.push_back(current_node_ind);
           else node_group2.push_back(current_node_ind);
           i = -1;
         }
      }
      if((current_bend_ind == bend_ind[max_indeces[1]]) && fill_group1){
        fill_group1 = false;
        node_group2.push_back(node_group1.back());
        node_group1.pop_back();
      }
      i++;
  }
  //if (node_group2.empty()) node_group2.push_back(sim.vbend[bend_ind[max_indeces[0]]].n3->index);
  bend_ind = sorted_bend_ind;
}

void changeTopology(sim &sim, int ind_broke){
  /*
    This function handles a topological change in the lattice when a node fails.
    It identifies affected bends, determines how to split the connections, and updates
    the node, bond, and bend data structures accordingly.
  */

  //Step 1: Find all bends that have a middle node that belongs to the above node. Store their indeces
  vector<int> bend_indeces;
  double dtheta;

  for (int i = 0; i<sim.n_bends; i++){ //Collect relevant bends to vectors corresponding to the bond's end-points.
    if (sim.vbend[i].n2->index == ind_broke)
      bend_indeces.push_back(i);
  }

  //Step 2: Calculate angular deviation (dtheta) for each affected bend
  vector<double> bend_dthetas;
  for (int i = 0; i<bend_indeces.size(); i++){
    dtheta = sim.vbend[bend_indeces[i]].a0 - thetaOfBend(sim.vbend[bend_indeces[i]]);
    bend_dthetas.push_back(dtheta);
    #ifdef VRBS
      printf("\tBend %d (%d %d %d) has Dtheta = %f\n", bend_indeces[i], sim.vbend[bend_indeces[i]].n1->index, sim.vbend[bend_indeces[i]].n2->index, sim.vbend[bend_indeces[i]].n3->index, dtheta);
    #endif
  }

  //Step 3: Create new node in sim.vnode and change already existing bonds that contain the two breakout nodes to point to it.
  node newn;
  for (int d = 0; d<DIM; d++){
    newn.r[d] = (sim.vbend[bend_indeces[0]].n2)->r[d];
    newn.f[d] = (sim.vbend[bend_indeces[0]].n2)->f[d];
    newn.v[d] = (sim.vbend[bend_indeces[0]].n2)->v[d];
  }
  //newn.boundary = (sim.vbend[bend_indeces[0]].n2)->boundary;
  newn.boundary[0] = 0; newn.boundary[0] = 0; //Assumed that new node is not a boundary node. THIS IS NOT GENERALLY TRUE AND DEPENDS ON HOW OTHER NODES ATTACH TO IT.
  newn.index = sim.n_nodes;
  sim.vnode.push_back(newn);
  sim.n_nodes++;




  //Step 4: Get indeces of bends to be "removed" - depending on the coordination number of the node, these are found differently
  vector<int> max_indeces;
  if (bend_indeces.size() <= 4){
    //IN CASE LATTICE IS NOT OVERCOORDINATED, FIND max_indeces ACCORDING TO bend_dthetas
    find_max_indeces(bend_dthetas, max_indeces); //Ideces in bend_theta vector of two maximum angles, fill up max_indeces
  }
  else{
    //FUNCTIONALITY FOR OVERCOORDINATED LATTICES - IF MORE THAN 4 BENDS IN BEND VECTOR, JUST CUT ONE NODE OFF.
    //USING EXISTING FUNCTIONALITY, THIS CAN BE DONE BY "FORCING" max_indeces TO BE ADJACENT TO BOND WITH MAX STRESS.
    vector<int> bond_indeces;
    vector<double> bond_stresses;
    //get all bonds that connect to node - stresses are already stored in the memory, so no need to re-calculate.
    //You do need to divide by the corresponding failure stress of each bond nonetheless - NEEDED IN CASE OF "NOISY" FAILURE STRESSES!!!
    for(int i = 0; i<sim.n_bonds; i++){
      if(ind_broke == sim.vbond[i].n1->index){
        bond_indeces.push_back(i);
        bond_stresses.push_back(abs(sim.vbond[i].fa / sim.vbond[i].famax) + abs(sim.vbond[i].fb[0] / sim.vbond[i].fbmax));
      }
      else if (ind_broke == sim.vbond[i].n2->index){
        bond_indeces.push_back(i);
        bond_stresses.push_back(abs(sim.vbond[i].fa / sim.vbond[i].famax) + abs(sim.vbond[i].fb[1] / sim.vbond[i].fbmax));
      }
    }
    //Find the index of the bond with maximum stress amongst the bonds connected to the breakout node.
    auto max_stress_iter = max_element(bond_stresses.begin(), bond_stresses.end());
    int max_bondind = distance(bond_stresses.begin(), max_stress_iter);
    int target_node_ind = (sim.vbond[bond_indeces[max_bondind]].n1->index == ind_broke)
                          ? sim.vbond[bond_indeces[max_bondind]].n2->index
                          : sim.vbond[bond_indeces[max_bondind]].n1->index;

    //Go over bends and find which bends have this target node as one of n1/n3. Their indeces shall be placed in the max_indeces array
    max_indeces.push_back(0); max_indeces.push_back(0);
    for(int i=0; i<bend_indeces.size(); i++){
      if((sim.vbend[bend_indeces[i]].n1->index == target_node_ind) || (sim.vbend[bend_indeces[i]].n3->index == target_node_ind)){
        max_indeces[1] = max_indeces[0];
        max_indeces[0] = i;
      }
    }
  }
  //printf("%d %d\n", max_indeces[0], max_indeces[1]);



  // Step 5: Split affected nodes into two groups
  // Find the nodes that should be broken off (connected to breakout node)
  // Go between the two max_indeces and get all of the bends that are to be broken.
  vector<int> node_group1, node_group2; //Indices of nodes to be broken off.
  broken_nodes(sim, node_group1, node_group2, bend_indeces, max_indeces);

  //for(auto i : node_group1) printf("Node group1: %d\n", i);
  //for(auto i : node_group2) printf("Node group2: %d\n", i); printf("Bend Indeces: ");
  //for(auto i : bend_indeces) printf("%d, ", i); printf("\n");


  //Step 6: Update bends and bonds of the nodes in node_group1 to point to the new node
  for (int j = 0; j < node_group1.size(); j++){
    //Change all relevant bends that have middle node as one of the breakout nodes.
    for (int k = 0; k < sim.n_bends; k++){
      if ((sim.vbend[k].n2->index == node_group1[j]) && (sim.vbend[k].n1 == sim.vbend[bend_indeces[0]].n2))
        sim.vbend[k].n1 = &(sim.vnode[sim.vnode.size()-1]);
      else if ((sim.vbend[k].n2->index == node_group1[j]) && (sim.vbend[k].n3 == sim.vbend[bend_indeces[0]].n2))
        sim.vbend[k].n3 = &(sim.vnode[sim.vnode.size()-1]);
    }

    //Change all bonds that end up on breakout nodes of group 1, to point to the new node.
    for(int i = 0; i < sim.n_bonds; i++){
      if ((sim.vbond[i].n1->index == node_group1[j]) && (sim.vbond[i].n2 == sim.vbend[bend_indeces[0]].n2)){
        #ifdef VRBS
          printf("\t\t-Changing bond %d (%d, %d), to point to last node.\n", i, sim.vbond[i].n1->index, sim.vbond[i].n2->index);
        #endif
        sim.vbond[i].n2 = &(sim.vnode[sim.vnode.size()-1]);
      }
      else if ((sim.vbond[i].n2->index == node_group1[j]) && (sim.vbond[i].n1 == sim.vbend[bend_indeces[0]].n2)){
        #ifdef VRBS
          printf("\t\t-Changing bond %d (%d, %d), to point to last node.\n", i, sim.vbond[i].n1->index, sim.vbond[i].n2->index);
        #endif
        sim.vbond[i].n1 = &(sim.vnode[sim.vnode.size()-1]);
      }
    }
  }


  /* ONE POSSIBLE SENARIO FOR ENDING TOPOLOGY CHANGE FUNCTION PREMATURELY:
  If there is only one bend assosciated with breakage, this is just an elbow to be broken.
  Reassign one node and assosciated bends (done above) and break void function prematurely.
  */
  if (bend_indeces.size() == 1){
    sim.vbend.erase(sim.vbend.begin() + bend_indeces[0]);
    #ifdef VRBS
      printf("\t\t-Removing Bend %d\n", bend_indeces[0]);
    #endif
    sim.n_bends--;
    return;
  }


  //Step 7: Change the bends that are not to be deleted, so that their middle nodes point to the new node.
  for (int j = 1; j < node_group1.size(); j++){ //If there are N nodes assosciated with new group, then there are N-1 bends
    sim.vbend[bend_indeces[j]].n2 = &(sim.vnode[sim.vnode.size()-1]);
    #ifdef VRBS
      printf("\t\t-Changing Bend %d to point to new middle node\n", bend_indeces[j]);
    #endif
  }


  //Step 8: if the fracture has left a group of 3 or more bonds together, a bend needs to be re-introduced.
  if (node_group1.size() >= 3){ //Create an extra bend wrapping around the exterior of the group of nodes.
    bend bnew;
    //Node pointers
    bnew.n1 = &(sim.vnode[node_group1[0]]);
    bnew.n2 = &(sim.vnode[sim.vnode.size()-1]);
    bnew.n3 = &(sim.vnode[node_group1.back()]);

    //Bend Angle
    bnew.a0 = 0;
    for (int j = 1; j < node_group1.size(); j++) bnew.a0 += (PI - abs(sim.vbend[bend_indeces[j]].a0));
    bnew.a0 = PI - bnew.a0;
    //printf("%f\n", 180*thetaOfBend(bnew)/PI);
    //bnew.a0 *= pow(-1, signbit(thetaOfBend(bnew)));
    bnew.a0 = abs(bnew.a0) * pow(-1, signbit(thetaOfBend(bnew)));

    //Stiffness of Bend
    double lengths[2];
    for(int i = 0; i<sim.n_bonds; i++){
      if (((sim.vbond[i].n1 == bnew.n1) && (sim.vbond[i].n2 == bnew.n2)) || ((sim.vbond[i].n2 == bnew.n1) && (sim.vbond[i].n1 == bnew.n2)))
        lengths[0] = sim.vbond[i].d0;
      else if (((sim.vbond[i].n1 == bnew.n3) && (sim.vbond[i].n2 == bnew.n2)) || ((sim.vbond[i].n2 == bnew.n3) && (sim.vbond[i].n1 == bnew.n2)))
        lengths[1] = sim.vbond[i].d0;
    }
    bnew.k = sim.kbe/min(lengths[0], lengths[1]);
    //printf("lengths = %f %f\n", lengths[0], lengths[1]);

    #ifdef VRBS
    printf("\t\t-1 Adding Bend (index %d), with angle %f and stiffness %f\n", sim.n_bends, bnew.a0*180/PI, bnew.k);
    printf("\t\t between nodes %d %d %d\n", bnew.n1->index, bnew.n2->index, bnew.n3->index);
    #endif
    sim.vbend.push_back(bnew);
    sim.n_bends++;
  }

  if (node_group2.size() >= 3){ //Same for group 2
    bend bnew;
    //Node pointers
    bnew.n1 = &(sim.vnode[node_group2[0]]);
    bnew.n2 = &(sim.vnode[ind_broke]);
    bnew.n3 = &(sim.vnode[node_group2.back()]);

    //Bend Angle
    bnew.a0 = 0;
    for (int j = node_group1.size() + 1; j < bend_indeces.size(); j++) bnew.a0 += (PI - abs(sim.vbend[bend_indeces[j]].a0));
    bnew.a0 = PI - bnew.a0;
    //printf("%f\n", 180*thetaOfBend(bnew)/PI);
    //bnew.a0 *= pow(-1, signbit(thetaOfBend(bnew)));
    bnew.a0 = abs(bnew.a0) * pow(-1, signbit(thetaOfBend(bnew)));

    //Stiffness of Bend
    double lengths[2];
    for(int i = 0; i<sim.n_bonds; i++){
      if (((sim.vbond[i].n1 == bnew.n1) && (sim.vbond[i].n2 == bnew.n2)) || ((sim.vbond[i].n2 == bnew.n1) && (sim.vbond[i].n1 == bnew.n2)))
        lengths[0] = sim.vbond[i].d0;
      else if (((sim.vbond[i].n1 == bnew.n3) && (sim.vbond[i].n2 == bnew.n2)) || ((sim.vbond[i].n2 == bnew.n3) && (sim.vbond[i].n1 == bnew.n2)))
        lengths[1] = sim.vbond[i].d0;
    }
    bnew.k = sim.kbe/min(lengths[0], lengths[1]);
    //printf("lengths = %f %f\n", lengths[0], lengths[1]);

    #ifdef VRBS
      printf("\t\t-2 Adding Bend (index %d), with angle %f and stiffness %f\n", sim.n_bends, bnew.a0*180/PI, bnew.k);
      printf("\t\t between nodes %d %d %d\n", bnew.n1->index, bnew.n2->index, bnew.n3->index);
    #endif
    sim.vbend.push_back(bnew);
    sim.n_bends++;
  }


  //Step 9: Delete the two bends that have been maximized.
  if (bend_indeces[0] < bend_indeces[node_group1.size()]){ //If remove lower element first, larger element index drops by one
    sim.vbend.erase(sim.vbend.begin() + bend_indeces[node_group1.size()]);
    sim.vbend.erase(sim.vbend.begin() + bend_indeces[0]);
    #ifdef VRBS
      printf("\t\t-Removing Bend %d and then Bend %d\n", bend_indeces[node_group1.size()], bend_indeces[0]);
    #endif
  }
  else{
    sim.vbend.erase(sim.vbend.begin() + bend_indeces[0]);
    sim.vbend.erase(sim.vbend.begin() + bend_indeces[node_group1.size()]);
    #ifdef VRBS
      printf("\t\t-Removing Bend %d and then Bend %d\n", bend_indeces[0], bend_indeces[node_group1.size()]);
    #endif
  }
  sim.n_bends -= 2;
}



/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ***************** MAIN ******************* */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
int main(int argc, char* argv[]){
  //Random Seed
  srand(time(NULL));

  sim sim;
  in inf;
  out outf;
  string strBuff;

  //READ ARGUMENTS
  //SET-UP INPUT AND OUTPUT FILES
  inf.inputDir = "./";
  for(int i=1; i<argc; i++) {
    // Check for a dash
    if(argv[i][0]=='-') {
      switch (argv[i][1]) {
        case 'i':
          i++;
          strBuff = argv[i];
          // Make sure that the directory ends with a "/"
          if( strBuff[strBuff.length()-1]!='/' ) strBuff += "/";
          printf("Setting up input from directory '%s'\n", strBuff.c_str());
          inf.inputDir = strBuff;
          break;
        case 'o':
          i++;
          strBuff = argv[i];
          //Make sure the directory ends with a '/'
          if( strBuff[strBuff.length()-1]!='/' ) strBuff += "/";
          printf("Setting up output in directory '%s'\n", strBuff.c_str());
          outf.outputDir = strBuff;
          break;
        case 'h':
          printf("\n\t\tSimulator for Lattice Damage (SimLaD).\n");
          printf("\tBy the people of MEGASLab, for the people of MEGASLab.\n\n");
          stateHelp();
          exit(EXIT_SUCCESS);
          break;
        default:
          printf("Error:\n\tInvalid arument: %s\n",argv[i]);
          stateHelp();
          exit(EXIT_FAILURE);
      }
    }
  }

  //INITIALIZE SIMULATOR/READ FROM FILES
  printf("Initializing System \n");
  if( access( inf.inputDir.c_str(), F_OK ) != -1 ) {
    // Read the simulation file
		readSim( &inf, &outf, &sim );
		// Read the nodes
		readNodes( &inf,&sim );
    //for(int i=0; i<sim.n_nodes; i++) printf( "\tNode %d: %f %f\n",i,sim.vnode[i].r[_x],sim.vnode[i].r[_y]);
    // Read the bonds
    readBonds( &inf,&sim );
    //printf("%p\n", (sim.vbond[0].n2));
    // Read the bends
    readBends( &inf,&sim );
  }
  printf("Setting up finished \n");


  //INITIALIZE OUTPUT INTO POSITION FILE AND PRINT
	openOutput(&outf.posF, outf.outputDir, "pos.txt");
  openOutput(&outf.bondF, outf.outputDir, "bondnew.txt");
  openOutput(&outf.systemE, outf.outputDir, "energy.txt");
  openOutput(&outf.bforce, outf.outputDir, "bforce.txt");
  //Don't need this line if printing before/after each avalance, because the first print is the one you need to get initial number of nodes
  fprintf(outf.posF, "%d InitNodes\n", sim.n_nodes);


  printf("Running... \n\n");
  struct timespec now, tmstart;
  clock_gettime(CLOCK_REALTIME, &tmstart);


  /*
  //Stretch to catastrophe
  int cnt = 0, ind_to_break;
  while ((checkCatastrophe(sim) == 0) || (cnt == 0)){
    cnt++;
    printf("\tFracture event no. %d.\n", cnt);
    ind_to_break = stretchToFracture(sim);
    changeTopology(sim, ind_to_break);
    FIRE2(sim);printf("\n");
    ind_to_break = stressCheck(sim, 1);
    while (ind_to_break != -1){
      printf("\tAvalance - Fracture event no. %d. Breaking node %d\n", cnt+1, ind_to_break);
      cnt++;
      changeTopology(sim, ind_to_break);
      FIRE2(sim);
      ind_to_break = stressCheck(sim, 1);
    }
  }
  printf("\n\n\t\t\tLattice broke all the way\n\n\n");
  */

  
  //Print out node positions and bond information
  outputPos(outf.posF, sim); //Output Positions
  outputBonds(outf.bondF, sim); //Output Bond stresses and assignments
  calcEnergy(sim); outputEnergy(outf.systemE, sim); //Calculate and output System axial/bending energy in file
  outputForceDisp(outf.bforce, sim); //Output Forces on Boundary

  //Stretch to n fracture events
  int n = 150, ind_to_break;
  for(int i = 0; i < n; i++){
    printf("\tFracture event no. %d.\n", i);
    ind_to_break = stretchToFracture(sim);

    //Print out node positions and bond information
    outputPos(outf.posF, sim); //Output Positions
    outputBonds(outf.bondF, sim); //Output Bond stresses and assignments
    calcEnergy(sim); outputEnergy(outf.systemE, sim); //Calculate and output System axial/bending energy in file
    outputForceDisp(outf.bforce, sim); //Output Forces on Boundary

    changeTopology(sim, ind_to_break);
    FIRE2(sim);printf("\n");
    ind_to_break = stressCheck(sim, 1);
    while ((ind_to_break != -1) && (i < n-1)){
      printf("\tAvalance - Fracture event no. %d. Breaking node %d\n", i+1, ind_to_break);
      i++;
      changeTopology(sim, ind_to_break);
      FIRE2(sim); printf("\n");
      ind_to_break = stressCheck(sim, 1);
    }

    //Print out node positions and bond information
    outputPos(outf.posF, sim); //Output Positions
    outputBonds(outf.bondF, sim); //Output Bond stresses and assignments
    calcEnergy(sim); outputEnergy(outf.systemE, sim); //Calculate and output System axial/bending energy in file
    outputForceDisp(outf.bforce, sim); //Output Forces on Boundary

    if ((checkCatastrophe(sim) == 1) && (i >= 20)){ //50 is to make sure the crack has advanced beyond initial diffuse phase, where it's likely lattice is unstressed after fracture.
      printf("\n\n\t\t\tLattice broke all the way\n\n\n");
      break;
    }
  }
  


  /*
  int ind_to_break;
  ind_to_break = stretchToFracture(sim);
  changeTopology(sim, ind_to_break);
  FIRE2(sim);
  stressCheck(sim, 0);
  outputPos(outf.posF, sim); //Output Positions
  outputBonds(outf.bondF, sim); //Output Bond stresses and assignments
  */


  /*
  //Stretch to failure, without changing topology, and output data at the end when you find failed node.
  //This is used for generating stress ratio and probability of failure estimates.
  int ind_to_break;
  ind_to_break = stretchToFracture(sim);
  printf("Node %d reached failure first.\n", ind_to_break);
  outputPos(outf.posF, sim); //Output Positions
  outputBonds(outf.bondF, sim); //Output Bond stresses and assignments
  */
  

  /*
  //FIRE LOOP only, no fracture
  int strain_increments = 10;
  //FIRE2(sim); //Initial equilibration should not be needed, if lattice is relaxed to begin with
  //stressCheck(sim, 0);

  //Print out node positions and bond information
  outputPos(outf.posF, sim); //Output Positions
  outputBonds(outf.bondF, sim); //Output Bond stresses and assignments
  calcEnergy(sim); outputEnergy(outf.systemE, sim); //Calculate and output System axial/bending energy in file
  outputForceDisp(outf.bforce, sim); //Output Forces on Boundary

  for (int iter = 0; iter < strain_increments; iter++){
    moveBoundary(sim, 1);
    FIRE2(sim);
    stressCheck(sim, 0);

    //Print out node positions and bond information
    outputPos(outf.posF, sim); //Output Positions
    outputBonds(outf.bondF, sim); //Output Bond stresses and assignments
    calcEnergy(sim); outputEnergy(outf.systemE, sim); //Calculate and output System axial/bending energy in file
    outputForceDisp(outf.bforce, sim); //Output Forces on Boundary

  }
  //stressCheck(sim, 0);
  */

  clock_gettime(CLOCK_REALTIME, &now);
  double seconds = (double)((now.tv_sec+now.tv_nsec*1e-9) - (double)(tmstart.tv_sec+tmstart.tv_nsec*1e-9));
  printf("Execution wall time: %fs\n", seconds);
  printf("Run successful - Terminating.\n");

  return 0;
}
