#ifndef VAR_H
#define VAR_H

#include <vector>
#include <limits>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <sstream>
#include <omp.h>

#include "Vectors.h"

#include "Third_Party/NanoFLANN/nanoflann.hpp"
#include "Third_Party/NanoFLANN/utils.h"
#include "Third_Party/NanoFLANN/KDTreeVectorOfVectorsAdaptor.h"

using std::vector;
using std::cout;
using std::cerr;
using std::ifstream;
using std::ofstream;
using std::endl;
using std::string; 
using std::setw;
using std::fabs;

// Define pi
#ifndef M_PI
#define M_PI (4.0*atan(1.0))
#endif

/* Define data type. */
#ifndef FOD
#define FOD 1 /*0 = float, 1 = double*/
#endif

#if FOD == 1
typedef double real;
#else
typedef float real;
#endif

typedef unsigned int uint;

/*Get machine bit precision for Simulation of Simplicity*/
#ifndef MEPSILON
#define MEPSILON std::numeric_limits<real>::epsilon() /*For  float, power is -24*/
#endif

#ifndef MERROR
#define MERROR (7*MEPSILON + 56*MEPSILON*MEPSILON)
#endif

#if UINTPTR_MAX == 0xffffffffffffffff
/* it's 64bits pointers */
typedef int64_t lint;
#elif UINTPTR_MAX == 0xffffffff
/* it's 32bits pointers */
typedef int32_t lint;
#endif



// #else
// 
// #endif

const std::string WHITESPACE = " \n\r\t\f\v";

#pragma omp declare reduction(+: std::vector<vec<real,3>> : \
        std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), \
        	[](vec<real,3> lhs, vec<real,3> rhs){return lhs + rhs;})) \
                    initializer(omp_priv = omp_orig)

#pragma omp declare reduction(+: std::vector<real> : \
        std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), \
        	[](real lhs, real rhs){return lhs + rhs;})) \
                    initializer(omp_priv = omp_orig)  

/* Simulation settings structure */
struct SETT
{
	SETT()
	{
		isBinary = 0;
		labelSize = -1;
		scalarSize = -1;
		TAUorFOAM = 0;

		scale = 1.0;
		discPos = vec<real,3>(-500000.0);
		discDiam = 0.0;
		discSpace = 0.0;
		rotate = matrix<real,3>(0.0);
		transp = matrix<real,3>(0.0);	

		gamma = 0.5;
		beta = 0.25;

		t = 0.0;
		dt = 0.0;
		framet = 0.0;
		frame = 0;
		max_frame = -1;
		max_iter = 30;
	}

	/* TAU input options*/
	string taumesh;
	string taubmap;
	string tausol;
	real scale;
	vector<int> markers;

	/* OpenFOAM input options */
	string foamdir;
	string foamsol;
	int isBinary;
	int labelSize;
	int scalarSize;

	int TAUorFOAM;

	/* Output filename */
	string outfile;
	string streakfile;

	/* Starting area conditions */
	vec<real,3> discPos;
	vec<real,3> discAngles;
	real discDiam;
	real discSpace;
	matrix<real,3> rotate;
	matrix<real,3> transp;

	real gamma, beta;

	/* Simulation settings */
	real t; /* Simulation time */
	real dt; /* Simulation timestep */
	real framet; /* Frame time */
	real nsteps;
	int frame; /* Current frame */
	int max_frame; /* Max output frames */
	int max_iter; /* Max Newmark Beta iterations */
};

struct FLUID
{	/* Initialise default values */
	FLUID()
	{
		g_rho = 1.29251;
		d_rho = 1000;
		mu_g = 1.716e-5;
		d_0 = -1;
		v_0 = vec<real,3>(-1);
	}

	real g_rho; /* Gas rho */
	real d_rho; /* dispersed rho */
	real mu_g; /* Gas viscosity */
	real d_0; /* Starting particle diameter */
	vec<real,3> v_0; /* Starting velocity */
};

/* Particle structure  */
struct part
{
	part()
	{
		partID = 0;
		going = 1;
		nNotFound = 0;
		mass = 0;
		d = 0;
		A = 0;
		xi = vec<real,3>(0.0);
		v = vec<real,3>(0.0);
		acc = vec<real,3>(0.0);

		cellID = 0;
		cellRho = 0.0;
		cellV = vec<real,3>(0.0);
	}

	part(vec<real,3> const xi_, vec<real,3> const v_, real const mass_, real const d_)
	{
		going = 1;
		partID = 0;
		nNotFound = 0;
		t = 0; dt = 0;

		/* Set ininial IDs to a nonsense value, so that they don't interfere */
		faceID = -1; 
		faceV = vec<real,3>(0.0);
		faceRho = 0.0;
		
		// faceIDm1 = std::numeric_limits<size_t>::max();
		// faceVm1 = vec<real,3>(0.0);
		// faceRhom1 = 0.0;

		cellID = -3;
		cellV = vec<real,3>(0.0); 
		cellRho = 0.0; 
		
		// cellIDm1 = std::numeric_limits<size_t>::max();
		// cellVm1 = vec<real,3>(0.0); 
		// cellRhom1 = 0.0;	
		
		acc = vec<real,3>(0.0);
		v = v_;
		xi = xi_;
		
		mass = mass_;
		d = d_;
		A = M_PI * d_*d_/4.0; 
	}

	size_t partID;
	size_t going; /* Is particle still being integrated */
	size_t nNotFound;

	/* Timestep properties */
	real t, dt;

	/* Face properties */
	lint faceID;
	vec<real,3> faceV;
	real faceRho;

	// size_t faceIDm1;
	// vec<real,3> faceVm1;
	// real faceRhom1;

	/* Containing cell properties */
	lint cellID;
	vec<real,3> cellV;
	real cellRho;

	// size_t cellIDm1;
	// vec<real,3> cellVm1;
	// real cellRhom1;

	vec<real,3> acc, v, xi; /* State variables */

	real mass; /* Particle mass */
	real d, A; /* Particle diameter and area */
};

struct MESH
{
	/*Standard contructor*/
	MESH(){}

	void SetCells()
	{
		uint nC = elems.size();
		cVel = std::vector<vec<real,3>>(nC, 0.0);		
	}

	size_t size() {return elems.size();}

	/*Zone info*/
	std::string zone;
	size_t nPnts, nElem, nFace, nSurf;
	real scale;

	/*Point based data*/
	vector<vec<real,3>> verts;

	/*Face based data*/
	vector<vector<size_t>> faces;
	vector<std::pair<lint,lint>> leftright;
	vector<std::pair<size_t,int>> smarkers;

	/*Cell based data*/
	vector<vector<size_t>> elems;
	vector<vector<size_t>> cFaces;
	vector<vec<real,3>> cCentre;
	
	/*Solution vectors*/
	vector<vec<real,3>> cVel;
	vector<real> cRho;
	
};

typedef vector<part> State;

typedef KDTreeVectorOfVectorsAdaptor<std::vector<vec<real,3>>,real,3,nanoflann::metric_L2_Simple,size_t> Vec_Tree;

#endif
