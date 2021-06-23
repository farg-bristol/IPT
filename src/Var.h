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
		partout = 0;
		streakout = 0;
		cellsout = 0;

		scale = 1.0;
		discPos = vec<real,3>(-500000.0);
		discDiam = 0.0;
		discSpace = 0.0;
		rotate = matrix<real,3>(0.0);
		transp = matrix<real,3>(0.0);	

		startType = -1;
		grid_verts = vector<vec<real,3>>(4,-500000);
		n_i = -1; n_j = -1;

		gamma = 0.5;
		beta = 0.25;

		integration_type = -1;
		t = 0.0;
		dt = 0.0;
		framet = 0.0;
		frame = 0;
		max_frame = 999999;
		max_iter = 30;
		eqOrder = 1;

		nSuccess = 0;
		nFailed = 0;
	}

	/* TAU input options*/
	string taumesh;
	string taubmap;
	string tausol;
	real scale;
	vector<int> markers;
	vector<string> bnames;
	vector<int> bwrite; /* If surface wants data writing */

	/* OpenFOAM input options */
	string foamdir;
	string foamsol;
	int isBinary;
	int labelSize;
	int scalarSize;

	int TAUorFOAM;

	/* Output folders */
	string outdir;
	string streakdir;
	string celldir;
	vector<ofstream> partfiles;
	string surfacefile;
	int partout, streakout, cellsout;

	/* Starting area conditions */
	string startName;
	int startType;
	vector<vec<real,3>> grid_verts;
	int n_i, n_j;
	vec<real,3> discPos;
	vec<real,3> discAngles;
	real discDiam;
	real discSpace;
	matrix<real,3> rotate;
	matrix<real,3> transp;

	real gamma, beta; /* Newmark beta parameters */

	/* Simulation settings */
	string integration_type; /* What kind of integration to use */
	int explicit_or_implicit;  /* Implicit or explicit integration */
	real t; /* Simulation time */
	real dt; /* Simulation timestep */
	real framet; /* Frame time */
	real nsteps;
	int frame; /* Current frame */
	int max_frame; /* Max output frames */
	int max_iter; /* Max Newmark Beta iterations */
	int eqOrder;

	size_t nSuccess, nFailed;
};

struct FLUID
{	/* Initialise default values */
	FLUID()
	{
		g_rho = 1.29251;
		d_rho = 1000;
		mu_g = 1.716e-5;
		d_0 = 1;
	}

	real g_rho; /* Gas rho */
	real d_rho; /* dispersed rho */
	real mu_g; /* Gas viscosity */
	real d_0; /* Starting particle diameter */
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

	part(vec<real,3> const xi_, real const mass_, real const d_)
	{
		going = 1;
		partID = 0;
		nNotFound = 0;
		t = 0; dt = 0;

		/* Set ininial IDs to a nonsense value, so that they don't interfere */
		faceID = -1; 
		faceV = vec<real,3>(0.0);
		faceRho = 0.0;
	
		cellID = -3;
		cellV = vec<real,3>(0.0); 
		cellRho = 0.0; 
			
		acc = vec<real,3>(0.0);
		v = vec<real,3>(0.0);
		xi = xi_;
		
		mass = mass_;
		d = d_;
		A = M_PI * d_*d_/4.0; 
	}

	part(vec<real,3> const xi_, part const& pi_)
	{
		going = 1;
		partID = 0;
		nNotFound = 0;
		t = 0; dt = 0;

		/* Set ininial IDs to a nonsense value, so that they don't interfere */
		faceID = -1; 
		faceV = vec<real,3>(0.0);
		faceRho = 0.0;
	
		cellID = pi_.cellID;
		cellV = pi_.cellV; 
		cellRho = pi_.cellRho; 
			
		acc = pi_.acc;
		v = pi_.v;
		xi = xi_;
		
		mass = pi_.mass;
		d = pi_.d;
		A = pi_.A; 
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

	/* Containing cell properties */
	lint cellID;
	vec<real,3> cellV;
	real cellRho;

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
		cRho = std::vector<real>(nC, 0.0);
	}

	size_t const size() const {return elems.size();}

	size_t nPnts, nElem, nFace, nSurf;
	real scale;

	/*Point based data*/
	vector<vec<real,3>> verts;

	/*Face based data*/
	vector<vector<size_t>> faces; /* Index of vertices defining the face */
	vector<std::pair<lint,lint>> leftright;
	vector<std::pair<size_t,int>> smarkers;

	/*Cell based data*/
	vector<vector<size_t>> elems; /* Index of vertices the cell uses (uniquely) */
	vector<vector<size_t>> cFaces; /* Index of faces defining the cell */
	vector<vec<real,3>> cCentre;
	
	/*Solution vectors*/
	vector<vec<real,3>> cVel;
	vector<real> cRho;
	
};

/* Structure of data for surface face impacts */
struct IMPACT
{
	IMPACT(){ length = 0;}

	void alloc(size_t n)
	{
		length = n;
		count = vector<uint>(n);
		pIDs = vector<vector<size_t>>(n,vector<size_t>());
		pos = vector<vector<vec<real,3>>>(n,vector<vec<real,3>>());
		mass = vector<real>(n);
		MVDs = vector<vector<real>>(n,vector<real>());
		colEff = vector<real>(n);
		area = vector<real>(n);
	}

	size_t const size() const { return length;}

	vector<size_t> faceIDs; /* Face IDs of the mesh */
	vector<unsigned> count; /* Number of collisions */
	vector<vector<size_t>> pIDs; /* Particle IDs (for collection efficiency)*/
	vector<vector<vec<real,3>>> pos; /* Particle positions (for collection efficiency) */
	vector<real> mass;   /* Total mass of particles hitting the face */
	vector<real> colEff; /* Collision efficiency (same as collection eff) */
	vector<real> area; /* area between points (To see how the calculation goes) */
	vector<vector<real>> MVDs; /* Record of all the MVDs that hit, for statistics */
	size_t length;
};

/* Another structure simply to identify characteristics for a marker */
struct SURF
{
	SURF()
	{ 
		count = 0;
		faceID = 0;
		marker = 0;
	}

	SURF(size_t const& face, int const& mark): faceID(face), marker(mark) {count = 0;}

	SURF(std::pair<size_t,int> const& mark)
	{
		faceID = mark.first;
		marker = mark.second;
		count = 0;
	}

	SURF(SURF const& bi_)
	{
		name = bi_.name;
		marker = bi_.marker;
		count = 0;
	}

	string name;
	unsigned count;
	size_t faceID;
	int marker;
	real colEff; /* Collision efficiency (same as collection eff) */
	real area;
	vector<size_t> pIDs;
	vector<vec<real,3>> end_pos;
	vector<vec<real,3>> start_pos;
	
};

typedef vector<part> State;

typedef KDTreeVectorOfVectorsAdaptor<std::vector<vec<real,3>>,real,3,nanoflann::metric_L2_Simple,size_t> Vec_Tree;

#endif
