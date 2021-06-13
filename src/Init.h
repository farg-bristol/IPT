#ifndef INIT_H
#define INIT_H

#include <random>
#include "Var.h"
#include "Containment.h"

real random(int const& interval)
{
	return real(rand() % interval - interval / 2) * MERROR;
}

void Init_Points(SETT const& svar, FLUID const& fvar,Vec_Tree const& TREE, MESH const& cells, State& pn)
{
    // cout << "Adding points..." << endl;
	
	real jetR = 0.5*(svar.discDiam);

	real d_disc = svar.discSpace;
    real d_0 = fvar.d_0;
    real vol_0 = 4.0/3.0 * M_PI * pow(d_0/2,3);
    real mass_0 = fvar.d_rho*vol_0;

    vec<real,3> v(fvar.v_0);

	// int const interval = 1;
	// vec<real,3> perturb;

	for(real rad = jetR; rad > 0.99*d_disc; rad -= d_disc)
	{
		// % Find spacing to have a well defined surface.
		real dtheta = atan(d_disc / rad);
		int ncirc = floor(abs(2.0 * M_PI / dtheta));
		dtheta = 2.0 * M_PI / real(ncirc);
		
		for(real theta = 0.0; theta < 2*M_PI; theta += dtheta)
		{	/* Create a ring of points */
			real x = rad * sin(theta);
			real z = rad * cos(theta);

			vec<real,3> xi(x,0,z);
			// perturb = vec<real,3>(random(interval), random(interval), random(interval));
			// xi += perturb;
			pn.emplace_back(part(xi,v,mass_0,d_0));
		}
	}
	
	/* Create centre point */
	vec<real,3> xi(0.0);
	// perturb = vec<real,3>(random(interval), random(interval), random(interval));
	// xi += perturb;
	pn.emplace_back(part(xi,v,mass_0,d_0));

    for(size_t ii = 0; ii < pn.size(); ++ii)
	{
        pn[ii].partID = ii;
		pn[ii].xi = (svar.rotate * pn[ii].xi) + svar.discPos;

		FindCell(svar,TREE, cells, pn[ii]);
		
		pn[ii].v = pn[ii].cellV;
		// pn[ii].cellIDm1 = pn[ii].cellID;
		// pn[ii].cellVm1 = pn[ii].cellV;
		// pn[ii].cellRhom1 = pn[ii].cellRho;
		// pn[ii].faceV = pn[ii].cellV; 
		// pn[ii].faceVm1 = pn[ii].cellV;
		// pn[ii].faceRho = pn[ii].cellRho; 
		// pn[ii].faceRhom1 = pn[ii].cellRho; 

	}	
}

#endif