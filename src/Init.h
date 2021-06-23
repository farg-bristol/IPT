#ifndef INIT_H
#define INIT_H

#include <random>
#include "Var.h"
#include "Containment.h"

real random(int const& interval)
{
	return real(rand() % interval - interval / 2) * MERROR;
}

void Init_Disk_Points(SETT const& svar, FLUID const& fvar,Vec_Tree const& TREE, MESH const& cells, State& pn)
{
    // cout << "Adding points..." << endl;
	
	real jetR = 0.5*(svar.discDiam);

	real d_disc = svar.discSpace;
    real d_0 = fvar.d_0;
    real vol_0 = 4.0/3.0 * M_PI * pow(d_0/2,3);
    real mass_0 = fvar.d_rho*vol_0;

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
			pn.emplace_back(part(xi,mass_0,d_0));
		}
	}
	
	/* Create centre point */
	vec<real,3> xi(0.0);
	// perturb = vec<real,3>(random(interval), random(interval), random(interval));
	// xi += perturb;
	pn.emplace_back(part(xi,mass_0,d_0));

    for(size_t ii = 0; ii < pn.size(); ++ii)
	{
        pn[ii].partID = ii;
		pn[ii].xi = (svar.rotate * pn[ii].xi) + svar.discPos;

		FindCell(svar,TREE, cells, pn[ii]);
		
		pn[ii].v = pn[ii].cellV;
	}	
}

void Init_Line_Points(SETT const& svar, FLUID const& fvar,Vec_Tree const& TREE, MESH const& cells, State& pn)
{
    // cout << "Adding points..." << endl;
	
    real d_0 = fvar.d_0;
    real vol_0 = 4.0/3.0 * M_PI * pow(d_0/2,3);
    real mass_0 = fvar.d_rho*vol_0;

	// int const interval = 1;
	// vec<real,3> perturb;
	/* Need to find deltas between each vertex, so create the grid upon */
	real n_i = real(svar.n_i - 1);
	vec<real,3> di1;

	di1 = (svar.grid_verts[1] - svar.grid_verts[0])/n_i;

	vec<real,3> const& xi_0 = svar.grid_verts[0];

	for(int ii = 0; ii < svar.n_i; ++ii)
	{
		vec<real,3> const dx_i = real(ii)* di1;

		vec<real,3> xi = xi_0 + dx_i; 
		pn.emplace_back(part(xi,mass_0,d_0));
	}

    for(size_t ii = 0; ii < pn.size(); ++ii)
	{
        pn[ii].partID = ii;

		FindCell(svar,TREE, cells, pn[ii]);
		
		pn[ii].v = pn[ii].cellV;
	}	
}

void Init_Grid_Points(SETT const& svar, FLUID const& fvar,Vec_Tree const& TREE, MESH const& cells, State& pn)
{
    // cout << "Adding points..." << endl;
	
    real d_0 = fvar.d_0;
    real vol_0 = 4.0/3.0 * M_PI * pow(d_0/2,3);
    real mass_0 = fvar.d_rho*vol_0;

	// int const interval = 1;
	// vec<real,3> perturb;
	/* Need to find deltas between each vertex, so create the grid upon */
	real n_i = real(svar.n_i - 1);
	real n_j = real(svar.n_j - 1);
	vec<real,3> di1, di2, dj1, dj2;

	di1 = (svar.grid_verts[1] - svar.grid_verts[0])/n_i;
	di2 = (svar.grid_verts[2] - svar.grid_verts[3])/n_i;

	dj1 = (svar.grid_verts[3] - svar.grid_verts[0])/n_j;
	dj2 = (svar.grid_verts[2] - svar.grid_verts[1])/n_j;

	vec<real,3> const& xi_0 = svar.grid_verts[0];

	for(int ii = 0; ii < svar.n_i; ++ii)
	{	
		for(int jj = 0; jj < svar.n_j; ++jj)
		{
			vec<real,3> const dx_i = real(ii)*(
			(n_j - real(jj))/n_j * di1 + real(jj)/n_j * di2);

			vec<real,3> const dx_j = real(jj)*(
			(n_i - real(ii))/n_i * di1 + real(ii)/n_i * di2);

			vec<real,3> xi = xi_0 + dx_i + dx_j; 
			pn.emplace_back(part(xi,mass_0,d_0));
		}
	}

    for(size_t ii = 0; ii < pn.size(); ++ii)
	{
        pn[ii].partID = ii;
		FindCell(svar,TREE, cells, pn[ii]);
		pn[ii].v = pn[ii].cellV;
	}	
}

void Init_Points(SETT const& svar, FLUID const& fvar,Vec_Tree const& TREE, MESH const& cells, State& pn)
{
	if(svar.startType == 0)
	{
		Init_Grid_Points(svar,fvar,TREE,cells,pn);
	}
	else if(svar.startType == 1)
	{
		Init_Line_Points(svar,fvar,TREE,cells,pn);
	}
	else if (svar.startType == 2)
	{
		Init_Disk_Points(svar,fvar,TREE,cells,pn);
	}
}

void Init_Surface(SETT const& svar, MESH const& cells, vector<vector<SURF>>& surf_faces, vector<SURF>& surf_marks)
{
	surf_faces = vector<vector<SURF>>(svar.markers.size(), vector<SURF>());
	surf_marks = vector<SURF>(svar.markers.size());

	vector<vector<size_t>> faceIDs(svar.markers.size());
	vector<vector<int>> markers(svar.markers.size());
	/* for each surface, find how many faces are in it. */
	for(std::pair<size_t,int> const& marker:cells.smarkers)
	{
		auto index = find(svar.markers.begin(),svar.markers.end(),marker.second);
		if(index != svar.markers.end())
		{
			size_t mark = index - svar.markers.begin();
			faceIDs[mark].emplace_back(marker.first);
			markers[mark].emplace_back(marker.second);

			// surf_faces[mark].back().faceID  = marker.first;
			// if()
			// cout << mark << "  " << markers.first << endl;
			// cout << surf_faces[mark].back().faceID << "  " << markers.first << endl;
			// surf_faces[mark].back().marker  = marker.second;
		}
		else
		{
			cout << "Couldn't find the marker in the index" << endl;
		}
		
		
	}

	for(size_t ii = 0; ii < svar.markers.size(); ii++)
	{
		surf_marks[ii].name = svar.bnames[ii];
		surf_marks[ii].marker = svar.markers[ii];

		surf_faces[ii] = vector<SURF>(faceIDs[ii].size());


		for(size_t jj = 0; jj < faceIDs[ii].size(); jj++)
		{
			surf_faces[ii][jj].faceID = faceIDs[ii][jj];
			surf_faces[ii][jj].marker = markers[ii][jj];
		}
	
	}

	// /* allocate the impacte vector */
	// for(size_t ii = 0; ii < nSurf; ++ii)
	// {
	// 	surfs[ii].alloc(nFaces[ii]);
	// }
}

#endif