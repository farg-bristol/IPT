#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "Var.h"
#include "Containment.h"

void Terminate_Particle(SETT const& svar, MESH const& cells, vector<State> const& time_record,
    part& pnp1, vector<vector<SURF>>& surface_faces, vector<SURF>& surface_data )
{
    pnp1.going = 0;

    // cout << "Particle " << pnp1.partID << " stopped iterating"; 
    if(pnp1.cellID < 0)
    {
        /* Find which surface it has impacted  */
        size_t const& face = pnp1.faceID;

        auto const index = std::find_if(cells.smarkers.begin(), cells.smarkers.end(), 
        [face](std::pair<size_t,int> const& p1){return p1.first == face;});

        if(index != cells.smarkers.end())
        {
            /* Find which surface to add to, and which face of that*/
            auto const surface = std::find(svar.markers.begin(), svar.markers.end(), index->second);

            if(surface != svar.markers.end())
            {
                size_t const surf = surface - svar.markers.begin();

                // cout << " after hitting surface: " << svar.bnames[surf] << endl;
                surface_data[surf].count++;
                surface_data[surf].pIDs.emplace_back(pnp1.partID);
                surface_data[surf].end_pos.emplace_back(pnp1.xi);
                surface_data[surf].start_pos.emplace_back(time_record[0][pnp1.partID].xi);

                auto const index2 =
                    std::find_if(surface_faces[surf].begin(), surface_faces[surf].end(), 
                    [face](SURF const& p){return p.faceID == face;});

                if(index2 != surface_faces[surf].end())
                {
                    size_t fIndex = index2 - surface_faces[surf].begin();
                    surface_faces[surf][fIndex].count++;
                    surface_faces[surf][fIndex].pIDs.emplace_back(pnp1.partID);
                    surface_faces[surf][fIndex].end_pos.emplace_back(pnp1.xi);
                    surface_faces[surf][fIndex].start_pos.emplace_back(time_record[0][pnp1.partID].xi);
                }
                else
                {
                    cout << "Couldn't find which face of the surface the particle impacted" << endl;
                }
            }
            else
            {
                cout << "Couldn't find which surface index for the marker." << endl;
            }
            
        }
        else
        {
            cout << "Couldn't find which surface the particle impacted" << endl;
        }

        
    }
    else
    {
        // cout << " prematurely" << endl;
    }
    
    

    if(svar.streakout == 1)
    {
        Write_ASCII_Streaks(svar,time_record,pnp1, pnp1.partID);
    }

    if(svar.cellsout == 1)
    {
        Write_ASCII_Cells(svar,cells,time_record, pnp1.partID);
    }
}


real GetCd(real const& Re)
{
	return (1.0+0.197*pow(Re,0.63)+2.6e-04*pow(Re,1.38))*(24.0/(Re+0.0000001));

	// if( Re <= 1000.0)
	// 	return (24.0/Re)*(1+(1.0/6.0)*pow(Re,2.0/3.0));
	// else 
	// 	return 0.424;
}

vec<real,3> Explicit_AeroForce(vec<real,3> const& Vdiff, FLUID const& fvar, part const& pi)
{
	real const Re = 2.0*pi.cellRho*Vdiff.norm()*pi.d/fvar.mu_g;
	real const Cdi = GetCd(Re);

	// std::cout << "Reynolds: " << Re  << " Cd: " << Cdi << std::endl;
	return (0.5*pi.cellRho*Vdiff.norm()*Vdiff*Cdi*pi.A);
}

void Explicit_Newmark_Beta(FLUID const& fvar, real const& a, real const& b, real const& c, real const& d, 
                  real const& dt, real const& dt2, part const& pn, part& pnp1)
{
    vec<real,3> const Vdiff = pnp1.cellV - pnp1.v ;
    vec<real,3> const g(0,0,-9.81);
    vec<real,3> res = Explicit_AeroForce(Vdiff,fvar,pnp1);
    res /= pnp1.mass;
    res += g;
    // cout << res[0] << "  " << res[1] << "  " << res[2] << endl;

    pnp1.xi = pn.xi+dt*pn.v+dt2*(c*pn.acc+d*res);
    pnp1.v =  pn.v +dt*(a*pn.acc+b*res);
    pnp1.acc = res;
}


void Explicit_Integrate(SETT& svar, FLUID const& fvar, Vec_Tree const& TREE, MESH const& cells, State& pn, State& pnp1,
                        vector<State>& time_record, vector<vector<SURF>>& surface_faces, vector<SURF>& marker_data)
{
    real const& dt = svar.dt;
    real const& dt2 = dt*dt;
    real& t = svar.t;

    const real a = 1 - svar.gamma;
	const real b = svar.gamma;
	const real c = 0.5*(1-2*svar.beta);
	const real d = svar.beta;

    uint nGoing = pnp1.size();

    while(nGoing != 0)
    {

        for (size_t step = 0; step < svar.nsteps; ++step)
        {
            nGoing = 0; /* Number of particle still integrating */
            /* Integrate each particle */
            for(size_t ii = 0; ii < pnp1.size(); ++ii)
            {
                if(pnp1[ii].going == 1)
                {
                    /* Check containing cell */
                    FindCell(svar,TREE, cells, pnp1[ii]);

                    /* Perform one step of NB, then freeze containment for stability. */
                    Explicit_Newmark_Beta(fvar, a,b,c,d,dt,dt2,pn[ii],pnp1[ii]);
                    FindCell(svar,TREE, cells, pnp1[ii]);

                    int iter = 0;
                    int min_iter = 3;
                    real error = 1.0;
                    while(error > -7.0 && iter < min_iter)
                    {
                        vec<real,3> temp = pnp1[ii].xi;
                        /* Get the aerodynamic force */
                        Explicit_Newmark_Beta(fvar, a,b,c,d,dt,dt2, pn[ii],pnp1[ii]);

                        error = log10((pnp1[ii].xi-temp).norm());

                        if(iter < svar.max_iter)
                            iter++;
                        else
                            break;
                    }

                    if(pnp1[ii].going == 0)
                    {
                        /* Use face ID to test if it's crossed a boundary */
                        if(pnp1[ii].faceID == -2 )
                        {
                            Terminate_Particle(svar,cells,time_record,pnp1[ii],surface_faces,marker_data);
                            svar.nSuccess++;
                        }
                        else
                        {
                            svar.nFailed++;
                        }

                    }
                    
                    /* Move forward in time */
                    pn[ii] = pnp1[ii];
                    nGoing += pnp1[ii].going;
                }
                
            }
            if(nGoing == 0)
            {
                cout << "All particles finished simulating" << endl;
                break;
            }

            t += dt;
        }

        for(size_t ii = 0; ii < pnp1.size(); ++ii)
        {
            if(pnp1[ii].going == 1)
            {
                if(svar.partout == 1)
                    Write_ASCII_Point(svar.partfiles[ii], svar.scale, pnp1[ii], ii);
            }
        }

        if(svar.frame < svar.max_frame)
            svar.frame++;
        else
        {
            /* Wrap up the simulation */
            cout << "Reached the maximum number of frames. Ending..." << endl;

            for(size_t ii = 0; ii < pnp1.size(); ++ii)
            {
                if(pnp1[ii].going == 1)
                {
                    Terminate_Particle(svar,cells,time_record,pnp1[ii],surface_faces,marker_data);
                    svar.nSuccess++;
                }
            }

            nGoing = 0;
        }
        
        time_record.emplace_back(pnp1);
    }
}


real Implicit_AeroForce(vec<real,3> const& Vdiff, FLUID const& fvar, part const& pi)
{
    real const Re = 2.0*pi.faceRho*Vdiff.norm()*pi.d/fvar.mu_g;
	real const Cd = GetCd(Re);


    return Vdiff.norm() * (3.0*Cd*pi.faceRho)/(4.0 * pi.d * fvar.d_rho);
}

void Implicit_BFD1(FLUID const& fvar, part const& pn, part& pnp1)
{
    vec<real,3> const Vdiff = pnp1.faceV - pnp1.v ;
    // vec<real,3> const g(0,0,-9.81);
    vec<real,3> const g(0.0,0.0,0.0);

    real res = Implicit_AeroForce(Vdiff, fvar, pnp1);
    
    pnp1.acc = res;
    pnp1.v = (pn.v + pnp1.dt*res*pnp1.faceV + pnp1.dt*g)/(1.0 + pnp1.dt*res);
    pnp1.xi = pn.xi+0.5*pnp1.dt*(pnp1.v + pn.v); 
}

void Implicit_BFD2(FLUID const& fvar, real const& dt, real const& dtm1, part const& pnm1, part const& pn, part& pnp1)
{
    vec<real,3> const Vdiff = pnp1.faceV - pnp1.v ;
    // vec<real,3> const g(0,0,-9.81);
    vec<real,3> const g(0.0,0.0,0.0);

    real res = Implicit_AeroForce(Vdiff, fvar, pnp1);

    pnp1.acc = res;
    /* Second order velocity calculation */
    pnp1.v = ((dt+dtm1)*(dt*dtm1*(res*pnp1.faceV+g) + (dt+dtm1)* pn.v) - dt*dt*pnm1.v)
                /(dtm1*((2*dt+dtm1)+dt*(dt+dtm1)*res));
                
    // pnp1.xi = pn.xi-pnp1.dt*(pnp1.v + pn.v)+0.5*pn.dt*(pn.v + pnm1.v); 
    pnp1.xi = pn.xi+0.5*pnp1.dt*(pnp1.v + pn.v); /* Just a first order space integration. */
}

void Implicit_Integrate(SETT& svar, FLUID const& fvar, MESH const& cells, 
                        State& pnm1, State& pn, State& pnp1, 
                        vector<State>& time_record, 
                        vector<vector<SURF>>& surface_faces, vector<SURF>& marker_data)
{
    uint nGoing = pnp1.size();

    /* Do the first step */
    for(size_t ii = 0; ii < pnp1.size(); ++ii)
    {
        if(pnp1[ii].going == 1)
        {
            /* Find intersecting face for np1 */
            // pnp1[ii].faceV = pnp1[ii].cellV; /* Start by guessing the face as the cell*/
            // pnp1[ii].faceRho = pnp1[ii].cellRho;
            FindFace(svar, cells, pn[ii], pnp1[ii]);
            /* Now c+1 is known, use the central difference face properties*/
            pnp1[ii].faceV = 0.5*(pnp1[ii].cellV + pn[ii].cellV);
            pnp1[ii].faceRho = 0.5*(pnp1[ii].cellRho + pn[ii].cellRho);

            int iter = 0;
            int min_iter = 3;
            real error = 1.0;
            while(error > -7.0 && iter < min_iter)
            {
                // cout << "Iteration: " << iter << endl; 
                vec<real,3> temp = pnp1[ii].xi;
                // Implicit_BFD1(fvar, pn[ii], pnp1[ii]);
                real const dt = pnp1[ii].dt;
                real const dtm1 = pnp1[ii].dt;
                
                if(svar.eqOrder == 2)
                    Implicit_BFD2(fvar, dt, dtm1, pnm1[ii], pn[ii], pnp1[ii]);
                else
                    Implicit_BFD1(fvar, pn[ii], pnp1[ii]);

                FindFace(svar, cells, pn[ii], pnp1[ii]);
                /* Now c+1 is known, use the central difference face properties*/
                pnp1[ii].faceV = 0.5*(pnp1[ii].cellV + pn[ii].cellV);
                pnp1[ii].faceRho = 0.5*(pnp1[ii].cellRho + pn[ii].cellRho);

                error = log10((pnp1[ii].xi-temp).norm());

                if(iter < svar.max_iter)
                    iter++;
                else
                    break;
            }

            if(svar.partout == 1)
                Write_ASCII_Point(svar.partfiles[ii], svar.scale, pnp1[ii], ii);

            if(pnp1[ii].going == 0)
            {
                if(pnp1[ii].cellID < 0 )
                {
                    Terminate_Particle(svar,cells,time_record,pnp1[ii],surface_faces,marker_data);
                    svar.nSuccess++;
                }
                else
                {
                    svar.nFailed++;
                }

            } 

            
            /* Move forward in time */
            // cout << "Time: " << pnp1[ii].t << " dt: " << pnp1[ii].dt  << " x: " <<
            // pnp1[ii].xi[0] << " y: " << pnp1[ii].xi[1] << " z: " << pnp1[ii].xi[2] << 
            // " old cell: " << pn[ii].cellID << " new cell: " << pnp1[ii].cellID << endl;

            /* March forward in time */
            if(svar.eqOrder == 2)
                pnm1[ii] = pn[ii];

            pn[ii] = pnp1[ii];

        }

        nGoing += pnp1[ii].going;
        pnp1[ii].t += pnp1[ii].dt;
    }
    time_record.emplace_back(pnp1);

    while(nGoing != 0)
    {
        nGoing = 0; /* Number of particle still integrating */
        /* Integrate each particle */
        for(size_t ii = 0; ii < pnp1.size(); ++ii)
        {
            if(pnp1[ii].going == 1)
            {
                /* Find intersecting face for np1 */
                // pnp1[ii].faceV = pnp1[ii].cellV; /* Start by guessing the face as the cell*/
                // pnp1[ii].faceRho = pnp1[ii].cellRho;
                FindFace(svar, cells, pn[ii], pnp1[ii]);
                /* Now c+1 is known, use the central difference face properties*/
                pnp1[ii].faceV = 0.5*(pnp1[ii].cellV + pn[ii].cellV);
                pnp1[ii].faceRho = 0.5*(pnp1[ii].cellRho + pn[ii].cellRho);

                int iter = 0;
                int min_iter = 3;
                real error = 1.0;
                while(error > -7.0 && iter < min_iter)
                {
                    // cout << "Iteration: " << iter << endl; 
                    vec<real,3> temp = pnp1[ii].xi;
                    
                    real const dt = pnp1[ii].dt;
                    real dtm1 = pnp1[ii].dt;
                    if(iter > 0)
                        dtm1 = pn[ii].dt;

                    if(svar.eqOrder == 2)
                        Implicit_BFD2(fvar, dt, dtm1, pnm1[ii], pn[ii], pnp1[ii]);
                    else
                        Implicit_BFD1(fvar, pn[ii], pnp1[ii]);


                    FindFace(svar, cells, pn[ii], pnp1[ii]);
                    /* Now c+1 is known, use the central difference face properties*/
                    pnp1[ii].faceV = 0.5*(pnp1[ii].cellV + pn[ii].cellV);
                    pnp1[ii].faceRho = 0.5*(pnp1[ii].cellRho + pn[ii].cellRho);

                    error = log10((pnp1[ii].xi-temp).norm());

                    if(iter < svar.max_iter)
                        iter++;
                    else
                        break;
                }

                if(svar.partout == 1)
                    Write_ASCII_Point(svar.partfiles[ii], svar.scale, pnp1[ii], ii);

                if(pnp1[ii].cellID < 0 )
                {
                    Terminate_Particle(svar,cells,time_record,pnp1[ii],surface_faces,marker_data);
                } 
                    
                // /* Move forward in time */
                // cout << "Time: " << pnp1[ii].t << " dt: " << pnp1[ii].dt  << " x: " <<
                // pnp1[ii].xi[0] << " y: " << pnp1[ii].xi[1] << " z: " << pnp1[ii].xi[2] << 
                // " old cell: " << pn[ii].cellID << " new cell: " << pnp1[ii].cellID << endl;

                /* March forward in time */
                if(svar.eqOrder == 2)
                    pnm1[ii] = pn[ii];
                
                pn[ii] = pnp1[ii];

            }

            nGoing += pnp1[ii].going;
            pnp1[ii].t += pnp1[ii].dt;
        }
        time_record.emplace_back(pnp1);
    }
    
}


#endif