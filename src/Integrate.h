#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "Var.h"
#include "Containment.h"

void Terminate_Particle(SETT const& svar, MESH const& cells, vector<part> const& time_record, size_t const& ii,
    part& pnp1, vector<SURF>& surface_data )
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
                #pragma omp critical
                {
                    surface_data[surf].marker_count++;
                    surface_data[surf].marker_pIDs.emplace_back(ii);
                    surface_data[surf].end_pos.emplace_back(pnp1.xi);
                    surface_data[surf].start_pos.emplace_back(time_record[0].xi);
                }
                auto const index2 =
                    std::find(surface_data[surf].faceIDs.begin(), surface_data[surf].faceIDs.end(), face);

                if(index2 != surface_data[surf].faceIDs.end())
                {
                    size_t fIndex = index2 - surface_data[surf].faceIDs.begin();
                    #pragma omp critical
                    {
                        surface_data[surf].face_count[fIndex]++;
                        surface_data[surf].impacted_face.emplace_back(fIndex);
                    }
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
        Write_ASCII_Streaks(svar,time_record,pnp1);
    }

    if(svar.cellsout == 1)
    {
        Write_ASCII_Cells(svar,cells,time_record);
    }
}


inline real GetCd(real const& Re)
{
	return (1.0+0.197*pow(Re,0.63)+2.6e-04*pow(Re,1.38))*(24.0/(Re+0.0000001));

	// if( Re <= 1000.0)
	// 	return (24.0/Re)*(1+(1.0/6.0)*pow(Re,2.0/3.0));
	// else 
	// 	return 0.424;
}

inline vec<real,3> Explicit_AeroForce(vec<real,3> const& Vdiff, FLUID const& fvar, part const& pi)
{
	real const Re = 2.0*pi.cellRho*Vdiff.norm()*pi.d/fvar.mu_g;
	real const Cdi = GetCd(Re);

	// std::cout << "Reynolds: " << Re  << " Cd: " << Cdi << std::endl;
	return (0.5*pi.cellRho*Vdiff.norm()*Vdiff*Cdi*pi.A);
}

inline void Explicit_Newmark_Beta(FLUID const& fvar, real const& a, real const& b, real const& c, real const& d, 
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


void Explicit_Integrate(SETT& svar, FLUID const& fvar, Vec_Tree const& TREE, MESH const& cells,
                        size_t const& ii, part& pn, part& pnp1, 
                        vector<part>& time_record, vector<SURF>& marker_data)
{
    real const& dt = svar.dt;
    real const& dt2 = dt*dt;
    real& t = svar.t;

    const real a = 1 - svar.gamma;
	const real b = svar.gamma;
	const real c = 0.5*(1-2*svar.beta);
	const real d = svar.beta;

    while(pnp1.going != 0)
    {
        for (size_t step = 0; step < svar.nsteps; ++step)
        {
            if(pnp1.going == 1)
            {
                /* Check containing cell */
                FindCell(svar,TREE, cells, pnp1);

                /* Perform one step of NB, then freeze containment for stability. */
                Explicit_Newmark_Beta(fvar, a,b,c,d,dt,dt2,pn,pnp1);
                FindCell(svar,TREE, cells, pnp1);

                uint iter = 0;
                uint min_iter = 3;
                real error = 1.0;
                while(error > -7.0 && iter < min_iter)
                {
                    vec<real,3> temp = pnp1.xi;
                    /* Get the aerodynamic force */
                    Explicit_Newmark_Beta(fvar, a,b,c,d,dt,dt2, pn,pnp1);

                    error = log10((pnp1.xi-temp).norm());

                    if(iter < svar.max_iter)
                        iter++;
                    else
                        break;
                }

                if(pnp1.going == 0)
                {
                    /* Use face ID to test if it's crossed a boundary */
                    if( pnp1.faceID == -2 )
                    {
                        Terminate_Particle(svar,cells,time_record,ii,pnp1,marker_data);
                        svar.nSuccess++;
                    }
                    else
                    {
                        svar.nFailed++;
                    }
                }
                else if (pnp1.xi[0] > svar.max_x)
                { /* If trajectory is past the x-limit, then stop the trajectory */
                    Terminate_Particle(svar, cells, time_record, ii, pnp1, marker_data);
                    svar.nSuccess++;
                }

                /* Move forward in time */
                pn = pnp1;
                t += dt;
            }
            else
            {
                break;
            }
        }

       
        if(pnp1.going == 1)
        {
            if(svar.partout == 1)
                Write_ASCII_Point(svar.partfiles[pnp1.partID], svar.scale, pnp1);
        }
        else
        {
            break;
        }
        

        if(svar.frame < svar.max_frame)
            svar.frame++;
        else
        {
            /* Wrap up the simulation */      
            if(pnp1.going == 1)
            {
                Terminate_Particle(svar,cells,time_record,ii,pnp1,marker_data);
                svar.nSuccess++;
            }
        }
        
        /* Only need the information if streaks or cells are output */
        if(svar.cellsout == 1 || svar.streakout == 1)
            time_record.emplace_back(pnp1);
    }
}


inline real Implicit_AeroForce(vec<real,3> const& Vdiff, FLUID const& fvar, part const& pi)
{
    real const Re = 2.0*pi.faceRho*Vdiff.norm()*pi.d/fvar.mu_g;
	real const Cd = GetCd(Re);


    return Vdiff.norm() * (3.0*Cd*pi.faceRho)/(4.0 * pi.d * fvar.d_rho);
}

inline void Implicit_BFD1(FLUID const& fvar, part const& pn, part& pnp1)
{
    vec<real,3> const Vdiff = pnp1.faceV - pnp1.v ;
    // vec<real,3> const g(0,0,-9.81);
    vec<real,3> const g(0.0,0.0,0.0);

    real res = Implicit_AeroForce(Vdiff, fvar, pnp1);
    
    pnp1.acc = res;
    pnp1.v = (pn.v + pnp1.dt*res*pnp1.faceV + pnp1.dt*g)/(1.0 + pnp1.dt*res);
    pnp1.xi = pn.xi+0.5*pnp1.dt*(pnp1.v + pn.v); 
}

inline void Implicit_BFD2(FLUID const& fvar, real const& relax, real const& dt, real const& dtm1, 
                        part const& pnm1, part const& pn, part& pnp1)
{
    vec<real,3> const Vdiff = pnp1.faceV - pnp1.v ;
    // vec<real,3> const g(0,0,-9.81);
    vec<real,3> const g(0.0,0.0,0.0);

    real res = Implicit_AeroForce(Vdiff, fvar, pnp1);

    pnp1.acc = res;
    /* Second order velocity calculation */
    pnp1.v = ((dt+dtm1)*(dt*dtm1*(relax*res*pnp1.faceV+g) + (dt+dtm1)* pn.v) - dt*dt*pnm1.v)
                /(dtm1*((2*dt+dtm1)+dt*(dt+dtm1)*res));

                
    // pnp1.xi = pn.xi-pnp1.dt*(pnp1.v + pn.v)+0.5*pn.dt*(pn.v + pnm1.v); 
    pnp1.xi = pn.xi+0.5*pnp1.dt*(pnp1.v + pn.v); /* Just a first order space integration. */
}

void Implicit_Integrate(SETT& svar, FLUID const& fvar, MESH const& cells, size_t const& ii, 
                        part& pnm1, part& pn, part& pnp1, 
                        vector<part>& time_record, vector<SURF>& marker_data)
{
    /* Do the first step */

    /* Find intersecting face for np1 */
    // pnp1.faceV = pnp1.cellV; /* Start by guessing the face as the cell*/
    // pnp1.faceRho = pnp1.cellRho;
    FindFace(svar, cells, pn, pnp1);
    /* Now c+1 is known, use the central difference face properties*/
    pnp1.faceV = 0.5*(pnp1.cellV + pn.cellV);
    pnp1.faceRho = 0.5*(pnp1.cellRho + pn.cellRho);

    // int iter = 0;
    // int min_iter = 3;
    // real error = 1.0;
    real relaxation = 1.0; /* Relaxation factor for 2nd order */
    for(size_t iter = 0; iter < 10; iter++)
    // while(error > -7.0 && iter < min_iter)
    {
        relaxation = (real(iter)+1.0)/6.0;
        if (relaxation > 1.0)
            relaxation = 1.0;
        // cout << "Iteration: " << iter << endl; 
        // vec<real,3> temp = pnp1.xi;

        if (svar.eqOrder == 2)
        {
            real const dt = pnp1.dt;
            real const dtm1 = pnp1.dt;
            Implicit_BFD2(fvar, relaxation, dt, dtm1, pnm1, pn, pnp1);
        }
        else
            Implicit_BFD1(fvar, pn, pnp1);

        FindFace(svar, cells, pn, pnp1);
        /* Now c+1 is known, use the central difference face properties*/
        pnp1.faceV = 0.5*(pnp1.cellV + pn.cellV);
        pnp1.faceRho = 0.5*(pnp1.cellRho + pn.cellRho);

        // error = log10((pnp1.xi-temp).norm());

        // if(iter < svar.max_iter)
        //     iter++;
        // else
        //     break;
    }

    if(svar.partout == 1)
        Write_ASCII_Point(svar.partfiles[pnp1.partID], svar.scale, pnp1);

    if (pnp1.going == 0 || pnp1.v.norm() > 1e3
        || (pnp1.xi - pn.xi).norm() > cells.maxlength )
    {
        pnp1.failed = 1;
        Terminate_Particle(svar, cells, time_record, ii, pnp1, marker_data);
        svar.nFailed++;
    }
    else if (pnp1.cellID < 0 || pnp1.xi[0] > svar.max_x)
    {      
        Terminate_Particle(svar, cells, time_record, ii, pnp1, marker_data);
        svar.nSuccess++;
    }

    /* Move forward in time */
    // cout << "Time: " << pnp1.t << " dt: " << pnp1.dt  << " x: " <<
    // pnp1.xi[0] << " y: " << pnp1.xi[1] << " z: " << pnp1.xi[2] << 
    // " old cell: " << pn.cellID << " new cell: " << pnp1.cellID << endl;

    /* March forward in time */
    if(svar.eqOrder == 2)
        pnm1 = pn;

    pn = pnp1;
    pnp1.t += pnp1.dt;
    

    /* Only need the information if streaks or cells are output */
    if(svar.cellsout == 1 || svar.streakout == 1)
        time_record.emplace_back(pnp1);

    while(pnp1.going != 0)
    {
       
        /* Find intersecting face for np1 */
        // pnp1.faceV = pnp1.cellV; /* Start by guessing the face as the cell*/
        // pnp1.faceRho = pnp1.cellRho;
        
        if(pnp1.d < 1.0e-5)
        {   /* If particle is small, use cell velocity to start */
            pnp1.v = pnp1.cellV;
        }
        else if (pnp1.d < 100.0e-6)
        {   /* If particle is heavier, use sliding scale up to 100 microns */
            pnp1.v = (1.0 - (pnp1.d-1e-05)/(9e-5)) * pnp1.cellV + 
            (pnp1.d-1e-05)/(9e-5) * pnp1.v;
        }


        FindFace(svar, cells, pn, pnp1);
        /* Now c+1 is known, use the central difference face properties*/
        pnp1.faceV = 0.5*(pnp1.cellV + pn.cellV);
        pnp1.faceRho = 0.5*(pnp1.cellRho + pn.cellRho);

        // int iter = 0;
        // int min_iter = 3;
        // real error = 1.0;
        real relaxation = 1.0;
        for(size_t iter = 0; iter < 10; iter++)
        // while(error > -7.0 && iter < min_iter)
        {
            relaxation = (real(iter)+1.0)/6.0;
            if (relaxation > 1.0)
                relaxation = 1.0;

            // cout << "Iteration: " << iter << endl; 
            // vec<real,3> temp = pnp1.xi;
            
            if(svar.eqOrder == 2)
            {
                real const dt = pnp1.dt;
                real dtm1 = pnp1.dt;
                if(iter > 0)
                    dtm1 = pn.dt;
                Implicit_BFD2(fvar, relaxation, dt, dtm1, pnm1, pn, pnp1);
            }
            else
                Implicit_BFD1(fvar, pn, pnp1);


            FindFace(svar, cells, pn, pnp1);
            /* Now c+1 is known, use the central difference face properties*/
            pnp1.faceV = 0.5*(pnp1.cellV + pn.cellV);
            pnp1.faceRho = 0.5*(pnp1.cellRho + pn.cellRho);

            // error = log10((pnp1.xi-temp).norm());

            // if(iter < svar.max_iter)
            //     iter++;
            // else
            //     break;
        }

        if(svar.partout == 1)
            Write_ASCII_Point(svar.partfiles[pnp1.partID], svar.scale, pnp1);

        if ( pnp1.going == 0 || pnp1.v.norm() > 1e3
                || (pnp1.xi - pn.xi).norm() > cells.maxlength )
        {
            pnp1.failed = 1;
            Terminate_Particle(svar, cells, time_record, ii, pnp1, marker_data);
            svar.nFailed++;
            break;
        }
        else if ( pnp1.cellID < 0 || pnp1.xi[0] > svar.max_x ||
                pnp1.nIters > svar.max_frame)
        {      
            Terminate_Particle(svar, cells, time_record, ii, pnp1, marker_data);
            svar.nSuccess++;
            break;
        }

        // /* Move forward in time */
        // cout << "Time: " << pnp1.t << " dt: " << pnp1.dt  << " x: " <<
        // pnp1.xi[0] << " y: " << pnp1.xi[1] << " z: " << pnp1.xi[2] << 
        // " old cell: " << pn.cellID << " new cell: " << pnp1.cellID << endl;

        /* March forward in time */
        if(svar.eqOrder == 2)
            pnm1 = pn;
        
        pn = pnp1;
        pnp1.t += pnp1.dt;
        pnp1.nIters++;
        
        /* Only need the information if streaks or cells are output */
        if(svar.cellsout == 1 || svar.streakout == 1)
            time_record.emplace_back(pnp1);

    }
    
}


#endif