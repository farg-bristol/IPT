#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "Var.h"
#include "Containment.h"

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
    vec<real,3> Vdiff = pnp1.cellV - pnp1.v ;

    vec<real,3> res = Explicit_AeroForce(Vdiff,fvar,pnp1);
    res /= pnp1.mass;
    // cout << res[0] << "  " << res[1] << "  " << res[2] << endl;

    pnp1.xi = pn.xi+dt*pn.v+dt2*(c*pn.acc+d*res);
    pnp1.v =  pn.v +dt*(a*pn.acc+b*res);
    pnp1.acc = res;
}


int Explicit_Integrate(SETT& svar, FLUID const& fvar, Vec_Tree const& TREE, MESH const& cells, State& pn, State& pnp1)
{
    real const& dt = svar.dt;
    real const& dt2 = dt*dt;
    real& t = svar.t;

    const real a = 1 - svar.gamma;
	const real b = svar.gamma;
	const real c = 0.5*(1-2*svar.beta);
	const real d = svar.beta;

    uint nGoing = 0;

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
                
                /* Move forward in time */
                pn[ii] = pnp1[ii];
            }

            nGoing += pnp1[ii].going;
        }
        if(nGoing == 0)
        {
            cout << "All particles finished simulating" << endl;
            break;
        }

        t += dt;
    }

    return nGoing;
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
    vec<real,3> const g(0,0,0.0);

    real res = Implicit_AeroForce(Vdiff, fvar, pnp1);
    
    pnp1.acc = res;
    pnp1.v = (pn.v + pnp1.dt*res*pnp1.faceV + pnp1.dt*g)/(1 + pnp1.dt*res);
    pnp1.xi = pn.xi+0.5*pnp1.dt*(pnp1.v + pn.v); 
}

void Implicit_BFD2(FLUID const& fvar, real const& dt, real const& dtm1, part const& pnm1, part const& pn, part& pnp1)
{
    vec<real,3> const Vdiff = pnp1.faceV - pnp1.v ;
    vec<real,3> const g(0,0,0.0);

    real res = Implicit_AeroForce(Vdiff, fvar, pnp1);

    pnp1.acc = res;
    /* Second order velocity calculation */
    pnp1.v = ((dt+dtm1)*(dt*dtm1*(res*pnp1.faceV+g) + (dt+dtm1)* pn.v) - dt*dt*pnm1.v)
                /(dtm1*((2*dt+dtm1)+dt*(dt+dtm1)*res));
                
    // pnp1.xi = pn.xi-pnp1.dt*(pnp1.v + pn.v)+0.5*pn.dt*(pn.v + pnm1.v); 
    pnp1.xi = pn.xi+0.5*pnp1.dt*(pnp1.v + pn.v); /* Just a first order space integration. */
}

int Implicit_Integrate(SETT& svar, FLUID const& fvar, Vec_Tree const& TREE, MESH const& cells, 
                        State& pnm1, State& pn, State& pnp1, vector<State>& time_record, ofstream& fout)
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
            int min_iter = 2;
            real error = 1.0;
            while(error > -7.0 && iter < min_iter)
            {
                // cout << "Iteration: " << iter << endl; 
                vec<real,3> temp = pnp1[ii].xi;
                // Implicit_BFD1(fvar, pn[ii], pnp1[ii]);
                real const dt = pnp1[ii].dt;
                real const dtm1 = pnp1[ii].dt;
                
                Implicit_BFD2(fvar, dt, dtm1, pnm1[ii], pn[ii], pnp1[ii]);

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

            if(pnp1[ii].cellID < 0 )
            {
                pnp1[ii].going = 0;
                cout << "Particle " << ii << " stopped iterating" << endl;
            } 
            
            /* Move forward in time */
            cout << "Time: " << pnp1[ii].t << " dt: " << pnp1[ii].dt  << " x: " <<
            pnp1[ii].xi[0] << " y: " << pnp1[ii].xi[1] << " z: " << pnp1[ii].xi[2] << 
            " old cell: " << pn[ii].cellID << " new cell: " << pnp1[ii].cellID << endl;

            /* March forward in time */
            pnm1[ii] = pn[ii];
            pn[ii] = pnp1[ii];

        }

        nGoing += pnp1[ii].going;
        pnp1[ii].t += pnp1[ii].dt;
    }

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
                int min_iter = 5;
                real error = 1.0;
                while(error > -7.0 && iter < min_iter)
                {
                    // cout << "Iteration: " << iter << endl; 
                    vec<real,3> temp = pnp1[ii].xi;
                    // Implicit_BFD1(fvar, pn[ii], pnp1[ii]);
                    real const dt = pnp1[ii].dt;
                    real dtm1 = pnp1[ii].dt;
                    if(iter > 0)
                        dtm1 = pn[ii].dt;

                    Implicit_BFD2(fvar, dt, dtm1, pnm1[ii], pn[ii], pnp1[ii]);

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

                if(pnp1[ii].cellID < 0 )
                {
                    pnp1[ii].going = 0;
                    cout << "Particle " << ii << " stopped iterating" << endl;
                } 
                
                /* Move forward in time */
                cout << "Time: " << pnp1[ii].t << " dt: " << pnp1[ii].dt  << " x: " <<
                pnp1[ii].xi[0] << " y: " << pnp1[ii].xi[1] << " z: " << pnp1[ii].xi[2] << 
                " old cell: " << pn[ii].cellID << " new cell: " << pnp1[ii].cellID << endl;

                /* March forward in time */
                pnm1[ii] = pn[ii];
                pn[ii] = pnp1[ii];

            }

            nGoing += pnp1[ii].going;
            pnp1[ii].t += pnp1[ii].dt;
        }
        time_record.emplace_back(pnp1);
        Write_ASCII_Timestep(fout,svar,pnp1);
    }
    cout << "All particles finished simulating" << endl;
    return nGoing;
}


#endif