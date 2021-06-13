
#include "Var.h"
#include "Init.h"
#include "IO.h"
#include "Containment.h"
#include "CDFIO.h"
#include "Integrate.h"



int main(int argc, char** argv)
{
    SETT svar;
    FLUID fvar;

    Read_Settings(argv,svar,fvar);

    State pnm1,pn,pnp1;

    vector<State> part_time;

    MESH cells;

    Read_TAUMESH_FACE(svar,cells);

    if(cells.cCentre.empty())
    {
        cout << "Mesh not initialised. Cannot continue." << endl;
        exit(-1);
    }

    Vec_Tree TREE(3,cells.cCentre,20);

    Init_Points(svar,fvar,TREE,cells,pn);
    pnm1 = pn; 
    pnp1 = pn;

    part_time.emplace_back(pnp1);

    ofstream fout(svar.outfile);
    ofstream fstreak(svar.streakfile);

    if(!fout.is_open())
    {
        cout << "Failed to open output file. Path: " << svar.outfile << endl;
        exit(-1);
    }

    if(!fstreak.is_open())
    {
        cout << "Failed to open output streak file. Path: " << svar.streakfile << endl;
        exit(-1);
    }

    cout << "Opened output file, writing header..." << endl;
    Write_ASCII_header(fout);
    
    Write_ASCII_Timestep(fout, svar, pnp1);

    cout << "Starting simulation..." << endl;
    uint nGoing = pnp1.size();
    if(svar.max_frame > 0)
    {
        while (svar.frame < svar.max_frame && nGoing != 0)
        {
            nGoing = Implicit_Integrate(svar,fvar,TREE,cells,pnm1,pn,pnp1,part_time, fout); 
            // cout << "Time: " << svar.t << " Particles going: " << nGoing << endl; 
            // Write_ASCII_Timestep(fout,svar,pnp1);
            // part_time.emplace_back(pnp1);
        }
    }
    else
    {
        while (nGoing != 0)
        {
           nGoing = Implicit_Integrate(svar,fvar,TREE,cells,pnm1,pn,pnp1,part_time,fout); 
        //    cout << "Time: " << svar.t << " Particles going: " << nGoing << endl; 
        //    Write_ASCII_Timestep(fout,svar,pnp1);
        //    part_time.emplace_back(pnp1);
        }
    }

    
    if(fstreak.is_open())
        Write_ASCII_Streaks(fstreak,svar,part_time);

    fstreak.close();
    fout.close();

    cout << "Finished!" << endl;
    return 0;
}