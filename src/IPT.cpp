
#include "Var.h"
#include "Containment.h"
#include "Init.h"
#include "IO.h"
#include "CDFIO.h"
#include "Post.h"
#include "Integrate.h"

#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>  /* defines FILENAME_MAX */
#include <dirent.h>

#ifdef WINDOWS
    #include <direct.h>
    #define GetCurrentDir _getcwd
	#define stat _stat
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
#endif

inline bool file_exists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

int main(int argc, char** argv)
{
    SETT svar;
    FLUID fvar;

    Read_Settings(argv,svar,fvar);

    State pnm1,pn,pnp1;

    vector<State> part_time;

    MESH cells;
    
    vector<vector<SURF>> surface_faces;
    vector<SURF> surface_data;

    Read_BMAP(svar);

    Read_TAUMESH_FACE(svar,cells);

    if(cells.cCentre.empty())
    {
        cout << "Mesh not initialised. Cannot continue." << endl;
        exit(-1);
    }



    Vec_Tree TREE(3,cells.cCentre,20);

    Init_Points(svar,fvar,TREE,cells,pn);
    if(svar.eqOrder == 2)
        pnm1 = pn; 

    pnp1 = pn;

    /* Initialise the surface tracking structure. */
    Init_Surface(svar,cells,surface_faces, surface_data);

    part_time.emplace_back(pnp1);

    if(svar.partout == 1)
    {
        /* Remove previous particle files */
        string file = svar.outdir;
        file.append("/particle_0_scatter.dat");
        string cmd = "exec rm ";
        cmd.append(svar.outdir);
        cmd.append("/particle_*_scatter.dat");

  		if(file_exists(file))
  		{
            if (system(cmd.c_str()))
            {
                cout << "System command failed to execute." << endl;
                cout << "Command: " << cmd << endl;
                exit(-1);
            }
        }

        for(size_t ii = 0; ii < pnp1.size(); ++ii)
        {
            string pfile = svar.outdir;
            pfile.append("/particle_");
            pfile.append(std::to_string(ii));
            pfile.append("_scatter.dat");
            // ofstream fout(pfile);
            
            svar.partfiles.emplace_back(ofstream{pfile});
            
            if(!svar.partfiles[ii].is_open())
            {
                cout << "Failed to open output file. Path: " << pfile << endl;
                exit(-1);
            }
            Write_ASCII_variables(svar.partfiles[ii]);
            Write_ASCII_Scatter_Header(svar.partfiles[ii],ii);
            Write_ASCII_Point(svar.partfiles[ii], svar.scale, pnp1[ii], ii);
        }
    }

    if(svar.streakout == 1)
    {
        string cmd = "exec rm ";
        cmd.append(svar.streakdir);
        cmd.append("/particle_*_streaks.dat");

        string file = svar.streakdir;
        file.append("/particle_0_streaks.dat");

        if(file_exists(file))
  		{
            if (system(cmd.c_str()))
            {
                cout << "System command failed to execute." << endl;
                cout << "Command: " << cmd << endl;
                exit(-1);
            }
        }
    }

    if (svar.cellsout == 1)
    {
        string cmd = "exec rm ";
        cmd.append(svar.celldir);
        cmd.append("/particle_*_cells.dat");
        
        string file = svar.celldir;
        file.append("/particle_0_cells.dat");

        if(file_exists(file))
  		{
            if (system(cmd.c_str()))
            {
                cout << "System command failed to execute." << endl;
                cout << "Command: " << cmd << endl;
                exit(-1);
            }
        }
    }
    
    /* Do an initial sweep of trajectories */
    if(svar.explicit_or_implicit == 0)
    {
        cout << "Starting Newmark-Beta simulation..." << endl;
        Explicit_Integrate(svar,fvar,TREE,cells,pn,pnp1,part_time,surface_faces,surface_data);
    }
    else
    {
        cout << "Starting Implicit simulation..." << endl;
        Implicit_Integrate(svar,fvar,cells,pnm1,pn,pnp1,part_time,surface_faces,surface_data);
    }

    cout << "All particles finished simulating" << endl;

    cout << "Doing more particles for catch data" << endl;
    /* Now do additional trajectories to get the catch efficiency */
    vector<int> marks;
    for(size_t ii = 0; ii < svar.bwrite.size(); ii++)
    {
        if(svar.bwrite[ii] == 1)
        {
            marks.emplace_back(ii);
        }
    }

    SURF const& test_surface = surface_data[marks[0]];
    
    // test_surface.name = svar.bnames[marks[0]];
    // test_surface.marker = svar.markers[marks[0]];

    // test_surface.count = surface_data[marks[0]].count;
    // test_

    /* Prevent the writing of data for particles used to get catch efficiency */
    svar.streakout = 0;
    svar.cellsout = 0;
    svar.partout = 0;


    vector<real> beta_data;
    vector<real> area_data;
    vector<vec<real,3>> pos_data;

    for(size_t ii = 0; ii < test_surface.count; ++ii)
    {
        /* For each impacting point, do a trajectory above and below, to get the catch efficiency */
        real delta = 1e-3;
        vec<real,3> const& xi = test_surface.start_pos[ii];

        vec<real,3> xi_1 = xi;
        vec<real,3> xi_2 = xi;

        xi_1[2] += delta;
        xi_2[2] -= delta;

        State beta_test, beta_testnp1, beta_testnm1;
        beta_test.emplace_back(part(xi_1,part_time[0][test_surface.pIDs[ii]]));
        beta_test.emplace_back(part(xi_2,part_time[0][test_surface.pIDs[ii]]));

        beta_testnp1 = beta_test;

        vector<State> beta_time;
        beta_time.emplace_back(beta_test);

        vector<SURF> beta_marker;
        // vector<vector<SURF>> beta_faces = surface_faces;

        for(size_t jj = 0; jj < svar.markers.size(); jj++)
        {
            beta_marker.emplace_back(SURF(surface_data[jj]));
            // beta_faces.emplace_back(surface_faces[jj]);
        }

        Implicit_Integrate(svar,fvar,cells,beta_testnm1,beta_test,beta_testnp1,part_time,surface_faces,beta_marker);

        /* Make sure they have both hit the same surface */
        if (beta_marker[marks[0]].count == 2)
        {
            /* Get the difference in end area. */
            vec<real,3> diff = beta_marker[marks[0]].end_pos[0] - beta_marker[marks[0]].end_pos[1];

            real dist = diff.norm();

            real beta_i = 2*delta/dist;

            beta_data.emplace_back(beta_i);
            area_data.emplace_back(dist);
            vec<real,3> point = 0.5*(beta_marker[marks[0]].end_pos[0] + beta_marker[marks[0]].end_pos[1]);
            pos_data.emplace_back(point);

        }
        else if (beta_marker[marks[0]].count == 1)
        {
            /* Only one has hit the same marker. So use the original point and the successful one */
            vec<real,3> diff = beta_marker[marks[0]].end_pos[0] - test_surface.end_pos[ii];

            real dist = diff.norm();

            real beta_i = delta/dist;

            beta_data.emplace_back(beta_i);
            area_data.emplace_back(dist);

            vec<real,3> point = 0.5*(beta_marker[marks[0]].end_pos[0] + test_surface.end_pos[ii]);
            pos_data.emplace_back(point);
            
        }
        /* One will have to have made the same marker, surely... */
        
    }


    /* Now get the arc d */





    /* Write failure statistics */
    cout << "Number of failed particles: " << svar.nFailed << 
    "  \% of particles: " << real(svar.nFailed)/real(pnp1.size()) * 100.0 << "\%" << endl;
    //    cout << "Time: " << svar.t << " Particles going: " << nGoing << endl; 
    //    Write_ASCII_Timestep(fout,svar,pnp1);
    //    part_time.emplace_back(pnp1);
     
    // if(fstreak.is_open())
    //     Write_ASCII_Streaks(fstreak,svar,part_time);

    // fstreak.close();
    // fout.close();

    // Write_ASCII_Cells(svar,cells,part_time);
    cout << "Post processing" << endl;
    Post_Process(svar,part_time,surface_faces,surface_data);

    cout << "Writing catch efficiency data" << endl;
    Write_Aerofoil_Catch_Efficiency(svar,cells,marks[0],beta_data,area_data,pos_data);

    cout << "Writing surface impact data" << endl;
    Write_ASCII_Impacts(svar,cells,surface_faces);


    cout << "Finished!" << endl;
    return 0;
}