
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

int isDirectoryEmpty(string const& dirname) {
  int n = 0;
  struct dirent *d;
  DIR *dir = opendir(dirname.c_str());
  if (dir == NULL) //Not a directory or doesn't exist
    return 1;
  while ((d = readdir(dir)) != NULL) {
    if(++n > 2)
      break;
  }
  closedir(dir);
  if (n <= 2) //Directory Empty
    return 1;
  else
    return 0;
}

void Run_Simulation(SETT& svar, FLUID& fvar, MESH const& cells, Vec_Tree const& TREE,
             vector<int> const& marks, vector<real>& avg_beta, vector<real>& avg_area, vector<int>& avg_count)
{
    State pnm1,pn,pnp1;

    vector<vector<part>> part_time;

    vector<SURF> surface_data;

    Init_Points(svar,fvar,TREE,cells,pn);

    if(svar.eqOrder == 2)
        pnm1 = pn; 

    pnp1 = pn;

    part_time = vector<vector<part>>(pnp1.size(),vector<part>());

    /* Initialise the surface tracking structure. */
    Init_Surface(svar,cells, surface_data);


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
            Write_ASCII_Point(svar.partfiles[ii], svar.scale, pnp1[ii]);
        }
    }

    if(svar.streakout == 1)
    {
        string cmd = "exec rm ";
        cmd.append(svar.streakdir);
        cmd.append("/particle_*_streaks.dat");

        if(!isDirectoryEmpty(svar.streakdir))
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
        
        if(!isDirectoryEmpty(svar.celldir))
  		{
            if (system(cmd.c_str()))
            {
                cout << "System command failed to execute." << endl;
                cout << "Command: " << cmd << endl;
                exit(-1);
            }
        }
    }

    for(size_t ii = 0; ii < pnp1.size(); ++ii)
        part_time[ii].emplace_back(pnp1[ii]);

    svar.nSuccess = 0;
    svar.nFailed = 0;
    
    /* Do an initial sweep of trajectories */
    if(svar.explicit_or_implicit == 0)
    {
        cout << "Starting Newmark-Beta simulation. Particle diameter: " << fvar.d_0*1e6 << endl;
        cout << "Number of particles: " << pnp1.size() << endl;
        #pragma omp parallel for schedule(dynamic) shared(pn,pnp1,part_time)
        for(size_t ii = 0; ii < pnp1.size(); ++ii)
            Explicit_Integrate(svar,fvar,TREE,cells,ii,pn[ii],pnp1[ii],part_time[ii],surface_data);
    }
    else
    {
        cout << "Starting implicit simulation. Particle diameter: " << fvar.d_0*1e6 << endl;
        cout << "Number of particles: " << pnp1.size() << endl;
        #pragma omp parallel for schedule(dynamic) shared(pnm1,pn,pnp1,part_time)
        for(size_t ii = 0; ii < pnp1.size(); ++ii)
            Implicit_Integrate(svar,fvar,cells,ii,pnm1[ii],pn[ii],pnp1[ii],part_time[ii],surface_data);
    }

    /* Prevent the writing of data for particles used to get catch efficiency */
    // svar.streakout = 0;
    svar.cellsout = 0;
    svar.partout = 0;

    cout << "All particles finished simulating" << endl;

    /* Now do additional trajectories to get the catch efficiency */
    cout << "Doing more particles for catch data" << endl;

    vector<vector<real>> beta_data;
    vector<vector<real>> area_data;
    vector<vector<vec<real,3>>> pos_data;

    Get_Catch_Data(svar,fvar,TREE,cells,surface_data,marks,part_time,
                    beta_data,area_data,pos_data);

    svar.ntot_success += svar.nSuccess;
    svar.ntot_fail += svar.nFailed;

    /* Write failure statistics */
    cout << "Number of failed particles: " << svar.nFailed << 
    "  \% of particles: " << real(svar.nFailed)/real(svar.nFailed+svar.nSuccess) * 100.0 << "\%" << endl;
    cout << "Total number of failed particles: " << svar.ntot_fail << 
    " total \% of particles: " << real(svar.ntot_fail)/real(svar.ntot_fail+svar.ntot_success) * 100.0 << "\%" << endl;

    
    // auto ptr = std::max_element(pnp1.begin(), pnp1.end(), 
    //         [](part const& a, part const& b){return a.nIters > b.nIters;});
    // uint max_integ = ptr->nIters;

    // cout << "Maximum number of cells passed through: " << max_integ << endl;

    cout << "Post processing" << endl;
    Post_Process(svar, marks, beta_data, area_data, surface_data);

    // cout << "Writing catch efficiency data" << endl;
    // Write_Aerofoil_Catch_Efficiency(svar,fvar,cells,marks,beta_data,area_data,pos_data);

    cout << "Writing surface impact data" << endl;
    Write_ASCII_Impacts(svar, fvar, cells, surface_data, beta_data);

    for(size_t ii = 0; ii < marks.size(); ++ii)
    {
        avg_beta[ii] = surface_data[marks[ii]].marker_beta;
        avg_area[ii] = surface_data[marks[ii]].marker_area;
        avg_count[ii] = surface_data[marks[ii]].marker_count;
    }
}

int main(int argc, char** argv)
{
    SETT svar;
    FLUID fvar;

    Read_Settings(argv,svar,fvar);

    MESH cells;

    Read_BMAP(svar);

    Read_TAUMESH_FACE(svar,cells);

    if(cells.cCentre.empty())
    {
        cout << "Mesh not initialised. Cannot continue." << endl;
        exit(-1);
    }

    Vec_Tree TREE(3,cells.cCentre,10);

    vector<int> marks;
    for(size_t ii = 0; ii < svar.bwrite.size(); ii++)
    {
        if(svar.bwrite[ii] == 1)
        {
            marks.emplace_back(ii);
        }
    }

    if(svar.particle_sweep == 0)
    {
        vector<real> avg_beta(marks.size(),0.0);
        vector<real> avg_area(marks.size(),0.0);
        vector<int> avg_count(marks.size(),0.0);
        Run_Simulation(svar,fvar,cells,TREE, marks, avg_beta, avg_area, avg_count);
    }
    else
    {
        vector<vector<real>> beta_data(marks.size(),vector<real>(svar.nParticle_diams,0.0));
        vector<vector<real>> area_data(marks.size(),vector<real>(svar.nParticle_diams,0.0));
        vector<vector<real>> particle_data(marks.size(),vector<real>(svar.nParticle_diams,0.0));
        vector<vector<int>>  count_data(marks.size(),vector<int>(svar.nParticle_diams,0.0));

        vector<int> total_count(marks.size(),0);

        /* Do a sweep of particle diameters */
        for(int jj = 0; jj < svar.nParticle_diams; ++jj)
        {
            fvar.d_0 = fvar.diams[jj];
            vector<real> avg_beta(marks.size());
            vector<real> avg_area(marks.size());
            vector<int> avg_count(marks.size());
            
            cout << endl << endl;
            Run_Simulation(svar,fvar,cells,TREE, marks, avg_beta, avg_area, avg_count);

            for(size_t ii = 0; ii < marks.size(); ++ii)
            {
                beta_data[ii][jj] = avg_beta[ii];
                area_data[ii][jj] = avg_area[ii];
                particle_data[ii][jj] = fvar.d_0*1e6;
                count_data[ii][jj] = avg_count[ii];
                total_count[ii] += avg_count[ii];
            }
        }

        /* Output graph of beta vs particle size */
        string file = svar.surfacefile;
        file.append("_");
        file.append("Beta_Sweep.dat");
        
        ofstream fout(file);

        if(!fout.is_open())
        {
            cout << "Couldn't open the particle sweep file." << endl;
            exit(-1);
        }

        fout << "TITLE=\"Catch Efficiency vs Particle size\"\n";
        fout << "VARIABLES= \"Particle diameter (microns)\" \"beta\" \"area\""
        << " \"count\" \"normalized count\" \"cumulative count\" \"normalized cumulative count\"\n";

        

        for(size_t ii = 0; ii < marks.size(); ++ii)
        {
            fout << "ZONE T=\"" << svar.bnames[marks[ii]] << "\"\n";

            size_t w = 15;
            size_t w2 = 9;
            size_t preci = 6;
            fout << std::left << std::scientific << std::setprecision(preci);
            fout << std::setw(1);

            int cum_count = 0;
            real norm_cum_count = 0;

            for(int jj = 0; jj < svar.nParticle_diams; ++jj)
            {
                real norm_count = 0;
                if(total_count[ii] > 0)
                    norm_count = real(count_data[ii][jj])/real(total_count[ii]);

                cum_count += count_data[ii][jj];
                
                if(total_count[ii] > 0)
                    norm_cum_count += real(count_data[ii][jj])/real(total_count[ii]);
                
                fout << std::setw(w) << particle_data[ii][jj] << std::setw(w) << 
                beta_data[ii][jj] << std::setw(w) << area_data[ii][jj]
                 << std::setw(w2) << count_data[ii][jj] << std::setw(w) << norm_count
                 << std::setw(w2) << cum_count << std::setw(w2) << norm_cum_count << endl;

                
            }

            // fout << std::setw(w) << particle_data[ii].back() << std::setw(w) << 
            //     beta_data[ii].back() << std::setw(w) << area_data[ii].back()
            //         << std::setw(w2) << count_data[ii].back() << std::setw(w) << real(count_data[ii].back())/real(total_count[ii])
            //         << std::setw(w2) << cum_count << std::setw(w2) << norm_cum_count << endl;

        }

        fout.close();

    }


    cout << "Finished!" << endl;
    return 0;
}