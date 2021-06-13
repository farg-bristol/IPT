#ifndef IO_H
#define IO_H

#include "Var.h"


std::string ltrim(const std::string &s)
{
    size_t start = s.find_first_not_of(WHITESPACE);
    return (start == std::string::npos) ? "" : s.substr(start);
}
 
std::string rtrim(const std::string &s)
{
    size_t end = s.find_last_not_of(WHITESPACE);
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

string Get_Parameter_Value(string const& line)
{
    size_t pos = line.find(":");
    size_t end = line.find("#",pos+1); /* Check if a comment exists on the line */

    if (end != string::npos)
    {
        string value = line.substr(pos + 1, (end-pos+2) );
        return ltrim(rtrim(value));
    }

    string value = line.substr(pos + 1);
    return ltrim(rtrim(value));
}

void Get_String(string const& line, string const& param, string &value)
{
    if(line.find(param) != string::npos)
    {
        value = Get_Parameter_Value(line);
    }
}

template<typename T>
void Get_Number(string const& line, string const& param, T &value)
{
    if(line.find(param) != string::npos)
    {
        string temp = Get_Parameter_Value(line);
        std::istringstream iss(temp);
        iss >> value;
    }
}

void Get_Vector(string const& line, string const& param, vec<real,3> &value)
{
    if(line.find(param) != string::npos)
    {
        string temp = Get_Parameter_Value(line);
        std::istringstream iss(temp);
        real a, b, c;
        iss >> a; iss >> b; iss >> c;
        
        value = vec<real,3>(a,b,c);
    }
}

matrix<real,3> GetRotationMat(vec<real,3>& angles)
{

    matrix<real,3> rotx, roty, rotz;
    rotx = matrix<real,3>(1.0, 0.0            , 0.0           ,
                0.0, cos(angles[0]) , sin(angles[0]),
                0.0, -sin(angles[0]), cos(angles[0]));

    roty = matrix<real,3>(cos(angles[1]) , 0.0 , -sin(angles[1]),
                0.0            , 1.0 , 0.0            ,
                sin(angles[1]) , 0.0 , cos(angles[1]));

    rotz = matrix<real,3>(cos(angles[2]) , sin(angles[2]) , 0.0 ,
                -sin(angles[2]), cos(angles[2]) , 0.0 ,
                0.0            , 0.0            , 1.0 );

    return rotx*roty*rotz;
}


void Read_Settings(char** argv, SETT &svar, FLUID& fvar)
{
    ifstream fin(argv[1]);

    string line;
    while (getline(fin,line))
    {
        line = ltrim(line);
        if(line[0] == '#') /* Skip commented line */
            continue;

        /* Get the file IO */
        Get_String(line, "Primary grid filename", svar.taumesh);
        Get_String(line, "Boundary mapping filename", svar.taubmap);
        Get_String(line, "Restart-data prefix", svar.tausol);
        Get_String(line, "OpenFOAM folder", svar.foamdir);
        Get_String(line, "OpenFOAM solution folder", svar.foamsol);
        Get_String(line, "Particle output filename", svar.outfile);
        Get_String(line, "Particle streak output filename", svar.streakfile);
        Get_Number(line, "OpenFOAM binary (0/1)", svar.isBinary);
        Get_Number(line, "Label size (32/64)", svar.labelSize);
        Get_Number(line, "Scalar size (32/64)", svar.scalarSize);

        /* Get the fluid data */
        Get_Number(line, "Grid scale", svar.scale);
        Get_Number(line, "Reference density", fvar.g_rho);
        Get_Number(line, "Reference dispersed density", fvar.d_rho);
        Get_Number(line, "Sutherland reference viscosity", fvar.mu_g);

        /* Get Simulation settings */
        Get_Number(line, "Particle timestep", svar.dt);
        Get_Number(line, "Particle frame time", svar.framet);
        Get_Number(line, "Newmark-Beta iteration limit", svar.max_iter);
        Get_Number(line, "Particle maximum frame time", svar.max_frame);

        /* Starting area conditions */
        Get_Vector(line, "Seed disk center", svar.discPos);
        Get_Vector(line, "Seed disk angles", svar.discAngles);
        Get_Number(line, "Seed disk diameter", svar.discDiam);
        Get_Number(line, "Seed disk spacing", svar.discSpace);
        Get_Number(line, "Particle diameter", fvar.d_0);
        Get_Vector(line, "Particle velocity", fvar.v_0);
    }

    fin.close();

    if(svar.taumesh.empty())
    {
        cout << "Input TAU file not defined." << endl;
        exit(-1);
    }

    if(svar.taubmap.empty())
    {
        cout << "Input TAU bmap file not defined." << endl;
        exit(-1);
    }

    if(svar.tausol.empty())
    {
        cout << "Input TAU solution file not defined." << endl;
        exit(-1);
    }

    if(svar.outfile.empty())
    {
        cout << "Output particle filename not defined." << endl;
        exit(-1);
    }

    svar.streakfile = svar.outfile;
    size_t pos = svar.streakfile.find_last_of(".");
    svar.streakfile.insert(pos,"_streak");

    if(svar.discPos[0] == -500000)
    {
        cout << "Seed disk position has not been initialised." << endl;
        exit(-1);
    }

    svar.rotate = GetRotationMat(svar.discAngles);
    svar.transp = svar.rotate.transpose();


    /* Find a timestep that gets a perfect frame interval */
    int nsteps = ceil(svar.framet/svar.dt);
    svar.nsteps = nsteps;
    svar.dt = svar.framet / real(nsteps);

    // if(svar.foamdir.empty())
    // {
    //     cout << "Output OpenFOAM directory not defined." << endl;
    //     exit(-1);
    // }

} 

void Write_ASCII_header(ofstream& fp)
{
	string variables;

	variables = 
"\"X\", \"Y\", \"Z\", \"v\", \"a\", \"ptID\", \"Cell_V\", \"Cell_Rho\", \"Cell_ID\"";

	fp << "VARIABLES = " << variables << "\n";
}

void Write_ASCII_Point(ofstream& fp, real const& scale, part const& pnp1, size_t ptID)
{
    const static uint width = 15;

    for(uint dim = 0; dim < 3; ++dim)
        fp << setw(width) << pnp1.xi[dim]/scale;
    
    fp << setw(width) << pnp1.v.norm(); 

    fp << setw(width) << pnp1.acc.norm(); 

    fp << setw(width) << ptID;

    fp << setw(width) << pnp1.cellV.norm();

    fp << setw(width) << pnp1.cellRho;
    fp << setw(width) << pnp1.cellID << "\n"; 
}

void Write_ASCII_Timestep(ofstream& fp, SETT& svar, State const& pnp1)
{
    fp <<  "ZONE T=\"" << "Ash particles" << "\"";
    fp <<", I=" << pnp1.size() << ", F=POINT" <<", STRANDID=1, SOLUTIONTIME=" << svar.t  << "\n";
    fp << std::left << std::scientific << std::setprecision(6);
    
    for (size_t ii = 0; ii < pnp1.size(); ++ii)
    {
        Write_ASCII_Point(fp, svar.scale, pnp1[ii], ii);
    }
    fp << std::flush;
}

void Write_ASCII_Streaks(ofstream& fp, SETT& svar, vector<State> const& t_pnp1)
{
    fp << "TITLE = \"Ash Particle Streaks\"\n";
	Write_ASCII_header(fp);
     
    size_t nTimes = t_pnp1.size();
    size_t nParts = t_pnp1[0].size();
    size_t time = 0;
    size_t part = 0;

    /* Outer loop is for each particle, writing the time data */
    for(part = 0; part < nParts; ++part)
    {
        size_t nGoing = 0;
        for (time = 0; time < nTimes; ++time)
        {   /* Find how many times to write */
            nGoing += t_pnp1[time][part].going;
        }

        fp <<  "ZONE T=\"" << "Ash particle " << part << "\"";
        fp <<", I=" << nGoing+1 << ", J=1, K=1, DATAPACKING=POINT"<< "\n";
        fp << std::left << std::scientific << std::setprecision(6);
        for (time = 0; time <= nGoing; ++time)
        { /* Inner loop to write the times of the particle */
            Write_ASCII_Point(fp, svar.scale, t_pnp1[time][part], part);    
        }
    
    }
    fp << std::flush;
}



#endif