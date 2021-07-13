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
        string temp2;
        
        std::getline(iss,temp2,',');
        std::istringstream iss2(temp2);
        iss2 >> a;

        std::getline(iss,temp2,',');
        iss2 = std::istringstream(temp2);
        iss2 >> b;

        std::getline(iss,temp2,',');
        iss2 = std::istringstream(temp2);
        iss2 >> c;
        
        value = vec<real,3>(a,b,c);
    }
}

matrix<real,3> GetRotationMat(vec<real,3> const& angles)
{

    matrix<real,3> rotx, roty, rotz;
    rotx = matrix<real,3>(1.0, 0.0            , 0.0             ,
                          0.0, cos(angles[0]) , -sin(angles[0]) ,
                          0.0, sin(angles[0]) , cos(angles[0])  );

    roty = matrix<real,3>(cos(angles[1])  , 0.0 , sin(angles[1])  ,
                          0.0             , 1.0 , 0.0             ,
                          -sin(angles[1]) , 0.0 , cos(angles[1]) );

    rotz = matrix<real,3>(cos(angles[2]) , -sin(angles[2]) , 0.0 ,
                          sin(angles[2]) , cos(angles[2])  , 0.0 ,
                          0.0            , 0.0             , 1.0 );

    return rotz*roty*rotx;
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

        /* File Inputs */
        Get_String(line, "Primary grid face filename", svar.taumesh);
        Get_String(line, "Boundary mapping filename", svar.taubmap);
        Get_String(line, "Restart-data prefix", svar.tausol);
        Get_String(line, "OpenFOAM folder", svar.foamdir);
        Get_String(line, "OpenFOAM solution folder", svar.foamsol);
        Get_Number(line, "OpenFOAM binary (0/1)", svar.isBinary);
        Get_Number(line, "Label size (32/64)", svar.labelSize);
        Get_Number(line, "Scalar size (32/64)", svar.scalarSize);

        /* File outputs */
        Get_Number(line, "Particle scatter output (0/1)", svar.partout);
        Get_Number(line, "Particle streak output (0/1)", svar.streakout);
        Get_Number(line, "Particle cell intersection output (0/1)", svar.cellsout);
        Get_String(line, "Particle scatter output directory", svar.outdir);
        Get_String(line, "Particle streak output directory", svar.streakdir);
        Get_String(line, "Particle cell intersection output directory", svar.celldir);
        Get_String(line, "Particle surface impact filename", svar.surfacefile);

        /* Fluid data */
        Get_Number(line, "Grid scale", svar.scale);
        Get_Number(line, "Reference density", fvar.g_rho);
        Get_Number(line, "Reference dispersed density", fvar.d_rho);
        Get_Number(line, "Sutherland reference viscosity", fvar.mu_g);

        /* Simulation settings */
        Get_Number(line, "2D offset vector (0 / x=1,y=2,z=3)",svar.offset_axis);
        Get_String(line, "Trajectory ODE integrator", svar.integration_type);
        Get_Number(line, "Particle timestep", svar.dt);
        Get_Number(line, "Particle frame time", svar.framet);
        Get_Number(line, "Newmark-Beta iteration limit", svar.max_iter);
        Get_Number(line, "Particle maximum frame number", svar.max_frame);
        Get_Number(line, "Velocity equation order (1/2)", svar.eqOrder);
        Get_Number(line, "Maximum x trajectory coordinate", svar.max_x);
        Get_Number(line, "Catch efficiency delta (micron)", svar.delta);
        Get_Number(line, "Catch efficiency delta factor", svar.delta_fac);
        Get_Number(line, "Catch efficiency max delta (micron)", svar.max_delta);        

        /* Starting area conditions */
        Get_String(line, "Start grid type", svar.startName);

        /* Disk */
        Get_Vector(line, "Seed disk center", svar.discPos);
        Get_Vector(line, "Seed disk angles", svar.discAngles);
        Get_Number(line, "Seed disk diameter", svar.discDiam);
        Get_Number(line, "Seed disk spacing", svar.discSpace);

        /* Grid */
        Get_Vector(line, "First corner (x,y,z)", svar.grid_verts[0]);
        Get_Vector(line, "Second corner (x,y,z)", svar.grid_verts[1]);
        Get_Vector(line, "Third corner (x,y,z)", svar.grid_verts[2]);
        Get_Vector(line, "Fourth corner (x,y,z)", svar.grid_verts[3]);
        Get_Number(line, "Points in i-direction", svar.n_i);
        Get_Number(line, "Points in j-direction", svar.n_j);
        
        /* Line */
        Get_Vector(line, "First point (x,y,z)", svar.line_verts[0]);
        Get_Vector(line, "Second point (x,y,z)", svar.line_verts[1]);
        Get_Number(line, "Resolution of line", svar.n_i);

        /* Particle diameter */
        Get_Number(line, "Particle diameter sweep (0/1)", svar.particle_sweep);
        Get_Number(line, "Particle diameter lower bound",fvar.d_lower);
        Get_Number(line, "Particle diameter upper bound",fvar.d_upper);
        Get_Number(line, "Particle diameter interval",fvar.d_interval);
        Get_Number(line, "Particle diameter (micron)", fvar.d_0);
        Get_Number(line, "Particle sweep resolution",svar.nParticle_diams);
        Get_String(line, "Particle sweep distribution (linear/log)",svar.distr_type);

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

    if(svar.partout == 1)
    {
        if(svar.outdir.empty())
        {
            cout << "Output particle directory not defined." << endl;
            exit(-1);
        }
    }

    if(svar.streakout == 1)
    {
        if(svar.streakdir.empty())
        {
            cout << "Output particle streaks directory not defined." << endl;
            exit(-1);
        }
    }

    if(svar.cellsout == 1)
    {
        if(svar.celldir.empty())
        {
            cout << "Output cell intersections directory not defined." << endl;
            exit(-1);
        }
    }

    if(svar.eqOrder > 2 || svar.eqOrder < 1)
    {
        cout << "Equation order not 1 or 2. Please choose between these." << endl;
        exit(-1);
    }

    // svar.streakdir = svar.outfile;
    // size_t pos = svar.streakfile.find_last_of(".");
    // svar.streakfile.insert(pos,"_streak");
    
    if(svar.offset_axis != 0)
    {
        svar.dimension = 2;
    }
    else if (svar.offset_axis > 3)
    {
        cout << "2D offset axis option out of bounds" << endl;
        exit(-1);
    }

    /* Find a timestep that gets a perfect frame interval */
    if(svar.integration_type.empty())
    {
        cout << "Integration method not defined" << endl;
        exit(-1);
    }
    else if(svar.integration_type == "Newmark-Beta")
    {
        svar.explicit_or_implicit = 0;
        /* Newmark-Beta only for now */
        if (svar.dt == 0.0)
        {
            cout << "timestep has not been initialised." << endl;
            exit(-1);
        }

        if(svar.framet == 0.0)
        {
            cout << "frame timestep has not been initialised." << endl;
            exit(-1);
        }

        int nsteps = ceil(svar.framet/svar.dt);
        svar.nsteps = nsteps;
        svar.dt = svar.framet / real(nsteps);
    }
    else if (svar.integration_type == "Implicit")
    {
        svar.explicit_or_implicit = 1;
    }
    
    if(svar.delta > 0)
    {
        svar.fac_or_fixed_delta = 1;
        svar.delta /=1e6;
    }
    svar.max_delta/=1e6;


    if(svar.startName == "Disk")
    { 
        svar.startType = 2;

        if(svar.discPos[0] == -500000)
        {
            cout << "Seed disk position has not been initialised." << endl;
            exit(-1);
        }

        if(svar.discDiam == 0.0)
        {
            cout << "Seed disk diameter has not been initialised." << endl;
            exit(-1);
        }

        if(svar.discSpace == 0.0)
        {
            cout << "Seed disk spacing has not been initialised." << endl;
            exit(-1);
        }
        svar.discAngles *= M_PI/180.0;
        svar.rotate = GetRotationMat(svar.discAngles);
        svar.transp = svar.rotate.transpose();
    }
    else if(svar.startName == "Grid")
    {
        svar.startType = 0;

        if(svar.n_i == -1)
        {
            cout << "Points in i-direction has not been initialised." << endl;
            exit(-1);
        }

        if(svar.n_j == -1)
        {
            cout << "Points in j-direction has not been initialised." << endl;
            exit(-1);
        }

        /* Check all corners are defined */
        if(svar.grid_verts[0][0] == -500000 || 
            svar.grid_verts[0][1] == -500000 ||
            svar.grid_verts[0][2] == -500000)
        {
            cout << "Some or all of the first corner is not initialised." << endl;
            exit(-1);
        }

        if(svar.grid_verts[1][0] == -500000 || 
            svar.grid_verts[1][1] == -500000 ||
            svar.grid_verts[1][2] == -500000)
        {
            cout << "Some or all of the second corner is not initialised." << endl;
            exit(-1);
        }

        if(svar.grid_verts[2][0] == -500000 || 
            svar.grid_verts[2][1] == -500000 ||
            svar.grid_verts[2][2] == -500000)
        {
            cout << "Some or all of the third corner is not initialised." << endl;
            exit(-1);
        }

        if(svar.grid_verts[3][0] == -500000 || 
            svar.grid_verts[3][1] == -500000 ||
            svar.grid_verts[3][2] == -500000)
        {
            cout << "Some or all of the fourth corner is not initialised." << endl;
            exit(-1);
        }

    }
    else if (svar.startName == "Line")
    {
        svar.startType = 1;

        if(svar.n_i == -1)
        {
            cout << "Resolution of line has not been initialised." << endl;
            exit(-1);
        }

        /* Check both points of the line are defined */
        if(svar.line_verts[0][0] == -500000 || 
            svar.line_verts[0][1] == -500000 ||
            svar.line_verts[0][2] == -500000)
        {
            cout << "Some or all of the first point is not initialised." << endl;
            exit(-1);
        }

        if(svar.line_verts[1][0] == -500000 || 
            svar.line_verts[1][1] == -500000 ||
            svar.line_verts[1][2] == -500000)
        {
            cout << "Some or all of the second point is not initialised." << endl;
            exit(-1);
        }
    }
    else
    {
        cout << "Starting grid type has not been initialised with a correct value." << endl;
        exit(-1);
    }

    if(svar.particle_sweep == 1)
    {
        if(fvar.d_lower == -1)
        {
            cout << "Lower particle boundary for sweep has not been defined." << endl;
            exit(-1);
        }

        if(fvar.d_upper == -1)
        {
            cout << "Upper particle boundary for sweep has not been defined." << endl;
            exit(-1);
        }


        if(svar.nParticle_diams != 0)
        {
            if(svar.distr_type == "log")
            {
                svar.particle_distr = 1;
                real log_interval = (log10(fvar.d_upper) - log10(fvar.d_lower))/real(svar.nParticle_diams-1);
                fvar.diams = vector<real>(svar.nParticle_diams,0.0);
                fvar.diams[0] = 1e-6*fvar.d_lower;
                for(int ii = 1; ii < svar.nParticle_diams; ++ii)
                {
                    real log_diam = log10(fvar.d_lower) + ii*log_interval;
                    fvar.diams[ii] = 1e-6 * (pow(10,log_diam));
                }

                fvar.d_interval = log_interval; 
            }
            else if (svar.distr_type == "linear")
            {
                svar.particle_distr = 0;
                fvar.d_interval = (fvar.d_upper - fvar.d_lower)/real(svar.nParticle_diams-1); 
                fvar.diams = vector<real>(svar.nParticle_diams,0.0);
                for(int ii = 0; ii < svar.nParticle_diams; ++ii)
                {
                    fvar.diams[ii] = 1e-6 * (fvar.d_lower + real(ii)*fvar.d_interval);
                }
            }
            
        }
        else if (fvar.d_interval != -1)
        {
            if(svar.distr_type == "log")
            {
                cout << "A logarithmic distribution is only available with a resolution definition" << endl;
                exit(-1);
            }

            svar.nParticle_diams = ceil((fvar.d_upper - fvar.d_lower)/fvar.d_interval)+1;

            fvar.d_interval = (fvar.d_upper - fvar.d_lower)/real(svar.nParticle_diams-1);
            fvar.diams = vector<real>(svar.nParticle_diams,0.0);
            for(int ii = 0; ii < svar.nParticle_diams; ++ii)
            {
                fvar.diams[ii] = 1e-6 * (fvar.d_lower + real(ii)*fvar.d_interval);
            }
        }
        else
        {
            cout << "Particle interval metric mot defined. Either specify resolution or spacing" << endl;
            exit(-1);
        }
    }
 
    fvar.d_0 /= 1e6;
    fvar.d_lower /= 1e6;
    fvar.d_upper /= 1e6;
    fvar.d_interval /= 1e6;

    // if(svar.foamdir.empty())
    // {
    //     cout << "Output OpenFOAM directory not defined." << endl;
    //     exit(-1);
    // }

} 

void Write_ASCII_variables(ofstream& fp)
{
	string variables;

	variables = 
"\"X\", \"Y\", \"Z\", \"t\", \"dt\", \"v\", \"a\", \"ptID\", \"Cell_V\", \"Cell_Rho\", \"Cell_ID\"";

	fp << "VARIABLES = " << variables << "\n";
}

void Write_ASCII_Scatter_Header(ofstream& fp, size_t const& pID)
{
    fp <<  "TITLE =\"" << "Ash particle " << pID << "\" \n";
    fp << "F=POINT\n";
    fp << std::left << std::scientific << std::setprecision(6);
}

void Write_ASCII_Point(ofstream& fp, real const& scale, part const& pnp1)
{
    const static uint width = 15;

    for(uint dim = 0; dim < 3; ++dim)
        fp << setw(width) << pnp1.xi[dim]/scale;

    fp << setw(width) << pnp1.t; 

    fp << setw(width) << pnp1.dt; 
    
    fp << setw(width) << pnp1.v.norm(); 

    fp << setw(width) << pnp1.acc.norm(); 

    fp << setw(width) << pnp1.partID;

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
        Write_ASCII_Point(fp, svar.scale, pnp1[ii]);
    }
    fp << std::flush;
}

void Write_ASCII_Streaks(SETT const& svar, vector<part> const& t_pnp1, part const& pnp1)
{
    /* Write file for each particle. Gonna get way too confusing otherwise */

    string file = svar.streakdir;
    file.append("/particle_");
    file.append(std::to_string(pnp1.partID));
    file.append("_streaks.dat");
    ofstream fp(file);

    if(!fp.is_open())
    {
        cout << "Couldn't open the streak output file." << endl;
        exit(-1);
    }

    fp << "TITLE = \"Ash Particle " << pnp1.partID << " Streaks\"\n";
	Write_ASCII_variables(fp);
     
    size_t nTimes = t_pnp1.size();
    size_t time = 0;
    
    

    fp <<  "ZONE T=\"" << "Ash particle " << pnp1.partID << "\"";
    if(pnp1.failed == 0)
        fp <<", I=" << nTimes+1 << ", J=1, K=1, DATAPACKING=POINT"<< "\n";
    else
        fp <<", I=" << nTimes << ", J=1, K=1, DATAPACKING=POINT"<< "\n";

    fp << std::left << std::scientific << std::setprecision(6);
    for (time = 0; time < nTimes; ++time)
    { /* Inner loop to write the times of the particle */
        Write_ASCII_Point(fp, svar.scale, t_pnp1[time]);
    }

    if(pnp1.failed == 0)
        Write_ASCII_Point(fp, svar.scale, pnp1);

    fp.close();

}

void Write_ASCII_Cells(SETT const& svar, MESH const& cells, vector<part> const& t_pnp1)
{
    size_t nTimes = t_pnp1.size();
    size_t time = 0;

    /* Write file for each particle. Gonna get way too confusing otherwise */
    string file = svar.celldir;
    file.append("/particle_");
    file.append(std::to_string(t_pnp1[0].partID));
    file.append("_cells.dat");
    ofstream fout(file);

    if(!fout.is_open())
    {
        cout << "Couldn't open the cell intersection output file." << endl;
        exit(-1);
    }

    fout << "TITLE = \"Ash particle " << t_pnp1[0].partID << " intersecting cells\"\n";

    int TotalNumFaceNodes = 0;

    vector<size_t> vertIndexes;
    vector<size_t> faceIndexes;
    vector<lint> cellIndexes;

    vector<vec<real,3>> usedVerts;
    vector<vector<size_t>> faces;
    vector<lint> left;
    vector<lint> right;


    /* Get cell properties for the intersected cells (delete duplicates later)*/
    for (time = 0; time < nTimes; ++time)
    {
        lint cellID = t_pnp1[time].cellID;
        
        /* Find how many vertices the cell has */
        vertIndexes.insert(vertIndexes.end(), cells.elems[cellID].begin(), cells.elems[cellID].end());
        
        faceIndexes.insert(faceIndexes.end(), cells.cFaces[cellID].begin(), cells.cFaces[cellID].end());

        cellIndexes.emplace_back(cellID);
    }

    /* Delete repeats of vertex mentions*/
    std::sort(vertIndexes.begin(),vertIndexes.end());
    vertIndexes.erase( std::unique( vertIndexes.begin(), vertIndexes.end()), vertIndexes.end() );

    std::sort(faceIndexes.begin(),faceIndexes.end());
    faceIndexes.erase( std::unique( faceIndexes.begin(), faceIndexes.end()), faceIndexes.end() );

    for(size_t const& vert : vertIndexes)
    {
        usedVerts.emplace_back(cells.verts[vert]);
    }

    /* Get faces properties used in the cell */
    for(size_t const& face: faceIndexes )
    {
        faces.emplace_back(cells.faces[face]);
        left.emplace_back(cells.leftright[face].first);
        right.emplace_back(cells.leftright[face].second);
    }

    /* Go through the faces indexes, finding the appropriate index to change it to */

    for(vector<size_t>& face : faces)
    {
        for(size_t& vert : face)
        {
            vector<size_t>::iterator index = std::find(vertIndexes.begin(), vertIndexes.end(), vert);
            if(index != vertIndexes.end())
            {
                vert = index - vertIndexes.begin() + 1 ;
            }
            else
            {
                cout << "Couldn't find the vertex used in the prepared list." << endl;
                exit(-1);
            }                
        }
        TotalNumFaceNodes += face.size();
    }

    /* Find the cell used and change it's index */
    for(lint& leftCell:left)
    {
        vector<lint>::iterator index = std::find(cellIndexes.begin(), cellIndexes.end(), leftCell);
        if(index != cellIndexes.end())
        {
            leftCell = (index - cellIndexes.begin())+1;
        }
        else
        {   /* Cell isn't intersected, so it won't be drawn. Boundary face. */
            leftCell = 0;
        } 
    }

    for(lint& rightCell:right)
    {
        vector<lint>::iterator index = std::find(cellIndexes.begin(), cellIndexes.end(), rightCell);
        if(index != cellIndexes.end())
        {
            rightCell = index - cellIndexes.begin() + 1;
        }
        else
        {   /* Cell isn't intersected, so it won't be drawn. Boundary face. */
            rightCell = 0;
        }
    }

    
    fout << "VARIABLES= \"X\" \"Y\" \"Z\" " << endl;
    fout << "ZONE T=\"particle " << t_pnp1[0].partID << " intersecting cells\"" << endl;
    fout << "ZONETYPE=FEPOLYHEDRON" << endl;
    fout << "NODES=" << usedVerts.size() << " ELEMENTS=" << cellIndexes.size() << 
            " FACES=" << faces.size() << endl;
    fout << "TotalNumFaceNodes=" << TotalNumFaceNodes << endl;
    fout << "NumConnectedBoundaryFaces=0 TotalNumBoundaryConnections=0" << endl;

    size_t w = 15;
    size_t preci = 6;
    fout << std::left << std::scientific << std::setprecision(preci);

    size_t newl = 0;
    fout << std::setw(1);
    for(size_t DIM = 0; DIM < 3; ++DIM)
    {
        for(size_t ii = 0; ii < usedVerts.size(); ++ii)
        {
            fout << std::setw(w) << usedVerts[ii][DIM];
            newl++;

            if(newl>4)
            {
                fout << endl;
                fout << " ";
                newl=0;
            }
        }
    }
    fout << endl;

    fout << std::left << std::fixed;
    w = 9;
    /*Inform of how many vertices in each face*/
    fout << "#node count per face" << endl;
    newl = 0;
    for (size_t ii = 0; ii < faces.size(); ++ii)
    {
        fout << std::setw(w) << faces[ii].size();
        newl++;

        if(newl>4)
        {
            fout << endl;
            newl=0;
        }
    }
    fout << endl;
    /*Write the face data*/
    fout << "#face nodes" << endl;
    for (size_t ii = 0; ii < faces.size(); ++ii)
    {
        for(auto const& vertex:faces[ii])
        {	/*Write face vertex indexes*/
            fout << std::setw(w) << vertex;
            // if (vertex > fdata.nPnts)
            // {
            // 	cout << "Trying to write a vertex outside of the number of points." << endl;
            // }
        }
        fout << endl;
    }

    /*Write face left and right*/
    newl = 0;
    fout << "#left elements" << endl;
    for (size_t ii = 0; ii < left.size(); ++ii)
    {
        fout << std::setw(w) << left[ii] ;
        newl++;

        if(newl>4)
        {
            fout << endl;
            newl=0;
        }
    }
    fout << endl;

    newl = 0;
    fout << "#right elements" << endl;
    for (size_t ii = 0; ii < right.size(); ++ii)
    {
        fout << std::setw(w) << right[ii] ;
        newl++;

        if(newl>4)
        {
            fout << endl;
            newl=0;
        }
    }

    fout.close();
}


void Write_ASCII_Impacts(SETT const& svar, FLUID const& fvar, MESH const& cells, 
            vector<SURF> const& surfs, vector<vector<real>> const& beta_data)
{
    // uint nSurf = svar.markers.size();
    uint nSurf = 0;
    vector<int> marks;
    vector<string> names;
    vector<SURF> surfaces_to_write;
    for(size_t ii = 0; ii < surfs.size(); ii++)
    {
        if(surfs[ii].output == 1)
        {
            marks.emplace_back(svar.markers[ii]);
            names.emplace_back(svar.bnames[ii]);
            surfaces_to_write.emplace_back(surfs[ii]);
        }
        nSurf += surfs[ii].output;
    }

    nSurf = 0;
    for(size_t ii = 0; ii < surfs.size(); ii++)
    {
        if(surfs[ii].output == 1)
        {
            surfaces_to_write[nSurf].face_count = surfs[ii].face_count;
            surfaces_to_write[nSurf].face_beta = surfs[ii].face_beta;
            surfaces_to_write[nSurf].face_area = surfs[ii].face_area;
        }
        nSurf += surfs[ii].output;
    }
    

    vector<vector<vector<size_t>>> faces(nSurf);
    vector<std::pair<size_t,int>> smarkers = cells.smarkers;

    vector<vector<size_t>> vertIndexes(nSurf);
    vector<vector<vec<real,3>>> usedVerts(nSurf);
    
    /* Do I need to sort the markers? */
    std::sort(smarkers.begin(), smarkers.end(),
    [](std::pair<size_t,int> const& p1, std::pair<size_t,int> const& p2){return p1.second > p2.second;});

    
    for(size_t ii = 0; ii < surfaces_to_write.size(); ++ii)
    {

        for(size_t jj = 0; jj < surfaces_to_write[ii].faceIDs.size(); ++jj )
        {
            size_t faceID = surfaces_to_write[ii].faceIDs[jj];
            faces[ii].emplace_back(cells.faces[faceID]);
        }
    
        /* Get the vertex indexes, delete duplicates, reorder, and recast the face indexes */
        for (vector<size_t> const& face : faces[ii])
            vertIndexes[ii].insert(vertIndexes[ii].end(), face.begin(), face.end());
        

        std::sort(vertIndexes[ii].begin(),vertIndexes[ii].end());
        vertIndexes[ii].erase(std::unique( vertIndexes[ii].begin(), vertIndexes[ii].end()), vertIndexes[ii].end() );

        for(size_t const& vert : vertIndexes[ii])
            usedVerts[ii].emplace_back(cells.verts[vert]);
    

        for(vector<size_t>& face : faces[ii])
        {
            for(size_t& vert : face)
            {
                vector<size_t>::iterator index = std::find(vertIndexes[ii].begin(), vertIndexes[ii].end(), vert);
                if(index != vertIndexes[ii].end())
                {
                    vert = index - vertIndexes[ii].begin() + 1 ;
                }
                else
                {
                    cout << "Couldn't find the vertex used in the prepared list." << endl;
                    exit(-1);
                }                
            }
        }
    }

    
    string file = svar.surfacefile;
    file.append("_");

    std::stringstream stream;
    stream << std::fixed << std::setprecision(2) << fvar.d_0*1e6 ;
    file.append(stream.str());
    file.append("um.dat");
    ofstream fout(file);

    if(!fout.is_open())
    {
        cout << "Couldn't open the cell intersection output file." << endl;
        exit(-1);
    }

    fout << "TITLE=\"Surface collision metrics\"\n";
    fout << "VARIABLES= \"X\" \"Y\" \"Z\" \"Number of Impacts\" \"beta\" \"average area\"\n";
    
    for(size_t ii = 0; ii < faces.size(); ++ii)
    {
        fout << "ZONE T=\"" << names[ii] << "\"\n";
        fout << "N=" << usedVerts[ii].size() << ", E=" << faces[ii].size();
        fout << ", F=FEBLOCK ET=QUADRILATERAL, VARLOCATION=([1-3]=NODAL,[4-7]=CELLCENTERED)\n\n";
        
        
        size_t w = 15;
        size_t preci = 6;
        fout << std::left << std::scientific << std::setprecision(preci);

        size_t newl = 0;
        fout << std::setw(1);
        for(size_t DIM = 0; DIM < 3; ++DIM)
        {
            for(size_t jj = 0; jj < usedVerts[ii].size(); ++jj)
            {
                fout << std::setw(w) << usedVerts[ii][jj][DIM];
                newl++;

                if(newl>4)
                {
                    fout << endl;
                    fout << " ";
                    newl=0;
                }
            }
        }
        fout << endl << endl;

        /* Variable data goes here */
        fout << "#face impact count" << endl;
        for(size_t jj = 0; jj < surfaces_to_write[ii].face_count.size(); ++jj)
        {
            fout << std::setw(w) << surfaces_to_write[ii].face_count[jj];
            newl++;

            if(newl>4)
            {
                fout << endl;
                fout << " ";
                newl=0;
            }
        }
        fout << endl << endl;

        // for(size_t jj = 0; jj < surfs[surf].mass.size(); ++jj)
        // {
        //     fout << std::setw(w) << surfs[surf].mass[jj];
        //     newl++;

        //     if(newl>4)
        //     {
        //         fout << endl;
        //         fout << " ";
        //         newl=0;
        //     }
        // }
        // fout << endl;
        fout << "#face beta value" << endl; 
        for(size_t jj = 0; jj < surfaces_to_write[ii].face_beta.size(); ++jj)
        {
            fout << std::setw(w) << surfaces_to_write[ii].face_beta[jj];
            newl++;

            if(newl>4)
            {
                fout << endl;
                fout << " ";
                newl=0;
            }
        }
        fout << endl << endl;

        fout << "#face area value" << endl;
        for(size_t jj = 0; jj < surfaces_to_write[ii].face_area.size(); ++jj)
        {
            fout << std::setw(w) << surfaces_to_write[ii].face_area[jj];
            newl++;

            if(newl>4)
            {
                fout << endl;
                fout << " ";
                newl=0;
            }
        }
        fout << endl << endl;


        /*Write the face data*/
        fout << "#face nodes" << endl;
        for (size_t jj = 0; jj < faces[ii].size(); ++jj)
        {
            for(auto const& vertex:faces[ii][jj])
            {	/*Write face vertex indexes*/
                fout << std::setw(w) << vertex;
                // if (vertex > fdata.nPnts)
                // {
                // 	cout << "Trying to write a vertex outside of the number of points." << endl;
                // }
            }

            if(faces[ii][jj].size() == 3)
                fout << std::setw(w) << faces[ii][jj].back();

            fout << endl;
        }
    }
    fout.close();
}   

/* Write catch efficiency as a plot of efficiency vs arclength chord fraction. */
void Write_Aerofoil_Catch_Efficiency(SETT const& svar, FLUID const& fvar, MESH const& cells, vector<int> const& marks, 
            vector<vector<real>> const& beta, vector<vector<real>> const& area, vector<vector<vec<real,3>>> const& pos )
{
    /* Open file */
    string file = svar.surfacefile;
    file.append("_");
    std::stringstream stream;
    stream << std::fixed << std::setprecision(2) << fvar.d_0*1e6 ;
    file.append(stream.str());
    file.append("um_Arc.dat");

    ofstream fout(file);

    if(!fout.is_open())
    {
        cout << "Couldn't open the catch efficiency output file." << endl;
        exit(-1);
    }

    fout << "TITLE=\"Catch efficiency vs arclength chord fraction\"\n";
    fout << "VARIABLES= \"Arclength chord fraction\" \"beta\" \"average area\"\n";


    for(size_t jj = 0; jj < marks.size(); ++jj)
    {
        string name = svar.bnames[marks[jj]];

        /* Find the face centers, and their distance along the surface from the leading edge  */
        vector<real> fArc(pos[jj].size(),0.0);
        
        // real smallest = 10;
        real dx = 1e-5;
        for(size_t ii = 0; ii < pos[jj].size(); ++ii)
        {
            /* Ignore points where the calculation failed for whatever reason */
            if(beta[jj][ii] == 0)
                continue;

            /* Get the distance */
            real length = 0.0;

            /* Calculate the camber line */
            // real xi = pos[ii][0];
            real yi = pos[jj][ii][2];

            real r = 0.2025;
            real k1 = 15.957;
            real a0 = 0.2969;
            real a1 = -0.126;
            real a2 = -0.3516;
            real a3 = 0.2843;
            real a4 = -0.1036;
            
            // real dx = xi/100.0;

            vec<real,2> posm1(0.0);
            vec<real,2> posp1(0.0);
            real x = 0.0;
            while(fabs(posp1[1]) < fabs(yi))
            {   
                /* Get the y of the camber line */         
                real y_camber = 0.0;
                real dy_camber = 0.0;

                if(x < r)
                {
                    y_camber = k1/6.0 * (pow(x,3)-3*r*x*x + r*r*(3-r)*x);
                    dy_camber = k1/6.0 * (3*pow(x,2)-6*r*x + r*r*(3-r));
                }
                else
                {
                    y_camber = k1*r*r*r/6 * (1-x);
                    dy_camber = -k1*r*r*r/6;
                }

                /* Get the thickness */
                real y_thick = 0.12/0.2 * (a0*pow(x,0.5) + a1*x + a2*x*x + a3*x*x*x + a4*x*x*x*x);
                
                real theta = atan(dy_camber);

                if(yi > 0)
                {
                    posp1 = vec<real,2>(x-y_thick*sin(theta), y_camber + y_thick*cos(theta));
                    length += (posp1-posm1).norm();   
                }
                else
                {
                    posp1 = vec<real,2>(x+y_thick*sin(theta), y_camber - y_thick*cos(theta));   
                    length -= (posp1-posm1).norm();
                }

                /* Get the distance */
                
                posm1 = posp1;
                x += dx;
            }

            fArc[ii] = (length);
            // if(fabs(arc) < smallest)
            //     smallest = arc;
        }

        vector<vector<real>> data;

        for(size_t ii = 0; ii < beta[jj].size(); ++ii)
        {
            if(beta[jj][ii] == 0)
                continue;

            vector<real> temp = {fArc[ii],beta[jj][ii],area[jj][ii]};
            data.emplace_back(temp);
        } 

        /* Sort by the first index */
        if(!data.empty())
        {
            std::sort(data.begin(),data.end(),[](vector<real> const& p1, vector<real> const& p2){return p1[0] < p2[0];});

            /* Add a 0 value at the front and end. */
            vector<real> temp(3,0.0);
            temp[0] = data[0][0];

            data.insert(data.begin(),temp);
            temp[0] = data.back()[0];
            data.emplace_back(temp);

            fout << "ZONE T=\"" << name << "\"\n";

            size_t w = 15;
            size_t preci = 6;
            fout << std::left << std::scientific << std::setprecision(preci);
            fout << std::setw(1);
            
            for(size_t kk = 0; kk < data.size(); ++kk)
            {
                fout << std::setw(w) << data[kk][0] << std::setw(w) << 
                data[kk][1] << std::setw(w) << data[kk][2] << endl;
            }
        }
    }
    
    fout.close();

}



#endif