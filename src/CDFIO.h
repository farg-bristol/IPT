#ifndef CDFIO_H
#define CDFIO_H

#include "Var.h"
#include <netcdf>
// using namespace netCDF;
// using namespace netCDF::exceptions;
#define NC_ERR 2
#define ERR(e)                                 \
	{                                          \
		printf("Error: %s\n", nc_strerror(e)); \
		exit(-1);                              \
	}

uint index(uint ii, uint jj, uint nPts)
{
	return(ii*nPts + jj);
}

void Average_Point_to_Cell(vector<vec<real,3>> const &pData, vector<vec<real,3>> &cData,
						   vector<vector<size_t>> const &elems)
{
	vector<vec<real,3>> sum(elems.size(), vec<real,3>(0.0));

	#pragma omp parallel for reduction(+ : sum)
	for (uint ii = 0; ii < elems.size(); ++ii)
	{

		uint const nVerts = elems[ii].size();
		for (auto jj : elems[ii])
		{
			sum[ii] += pData[jj];
		}
		sum[ii] /= nVerts;
	}

	cData = sum;
}

void Average_Point_to_Cell(vector<real> const &pData, vector<real> &cData,
						   vector<vector<size_t>> const &elems)
{
	vector<real> sum(elems.size(), 0.0);

	#pragma omp parallel for reduction(+ : sum)
	for (uint ii = 0; ii < elems.size(); ++ii)
	{

		uint const nVerts = elems[ii].size();
		for (auto jj : elems[ii])
		{
			sum[ii] += pData[jj];
		}
		sum[ii] /= nVerts;
	}

	cData = sum;
}


/*****************************************************************************/
/*************** READING NETCDF CELL BASED DATA FUNCTIONS ********************/
/*****************************************************************************/
void Read_BMAP(SETT& svar)
{
	#ifdef DEBUG
	dbout << "Entering Find_Bmap_Markers..." << endl;
	#endif
	string const& bmapIn = svar.taubmap; 
	std::ifstream fin(bmapIn, std::ios::in);

	if (!fin.is_open())
	{
		cout << "Couldn't open the boundary map file." << endl;
		cout << "Attempted path: " << bmapIn << endl;
		exit(-1);
	}

	/* Find out how many blocks, i.e. how many boundaries to expect. */
	uint nBlocks = 0;
	string line;
	while (getline(fin, line))
	{
		if (line.find("block end") != string::npos)
		{ /*We're on a new block*/
			nBlocks++;
		}
	}
	fin.clear();
	fin.seekg(0);
	
	uint blockno = 0;

	svar.markers = vector<int>(nBlocks,0);
	svar.bnames = vector<string>(nBlocks);
	svar.bwrite = vector<int>(nBlocks,0);

	vector<string> types(nBlocks);

	while (getline(fin, line))
	{
		/* Remove whitespace */
		line = ltrim(line);
		/* Check if line is commented out */
		if (line[0] == '#')
			continue;

		if (line.find("Markers") != string::npos)
		{
			string value = Get_Parameter_Value(line);
			std::istringstream sstr(value);
			sstr >> svar.markers[blockno];
		}

		if (line.find("Name") != string::npos)
		{ /* Get the name */
			svar.bnames[blockno] = Get_Parameter_Value(line);
		}

		if (line.find("Type") != string::npos)
		{ /* Get the name */
			types[blockno] = Get_Parameter_Value(line);
		}

		if (line.find("Write surface data (0/1)") != string::npos)
		{
			string value = Get_Parameter_Value(line);
			std::istringstream sstr(value);
			sstr >> svar.bwrite[blockno];
		}

		if (line.find("block end") != string::npos)
		{ /*We're on a new block*/
			blockno++;
		}
	}

	fin.close();

	/* Check if the boundary has a name. If not, then give it the type name*/
	for(size_t ii = 0; ii < nBlocks;  ++ii)
	{
		if(svar.bnames[ii].empty())
		{
			svar.bnames[ii] = types[ii];
		}
	}

	#ifdef DEBUG
	dbout << "Exiting Find_Bmap_Markers..." << endl;
	#endif
}


/*To run on the solution file*/
vector<real> Get_Scalar_Property_real(int &fin, string const &variable, size_t const &nPts)
{
	#ifdef DEBUG
		dbout << "Reading variable: " << variable << endl;
	#endif
	int retval;
	int varID;

	if ((retval = nc_inq_varid(fin, variable.c_str(), &varID)))
	{
		cout << "Failed to get variable id for: " << variable << endl;
		ERR(retval);
	}

	#ifdef DEBUG
		dbout << "Allocating array of: " << nPts << endl;
	#endif

	double *array = new double[nPts];

	#ifdef DEBUG
		dbout << "Attempting to read NetCDF variable:  " << variable << endl;
	#endif

	if ((retval = nc_get_var_double(fin, varID, &array[0])))
	{
		cout << "Failed to get variable id for: " << variable << endl;
		ERR(retval);
	}

	/*Convert it to a vector to store*/
	vector<double> propVec;
	propVec.insert(propVec.end(), &array[0], &array[nPts]);
	vector<real> var = propVec;

	#ifdef DEBUG
		dbout << "returning vector" << endl;
	#endif
	return var;
}

vector<int> Get_Scalar_Property_int(int &fin, string variable, int nPts)
{
	#ifdef DEBUG
		dbout << "Reading variable: " << variable << endl;
	#endif
	int retval;
	int varID;

	if ((retval = nc_inq_varid(fin, variable.c_str(), &varID)))
	{
		cout << "Failed to get variable id for: " << variable << endl;
		ERR(retval);
	}

	#ifdef DEBUG
		dbout << "Allocating array of: " << nPts << endl;
	#endif

	int *array = new int[nPts];

	#ifdef DEBUG
		dbout << "Attempting to read NetCDF variable:  " << variable << endl;
	#endif

	if ((retval = nc_get_var_int(fin, varID, &array[0])))
	{
		cout << "Failed to get variable data for: " << variable << endl;
		ERR(retval);
	}

	/*Convert it to a vector to store*/
	vector<int> propVec;
	propVec.insert(propVec.end(), &array[0], &array[nPts]);

	#ifdef DEBUG
		dbout << "returning vector" << endl;
	#endif
	return propVec;
}

/*To run on the mesh file*/
vector<vector<size_t>> Get_Element(int &fin, string const &variable, size_t const &nElem, size_t const &nPpEd)
{
	#ifdef DEBUG
		dbout << "Reading \"" << variable << "\"" << endl;
	#endif
	int retval;
	int varID;

	if ((retval = nc_inq_varid(fin, variable.c_str(), &varID)))
	{
		cout << "Failed to get the variable ID of \"" << variable << "\"" << endl;
		ERR(retval);
		exit(-1);
	}

	// cout << nElem << "  " << nPoints << endl;
	#ifdef DEBUG
		dbout << "Allocating array of: " << nElem << " by " << nPpEd << endl;
	#endif

	vector<vector<size_t>> elemVec(nElem, vector<size_t>(nPpEd));
	
	try
	{
		/*Allocate on the heap (can be big datasets)*/
		int *elemArray = new int [nElem*nPpEd];
	
		size_t start[] = {0,0};
		size_t end[] = {nElem,nPpEd};

		/*Get the actual data from the file*/
		if ((retval = nc_get_vara_int(fin, varID, start,end, &elemArray[0])))
		{
			cout << "Failed to get the variable data of \"" << variable << "\"" << endl;
			ERR(retval);
			exit(-1);
		}

		#ifdef DEBUG
			dbout << "Attempting to read NetCDF elements." << endl;
		#endif

		cout << "Successfully read \"" << variable << "\"" << endl;
		cout << "Number of faces: " << nElem << endl;

		#ifdef DEBUG
			dbout << "Successfully read \"" << variable << "\"" << endl;
		#endif

		/*Convert it to a vector to store*/
		for (size_t ii = 0; ii < nElem; ++ii)
			for (size_t jj = 0; jj < nPpEd; ++jj)
				elemVec[ii][jj] = static_cast<size_t>(elemArray[index(ii,jj,nPpEd)]);

	}
	catch (std::bad_alloc &ba)
	{

		std::cerr << "Bad alloc caught. Failed to allocate \"" << variable << "\"" << endl;
		exit(-1);
	}

	#ifdef DEBUG
		dbout << "Returning vector" << endl;
	#endif
	return elemVec;
}

/*To run on the mesh file*/
void Get_Coordinates(int &fin, size_t const &nPnts, vector<vec<real,3>> &coordVec)
{
	#ifdef DEBUG
		dbout << "Reading coordinates." << endl;
	#endif
		
	int retval;
	int xcID, ycID, zcID;

	if ((retval = nc_inq_varid(fin, "points_xc", &xcID)))
	{
		cout << "Failed to get the variable ID of: "
			<< "points_xc" << endl;
		cout << "Stopping" << endl;
		exit(-1);
	}

	if ((retval = nc_inq_varid(fin, "points_yc", &ycID)))
	{
		cout << "Failed to get the variable ID of: "
			<< "points_yc" << endl;
		cout << "Stopping" << endl;
		exit(-1);
	}

	if ((retval = nc_inq_varid(fin, "points_zc", &zcID)))
	{
		cout << "Failed to get the variable ID of: "
			<< "points_zc" << endl;
		cout << "Stopping" << endl;
		exit(-1);
	}


	#ifdef DEBUG
		dbout << "Number of points: " << nPnts << endl;
	#endif
	/*Allocate on the heap (can be big datasets)*/
	real *coordX = new real[nPnts];
	real *coordY = new real[nPnts];

	coordVec = vector<vec<real,3>>(nPnts);

	real *coordZ = new real[nPnts];
	if ((retval = nc_get_var_double(fin, xcID, &coordX[0])))
		ERR(retval);
	if ((retval = nc_get_var_double(fin, ycID, &coordY[0])))
		ERR(retval);
	if ((retval = nc_get_var_double(fin, zcID, &coordZ[0])))
		ERR(retval);

	for (uint ii = 0; ii < nPnts; ++ii)
	{ /*Convert it to a vector to store*/
		coordVec[ii] = vec<real,3>(coordX[ii], coordY[ii], coordZ[ii]);
	}

	#ifdef DEBUG
		dbout << "Returning coordinates." << endl;
	#endif
}

void Get_Cell_Faces(vector<vector<uint>> const &cell,
					vector<vector<uint>> const &facenum, std::vector<std::vector<std::vector<uint>>> &cFaces)
{
	for (uint ii = 0; ii < cell.size(); ++ii)
	{
		uint jj = 0;
		cFaces.emplace_back();
		for (auto faces : facenum)
		{
			cFaces[ii].emplace_back();
			for (auto vert : faces)
			{
				cFaces[ii][jj].emplace_back(cell[ii][vert]);
			}
			++jj;
		}
	}
}

/*****************************************************************************/
/***************** READING NETCDF SOLUTION DATA FUNCTIONS ********************/
/*****************************************************************************/

void Read_SOLUTION(string const &solIn, MESH &cells)
{
	cout << "Reading solution file..." << endl;

	int retval;
	int solID;

	if ((retval = nc_open(solIn.c_str(), NC_NOWRITE, &solID)))
	{
		cout << "A netCDF error occured whilst trying to open the solution file:" << endl;
		cout << "\t" << solIn << endl
			 << endl;
		ERR(retval);
		exit(-1);
	}

	cout << "Solution file opened, reading contents..." << endl;

	int ptDimID;
	size_t solPts;

	if ((retval = nc_inq_dimid(solID, "no_of_points", &ptDimID)))
	{
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_inq_dimlen(solID, ptDimID, &solPts)))
	{
		ERR(retval);
		exit(-1);
	}

	#ifdef DEBUG
		dbout << "Solution points: " << solPts << endl;
	#endif


	if (solPts != cells.nPnts)
	{
		cout << "Solution file does not have the same number of vertices as the mesh." << endl;
		cout << "Please check again." << endl;
		exit(-1);
	}


	vector<real> realDens = Get_Scalar_Property_real(solID, "density", solPts);

	/*Get the velocities*/
	vector<real> xvel, yvel, zvel;
	xvel = Get_Scalar_Property_real(solID, "x_velocity", solPts);
	yvel = Get_Scalar_Property_real(solID, "y_velocity", solPts);
	zvel = Get_Scalar_Property_real(solID, "z_velocity", solPts);

	vector<vec<real,3>> vel(solPts);
	/*Test for size*/
	if (xvel.size() == solPts)
	{ /*Turn the arrays into a state vector*/
		#pragma omp parallel for
		for (uint ii = 0; ii < solPts; ++ii)
		{
			vel[ii] = vec<real,3>(xvel[ii], yvel[ii], zvel[ii]);
		}
	}
	else
	{
		cout << "velocities do not have the same number of vertices as the mesh." << endl;
		cout << xvel.size() << "  " << solPts << endl;
		cout << "Please check again." << endl;
		exit(-1);
	}

	vector<real> dens(solPts);

	if (vel.size() == 0)
	{
		cout << "The solution data has not been correctly ingested. Please check the solution file." << endl;
	}

	/*Average the data to a the cell*/
	cells.SetCells();

	cout << "Averaging points to cell centres..." << endl;
	Average_Point_to_Cell(vel, cells.cVel, cells.elems);
	Average_Point_to_Cell(realDens, cells.cRho, cells.elems);

	/*Check the data*/
	// for(auto value:dens)
	// {
	// 	if(value == 0)
	// 	{
	// 		cout << "Cell Rho data has not been initialised" << endl;
	// 	}
	// }
}

/*****************************************************************************/
/*************** READING NETCDF EDGE BASED DATA FUNCTIONS ********************/
/*****************************************************************************/

#if SIMDIM == 2

void Place_Edges(int &fin, size_t const &nElem, size_t const &nPnts, size_t const &nEdge, MESH &cells)
{
#ifdef DEBUG
	dbout << "Reading element left/right and placing faces" << endl;
#endif

	vector<int> left = Get_Scalar_Property_int(fin,"left_element_of_edges",nEdge);
	vector<int> right = Get_Scalar_Property_int(fin,"right_element_of_edges",nEdge);

	vector<std::pair<int, int>> leftright(nEdge);

#pragma omp parallel shared(leftright)
	{
		/*Create local of */

#pragma omp for schedule(static) nowait
		for (size_t ii = 0; ii < nEdge; ++ii)
		{
			int lindex = left[ii];
			int rindex = right[ii];
			leftright[ii] = std::pair<int, int>(lindex, rindex);
#pragma omp critical
			{
				cells.cFaces[lindex].emplace_back(ii);
				if (rindex >= 0)
					cells.cFaces[rindex].emplace_back(ii);
			}
		}
	}

	cells.leftright = leftright;

#ifdef DEBUG
	dbout << "End of placing edges in elements." << endl;
#endif

	/*Now go through the faces and see which vertices are unique, to get element data*/
	cout << "Building elements..." << endl;
#pragma omp parallel
	{
// vector<vector<uint>> local;
#pragma omp for schedule(static) nowait
		for (uint ii = 0; ii < cells.cFaces.size(); ++ii)
		{
			for (auto const &index : cells.cFaces[ii])
			{
				vector<size_t> const face = cells.faces[index];
				for (auto const &vert : face)
				{
					if (std::find(cells.elems[ii].begin(), cells.elems[ii].end(), vert) == cells.elems[ii].end())
					{ /*Vertex doesn't exist in the elems vector yet.*/
#pragma omp critical
						{
							cells.elems[ii].emplace_back(vert);
						}
					}
				}
			}
		}

/*Find cell centres*/
#pragma omp single
		{
			cout << "Finding cell centres..." << endl;
		}

		Average_Point_to_Cell(cells.verts, cells.cCentre, cells.elems);

/*Find cell centres*/
#pragma omp single
		{
			cout << "Finding cell volumes..." << endl;
		}

// Find cell volumes
#pragma omp for schedule(static) nowait
		for (size_t ii = 0; ii < cells.elems.size(); ++ii)
		{
			cells.cVol[ii] = Cell_Volume(cells.verts, cells.faces, cells.elems[ii], cells.cFaces[ii], cells.cCentre[ii]);
		}
	}
	cout << "Elements built." << endl;

#ifdef DEBUG
	dbout << "All elements defined." << endl
		  << endl;
#endif
}

void Read_TAUMESH_EDGE(SIM &svar, MESH &cells, FLUID const &fvar, AERO const &avar)
{
	string meshIn = svar.infolder;
	string solIn = svar.infolder;
	meshIn.append(svar.meshfile);
	solIn.append(svar.solfile);

#ifdef DEBUG
	dbout << "Attempting read of NetCDF file." << endl;
	dbout << "Mesh file: " << meshIn << endl;
	dbout << "Solution file: " << solIn << endl;
#endif

	/*Read the mesh data*/
	int retval;
	int meshID;

	if ((retval = nc_open(meshIn.c_str(), NC_NOWRITE, &meshID)))
	{
		cout << "A netCDF error occured whilst trying to open the mesh file:" << endl;
		cout << "\t" << meshIn << endl
			 << endl;
		ERR(retval);
		exit(-1);
	}

	cout << "Mesh file open. Reading face data..." << endl;

	int ptDimID, elemDimID, edgeDimID, nPpEDimID;
	size_t nPnts, nElem, nEdge, nPpEd;

	

	// Retrieve how many elements there are.
	if ((retval = nc_inq_dimid(meshID, "no_of_elements", &elemDimID)))
	{
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_inq_dimlen(meshID, elemDimID, &nElem)))
	{
		ERR(retval);
		exit(-1);
	}

	// Retrieve edge number
	if ((retval = nc_inq_dimid(meshID, "no_of_edges", &edgeDimID)))
	{
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_inq_dimlen(meshID, edgeDimID, &nEdge)))
	{
		ERR(retval);
		exit(-1);
	}

	// Retrieve points per edge
	if ((retval = nc_inq_dimid(meshID, "points_per_edge", &nPpEDimID)))
	{
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_inq_dimlen(meshID, nPpEDimID, &nPpEd)))
	{
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_inq_dimid(meshID, "no_of_points", &ptDimID)))
	{
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_inq_dimlen(meshID, ptDimID, &nPnts)))
	{
		ERR(retval);
		exit(-1);
	}

#ifdef DEBUG
	dbout << "nElem : " << nElem << " nPnts: " << nPnts << " nEdge: " << nEdge << endl;
#endif

	cout << "nElem : " << nElem << " nPnts: " << nPnts << " nEdge: " << nEdge << endl;

	cells.numPoint = nPnts;
	cells.numElem = nElem;
	cells.numFace = nEdge;

	cells.elems = vector<vector<size_t>>(nElem);
	cells.cFaces = vector<vector<size_t>>(nElem);
	cells.verts = vector<vec<real,3>>(nPnts);
	cells.cVol = vector<real>(nElem);
	// cells.leftright = vector<std::pair<int,int>>(nFace);

	/*Get the faces of the mesh*/
	
	cells.faces = Get_Element(meshID, "points_of_element_edges", nEdge, nPpEd);

	/*Get the vertices that are in use to take the values from the solution file*/
	vector<int> usedVerts = Get_Scalar_Property_int(meshID, "vertices_in_use", nPnts);

	#ifdef DEBUG
		dbout << "Successfully read vertices_in_use data." << endl;
	#endif

	vector<uint> uVerts(nPnts);
	for (size_t ii = 0; ii < nPnts; ++ii)
		uVerts[ii] = static_cast<uint>(usedVerts[ii]);

	/*Get the coordinates of the mesh*/
	uint ignored = Get_Coordinates(meshID, nPnts, cells.verts);
	if (cells.verts.size() != nPnts)
	{
		cout << "Some data has been missed.\nPlease check how many points." << endl;
	}

	/*Adjust the scale*/
	cells.scale = svar.scale;
	for (auto &vert : cells.verts)
	{
		vert *= cells.scale;
	}
	svar.Start *= cells.scale;
	// svar.Jet/=cells.scale;

	/*Get face left and right, and put the faces in the elements*/
	Place_Edges(meshID, nElem, nPnts, nEdge, cells);

#ifdef DEBUG
	dbout << "End of interaction with mesh file and ingested data." << endl
		  << endl;
	dbout << "Opening solultion file." << endl;
#endif

	Read_SOLUTION(solIn, fvar, avar, ignored, cells, uVerts);

	for (size_t ii = 0; ii < cells.cMass.size(); ii++)
	{
		cells.cMass[ii] = cells.cRho[ii] * cells.cVol[ii];
	}
}

#endif

/*****************************************************************************/
/*************** READING NETCDF FACE BASED DATA FUNCTIONS ********************/
/*****************************************************************************/

void Place_Faces(int &fin, size_t const &nFace, MESH &cells)
{
	#ifdef DEBUG
		dbout << "Reading element left/right and placing faces" << endl;
	#endif

	vector<int> left = Get_Scalar_Property_int(fin,"left_element_of_faces",nFace);
	vector<int> right = Get_Scalar_Property_int(fin,"right_element_of_faces",nFace);

	vector<int> markers = Get_Scalar_Property_int(fin,"boundarymarker_of_surfaces",cells.nSurf);

	vector<std::pair<lint, lint>> leftright(nFace);
	cells.smarkers = vector<std::pair<size_t,int>>(cells.nSurf);
	
	vector<size_t> surf_IDs;
	#pragma omp parallel shared(leftright)
	{
		vector<size_t> local;
		#pragma omp for schedule(static) nowait
		for (size_t ii = 0; ii < nFace; ++ii)
		{
			int lindex = left[ii];
			int rindex = right[ii];
			leftright[ii] = std::pair<lint, lint>(lindex, rindex);
			#pragma omp critical
			{
				cells.cFaces[lindex].emplace_back(ii);
				if(rindex >= 0)
					cells.cFaces[rindex].emplace_back(ii);
				else
				{
					local.emplace_back(ii);
				}
			}
		}

		#pragma omp for schedule(static) ordered
		for(int ii = 0; ii < omp_get_num_threads(); ++ii)
		{
			#pragma omp ordered
			surf_IDs.insert(surf_IDs.end(),local.begin(),local.end());
		}

		#pragma omp single
		{
			if(surf_IDs.size() != cells.nSurf)
			{
				cout << "Mismatch of number of surface faces identified, and the number given" << endl;
				cout << "Identified: " << surf_IDs.size() << "  Given: " << cells.nSurf << endl;
				exit(-1);
			}
		}

		#pragma omp for schedule(static) nowait
		for(size_t ii = 0; ii < cells.nSurf; ++ii)
		{
			cells.smarkers[ii].first = surf_IDs[ii];
			cells.smarkers[ii].second = markers[ii];
		}
	}

	#ifdef DEBUG
		dbout << "End of placing faces in elements." << endl;
		dbout << "Building boundary indexes." << endl;
	#endif

	cells.leftright = leftright;

	/*Now go through the faces and see which vertices are unique, to get element data*/
	cout << "Building elements..." << endl;
	#pragma omp parallel
	{
		// Create list of unique vertex indices for each element
		#pragma omp for schedule(static) nowait
		for (size_t ii = 0; ii < cells.cFaces.size(); ++ii)
		{
			for (auto const &index : cells.cFaces[ii])
			{
				vector<size_t> const face = cells.faces[index];
				for (auto const& vert : face)
				{
					if (std::find(cells.elems[ii].begin(), cells.elems[ii].end(), vert) == cells.elems[ii].end())
					{ /*Vertex doesn't exist in the elems vector yet.*/
						#pragma omp critical
						{
							cells.elems[ii].emplace_back(vert);
						}
					}
				}
			}
		}

		/*Find cell centres*/
		#pragma omp single
		{
			cout << "Finding cell centres..." << endl;
		}

		Average_Point_to_Cell(cells.verts, cells.cCentre, cells.elems);

		/*Find cell centres*/
		// #pragma omp single
		// {
		// 	cout << "Finding cell volumes..." << endl;
		// }

		// Find cell volumes
		// #pragma omp for schedule(static) nowait
		// for (size_t ii = 0; ii < cells.elems.size(); ++ii)
		// {
		// 	cells.cVol[ii] = Cell_Volume(cells.verts, cells.faces, cells.elems[ii], cells.cFaces[ii], cells.cCentre[ii]);
		// }
	}

	cout << "Elements built." << endl;

	#ifdef DEBUG
		dbout << "All elements defined." << endl << endl;
	#endif
}

void Read_TAUMESH_FACE(SETT &svar, MESH &cells)
{
	string meshIn = svar.taumesh;
	string solIn = svar.tausol;

	#ifdef DEBUG
		dbout << "Attempting read of NetCDF file." << endl;
		dbout << "Mesh file: " << meshIn << endl;
		dbout << "Solution file: " << solIn << endl;
	#endif

	/*Read the mesh data*/
	int retval;
	int meshID;

	if ((retval = nc_open(meshIn.c_str(), NC_NOWRITE, &meshID)))
	{
		cout << "A netCDF error occured whilst trying to open the mesh file:" << endl;
		cout << "\t" << meshIn << endl
			 << endl;
		ERR(retval);
		exit(-1);
	}

	cout << "Mesh file open. Reading face data..." << endl;

	int hasTrig, hasQuad;
	int ptDimID, elemDimID, faceDimID, surfDimID, faceTDimID, faceQDimID, nPpTFDimID, nPpQFDimID;
	size_t nPnts, nElem, nFace, nSurf, nTFace, nQFace, nPpTFc, nPpQFc;

	
	// Retrieve how many elements there are.
	if ((retval = nc_inq_dimid(meshID, "no_of_elements", &elemDimID)))
	{
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_inq_dimlen(meshID, elemDimID, &nElem)))
	{
		ERR(retval);
		exit(-1);
	}

	// Retrieve face number
	if ((retval = nc_inq_dimid(meshID, "no_of_faces", &faceDimID)))
	{
		ERR(retval);
		exit(-1);
	}
	
	if ((retval = nc_inq_dimlen(meshID, faceDimID, &nFace)))
	{
		ERR(retval);
		exit(-1);
	}

	/* Retrieve how many points there are */
	if ((retval = nc_inq_dimid(meshID, "no_of_points", &ptDimID)))
	{
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_inq_dimlen(meshID, ptDimID, &nPnts)))
	{
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_inq_dimid(meshID, "no_of_surfaceelements", &surfDimID)))
	{
		ERR(retval);
		exit(-1);
	}

	if ((retval = nc_inq_dimlen(meshID, surfDimID, &nSurf)))
	{
		ERR(retval);
		exit(-1);
	}

	// Retrieve triangle face dimensions 
	if ((retval = nc_inq_dimid(meshID, "no_of_triangles", &faceTDimID)))
	{
		cout << "No triangle faces" << endl;
		hasTrig = 0;
	}
	else
	{
		if ((retval = nc_inq_dimlen(meshID, faceTDimID, &nTFace)))
		{
			ERR(retval);
			exit(-1);
		}	

		if ((retval = nc_inq_dimid(meshID, "points_per_triangle", &nPpTFDimID)))
		{
			ERR(retval);
			exit(-1);
		}
		
		if ((retval = nc_inq_dimlen(meshID, nPpTFDimID, &nPpTFc)))
		{
			ERR(retval);
			exit(-1);
		}

		hasTrig = 1;
	}

	// Retrieve quadrilateral face dimensions
	if ((retval = nc_inq_dimid(meshID, "no_of_quadrilaterals", &faceQDimID)))
	{
		cout << "No quadrilateral faces" << endl;
		hasQuad = 0;
	}
	else
	{
		if ((retval = nc_inq_dimlen(meshID, faceQDimID, &nQFace)))
		{
			ERR(retval);
			exit(-1);
		}

		if ((retval = nc_inq_dimid(meshID, "points_per_quadrilateral", &nPpQFDimID)))
		{
			ERR(retval);
			exit(-1);
		}
		
		if ((retval = nc_inq_dimlen(meshID, nPpQFDimID, &nPpQFc)))
		{
			ERR(retval);
			exit(-1);
		}

		hasQuad = 1;
	}

	

#ifdef DEBUG
	dbout << "nElem : " << nElem << " nPnts: " << nPnts << " nFace: " << nFace << endl;
#endif

	cout << "nElem : " << nElem << " nPnts: " << nPnts << " nFace: " << nFace << endl;

	cells.nPnts = nPnts;
	cells.nElem = nElem;
	cells.nFace = nFace;
	cells.nSurf = nSurf;

	cells.elems = vector<vector<size_t>>(nElem);
	cells.cFaces = vector<vector<size_t>>(nElem);
	// cells.cVol = vector<real>(nElem);

	cells.verts = vector<vec<real,3>>(nPnts);
	// cells.leftright = vector<std::pair<int,int>>(nFace);

	/*Get the faces of the mesh*/
	vector<vector<size_t>> trig, quad;
	if(hasTrig)
	{
		trig = vector<vector<size_t>>(nTFace,vector<size_t>(nPpTFc));
		trig = Get_Element(meshID, "points_of_triangles", nTFace, nPpTFc);
	}

	if(hasQuad)
	{
		quad = vector<vector<size_t>>(nQFace,vector<size_t>(nPpQFc));
		quad = Get_Element(meshID, "points_of_quadrilaterals", nQFace, nPpQFc);
	}

	cells.faces.insert(cells.faces.end(),trig.begin(),trig.end());
	cells.faces.insert(cells.faces.end(),quad.begin(),quad.end());

	if(cells.faces.size() != cells.nFace)
	{
		cout << "Mismatch of number of faces to that defined." << endl;
		cout << "number of faces: " << cells.nFace << " faces size: " << cells.faces.size() << endl;
	}


	/*Get the coordinates of the mesh*/
	Get_Coordinates(meshID, nPnts, cells.verts);
	if (cells.verts.size() != nPnts)
	{
		cout << "Some data has been missed.\nPlease check how many points." << endl;
	}

	/*Adjust the scale*/
	cells.scale = svar.scale;
	for (auto &vert : cells.verts)
	{
		vert *= cells.scale;
	}

	/*Get face left and right, and put the faces in the elements*/
	Place_Faces(meshID, nFace, cells);

	#ifdef DEBUG
		dbout << "End of interaction with mesh file and ingested data." << endl
			<< endl;
		dbout << "Opening solultion file." << endl;
	#endif
	Read_SOLUTION(solIn, cells);

	// for (size_t ii = 0; ii < cells.cMass.size(); ii++)
	// {
	// 	cells.cMass[ii] = cells.cRho[ii] * cells.cVol[ii];
	// }

	cout << "Cell processing complete!" << endl;
}

#endif