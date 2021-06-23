#ifndef CONTAINMENT_H
#define CONTAINMENT_H

#include "Var.h"

#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif

#define PERTURB(i,j) pow(MEPSILON,pow(2,real(i)*3-real(j)))

void FindCellNeighbours(Vec_Tree const& CELL_INDEX, size_t const& num_results, vec<real,3> const& testp, vector<size_t>& outlist)
{
    outlist = vector<size_t>(num_results);
    vector<real> out_dists_sqr(num_results);

    nanoflann::KNNResultSet<real> resultSet(num_results);
    resultSet.init(&outlist[0], &out_dists_sqr[0]);
    
    CELL_INDEX.index->findNeighbors(resultSet, &testp[0], nanoflann::SearchParams(10));

}


/*Crossing test for 3 dimensions.*/
int LessThanREError(matrix<real,4> const& A)
{
    real a1, a2, a3;

    /*Calculate components of the absolute*/
    a1 = fabs(A(0,2)-A(3,2))*(fabs((A(1,0)-A(3,0))*(A(2,1)-A(3,1)))+fabs((A(1,1)-A(3,1))*(A(2,0)-A(3,0))));
    a2 = fabs(A(1,2)-A(3,2))*(fabs((A(2,0)-A(3,0))*(A(0,1)-A(3,1)))+fabs((A(2,1)-A(3,1))*(A(0,0)-A(3,0))));
    a3 = fabs(A(2,2)-A(3,2))*(fabs((A(0,0)-A(3,0))*(A(1,1)-A(3,1)))+fabs((A(0,1)-A(3,1))*(A(1,0)-A(3,0))));

    if(fabs(A.determinant()) <=  MERROR*(a1+a2+a3))
    {
        #pragma omp critical
        {
        cout << "Volume is inside the tolerance of round-off error. Need to do some form of tiebreaking." << endl;
        cout << "Matrix: " << endl;
        cout << std::scientific << std::setprecision(5);
        cout << setw(13) << A(0,0) << setw(13) << A(0,1) << setw(13) << A(0,2) << setw(13) << A(0,3) << endl;
        cout << setw(13) << A(1,0) << setw(13) << A(1,1) << setw(13) << A(1,2) << setw(13) << A(1,3) << endl;
        cout << setw(13) << A(2,0) << setw(13) << A(2,1) << setw(13) << A(2,2) << setw(13) << A(2,3) << endl;
        cout << setw(13) << A(3,0) << setw(13) << A(3,1) << setw(13) << A(3,2) << setw(13) << A(3,3) << endl;
        // cout << "MERROR: " << MERROR << endl;
        // cout << "Components: " << a1 << "  " << a2 << "  " << a3 << endl;
        cout << A.determinant() <<  " <= " <<  MERROR*(a1+a2+a3) << endl;
        }
        return TRUE;
    }

    return FALSE;
}

// Returns 1 if the lines intersect, otherwise 0. In addition, if the lines 
// intersect the intersection point may be stored in the floats i_x and i_y.
bool get_line_intersection(vector<vec<real,2>> const& verts, vector<size_t> const& edge, 
    vec<real,2> const& p1, vec<real,2> const& cellC)
{
    vec<real,2> const& e1 = verts[edge[0]];
    vec<real,2> const& e2 = verts[edge[1]];
    vec<real,2> s,r;
    s = cellC - p1; r = e2 - e1;

    // Find the denominator of the calculation 
    real denom =  (-r[0] * s[1] + s[0] * r[1]);

    // If the value of this is nearly 0, then 
    if(denom < 2*MEPSILON)
    {   // Lines are colinear
        return 0;
    }

    real u, t;
    u = (-s[1] * (p1[0] - e1[0]) + s[0] * (p1[1] - e1[1])) / denom;
    t = ( r[0] * (p1[1] - e1[1]) - r[1] * (p1[0] - e1[0])) / denom;

    if (u > 0 && u < 1 && t > 0 && t < 1)
    {   // Collision determinantected
        return 1;
    }

    return 0; // No collision
}


/* ======= Crossings algorithm ============================================  */
/* By Eric Haines, 3D/Eye Inc, erich@eye.com                                 */
/* Shoot a test ray along +X axis.  The strategy, from MacMartin, is to      */
/* compare vertex Y values to the testing point's Y and quickly discard      */
/* edges which are entirely to one side of the test ray.                     */
/*                                                                           */
/* Input 2D polygon _pgon_ with _numverts_ number of vertices and test point */
/* _point_, returns 1 if inside, 0 if outside.  WINDING and CONVEX can be    */
/* defined for this test.                                                    */
int Crossings2D(vector<vec<real,2>> const& verts, vector<size_t> const& edge, vec<real,2> const& point)
    {
        int  yflag0, yflag1, inside_flag;
        real  ty, tx;
        vec<real,2> vtx0, vtx1;

        tx = point[0];
        ty = point[1];

        inside_flag = 0;

        vtx0 = verts[edge[0]];
        vtx1 = verts[edge[1]];
        /* Move to the next pair of vertices, retaining info as possible. */
        yflag0 = ( vtx0[1] >= ty );
        yflag1 = ( vtx1[1] >= ty );

        /* Check if endpoints straddle (are on opposite sides) of X axis
         * (i.e. the Y's differ); if so, +X ray could intersect this edge.
         * Credit to Joseph Samosky to try dropping
         * the "both left or both right" part of my code.
         */
        if ( yflag0 != yflag1 ) 
        {
            /* Check intersection of pgon segment with +X ray.
             * Note if >= point's X; if so, the ray hits it.
             * The division operation is avoided for the ">=" test by checking
             * the sign of the first vertex wrto the test point; idea inspired
             * by Joseph Samosky's and Mark Haigh-Hutchinson's different
             * polygon inclusion tests.
             */
            if ( ((vtx1[1]-ty) * (vtx1[0]-vtx0[0]) >=
              (vtx1[0]-tx) * (vtx1[1]-vtx0[1])) == yflag1 )
            {
              inside_flag = !inside_flag;
            }

            /* For convex cells, further optimisation can be done: */
            /* A ray can only pass through a maximum of two faces.*/
            /* If this is second edge hit, then done testing. */ 
        }            
        return( inside_flag );
    }

int Crossings3D(vector<vec<real,3>> const& verts, vector<size_t> const& face, 
    vec<real,3> const& point, vec<real,3> const& point2, bool& perturb)
{   /*Using Signed volumes of tetrahedra*/
    /*Shewchuk J.R 1996, & 
    Robust Adaptive Floating-Point Geometric Predicates
    Michael Aftosmis, Cart3D Software*/
    /*https://www.nas.nasa.gov/publications/software/docs/cart3d/pages/boolintersection.html*/
    vec<real,3> const testp = point; 
    vec<real,3> const rayp = point2;
    matrix<real,4> vol1;
    int flag1, flag2;
    vol1  = matrix<real,4>(testp[0], testp[1], testp[2] , 1.0,
            verts[face[0]][0], verts[face[0]][1], verts[face[0]][2], 1.0,
            verts[face[1]][0], verts[face[1]][1], verts[face[1]][2], 1.0,
            verts[face[2]][0], verts[face[2]][1], verts[face[2]][2], 1.0);
     

    if(LessThanREError(vol1))
    {   // Perturb the test point so that it doesn't go into roundoff error  
        perturb = TRUE;
        return 0;  
    }
    
    flag1 = (vol1.determinant() < 0.0);

    vol1.row(0,vec<real,4>(rayp[0], rayp[1], rayp[2],1.0));

    // ray point is very far so I don't see this falling into roundoff error 
    if(LessThanREError(vol1))
    {     
        perturb = TRUE;
        return 0;  
    }

    flag2 = (vol1.determinant() < 0.0); 

    /*If signs of the volumes alternate, then the points lie either side of the plane*/
    /*Now check if the line drawn by the two points intersects inside the bounds of the triangle plane*/
    if(flag1 != flag2)
    {   
        matrix<real,4> vol;     
        int flag3, flag4;

        vec<real,3> vtx0, vtx1;
        vtx0 = verts[face.back()]; /*Start on the last - first point edge*/
        vtx1 = verts[face[0]];
        
        /*Find initial volume size*/
        
        vol = matrix<real,4>(testp[0], testp[1], testp[2], 1.0,
                  vtx0[0], vtx0[1], vtx0[2], 1.0,
                  vtx1[0], vtx1[1], vtx1[2], 1.0,
                  rayp[0], rayp[1], rayp[2],1.0);

        if(LessThanREError(vol))
        {   /*Perturb the test point, since all volume calculations need to be done with the new location*/
            perturb = TRUE;
            return 0; 
        }
        
        flag3 = (vol.determinant() < 0.0);

        /*Check for each face, if the signs of all the tets are the same.*/
        for (size_t ii = 1; ii < face.size(); ++ii)
        {   /*Change the face vertices used*/
            // cout << face.size() << "  " << ii << endl;
            vtx0 = vtx1;
            vtx1 = verts[face[ii]];

            vol.row(1, vec<real,4>(vtx0[0], vtx0[1], vtx0[2], 1.0));
            vol.row(2, vec<real,4>(vtx1[0], vtx1[1], vtx1[2], 1.0));

            if(LessThanREError(vol))
            {
                perturb = TRUE;
                return 0; 
            }
        
            flag4 = (vol.determinant() < 0.0);

            /*If the sign of the tet is different, this face isn't intersected.*/
            if (flag4 != flag3)
                return 0;   
        }  
        return 1;
    }  
    return 0;    
}

int Crossings3D_P(vector<vec<real,3>> const& verts, vector<size_t> const& face, 
    vec<real,3> const& point, vec<real,3> const& point2, bool& perturb)
{   /* Using Signed volumes of tetrahedra, Shewchuk J.R 1996, */
    /* Robust Adaptive Floating-Point Geometric Predicates, Michael Aftosmis, Cart3D Software*/
    /* https://www.nas.nasa.gov/publications/software/docs/cart3d/pages/boolintersection.html*/
    vec<real,3> const testp = point; 
    vec<real,3> const rayp = point2;
    matrix<real,4> vol1;
    int flag1, flag2;

    vec<real,3> v1(verts[face[0]][0], verts[face[0]][1], verts[face[0]][2]);
    vec<real,3> v2(verts[face[1]][0], verts[face[1]][1], verts[face[1]][2]);
    vec<real,3> v3(verts[face[2]][0], verts[face[2]][1], verts[face[2]][2]);

    /* Perturb based on the index in the array, and coordinate direction */
    /* Ensures the perturbation is unique for each point, and reproducible */
    vec<real,3> p1 = vec<real,3>(PERTURB(0,0),PERTURB(0,1),PERTURB(0,2));
    vec<real,3> p2 = vec<real,3>(PERTURB(1,0),PERTURB(1,1),PERTURB(1,2));
    vec<real,3> p3 = vec<real,3>(PERTURB(2,0),PERTURB(2,1),PERTURB(2,2));

    v1 += p1;
    v2 += p2;
    v3 += p3;


    vol1  = matrix<real,4>(testp[0], testp[1], testp[2], 1.0,
                 v1[0]   , v1[1]   , v1[2]   , 1.0,
                 v2[0]   , v2[1]   , v2[2]   , 1.0,
                 v3[0]   , v3[1]   , v3[2]   , 1.0);
     
    if(LessThanREError(vol1))
    {   // Perturb the test point so that it doesn't go into roundoff error  
        perturb = TRUE;
        return 0;  
    }
    
    flag1 = (vol1.determinant() < 0.0);

    vol1.row(0,vec<real,4>(rayp,1.0));

    flag2 = (vol1.determinant() < 0.0); 

    /*If signs of the volumes alternate, then the points lie either side of the plane*/
    /*Now check if the line drawn by the two points intersects inside the bounds of the triangle plane*/
    if(flag1 != flag2)
    {   
        matrix<real,4> vol;     
        int flag3, flag4;
        
        /*Find initial volume size of back edge*/
        vol = matrix<real,4>(testp[0], testp[1], testp[2], 1.0,
                  v3[0], v3[1], v3[2], 1.0,
                  v1[0], v1[1], v1[2], 1.0,
                  rayp[0], rayp[1], rayp[2],1.0);

        if(LessThanREError(vol))
        {   /*Perturb the test point, since all volume calculations need to be done with the new location*/
            perturb = TRUE;
            return 0; 
        }

        flag3 = (vol.determinant() < 0.0);

        /* Edge 1 */
        vol.row(1, vec<real,4>(v1, 1.0));
        vol.row(2, vec<real,4>(v2, 1.0));

        if(LessThanREError(vol))
        {
            perturb = TRUE;
            return 0; 
        }
    
        flag4 = (vol.determinant() < 0.0);
        if (flag4 != flag3)
            return 0;   

        /* Edge 2 */
        vol.row(1, vec<real,4>(v2, 1.0));
        vol.row(2, vec<real,4>(v3, 1.0));

        if(LessThanREError(vol))
        {
            perturb = TRUE;
            return 0; 
        }
    
        flag4 = (vol.determinant() < 0.0);
        if (flag4 != flag3)
            return 0;   
         
        return 1;
    }  
    return 0;    
}


/* Function to check if the cell contains the point being tested, and if it needs perturbing */
bool CheckCell(size_t const& cell, MESH const& cells, vec<real,3> const& testp, vec<real,3> const& rayp, uint& pert)
{
    bool line_flag = 0;
    bool inside_flag = 0;
    bool perturb = 0;
    for (auto const& cFaces:cells.cFaces[cell]) 
    {   // Step through cell faces
        vector<size_t> const& face = cells.faces[cFaces];
        // if(Crossings2D(cells.verts,face,testp))  

        if(Crossings3D(cells.verts,face,testp,rayp,perturb))             
        {   
            /* Intersects a face */
            inside_flag=!inside_flag;
            if ( line_flag ) //Convex assumption: can hit a max of two faces
                return FALSE; 
            //  note that one edge has been hit by the ray's line 
            line_flag = TRUE;
        }

        if(perturb == TRUE)
        {   /* If it needed perturbing, check the face with perturbation */
            cout << "face needs perturbing" << endl;
            // pert++;
            if(Crossings3D_P(cells.verts,face,testp,rayp,perturb))
            {
                /* Intersects a face */
                inside_flag=!inside_flag;
                if ( line_flag ) //Convex assumption: can hit a max of two faces
                    return FALSE; 
                //  note that one edge has been hit by the ray's line 
                line_flag = TRUE;
            }
            // return FALSE;
        }

    }
    return inside_flag;
}

/* Find the cells for all non-boundary particles. Checks for whether a paricle is free.  */
/* It is assumed that the particle does have a previous cell defined, and checks this cell */
/* before going to the KD Tree to find the nearest cells to iterate through. </summary> */
void FindCell(SETT const& svar, Vec_Tree const& TREE, MESH const& cells, part& pi)
{
    /*Find which cell the particle is in*/  
    vec<real,3> testp = pi.xi;  
    vec<real,3> const rayp(testp[0]+1e5,testp[1],testp[2]);
    uint inside_flag = 0;
    uint pert = 0;
    
    /*Perform a small search, since most cells are found within 1 or 2 indexes.*/
    vector<size_t> outlist;
    FindCellNeighbours(TREE, 5, testp, outlist);

    uint count = 0;
    for(auto const& cell:outlist)
    {   
        if(CheckCell(cell,cells,testp,rayp,pert))
        {
            #ifdef DEBUG
            // cout << "Moved to a neighbour cell." << endl;
            #endif
            pi.faceID = -1; /* We're inside a cell */
            pi.cellID = cell;
            pi.cellV = cells.cVel[cell];
            pi.cellRho = cells.cRho[cell];
            inside_flag = 1;
            pi.nNotFound = 0;
            break;
            // cout << "Found new cell at " << count << " in list" << endl;
        }

        // if(pert == 1)
        // {    
        //     testp = pi.xi + vec<real,3>(PERTURB(0,1),PERTURB(0,2),PERTURB(0,3));
        
        //     if(CheckCell(pi.cellID,cells,testp,rayp,pert))
        //     {
        //         inside_flag = 1;
        //         pi.nNotFound = 0;
        //         pi.cellID = cell;
        //         pi.cellV = cells.cVel[cell];
        //         pi.cellRho = cells.cRho[cell];
        //         break;
        //     }
        // }

        // if(pert > 1)
        // {
        //     cout << "Still can't identify which cell point is in even after perturbation" << endl; 
        // }

        count++;
    }
    

    
    if(inside_flag == 0)
    {   /*If first search fails to find the cell, perform a large search*/
        // In tets, the cell centre can be very far away, making the index large
        vector<size_t> outlist;
        FindCellNeighbours(TREE, 100, testp, outlist);

        uint count = 0;
        for(auto const& cell:outlist)
        {   
            if(CheckCell(cell,cells,testp,rayp,pert))
            {
                #ifdef DEBUG
                    // cout << "Moved to a neighbour cell." << endl;
                #endif
                pi.faceID = -1; /* We're inside a cell */
                pi.cellID = cell;
                pi.cellV = cells.cVel[cell];
                pi.cellRho = cells.cRho[cell];
                inside_flag = 1;
                pi.nNotFound = 0;

                // cout << "Found new cell at " << count << " in list" << endl;
            }

            if(pert == 1)
            {
                
                testp = pi.xi + vec<real,3>(PERTURB(0,1),PERTURB(0,2),PERTURB(0,3));

                if(CheckCell(pi.cellID,cells,testp,rayp,pert))
                {
                    inside_flag = 1;
                    pi.nNotFound = 0;
                    pi.faceID = -1; /* We're inside a cell */
                    pi.cellID = cell;
                    pi.cellV = cells.cVel[cell];
                    pi.cellRho = cells.cRho[cell];
                    break;
                }
            }

            if(pert > 1)
            {
                cout << "Still can't identify which cell point is in even after perturbation" << endl; 
            }

            count++;
        }

        if(inside_flag == 0)
        {
            // cout << "Checking if particle " << pi.partID << " has crossed a boundary." << endl;
            // If still not found, then the point could be across a boundary.
            // Check if the ray from the point to the cell centre crosses a boundary face.
            uint cross = 0;
            for(auto const& index:outlist)
            {
                vec<real,3> rayp = cells.cCentre[index];            
                for (size_t const& findex:cells.cFaces[index] ) 
                {   /*If its a boundary face, check if the point crosses it*/
                    if(cells.leftright[findex].second < 0)
                    {
                        vector<size_t> const& face = cells.faces[findex];
                        int ints;
                                        
                        bool perturb = FALSE;              
                        ints = Crossings3D(cells.verts,face,testp,rayp,perturb);
                        if(perturb == TRUE)
                        {
                            ints = Crossings3D_P(cells.verts,face,testp,rayp,perturb);
                        }
                        
                        if(ints)
                        {
                            cross=1;
                            if(cells.leftright[findex].second == -1)
                            {
                                // cout << "Particle " << pi.partID << " has crossed a boundary!" << endl;
                                /* Use the faceID (Not used otherwise) to identify if particle is outside the cell */
                                pi.faceID = -2;
                                goto matchfound;
                            }
                            
                        }  
                    }
                }

            }
            matchfound:	
            if(cross == 0)
            {
                if(pi.nNotFound > 10)
                {
                    cout << "Particle " << pi.partID << " containing cell not found. Something is wrong." << endl;
                    cout << "Position: " << pi.xi[0] << "  " << pi.xi[1] << "  " << pi.xi[2] << endl;
                    // Write_Containment(svar,ret_indexes,cells,testp);
                    // exit(-1);
                    pi.going = 0;
                }
                else
                {
                    cout << "Couldn't locate particle " << pi.partID << " cell this time." << endl;
                    cout << "Count: " << pi.nNotFound << endl;
                    pi.nNotFound++;
                }
            }

        }
    }

       
}


/********************** IMPLICIT TRACKING FUNCTIONS **************************/


/* Check for intersection with the infinite plane (I.e. just do two volumes) */
int Cross_Plane(vector<vec<real,3>> const& verts, vector<size_t> const& face, 
    vec<real,3> const& point, vec<real,3> const& point2, bool& perturb)
{   /*Using Signed volumes of tetrahedra*/
    /*Shewchuk J.R 1996 */
    vec<real,3> const testp = point; 
    vec<real,3> const rayp = point2;
    matrix<real,4> vol1;
    int flag1, flag2;
    vol1  = matrix<real,4>(testp[0], testp[1], testp[2] , 1.0,
            verts[face[0]][0], verts[face[0]][1], verts[face[0]][2], 1.0,
            verts[face[1]][0], verts[face[1]][1], verts[face[1]][2], 1.0,
            verts[face[2]][0], verts[face[2]][1], verts[face[2]][2], 1.0);
     

    if(LessThanREError(vol1))
    {   // Perturb the test point so that it doesn't go into roundoff error  
        perturb = TRUE;
        return 0;  
    }
    
    flag1 = (vol1.determinant() < 0.0);

    vol1.row(0,vec<real,4>(rayp[0], rayp[1], rayp[2],1.0));

    // ray point is very far so I don't see this falling into roundoff error 
    if(LessThanREError(vol1))
    {     
        perturb = TRUE;
        return 0;  
    }

    flag2 = (vol1.determinant() < 0.0); 

    /* If signs of the volumes alternate, */
    /* then the ray intersects the infinite plane*/
    if(flag1 != flag2)
    {   
        return 1;
    }  
    return 0;    
}


/* Check for intersection with the infinite plane (I.e. just do two volumes) */
int Cross_Plane_P(vector<vec<real,3>> const& verts, vector<size_t> const& face, 
    vec<real,3> const& point, vec<real,3> const& point2, bool& perturb)
{   /*Using Signed volumes of tetrahedra*/
    /*Shewchuk J.R 1996 */
    vec<real,3> const testp = point; 
    vec<real,3> const rayp = point2;
    matrix<real,4> vol1;
    int flag1, flag2;
    vol1  = matrix<real,4>(testp[0], testp[1], testp[2] , 1.0,
            verts[face[0]][0] + PERTURB(0,0), verts[face[0]][1] + PERTURB(0,1), verts[face[0]][2] + PERTURB(0,2), 1.0,
            verts[face[1]][0] + PERTURB(1,0), verts[face[1]][1] + PERTURB(1,1), verts[face[1]][2] + PERTURB(1,2), 1.0,
            verts[face[2]][0] + PERTURB(2,0), verts[face[2]][1] + PERTURB(2,1), verts[face[2]][2] + PERTURB(2,2), 1.0);
     

    if(LessThanREError(vol1))
    {   // Perturb the test point so that it doesn't go into roundoff error  
        perturb = TRUE;
        return 0;  
    }
    
    flag1 = (vol1.determinant() < 0.0);

    vol1.row(0,vec<real,4>(rayp[0], rayp[1], rayp[2],1.0));

    // ray point is very far so I don't see this falling into roundoff error 
    if(LessThanREError(vol1))
    {     
        perturb = TRUE;
        return 0;  
    }

    flag2 = (vol1.determinant() < 0.0); 

    /* If signs of the volumes alternate, */
    /* then the ray intersects the infinite plane*/
    if(flag1 != flag2)
    {   
        return 1;
    }  
    return 0;    
}

/* Function to check if the cell contains the point being tested, and if it needs perturbing */
/* Change function to find the shortest intersection with the infinite planes of each face */
vector<lint> CheckCellFace(size_t const& cell, MESH const& cells, size_t const& curr_face, 
                    vec<real,3> const& testp, vec<real,3> const& rayp)
{
    vector<lint> intersects;
    bool perturb = 0;
    // size_t face_intersect = std::numeric_limits<size_t>::max(); 
    for (auto const& cFace:cells.cFaces[cell]) 
    {   // Step through cell faces

        vector<size_t> const& face = cells.faces[cFace];
        // if(Crossings2D(cells.verts,face,testp))  

        if(Cross_Plane(cells.verts,face,testp,rayp,perturb))             
        {   /* Intersects a face */
            intersects.emplace_back(static_cast<lint>(cFace));
        }

        if(perturb)
        {
            perturb = FALSE;
            if(Cross_Plane_P(cells.verts,face,testp,rayp,perturb))
            {
                intersects.emplace_back(static_cast<lint>(cFace));
            }

            if(perturb)
            {
                cout << "Still can't tell if face was crossed after perturbing" << endl;
            }
            // cout << "point needs perturbing" << endl;
            // pert++;
            // return FALSE;
        }

        /* Assuming convex cells, and point is residing on a face, */
        /* stop at the first interecting face */
    }
    return intersects;
}

/* Möller–Trumbore intersection algorithm */
/* It has already been confirmed that the line intersects the triangle */
/* Only need to find the timestep */
void RayIntersectsTriangle(MESH const& cells, vec<real,3> const& rayOrigin, vec<real,3> const& rayVector,
                           vector<size_t> const& face, real& dt)
{
    vec<real,3> vertex0 = cells.verts[face[0]];
    vec<real,3> vertex1 = cells.verts[face[1]];  
    vec<real,3> vertex2 = cells.verts[face[2]];
    vec<real,3> edge1, edge2, h, s, q;
    real f;
    edge1 = vertex1 - vertex0;
    edge2 = vertex2 - vertex0;
    h = rayVector.cross(edge2);
    /* Check for sign to see if normal is pointing in or out */

    f = 1.0/(edge1.dot(h));
    
    s = rayOrigin - vertex0;
    q = s.cross(edge1);

    // At this stage we can compute t to find out where the intersection point is on the line.
    dt = f * edge2.dot(q);


    // if (t > EPSILON) // ray intersection
    // {
    //     outIntersectionPoint = rayOrigin + rayVector * t;
    //     return true;
    // }
    // else // This means that there is a line intersection but not a ray intersection.
    //     return false;
}

void RayNormalIntersection(MESH const& cells, vec<real,3> const& rayOrigin, vec<real,3> const& rayVector,
                           vector<size_t> const& face, lint const& cellID, real& dt, real& denom)
{
    /* find the normal vector of the face */
    vec<real,3> norm;
    if(face.size() == 3)
    {   /* Cross two edges of the triangle */
        norm = (cells.verts[face[1]]-cells.verts[face[0]]).cross(cells.verts[face[2]]-cells.verts[face[0]]);
    }   
    else
    {   /* Cross the two diagonals of the square */
        norm = (cells.verts[face[3]]-cells.verts[face[1]]).cross(cells.verts[face[2]]-cells.verts[face[0]]);
    }

    /* Check for normal orientation */
    vec<real,3> celldir;

    
    /* Use the face centre as the second part for the direction */
    vec<real,3> face_c(0.0);
    for(size_t ii = 0; ii < face.size(); ii++)
    {
        face_c += cells.verts[face[ii]];
    }

    face_c /= real(face.size());

    celldir = face_c - cells.cCentre[cellID];
    

    if(norm.dot(celldir) < 0)
    {
        /* Normal points inwards to the cell , so flip it*/
        norm = -1.0 * norm;
    }

    vec<real,3> temp_p(0.0);

    /* Find the most distant point from the current point */
    for(size_t ii = 0; ii < face.size(); ii++)
    {
        if( (cells.verts[face[ii]] - rayOrigin).squaredNorm() > temp_p.squaredNorm())
            temp_p = (cells.verts[face[ii]] - rayOrigin);
    }

    /* Find numerator */
    denom = rayVector.dot(norm);
    dt = temp_p.dot(norm)/denom; 
}


void FindFace(SETT const& svar, MESH const& cells, part const& pn, part& pnp1)
{
    /* Want to find the intersection of the droplet vector with a cell face. */

    /* Velocity vector = average of start and end */
    vec<real,3> testv = 0.5*(pnp1.v + pn.v);

    /* Test point: starting position*/
    vec<real,3> testp = pn.xi - 1e5*testv;
    vec<real,3> rayp = pn.xi + 1e5*testv;

    vector<lint> intersects;

    intersects = CheckCellFace(pn.cellID,cells,pn.faceID,testp,rayp);
    
    // {
    //     // if(pert == 1)
    //     // {
    //     //     testp = pn.xi + vec<real,3>(PERTURB(0,1),PERTURB(0,2),PERTURB(0,3));
            
    //     //     if(!CheckCellFace(pn.cellID,cells,pn.faceID,testp,rayp,pnp1.faceID,pert))
    //     //     {
    //     //         cout << "Still can't identify which cell point is in even after perturbation" << endl;
    //     //         pnp1.going = 0; 
    //     //     }
    //     // }
    //     // else
    //     // {
    //         cout << "Couldn't tell which face was crossed." << endl;
    //         pnp1.going = 0;
    //     // }
    // }

    /* Found the faces that were crossed. Now find the time spent in the cell */
    
    lint nextface = pn.faceID; /* At least use a reasonable face */
    real mindist = 1e6;
    int hasintersect = 0;
    for(size_t ii = 0; ii < intersects.size(); ++ii)
    {
        // /* Ignore the current face. It's not an option. */
        if(intersects[ii] == pn.faceID)
            continue;

        vector<size_t> const& face = cells.faces[intersects[ii]];
        
        real dt = 0.0;
        real denom = 0.0;
        // RayIntersectsTriangle(cells, pn.xi, testv, face, dt);

        RayNormalIntersection(cells, pn.xi, testv, face, pn.cellID, dt, denom);

        /* Check for if the face normal points in or out of the cell */
        // if(cells.leftright[intersects[ii]].first != pn.cellID)
        // {
        //     /* face normal is pointing into the cell, so flip the sign*/
        //     dt = -dt;
        // }

        if(denom > 0)
        {
            /* face normal product with the velocity vector is positive */
            /* so it isn't behind the particle */
            
            if ( dt < mindist)
            {
                nextface = intersects[ii];
                if(dt > MEPSILON)
                {
                    mindist = dt;
                    // cout << "Updated dt: " << dt  << " Crossing face: " << nextface << endl;
                }
                else
                {
                    /* Calculation has resulted in a roundoff suggesting a */
                    /* negative intersection distance, so set distance to zero */
                    /* so particle containment is updated, but not position */
                    // cout << "distance is near 0" << endl;
                    mindist = MEPSILON;
                }
            }

            hasintersect = 1;
        }        
    }

    if (hasintersect == 0)
    {
        cout << "Could not find a positive intersection for some reason..." << endl;
        cout << "No of intersections with planes: " << intersects.size() << endl;
        cout << "ParticleID: " << pnp1.partID << " position: " << pn.xi[0] << "  " << pn.xi[1]
        << "  " << pn.xi[2] << " velocity n: " << pn.v[0] << "  " << pn.v[1] << "  " << pn.v[2] <<
         " velocity np1: " << pnp1.v[0] << "  " << pnp1.v[1] << "  " << pnp1.v[2] << endl;
        pnp1.xi = pn.xi;
        pnp1.v = pn.v;
        pnp1.dt = 0.0;

        pnp1.going = 0;
    }
    else
    {
        // auto min_t = std::min_element(dists.begin(),dists.end(),
        // [](std::pair<real,size_t> const& p1, std::pair<real,size_t> const& p2)
        // {return p1.first < p2.first;});

        /* This is the face that has actually been crossed */
        // pnp1.dt = min_t->first;
        // pnp1.faceID = min_t->second;

        pnp1.dt = mindist;
        pnp1.faceID = nextface;
            
        /* Update the face and cell data */
        lint newcell=0;
        
        if(cells.leftright[nextface].first == static_cast<lint>(pn.cellID))    
            newcell = cells.leftright[nextface].second;
        else
            newcell = cells.leftright[nextface].first;

        // cout << "first: " << cells.leftright[pnp1.faceID].first << " second: " <<
        //     cells.leftright[pnp1.faceID].second << " cellID: " << pn.cellID << " newcell: "
        //     << newcell << endl;
        
        if(newcell > -1)
        {
            pnp1.cellID = newcell; 
            pnp1.cellV = cells.cVel[newcell];
            pnp1.cellRho = cells.cRho[newcell];
        }
        else
        {
            pnp1.cellID = newcell; 
        }
    }  
}

#endif