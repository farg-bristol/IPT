#ifndef POST_H
#define POST_H

#include "Var.h"
#include "Integrate.h"

void Post_Process(SETT const& svar,  vector<int> const& marks,
    vector<vector<real>> const& beta_data, vector<vector<real>> const& area_data, 
            vector<SURF>& surface_data)
{

    /* Get impact data */ 
    for(size_t ii = 0; ii < marks.size(); ++ii)
    {
        SURF& surface = surface_data[marks[ii]];

        /* Create a KD Tree of the finished points */
        // Vec_Tree TREE(3,surface.pos,20);
        
        vector<vector<real>> beta_values(surface.faceIDs.size());
        vector<vector<real>> area_values(surface.faceIDs.size());
        

        for(size_t jj = 0; jj < surface.impacted_face.size(); ++jj)
        {
            /* Add to the face with the beta value */
            size_t faceID = surface.impacted_face[jj];
            beta_values[faceID].emplace_back(beta_data[ii][jj]);
            area_values[faceID].emplace_back(area_data[ii][jj]);
            
        }

        

        real avg_beta = 0.0;
        real avg_area = 0.0;
        for(size_t jj = 0; jj < surface.faceIDs.size(); ++jj)
        {
            /* Create a vector for the beta values */
            
            if(!beta_values[jj].empty())
            {
                /* Average the values */
                real beta_sum = 0.0;              
                real area_sum = 0.0;
                for(size_t kk = 0; kk < beta_values[jj].size(); kk++ )
                {
                    beta_sum += beta_values[jj][kk];
                    area_sum += area_values[jj][kk];
                }

                beta_sum /= real(beta_values[jj].size());
                area_sum /= real(beta_values[jj].size());

                surface.face_beta[jj] = beta_sum;
                surface.face_area[jj] = area_sum;
            }
            else
            {
                surface.face_beta[jj] = 0.0;
                surface.face_area[jj] = 0.0;
            }

            avg_beta += surface.face_beta[jj];
            avg_area += surface.face_area[jj];
            // real MVDavg = 0.0;
            // real sum = 0.0;

        
            // if(!surface.MVDs.empty())
            // {
            //     for(real const& MVD : surface.MVDs[jj])
            //     {
            //         MVDavg += MVD;
            //     }

            //     MVDavg /= real(surface.MVDs[jj].size());

            //     for(real const& MVD : surface.MVDs[jj])
            //     {
            //         sum += (MVD - MVDavg)*(MVD - MVDavg);
            //     }

            //     sum /= real(surface.MVDs[jj].size());

            //     sum = sqrt(sum);
                
            
                
            // }
            
            // MVD_means[surf].emplace_back(MVDavg);
            // MVD_devs[surf].emplace_back(sum);

            // counts[surf].emplace_back(surface.count[jj]);
            // mass[surf].emplace_back(surface.mass[jj]);
            
        }

        avg_beta /= surface.face_beta.size();
        avg_area /= surface.face_area.size();
        surface.marker_beta = avg_beta;
        surface.marker_area = avg_area;
    }
}


void Get_Catch_Data(SETT& svar, FLUID const& fvar, Vec_Tree const& TREE, MESH const& cells, 
    vector<SURF> const& surface_data, vector<int> const& marks, vector<State> const& part_time,
    vector<vector<real>>& beta_data, vector<vector<real>>& area_data, vector<vector<vec<real,3>>>& pos_data)
{
    size_t pID = part_time.size();
    /* For each face that wants a surface output, get the catch efficiency data */
    for(size_t jj = 0; jj < marks.size(); ++jj)
    {
        cout << "Getting catch data for marker: " << marks[jj] << endl;
        SURF const& test_surface = surface_data[marks[jj]];  
        
        cout << "Number of collisions: " << test_surface.marker_count << endl;

        beta_data.emplace_back();
        area_data.emplace_back();
        pos_data.emplace_back();

        if(svar.dimension == 2)
        {
            #pragma omp parallel 
            {
                vector<real> beta_local, area_local;
                vector<vec<real,3>> pos_local;

                #pragma omp for schedule(static) nowait
                for(size_t ii = 0; ii < test_surface.marker_count; ++ii)
                {
                    /* For each impacting point, do a trajectory above and below, to get the catch efficiency */
                    real delta;
                    if(svar.fac_or_fixed_delta == 0)
                    {
                        delta = svar.delta_fac*fvar.d_0;
                        if (delta > svar.max_delta)
                            delta = svar.max_delta;
                    }
                    else
                        delta = svar.delta;

                    vec<real,3> const& xi = test_surface.start_pos[ii];

                    vec<real,3> xi_1 = xi;
                    vec<real,3> xi_2 = xi;

                    xi_1[2] += delta;
                    xi_2[2] -= delta;

                    State beta_test, beta_testnp1, beta_testnm1;
                    beta_test.emplace_back(part(xi_1,part_time[test_surface.marker_pIDs[ii]][0],pID));
                    pID++;
                    beta_test.emplace_back(part(xi_2,part_time[test_surface.marker_pIDs[ii]][0],pID));
                    pID++;

                    /* Find the cell they reside it (not guaranteed to be the same as the original particle) */
                    FindCell(svar, TREE, cells, beta_test[0]);
                    FindCell(svar, TREE, cells, beta_test[1]);
                    beta_test[0].v = beta_test[0].cellV;
                    beta_test[1].v = beta_test[1].cellV;

                    beta_testnp1 = beta_test;
                    beta_testnm1 = beta_test;

                    vector<vector<part>> beta_time(beta_test.size());
                    
                    if(svar.partout == 1)
                    {
                        for(size_t part = 0; part < beta_test.size(); ++part)
                        {
                            string pfile = svar.outdir;
                            pfile.append("/particle_");
                            pfile.append(std::to_string(beta_test[ii].partID));
                            pfile.append("_scatter.dat");
                            // ofstream fout(pfile);
                            
                            svar.partfiles.emplace_back(ofstream{pfile});
                            
                            if(!svar.partfiles.back().is_open())
                            {
                                cout << "Failed to open output file. Path: " << pfile << endl;
                                exit(-1);
                            }
                            Write_ASCII_variables(svar.partfiles.back());
                            Write_ASCII_Scatter_Header(svar.partfiles.back(),ii);
                            Write_ASCII_Point(svar.partfiles.back(), svar.scale, beta_test[ii]);
                        }
                    }

                    vector<SURF> beta_marker;
                    // vector<vector<SURF>> beta_faces = surface_faces;

                    for(size_t jj = 0; jj < svar.markers.size(); jj++)
                    {
                        beta_marker.emplace_back(SURF(surface_data[jj]));
                        // beta_faces.emplace_back(surface_faces[jj]);
                    }
                    
                    if(svar.explicit_or_implicit == 0)
                    {
                        for(size_t part = 0; part < beta_test.size(); ++part)
                        {
                            beta_time[part].emplace_back(beta_test[part]);
                            Explicit_Integrate(svar, fvar, TREE, cells, part, beta_test[part], 
                                        beta_testnp1[part], beta_time[part], beta_marker);
                        }
                    }
                    else
                    {
                        for(size_t part = 0; part < beta_test.size(); ++part)
                        {
                            beta_time[part].emplace_back(beta_test[part]);
                            Implicit_Integrate(svar,fvar,cells, part, beta_testnm1[part], beta_test[part], 
                                        beta_testnp1[part], beta_time[part], beta_marker);
                        }
                    }

                    /* Make sure they have both hit the same surface */
                    if (beta_marker[marks[jj]].marker_count == 2)
                    {
                        /* Get the difference in end area. */
                        vec<real,3> diff = beta_marker[marks[jj]].end_pos[0] - beta_marker[marks[jj]].end_pos[1];

                        real dist = diff.norm();

                        real beta_i = 2*delta/dist;
                        
                        beta_local.emplace_back(beta_i);
                        area_local.emplace_back(dist);
                        pos_local.emplace_back(test_surface.end_pos[ii]);
                    }
                    // else if (beta_marker[marks[jj]].marker_count == 1)
                    // {
                    //     /* Only one has hit the same marker. So use the original point and the successful one */
                    //     vec<real,3> diff = beta_marker[marks[jj]].end_pos[0] - xi;

                    //     real dist = diff.norm();

                    //     real beta_i = delta/dist;

                    //     #pragma omp critical
                    //     {
                    //         beta_data[jj].emplace_back(beta_i);
                    //         area_data[jj].emplace_back(dist);
                    //         pos_data[jj].emplace_back(test_surface.end_pos[ii]);
                    //     }
                        
                    // }
                    else
                    {   /* Must have the beta structure be the same size as the count */
                        beta_local.emplace_back(0.0);
                        area_local.emplace_back(0.0);
                        pos_local.emplace_back(test_surface.end_pos[ii]);
                    }
                    
                    /* One will have to have made the same marker, surely... */
                    
                } /* End particle for */

                #pragma omp for schedule(static) ordered
                for(int i=0; i<omp_get_num_threads(); i++)
                {
                    #pragma omp ordered
                    {
                        beta_data[jj].insert(beta_data[jj].end(),beta_local.begin(),beta_local.end());
                        area_data[jj].insert(area_data[jj].end(),area_local.begin(),area_local.end());
                        pos_data[jj].insert(pos_data[jj].end(),pos_local.begin(),pos_local.end());
                    }
                }

            }/* End parallel section */
        }
        else
        {
            #pragma omp parallel 
            {
                vector<real> beta_local, area_local;
                vector<vec<real,3>> pos_local;

                #pragma omp for schedule(static) nowait
                for(size_t ii = 0; ii < test_surface.marker_count; ++ii)
                {
                    
                    /* For each impacting point, do a square, to get the catch efficiency */
                    real delta;
                    if(svar.fac_or_fixed_delta == 0)
                    {
                        delta = svar.delta_fac*fvar.d_0;
                        if (delta > svar.max_delta)
                            delta = svar.max_delta;
                    }
                    else
                        delta = svar.delta;

                    vec<real,3> const& xi = test_surface.start_pos[ii];

                    vec<real,3> xi_1 = xi;
                    vec<real,3> xi_2 = xi;
                    vec<real,3> xi_3 = xi;
                    vec<real,3> xi_4 = xi;

                    xi_1 += vec<real,3>(0.0,delta,delta);
                    xi_2 += vec<real,3>(0.0,-delta,delta);
                    xi_3 += vec<real,3>(0.0,-delta,-delta);
                    xi_4 += vec<real,3>(0.0,delta,-delta);
                    
                    State beta_test, beta_testnp1, beta_testnm1;
                    beta_test.emplace_back(part(xi_1,part_time[test_surface.marker_pIDs[ii]][0],pID));
                    pID++;
                    beta_test.emplace_back(part(xi_2,part_time[test_surface.marker_pIDs[ii]][0],pID));
                    pID++;
                    beta_test.emplace_back(part(xi_3,part_time[test_surface.marker_pIDs[ii]][0],pID));
                    pID++;
                    beta_test.emplace_back(part(xi_4,part_time[test_surface.marker_pIDs[ii]][0],pID));
                    pID++;

                    /* Find the cell they reside it (not guaranteed to be the same as the original particle) */
                    FindCell(svar, TREE, cells, beta_test[0]);
                    FindCell(svar, TREE, cells, beta_test[1]);
                    FindCell(svar, TREE, cells, beta_test[2]);
                    FindCell(svar, TREE, cells, beta_test[3]);
                    beta_test[0].v = beta_test[0].cellV;
                    beta_test[1].v = beta_test[1].cellV;
                    beta_test[2].v = beta_test[2].cellV;
                    beta_test[3].v = beta_test[3].cellV;

                    beta_testnp1 = beta_test;
                    beta_testnm1 = beta_test;

                    vector<vector<part>> beta_time(beta_test.size());

                    if(svar.partout == 1)
                    {
                        for(size_t part = 0; part < beta_test.size(); ++part)
                        {
                            string pfile = svar.outdir;
                            pfile.append("/particle_");
                            pfile.append(std::to_string(beta_test[ii].partID));
                            pfile.append("_scatter.dat");
                            // ofstream fout(pfile);
                            
                            svar.partfiles.emplace_back(ofstream{pfile});
                            
                            if(!svar.partfiles.back().is_open())
                            {
                                cout << "Failed to open output file. Path: " << pfile << endl;
                                exit(-1);
                            }
                            Write_ASCII_variables(svar.partfiles.back());
                            Write_ASCII_Scatter_Header(svar.partfiles.back(),ii);
                            Write_ASCII_Point(svar.partfiles.back(), svar.scale, beta_test[ii]);
                        }
                    }

                    vector<SURF> beta_marker;
                    // vector<vector<SURF>> beta_faces = surface_faces;

                    for(size_t jj = 0; jj < svar.markers.size(); jj++)
                    {
                        beta_marker.emplace_back(SURF(surface_data[jj]));
                        // beta_faces.emplace_back(surface_faces[jj]);
                    }

                    if(svar.explicit_or_implicit == 0)
                    {
                        for(size_t part = 0; part < beta_test.size(); ++part)
                        {
                            beta_time[part].emplace_back(beta_test[part]);
                            Explicit_Integrate(svar, fvar, TREE, cells, part, beta_test[part], 
                                        beta_testnp1[part], beta_time[part], beta_marker);
                        }
                    }
                    else
                    {
                        for(size_t part = 0; part < beta_test.size(); ++part)
                        {
                            beta_time[part].emplace_back(beta_test[part]);
                            Implicit_Integrate(svar,fvar,cells, part, beta_testnm1[part], beta_test[part], 
                                        beta_testnp1[part], beta_time[part], beta_marker);
                        }
                    }
                    
                    /* Make sure they have both hit the same surface */
                    if (beta_marker[marks[jj]].marker_count == 4)
                    {
                        /* Get the difference in end area. */
                        vec<real,3> ab = beta_marker[marks[jj]].end_pos[1] - beta_marker[marks[jj]].end_pos[0];
                        vec<real,3> ad = beta_marker[marks[jj]].end_pos[3] - beta_marker[marks[jj]].end_pos[0];

                        vec<real,3> cb = beta_marker[marks[jj]].end_pos[1] - beta_marker[marks[jj]].end_pos[2];
                        vec<real,3> cd = beta_marker[marks[jj]].end_pos[3] - beta_marker[marks[jj]].end_pos[2];

                        real area1 = 0.5 * (ad).cross(ab).norm();
                        real area2 = 0.5 * (cb).cross(cd).norm();

                        real beta_i = (2*delta)*(2*delta)/(area1+area2);

                        
                        beta_local.emplace_back(beta_i);
                        area_local.emplace_back(area1+area2);
                        pos_local.emplace_back(test_surface.end_pos[ii]);
                    }
                    else if (beta_marker[marks[jj]].marker_count == 3)
                    {
                        /* Only one has hit the same marker. So use the original point and the successful one */
                        vec<real,3> ab = beta_marker[marks[jj]].end_pos[1] - beta_marker[marks[jj]].end_pos[0];
                        vec<real,3> ac = beta_marker[marks[jj]].end_pos[2] - beta_marker[marks[jj]].end_pos[0];
                        real area1 = 0.5 * (ac).cross(ab).norm();

                        real beta_i = (2*delta*delta)/(area1);

                        beta_local.emplace_back(beta_i);
                        area_local.emplace_back(area1);
                        pos_local.emplace_back(test_surface.end_pos[ii]);
                    }
                    else
                    {   /* The size of the beta vector must be equal to that in the surface file */
                        beta_local.emplace_back(0.0);
                        area_local.emplace_back(0.0);
                        pos_local.emplace_back(test_surface.end_pos[ii]);
                    }
                    
                    // else if (beta_marker[marks[jj]].count == 2)
                    // {
                    //     /* Get the difference in end area. */
                    //     vec<real,3> diff = beta_marker[marks[jj]].end_pos[0] - beta_marker[marks[jj]].end_pos[1];

                    //     real dist = diff.norm();

                    //     real beta_i = 2*delta/dist;

                    //     beta_data.emplace_back(beta_i);
                    //     area_data.emplace_back(dist);
                    // }
                    // else if (beta_marker[marks[jj]].count == 1)
                    // {
                    //     /* Only one has hit the same marker. So use the original point and the successful one */
                    //     vec<real,3> diff = beta_marker[marks[jj]].end_pos[0] - test_surface.end_pos[ii];

                    //     real dist = diff.norm();

                    //     real beta_i = delta/dist;

                    //     beta_data.emplace_back(beta_i);
                    //     area_data.emplace_back(dist);

                    //     // vec<real,3> point = 0.5*(beta_marker[marks[jj]].end_pos[0] + test_surface.end_pos[ii]);
                    //     // pos_data.emplace_back(point);
                        
                    // }
                    // else
                    // {
                    //     beta_data.emplace_back(0.0);
                    //     area_data.emplace_back());
                    // }
                                    
                } /* End particle for */


                #pragma omp for schedule(static) ordered
                for(int i=0; i<omp_get_num_threads(); i++)
                {
                    #pragma omp ordered
                    {
                        beta_data[jj].insert(beta_data[jj].end(),beta_local.begin(),beta_local.end());
                        area_data[jj].insert(area_data[jj].end(),area_local.begin(),area_local.end());
                        pos_data[jj].insert(pos_data[jj].end(),pos_local.begin(),pos_local.end());
                    }
                }
            } /* End parallel section */
        } /* End if dimension */
    } /* End surface marker for */
}
#endif