#ifndef POST_H
#define POST_H

#include "Var.h"

void Post_Process(SETT const& svar, vector<State> const& trajectories, 
            vector<vector<SURF>>& surface_faces, vector<SURF>& marker_data)
{
    uint nSurf = 0;
    vector<int> marks;
    vector<string> names;
    for(size_t ii = 0; ii < svar.bwrite.size(); ii++)
    {
        nSurf += svar.bwrite[ii];
        if(svar.bwrite[ii] == 1)
        {
            marks.emplace_back(svar.markers[ii]);
            names.emplace_back(svar.bnames[ii]);
        }
    }

    /* Get impact data */
    /* Process the MVDs */    
    // vector<vector<real>> MVD_means(nSurf, vector<real>()); /* Average of Medial Volume Diameter */
    // vector<vector<real>> MVD_devs(nSurf, vector<real>()); /* Standard deviation of MVD */   
    // vector<vector<int>> counts(nSurf, vector<int>());
    // vector<vector<real>> mass(nSurf, vector<real>());

    vector<vector<real>> beta(nSurf, vector<real>()); /* Collection efficiency */
    vector<vector<real>> area(nSurf, vector<real>()); /* Area used to calculate */

    

    for(size_t ii = 0; ii < surface_faces.size(); ++ii)
    {
        vector<SURF>& surface = surface_faces[ii];

        /* Create a KD Tree of the finished points */
        // Vec_Tree TREE(3,surface.pos,20);
        

        for(size_t jj = 0; jj < surface.size(); ++jj)
        {
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
                
                if(surface[jj].count > 1)
                {
                    // Vec_Tree TREE(3,surface[jj].end_pos,2);

                    // cout << "Count is greater than 1" << endl;
                    vector<real> start_dists;
                    vector<real> end_dists;
                    if(svar.startType == 1)
                    {
                        /* Use 2D line distances */
                        /* Find the distances for each particle pair */                        
                        for(size_t kk = 0; kk < surface[jj].end_pos.size(); ++kk)
                        {
                            for(size_t ll = kk; ll < surface[jj].end_pos.size(); ++ll)
                            {
                                if(kk == ll)
                                    continue;

                                real start =  (surface[jj].start_pos[kk] - surface[jj].start_pos[ll]).norm();

                                real end = (surface[jj].end_pos[kk] - surface[jj].end_pos[ll]).norm();

                                start_dists.emplace_back(start);
                                end_dists.emplace_back(end);
                            }
                        }

                    }
                    else
                    {
                        /* Need to use triangle areas */
                        if(surface[jj].end_pos.size() > 2)
                        {                            
                            for(size_t kk = 0; kk < surface[jj].end_pos.size(); ++kk)
                            {
                                for(size_t ll = 0; ll  < surface[jj].end_pos.size(); ++ll)
                                {
                                    for(size_t mm = 0; mm < surface[jj].end_pos.size(); ++mm)
                                    {
                                        if(kk == ll || kk == mm || ll == mm)
                                            continue;

                                        vec<real,3> ab = surface[jj].start_pos[ll] - 
                                                    surface[jj].start_pos[kk];

                                        vec<real,3> ac = surface[jj].start_pos[mm] - 
                                                    surface[jj].start_pos[kk];        


                                        real start = 0.5 * (ab.cross(ac)).norm();

                                        ab = surface[jj].end_pos[ll] - surface[jj].end_pos[kk];
                                        ac = surface[jj].end_pos[mm] - surface[jj].end_pos[kk];


                                        real end = 0.5 * (ab.cross(ac)).norm();

                                        start_dists.emplace_back(start);
                                        end_dists.emplace_back(end);
                                    }
                                }
                            }
                        }
                    }

                    if(start_dists.empty())
                    {
                        surface[jj].colEff = 0.0;
                        surface[jj].area = 0.0;
                    }
                    else
                    {
                        /* get the ratios */
                        vector<real> betas(start_dists.size());
                        for(size_t kk = 0; kk < start_dists.size(); ++kk)
                        {
                            betas[kk] = start_dists[kk]/end_dists[kk];
                        }

                        /* Average the values */
                        real beta_ = 0.0;
                        real avg_area = 0.0;
                        for(size_t kk = 0; kk < start_dists.size(); ++kk)
                        {
                            beta_ += betas[kk];
                            avg_area += end_dists[kk];
                        }

                        beta_ /= real(betas.size());
                        avg_area /= real(betas.size());
                        surface[jj].colEff= beta_;
                        surface[jj].area = (avg_area);

                        // cout << beta_ << "  " << avg_area << endl;

                    }
                    
                }
                else
                {
                    surface[jj].colEff = 0.0;
                    surface[jj].area = 0.0;
                }
                
            // }
            
            // MVD_means[surf].emplace_back(MVDavg);
            // MVD_devs[surf].emplace_back(sum);

            // counts[surf].emplace_back(surface.count[jj]);
            // mass[surf].emplace_back(surface.mass[jj]);
            
        }
    }
}


void Get_Catch_Data()
{
    
}
#endif