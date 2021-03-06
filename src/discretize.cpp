#include "discretize.h"
#include "discretize.hu"
#include <iostream>
#include <math.h>
#include <stdexcept>
#include <stdlib.h>
#include <moab/CartVect.hpp>
#include <moab/Range.hpp>

#define CHECKERR(err) \
  if((err) != moab::MB_SUCCESS){ \
    std::cerr << "Error on line " << __LINE__ << "of discretize.cpp" \
              << std::endl; \
    return err;}

#define CHECKERR_HERE(err) \
  if((err) != moab::MB_SUCCESS) \
    std::cout << "Error on line " << __LINE__ << "of discretize.cpp" \
              << std::endl;

using moab::CartVect;
using moab::Core;
using moab::GeomTopoTool;
using moab::ErrorCode;
using moab::EntityHandle;
using moab::GeomQueryTool;

float VOL_FRAC_TOLERANCE = 1e-10;
moab::Interface* MBI;
GeomTopoTool* GTT;
GeomQueryTool* GQT;

std::vector<std::vector<double> > discretize_geom(
    std::vector<std::vector<double> > mesh,
    const char* filename,
    int num_rays,
    bool grid) {

  // This will store the information of the individual row and how the rays
  // are to be fired.
  ErrorCode rval;
  std::vector<EntityHandle> vol_handles;
  rval = load_geometry(filename, &vol_handles);
  CHECKERR_HERE(rval);
  // This holds information about this individual row.
  struct mesh_row row;
  row.num_rays = num_rays;
  row.grid = grid;
  if (grid && (pow((int)sqrt(num_rays), 2) != num_rays))
    throw std::runtime_error("For rays fired in a grid, num_rays must be "
                             "a perfect square.");

  // Initialize this here to avoid creating it multiple times during loops.
  std::vector<double> zeros(2,0.0);
  // This is needed for the function get_idx, in order to determine the
  // identities of the elements in the current row. It is also used often to
  // help with separating the work between cores, as it contains the number of
  // volume elements in each direction.
  int sizes[] = {mesh[0].size() - 1, mesh[1].size() - 1,
                  mesh[2].size() - 1};
  // This will be the end result.
  std::vector<std::vector<double> > result;
  // This stores the intermediate totals.
  std::vector<std::map<int, std::vector<double> > > row_totals;
  // Define the size now to use [] assignment.
  row_totals.resize(sizes[0]*sizes[1]*sizes[2]);

  // These for loops visit each individual row.
  for (int d1 = 0; d1 < 3; d1++) {

    // Set up the different direction indices.
    int d2 = (d1 + 1)%3;

    row.d3 = 3 - d1 - d2;
    row.d3divs = mesh[row.d3];

    // Iterate over each starting volume element in the d1/d2 directions.
    for (int i = 0; i < mesh[d1].size() - 1; i++) {
      row.d1div1 = mesh[d1][i];
      row.d1div2 = mesh[d1][i+1];

      for (int j = 0; j < mesh[d2].size() - 1; j++) {
        row.d2div1 = mesh[d2][j];
        row.d2div2 = mesh[d2][j+1];
        // int sizes[] = {mesh[0].size() - 1, mesh[1].size() - 1, mesh[2].size() - 1};
        std::vector<int> idx = get_idx(sizes, i, j, row.d3);

        // The rays are fired and totals collected here.
        std::vector<std::map<int, std::vector<double> > > ray_totals =
            fireRays(row, vol_handles);

        // Combine the totals collected from each direction by volume element.
        for (int k = 0; k < ray_totals.size(); k++) {
          for (std::map<int, std::vector<double> >::iterator ray_it =
              ray_totals[k].begin(); ray_it != ray_totals[k].end();
              ++ray_it) {
            if (ray_it->second[0] < VOL_FRAC_TOLERANCE)
              continue;
            std::map<int, std::vector<double> >::iterator row_it =
                row_totals[idx[k]].find(ray_it->first);
            if (row_it == row_totals[idx[k]].end()){
              row_totals[idx[k]].insert(row_it,
                                        std::make_pair(ray_it->first, zeros));
            }
            row_totals[idx[k]][ray_it->first][0] += ray_it->second[0];
            row_totals[idx[k]][ray_it->first][1] += ray_it->second[1];
          }
        }
      }
    }
  }

  // Evaluate the results and prepare them for return.
  int total_rays = num_rays*3;
  std::vector<double> ray_results(4,0.0);
  for (int i = 0; i < row_totals.size(); i++) {
    for (std::map<int, std::vector<double> >::iterator it =
         row_totals[i].begin(); it != row_totals[i].end(); ++it) {
      ray_results[0] = i;
      ray_results[1] = it->first;
      ray_results[2] = it->second[0]/total_rays;
      ray_results[3] = sqrt((it->second[1])/pow((it->second[0]),2.0)
                                    - 1.0/total_rays);
      result.push_back(ray_results);
    }
  }
  return result;
}

std::vector<std::map<int, std::vector<double> > > fireRays(mesh_row &row,
    std::vector<EntityHandle> vol_handles) {

  vec3 pt;
  pt[row.d3] = row.d3divs[0];
  int result = 0;
  vec3 dir = {0};
  dir[row.d3] = 1;
  EntityHandle eh;
  // This needs to come from the kernel; it doesn't.
  int num_intersections;
  // Same with this one.
  EntityHandle *volumes;

  // This vector would have to be huge. I don't think it would work to just 
  // figure out the distances on the GPU; if I were to redo this, I would assume
  // the calculations are also being done there. See comment at discretize.cu
  // line 44. Not that the memory required isn't already large.
  std::vector<double> holy_distances_batman = cuda_rayfire(MBI, GTT, GQT, row,
                                                           vol_handles); 

  for (int i = 0; i < row.num_rays; i++) {

    startPoints(row, i);
    pt[(row.d3+1)%3] = row.start_point_d1;
    pt[(row.d3+2)%3] = row.start_point_d2;
    // If the next point starts in the same volume as the last one, calling
    // this here can save calling find_volume, which gets expensive.
    if (i != row.init)
      GQT->point_in_volume(eh, pt, result, dir);
    if (!result){
      eh = find_volume(vol_handles, pt, dir);
    }
    std::vector<double> zeros(2,0.0);
    // This holds the numbers to add to the totals.
    double value;
    // Keep track of at which volume element we're currently looking.
    int count = 0;
    // Keep track of the "current" intersection in the geometry.
    int intersection = 0;
    bool complete = false;
    double curr_width = width[count];
    // Loops over the ray's intersections with the different volumes.
    while(!complete) {
      // while the distance to the next intersection is greater than the current
      // width, we calculate the value and we move to the next width, shortening
      // the distance.
      while (distances[i*row.num_rays + intersection] >= curr_width) {
        std::map<int, std::vector<double> >::iterator it =
            row_totals[count].find(eh);
        // If the current cell isn't in the totals yet, add it.
        if (it == row_totals[count].end()){
          row_totals[count].insert(it, std::make_pair(eh, zeros));
        }

        value = curr_width/width[count];
        row_totals[count][eh][0] += value;
        row_totals[count][eh][1] += value*value;
        distances[i*row.num_rays + intersection] -= curr_width;

        // If there are more volume elements, move to the next one.
        if (width.size() - 1 > count) {
          count++;
          curr_width = width[count];
        }
        // If you get here you've finished the mesh for this ray.
        else {
          complete = true;
          break;
        }
      }
      // If the distance to the next intersection is smaller than the distance
      // to the next element, move up to the intersection and move on to the
      // next 'distance'.
      if (distances[i*row.num_rays + intersection] < curr_width &&
          distances[i*row.num_rays + intersection] > VOL_FRAC_TOLERANCE &&
          !complete) {
        std::map<int, std::vector<double> >::iterator it =
            row_totals[count].find(eh);

        // If the current cell isn't in the totals yet, add it.
        if (it == row_totals[count].end()){
          row_totals[count].insert(it, std::make_pair(eh, zeros));
        }
        value = distances[i*row.num_rays + intersection]/width[count];
        row_totals[count][eh][0] += value;
        row_totals[count][eh][1] += value*value;
        curr_width -= distances[i*row.num_rays + intersection];
      }
      if (intersection < num_intersections) {
        // This would probably need to be volumes[i*row.num_rays + intersection]
        // if I knew how to return it.
        eh = volumes[intersection];
        intersection++;
      } else {
        complete = true;
        break;
      }
    }
  }
  return row_totals;
}

void startPoints(mesh_row &row, int iter) {
  // The dimensions of the current row.
  double dist_d1 = row.d1div2 - row.d1div1;
  double dist_d2 = row.d2div2 - row.d2div1;
  if (row.grid) {
    int points = sqrt(row.num_rays);
    // Iterate slowly from 1 to points once.
    int iter_d1 = iter/points + 1;
    // Iterate from 1 to points, repeat points times.
    int iter_d2 = iter%points + 1;

    row.start_point_d1 = row.d1div1 + dist_d1/(points + 1)*iter_d1;
    row.start_point_d2 = row.d2div1 + dist_d2/(points + 1)*iter_d2;
  }
  else {
    // Is there a better way to do random numbers between 0 and 1?
    row.start_point_d1 = row.d1div1 + dist_d1*(double)rand()/((double)RAND_MAX);
    row.start_point_d2 = row.d2div1 + dist_d2*(double)rand()/((double)RAND_MAX);
  }
}

EntityHandle find_volume(std::vector<EntityHandle> vol_handles,
                         vec3 pt, vec3 dir) {
  int result = 0;
  ErrorCode rval;

  // Check each volume in the list. This is why it can take so long.
  for (int i = 0; i < vol_handles.size(); i++) {
    void* ptr;
    rval = GQT->point_in_volume(vol_handles[i], pt, result, dir,
        static_cast<const GeomQueryTool::RayHistory*>(ptr));
    CHECKERR_HERE(rval);
    if (result)
      return vol_handles[i];
  }
  std::cerr << "(" << pt[0] << ", " << pt[1] <<", "<< pt[2] << ")" << std::endl;
  throw std::runtime_error("It appears that this point is not in any volume.");
}

std::vector<int> get_idx(int sizes[], int d1, int d2, int d3) {
  std::vector<int> idx;
  // The ids iterate over x, then y, then z. Which one of those is defined by
  // which variable depends on d3.
  if (d3 == 0){
    for (int i = 0; i < sizes[d3]; i++) {
      idx.push_back(d2 + sizes[2]*d1 + sizes[1]*sizes[2]*i);
    }
  }
  else if (d3 == 1){
    for (int i = 0; i < sizes[d3]; i++) {
      idx.push_back(d1 + sizes[2]*i + sizes[1]*sizes[2]*d2);
    }
  }
  else if (d3 == 2){
    for (int i = 0; i < sizes[d3]; i++) {
      idx.push_back(i + sizes[2]*d2 + sizes[1]*sizes[2]*d1);
    }
  }
  return idx;
}

// Dagmc bridge functions

ErrorCode dag_ray_follow(EntityHandle firstvol, vec3 ray_start, vec3 ray_dir,
                         double distance_limit, int* num_intersections,
                         EntityHandle** surfs, double** distances,
                         EntityHandle** volumes, ray_buffers* buf){

  ErrorCode err;

  EntityHandle vol = firstvol;
  double dlimit = distance_limit;
  CartVect ray_point(ray_start);
  EntityHandle next_surf;
  double next_surf_dist;

  CartVect uvw(ray_dir);

  // iterate over the ray until no more intersections are available
  while(vol) {
    err = GQT->ray_fire(vol, ray_point.array(), ray_dir, next_surf,
                        next_surf_dist, &(buf->history), dlimit);
    CHECKERR(err);

    if(next_surf) {
      ray_point += uvw * next_surf_dist;
      buf->surfs.push_back(next_surf);
      buf->dists.push_back(next_surf_dist);
      err = GQT->gttool()->next_vol(next_surf, vol, vol);
      CHECKERR(err);
      buf->vols.push_back(vol);
      if(dlimit != 0){
        dlimit -= next_surf_dist;
      }
    }
    else vol = 0;
  }

  // assign to the output variables
  *num_intersections = buf->surfs.size();
  *surfs = &(buf->surfs[0]);
  *distances = &(buf->dists[0]);
  *volumes = &(buf->vols[0]);

  return err;
}

ErrorCode load_geometry(const char* filename,
                        std::vector<EntityHandle>* vol_handles){
  std::vector<int> volList;

  // Initialize the tools taken from MOAB
  MBI = new Core();
  ErrorCode err = MBI->load_file(filename);
  CHECKERR(err);
  GTT = new GeomTopoTool(MBI);
  GQT = new GeomQueryTool(GTT);
  err = GQT->initialize();
  CHECKERR(err);

  int num_vols = GQT->gttool()->num_ents_of_dim(3);

  moab::Range vols;
  err = GQT->gttool()->get_gsets_by_dimension(3, vols);
  CHECKERR(err);
  vol_handles->resize(num_vols);
  int i = 0;
  for (moab::Range::iterator it = vols.begin(); it != vols.end(); ++it) {
    vol_handles->at(i++) = *it;
  }
  return err;
}
