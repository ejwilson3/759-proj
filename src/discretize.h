#ifndef DISCRETIZE_GEOM_H
#define DISCRETIZE_GEOM_H

#include <map>
#include <vector>
#include <moab/GeomQueryTool.hpp>
#include <moab/GeomTopoTool.hpp>
#include <moab/Interface.hpp>
#include <moab/Types.hpp>

// The maxiumum volume fraction to be considered valid

// From dagmc_bridge.h in Pyne
class ray_buffers {

public:
  moab::GeomQueryTool::RayHistory history;
  std::vector<moab::EntityHandle> surfs;
  std::vector<double> dists;
  std::vector<moab::EntityHandle> vols;

};

typedef double vec3[3];

moab::ErrorCode dag_pt_in_vol(moab::EntityHandle vol, vec3 pt, int* result,
                              vec3 dir, const void* history);

moab::ErrorCode dag_ray_follow(moab::EntityHandle firstvol,
                               vec3 ray_start,
                               vec3 ray_dir,
                               double distance_limit,
                               int* num_intersections,
                               moab::EntityHandle** surfs,
                               double** distances,
                               moab::EntityHandle** volumes,
                               ray_buffers* data_buffers);
// End of functions from dagmc_bridge

//Keep together those things which need to be passed to many different functions
struct mesh_row {
  // The number of rays fired each time
  int num_rays;
  // True if evenly spaced ray fires, false if random. If true, num_rays must be
  // a perfect square.
  bool grid;
  // The coordinates of the start of each ray fired.
  double start_point_d1, start_point_d2;
  // The different divisions that define the current row
  double d1div1, d1div2, d2div1, d2div2;
  // The divisions along the current row.
  std::vector<double> d3divs;
  // The index for the direction the row is facing (along which the rays will be
  // fired.
  int d3;
};

/*
 * This function discretizes the mesh attached to a geometry.
 *
 * Parameters:
 *    mesh:             a vector of three vectors, each containing the
 * coordinates of the divisions in one direction of the actual mesh.
 *    vol_handles:      a map of the EntityHandles of the various volumes of the
 * geometry.
 *    num_rays:         The number of rays to be fired.
 *    grid:             Whether the rays are to be fired evenly spaced apart or
 * at random. True if evenly spaced.
 *
 * Returns:
 *    A vector of the results for each row, consisting of the sum of the values
 * given by the fired rays as well as the sum of the values squared.
*/
std::vector<std::vector<double> > discretize_geom(
    std::vector<std::vector<double> > mesh,
    // std::vector<EntityHandle> vol_handles,
    const char* filename,
    int num_rays,
    bool grid);

/*
 * This function fires rays down a single row of the mesh.
 *
 * Parameters:
 *    row:          The various variables which define the current row, as well
 * as how the rays are to be fired.
 *    vol_handles:  a list of the EntityHandles of the various volumes of the
 * geometry.
 *
 * Output:
 *    row:  The "result" member of row will have been changed to reflect the
 * information obtained by firing the rays.
*/
std::vector<std::map<int, std::vector<double> > > fireRays(
    mesh_row &row,
    std::vector<moab::EntityHandle> vol_handles_ids);

/*
 * This function determines the starting coordinates for the next ray. It is
 * called from within fireRays.
 *
 * Parameters:
 *    row:  The various variables which define the current row, as well as how
 * the rays are to be fired.
 *    iter: The number of times this function has been called; only needed if
 * the rays are being fired in a grid; otherwise starting points are chosen
 * randomly.
 *
 * Output:
 *    row:  the "start_point_d1" and "start_point_d2" members of row will have
 * been changed to reflect the new start points.
*/
void startPoints(mesh_row &row, int iter);

/*
 * This function determines the current volume in which a point is.
 * Called from within fireRays.
 *
 * Parameters:
 *    vol_handles: The list of entity handles cooresponding to the volumes.
 *    pt:              The point for which the volume is being found.
 *    dir:             The direction in which we're firing rays.
 *
 * Returns:
 *    eh:              The entity handle for the desired volume.
*/
moab::EntityHandle find_volume(std::vector<moab::EntityHandle> vol_handles_ids, vec3 pt, vec3 dir);

/*
 * This function determines the ids of the volume elements in the current row.
 *
 * Parameters:
 *    sizes:  The number of volume elements in each direction.
 *    d1:     The current volume element count in direction d1.
 *    d2:     The current volume element count in direction d2.
 *    d3:     This value, between 0 and 2, determines in which direction the
 * ray is being fired.
 *
 * Returns:
 *    A vector with the ids of each volume element in the current row.
*/
std::vector<int> get_idx(int sizes[], int d1,
                         int d2, int d3);

moab::ErrorCode load_geometry(const char* filename,
                        std::vector<moab::EntityHandle>* vol_handles);


  /**\brief find the next surface crossing from a given point in a given direction
   *
   * This is the primary method to enable ray tracing through a geometry.
   * Given a volume and a ray, it determines the distance to the nearest intersection
   * with a bounding surface of that volume and returns that distance and the 
   * EntityHandle indicating on which surface that intersection occurs.
   * The caller can compute the location of the intersection by adding the
   * distance to the ray.
   *
   * When a series of calls to this function are made along the same ray (e.g. for
   * the purpose of tracking a ray through several volumes), the optional history
   * argument should be given.  The history prevents previously intersected facets
   * from being intersected again.  A single history should be used as long as a
   * ray is proceeding forward without changing direction.  This situation is
   * sometimes referred to as "streaming."
   *
   * If a ray changes direction at an intersection site, the caller should call
   * reset_to_last_intersection() on the history object before the next ray fire.
   *
   * @param volume The volume at which to fire the ray
   * @param ray_start An array of x,y,z coordinates from which to start the ray.
   * @param ray_dir An array of x,y,z coordinates indicating the direction of the ray.
   *                Must be of unit length.
   * @param next_surf Output parameter indicating the next surface intersected by the ray.
   *                If no intersection is found, will be set to 0.
   * @param next_surf_dist Output parameter indicating distance to next_surf.  If next_surf is
   *                0, this value is undefined and should not be used.
   * @param history Optional RayHistory object.  If provided, the facets in the history are
   *                assumed to not intersect with the given ray.  The facet intersected
   *                by this query will also be added to the history.
   * @param dist_limit Optional distance limit.  If provided and > 0, no intersections at a
   *                distance further than this value will be returned.
   * @param ray_orientation Optional ray orientation. If provided determines intersections
   *                along the normal provided, e.g. if -1 allows intersections back along the
   *                the ray direction, Default is 1, i.e. exit intersections
   * @param stats Optional TrvStats object used to measure performance of underlying OBB
   *              ray-firing query.  See OrientedBoxTreeTool.hpp for details.
   *
   */
   /*
  moab::ErrorCode moab::GeomQueryTool2::ray_fire(
      const moab::EntityHandle volume,
                     const double ray_start[3], const double ray_dir[3],
                     moab::EntityHandle& next_surf, double& next_surf_dist,
                     RayHistory* history, double dist_limit,
                     int ray_orientation,
                     moab::OrientedBoxTreeTool::TrvStats* stats);
                     */
      /*
                     const double ray_start[3], const double ray_dir[3],
                     moab::EntityHandle& next_surf, double& next_surf_dist,
                     RayHistory* history = NULL, double dist_limit = 0,
                     int ray_orientation = 1,
                     moab::OrientedBoxTreeTool::TrvStats* stats = NULL );
                     */
#endif // DISCRETIZE_GEOM_H
