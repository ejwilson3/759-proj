#include "discretize.h"
#include <cuda.h>
#include <map>
#include <vector>
#include <moab/ErrorCode>

#define CUDA_CHECK(value, label) {              \
   cudaError_t c = (value);                     \
   if (c != cudaSuccess) {                      \
   fprintf(stderr,                              \
     "Error: '%s' at line %d in %s\n",          \
     cudaGetErrorString(c),__LINE__,__FILE__);  \
   goto label;                                  \
   } }

std::vector<double> cuda_rayfire(moab::Instance MBI, moab::GeomTopoTool GTT,
                                 moab::GeomQueryTool GQT, mesh_row row,
                                 std::vector<moab::EntityHandle> vol_handles) {

  std::vector<std::map<int, std::vector<double> > > row_totals;
  std::vector<double> width;
  moab::ErrorCode rval; 
  for (int i = 0; i < row.d3divs.size() - 1; i++) {
    width.push_back(row.d3divs[i+1] - row.d3divs[i]);
  }

  row_totals.resize(width.size());
  moab::EntityHandle root;
  ErrorCode rval = geomTopoTool->get_root(vol_handles[0], root);
// somehow we need to get the OBBTree into an array of doubles.
// We'll call it obbs[]
// Unfortunately, we'll need to hope that the tree is balanced well, because
// we'll only be able to find things by having them ordered; root, level 1,
// level 1, level 2 ...
// Any leaves that aren't there will be left as gaps in the array.
// Each OBB will be taken up by
// center[0], center[1], center[2], length[0], length[1], length[2],
// axes[0], ..., axis[8], radius
  double obbs[size_of_OBBT] = {};
  // The size of this one is tricky; we don't know how many we'll end up hitting
  // so we give it enough space for each of the leaf entitites in the OBBT plus
  // each mesh boundary, assuming only two triangles per OBB.
  double distances[size_of_OBBT/16 + width.size()];
  // Oh shoot, that's how many we need PER RAY! and it's useless if we don't
  // fire a bunch of rays, as that's where the parallelization happens. There
  // is REALLY no memory for this.
  double *d_obbs, *d_distances, *d_width;

  CUDA_CHECK(cudaMalloc(d_obbs, size_of_OBBT*sizeof(double)), cuda_error)
  CUDA_CHECK(cudaMalloc(d_width, width.size()*sizeof(double)), cuda_error)
  CUDA_CHECK(cudaMalloc(d_distances, (size_of_OBBT/16 + width.size())*sizeof(double)), cuda_error)
  
  CUDA_CHECK(cudaMemcpy(d_obbs, obbs, size_of_OBBT*sizeof(double), cudaMemcpyHostToDevice), cuda_error)
  CUDA_CHECK(cudaMemcpy(d_width, &width[0], width.size()*sizeof(double), cudaMemcpyHostToDevice), cuda_error)
  
  N = row.num_rays;
  // We're not using shared memory, so we're going to go with the maximum
  // number of threads per block.
  double blocksPerGrid = (N+1023)/1024;
  double threadsPerBlock = 1024;
  // There also needs to be a way to get the starting points of each ray;
  // didn't get that far.
  bogus_kernel<<<blocksPerGrid,threadsperblock>>>(d_obbs,d_width,d_distances,N,row.grid);

  // Don't actually do this. This would mangle the data.
  CUDA_CHECK(cudaMemcpy(distances, d_distances, (size_of_OBBT/16 + width.size())*sizeof(double), cudaMemcpyDeviceToHost), cuda_error)
  std::vector<int> distances_vec(std::begin(distances), std::end(distances));
  cudaFree(d_distances);
  cudaFree(d_obbs);
  cudaFree(d_width);
  
  return distances_vec
}
