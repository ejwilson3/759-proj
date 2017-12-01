#include "src/discretize.h"
#include <iostream>

int main(int argc, char* argv[]){
  int num_rays = 10;
  bool grid = false;
  if (argc > 4) {
    throw std::runtime_error("discretize_geom called with too many argments.");
  } else if (argc < 2) {
    throw std::runtime_error("No geometry provided.");
  } else if (argc > 2) {
    if (isdigit(*argv[2]))
      num_rays = atoi(argv[2]);
    if (argc == 4){
      grid = (!strcmp(argv[3], "true") || !strcmp(argv[3], "1"));
    }
  }
  char* filename = argv[1];
  std::vector<std::vector<double> > mesh;
  mesh.resize(3);
  for (int i = 0; i < 3; i++) {
    mesh[i].push_back(0);
    mesh[i].push_back(1);
    // mesh[i].push_back(-4);
    // mesh[i].push_back(-1);
    // mesh[i].push_back(1);
    // mesh[i].push_back(4);
  }
  mesh[0][1] = 2;
  std::vector<std::vector<double> > results = discretize_geom(mesh, filename,
                                                              num_rays, grid);
  for (int i = 0; i < results.size(); i++) {
    std::cout << "[ ";
    for (int j = 0; j < results[i].size(); j++) {
      std::cout << results[i][j] << " ";
    }
    std::cout << "]" << std::endl;
  }
}
