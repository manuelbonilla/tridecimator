#include <mvbb_decimator.h>

int main(int argc ,char**argv)
{

  MyMesh mesh;

  MVBBDecimator mesh_decimator;
  mesh_decimator.decimateTriMesh(argv[1], 1000); 

  std::string file_name(argv[1]);
  std::string file_name_no_extension = file_name.substr(0, file_name.size()-4);
  
  
  Eigen::MatrixXd eigen_vertex = mesh_decimator.getEigenVertices();

  vcg::tri::io::ExporterPLY<MyMesh>::Save(mesh_decimator.mesh, (file_name_no_extension + std::string("_reduced.ply")).c_str());

  std::ofstream results_file;
  results_file.open(file_name_no_extension + std::string("_reduced.shl"));
  
  results_file << eigen_vertex.rows() << std::endl;
 

  for (int i=0; i < eigen_vertex.rows(); i++)
  {
      results_file << eigen_vertex(i,0) << " " << eigen_vertex(i,1) << " " << eigen_vertex(i,2) << std::endl;
  }
  
    return 0;

}
