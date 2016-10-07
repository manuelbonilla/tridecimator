#include <mvbb_decimator.h>
#include <algorithm>

MVBBDecimator::MVBBDecimator() : meshreductor(NULL) {}

MVBBDecimator::~MVBBDecimator()
{
    if(meshreductor != NULL)
        delete meshreductor;
}

void MVBBDecimator::decimateTriMesh(string filename,
                                    int targetFacesNumber)
{
    int err=vcg::tri::io::Importer<MyMesh>::Open(mesh, filename.c_str());
    if(err)
    {
      std::cout << "Not possible to open file: " << filename << std::endl;
      exit(-1);
    }

    this->decimate(targetFacesNumber);
}

void MVBBDecimator::decimateTriMesh(const Eigen::MatrixXd &vertices,
                                    const Eigen::MatrixXi &faces,
                                    int targetFacesNumber)
{
    MyMesh::VertexIterator vi = vcg::tri::Allocator<MyMesh>::AddVertices(mesh, vertices.rows());
    MyMesh::FaceIterator fi = vcg::tri::Allocator<MyMesh>::AddFaces(mesh, faces.rows());

    std::vector<MyMesh::VertexPointer> ivp(vertices.rows(), NULL);
    for(unsigned int i = 0; i < vertices.rows(); ++i)
    {
        ivp[i]=&*vi; vi->P()=MyMesh::CoordType ( vertices(i, 0),
                                                 vertices(i, 1),
                                                 vertices(i, 2));
        ++vi;
    }

    for(unsigned int i = 0; i < faces.rows(); ++i)
    {
        fi->V(0)=ivp[faces(i, 0)];
        fi->V(1)=ivp[faces(i, 1)];
        fi->V(2)=ivp[faces(i, 2)];
        ++fi;
    }

    this->decimate(targetFacesNumber);
}

Eigen::MatrixXd MVBBDecimator::getEigenVertices()
{
    return meshreductor->getEigenVertices();
}

Eigen::MatrixXi MVBBDecimator::getEigenFaces()
{
    return meshreductor->getEigenFaces();
}

void MVBBDecimator::decimate(int FinalSize)
{
    meshreductor = new MeshReductor(mesh);
    meshreductor->reduceMesh(FinalSize);
}




int MeshReductor::reduceMesh(int FinalNoFaces)
{

  TriEdgeCollapseQuadricParameter qparams;
  qparams.QualityThr  =.3;
  float TargetError=std::numeric_limits<float>::max();

  int dup = tri::Clean<MyMesh>::RemoveDuplicateVertex(mesh_);
  int unref =  tri::Clean<MyMesh>::RemoveUnreferencedVertex(mesh_);

  vcg::tri::UpdateBounding<MyMesh>::Box(mesh_);

  // decimator initialization
  vcg::LocalOptimization<MyMesh> DeciSession(mesh_,&qparams);


  DeciSession.Init<MyTriEdgeCollapse>();


  DeciSession.SetTargetSimplices(FinalNoFaces);
  DeciSession.SetTimeBudget(0.5f);
  if(TargetError< std::numeric_limits<float>::max() ) DeciSession.SetTargetMetric(TargetError);

  while(DeciSession.DoOptimization() && mesh_.FN() > FinalNoFaces && DeciSession.currMetric < TargetError)
  {
    continue;
  }
  
  tri::Clean<MyMesh>::RemoveUnreferencedVertex(mesh_);
  vcg::tri::Allocator<MyMesh>::CompactFaceVector(mesh_);
  vcg::tri::Allocator<MyMesh>::CompactVertexVector(mesh_);
  
  return 1;
}

Eigen::MatrixXd MeshReductor::getEigenVertices()
{

  int n_vertex= mesh_.VN();
  Eigen::MatrixXd eigen_vertices(n_vertex, 3);
  
   for (int i = 0; i < n_vertex; i++)
   {
     eigen_vertices.block<1,3>(i,0) << mesh_.vert[i].cP()[0] , mesh_.vert[i].cP()[1] , mesh_.vert[i].cP()[2];
   }
  
  return eigen_vertices;
}

Eigen::MatrixXi MeshReductor::getEigenFaces()
{

  int n_face= mesh_.FN();
  int n_vertex= mesh_.VN();

  Eigen::MatrixXi eigen_faces(n_face, 3);
  
  unsigned int i_f = 0;
  for(MyMesh::FaceIterator fi = mesh_.face.begin(); fi!=mesh_.face.end(); ++fi )
  {
    if(fi->IsD()) 
    {
        std::cout << "!" << std::endl;
        continue;
    }
    
    int p_i[3];
    for(unsigned int i_v = 0; i_v < 3; ++i_v)
    {
        p_i[i_v]  = -1;

        for(unsigned int j_v = 0; j_v < n_vertex; ++j_v)
        {
            if(&(mesh_.vert[j_v]) == fi->cV(i_v))
            {
                p_i[i_v] = j_v;
            }
        }

        if(p_i[i_v] < 0) // brute force
        {
            for(unsigned int j_v = 0; j_v < n_vertex; ++j_v)
            {
                    if(mesh_.vert[j_v].P() == fi->cV(i_v)->P())
                    {
                        p_i[i_v] = j_v;
                        std::cout << "!";
                    }
            }
        }
            
        if(p_i[i_v] < 0)
        {
            std::cerr << "Error Finding Vertex " << i_v << " of face " << i_f << " with coords: " << fi->cV(i_v) << std::endl;
            std::cerr << "coords:" << fi->cV(i_v)->P()[0] << " "
                                   << fi->cV(i_v)->P()[1] << " "
                                   << fi->cV(i_v)->P()[2] << std::endl;
            exit(-1);
        }
    }
    eigen_faces.block<1,3>(i_f,0) << p_i[0] , p_i[1] , p_i[2];
        
    ++i_f;
  }
  return eigen_faces;
}
