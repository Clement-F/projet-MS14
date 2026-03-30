#include "mesh.h"

int main(int argc, char* argv[])
{
  int    iTri, iVer;
  double to, ti;

  if (argc < 2) {
    printf(" usage : mesh file \n");
    return 0;
  }

  // --- read a mesh
  to        = clock();
  Mesh* Msh = msh_read(argv[1], 1);
  ti        = clock();

  double* color = calloc((Msh->NbrTri+1),sizeof(double));

  color = (connex_comp(Msh));
  // for(int i=0; i<Msh->NbrTri+1; i++){printf("color of element %d : %f \n",i,color[i]);}
  msh_write2dfield_Triangles("color.sol", Msh->NbrTri, color);

//   --- TODO: compute metric field
  double3d* Met = (double3d*)malloc(sizeof(double3d) * (Msh->NbrVer + 1));

  for (iVer = 1; iVer <= Msh->NbrVer; iVer++) {
    Met[iVer][0] = 1.;
    Met[iVer][1] = 0.;
    Met[iVer][2] = 1.;
  }
  
  if (color !=NULL){
    free(color);
    color = NULL;
  }

  return 0;
}