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

  Mesh* Msh_Q = msh_read(argv[1], 1);

  if (!Msh)
    return 0;

  printf("  Vertices   %10d \n", Msh->NbrVer);
  printf("  Triangles  %10d \n", Msh->NbrTri);
  printf("  time to read the mesh %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  //--- create neigbhors Q2 version
  to = clock();
  msh_neighborsQ2(Msh_Q);
  ti = clock();

  printf("  time q2 neigh.        %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  //--- create neigbhors with hash table
  to = clock();
  msh_neighbors(Msh,0);
  ti = clock();

  printf("  time hash tab neigh.  %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);



  printf("Quality evaluation :\n");
  to = clock();
  double* Qal = (double*)malloc(sizeof(double) * (Msh->NbrTri + 1));

  int i1,i2,i3;

  for (iTri = 1; iTri <= Msh->NbrTri; iTri++) {
  // printf("calc %d \n",iTri);
  i1 = Msh->Tri[iTri][0];        i2 = Msh->Tri[iTri][1];        i3 = Msh->Tri[iTri][2];
  double2d P1 = {Msh->Crd[i1][0], Msh->Crd[i1][1] };
  double2d P2 = {Msh->Crd[i2][0], Msh->Crd[i2][1] };
  double2d P3 = {Msh->Crd[i3][0], Msh->Crd[i3][1] };
  Qal[iTri] = quality(P1,P2,P3);
  }

  msh_write2dfield_Triangles("quality.sol", Msh->NbrTri, Qal);
  ti = clock();
  printf("quality evaluated in  %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  
  // msh_write2dmetric("metric.sol", Msh->NbrVer, Met);

  //--- Free memory
  if (Qal != NULL) {
    free(Qal);
    Qal = NULL;
  }
  // if (Met != NULL) {
  //   free(Met);
  //   Met = NULL;
  // }

  return 0;
}
