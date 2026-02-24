#include "mesh.h"

int main(int argc, char* argv[])
{
  int    iTri, iVer;
  double to, ti;

  if (argc < 2) {
    printf(" usage : mesh file \n");
    return 0;
  }

  //--- read a mesh
  to        = clock();
  Mesh* Msh = msh_read(argv[1], 0);
  ti        = clock();

  if (!Msh)
    return 0;

  printf("  Vertices   %10d \n", Msh->NbrVer);
  printf("  Triangles  %10d \n", Msh->NbrTri);
  printf("  time to read the mesh %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  // //--- create neigbhors Q2 version
  // to = clock();
  // msh_neighborsQ2(Msh);
  // ti = clock();
  // printf("  time q2 neigh.        %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  //--- create neigbhors with hash table
  to = clock();
  msh_neighbors(Msh);
  ti = clock();
  printf("  time hash tab neigh.  %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  //--- TODO: compute mesh quality

  printf("Quality evaluation :\n");
  double* Qal = (double*)malloc(sizeof(double) * (Msh->NbrTri + 1));

  double alpha = sqrt(3)/12;
  double K_surf;

  int i1,i2,i3; double a,b,c;

  for (iTri = 1; iTri <= Msh->NbrTri; iTri++) {
  printf("calc %d \n",iTri);
  i1 = Msh->Tri[iTri][0];        i2 = Msh->Tri[iTri][1];        i3 = Msh->Tri[iTri][2];
  double2d P1 = {Msh->Crd[i1][0], Msh->Crd[i1][1] };
  double2d P2 = {Msh->Crd[i2][0], Msh->Crd[i2][1] };
  double2d P3 = {Msh->Crd[i3][0], Msh->Crd[i3][1] };

  a = (P1[0]-P2[0])*(P1[0]-P2[0]) + (P1[1]-P2[1])*(P1[1]-P2[1]); 
  b = (P2[0]-P3[0])*(P2[0]-P3[0]) + (P2[1]-P3[1])*(P2[1]-P3[1]); 
  c = (P3[0]-P1[0])*(P3[0]-P1[0]) + (P3[1]-P1[1])*(P3[1]-P1[1]); 
  double K_surf = fabs(0.5*((P2[0]-P1[0])*(P3[1]-P1[1]) - (P2[1]-P1[1])*(P3[0]-P1[0])));
  printf("donnee : a= %lg, b= %lg, c= %lg et K= %lg \n", a,b,c,K_surf);

  Qal[iTri] = (a+b+c)/K_surf; printf("%d quality evaluated at : %lg \n", iTri, Qal[iTri]);
  }

  msh_write2dfield_Triangles("quality.sol", Msh->NbrTri, Qal);

  printf("quality met \n");

  //--- TODO: compute metric field
  double3d* Met = (double3d*)malloc(sizeof(double3d) * (Msh->NbrVer + 1));

  for (iVer = 1; iVer <= Msh->NbrVer; iVer++) {
    Met[iVer][0] = 1.;
    Met[iVer][1] = 0.;
    Met[iVer][2] = 1.;
  }

  msh_write2dmetric("metric.solb", Msh->NbrVer, Met);

  //--- Free memory
  if (Qal != NULL) {
    free(Qal);
    Qal = NULL;
  }
  if (Met != NULL) {
    free(Met);
    Met = NULL;
  }

  return 0;
}
