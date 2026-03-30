  #include "mesh.h"

int main(int argc, char* argv[])
{
  int    iTri;
  double to, ti;

  if (argc < 2) {
    printf(" usage : mesh file \n");
    return 0;
  }
  to        = clock();
  Mesh* Msh = msh_read(argv[1], 0);
  ti        = clock();

<<<<<<< HEAD
  Msh = Maillage_Delaunay(500,Msh);
=======
  Msh = Maillage_Delaunay(1000,Msh);
>>>>>>> tryout
  msh_neighbors(Msh);
  
  printf(" Mesh created \n");

  
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


    //--- Free memory
  if (Qal != NULL) {
    free(Qal);
    Qal = NULL;
  }
  if(Msh != NULL){
    free(Msh);
    Msh = NULL;
  }

  return 0;
}