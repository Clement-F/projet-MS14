#include "mesh.h"

int main(int argc, char* argv[])
{
  // int    iTri, iVer;
  // double to, ti;

  // if (argc < 2) {
  //   printf(" usage : mesh file \n");
  //   return 0;
  // }

  //--- read a mesh
  // to        = clock();
  // Mesh* Msh = msh_read(argv[1], 0);
  // ti        = clock();

  // Mesh* Msh_Q = msh_read(argv[1], 1);

  // if (!Msh)
  //   return 0;

  // printf("  Vertices   %10d \n", Msh->NbrVer);
  // printf("  Triangles  %10d \n", Msh->NbrTri);
  // printf("  time to read the mesh %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  // //--- create neigbhors Q2 version
  // // to = clock();
  // // msh_neighborsQ2(Msh_Q);
  // // ti = clock();

  // printf("  time q2 neigh.        %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  // //--- create neigbhors with hash table
  // to = clock();
  // msh_neighbors(Msh);
  // ti = clock();

  // printf("  time hash tab neigh.  %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  // --- check solution
  // for(int tri_k=0;tri_k<Msh->NbrTri+1; tri_k++)
  // {
  //   int Vois1 = Msh->TriVoi[tri_k][0]; int Vois2 = Msh->TriVoi[tri_k][1]; int Vois3 = Msh->TriVoi[tri_k][2];
  //   int Vois1_Q = Msh_Q->TriVoi[tri_k][0]; int Vois2_Q = Msh_Q->TriVoi[tri_k][1]; int Vois3_Q = Msh_Q->TriVoi[tri_k][2];

  //   if(!((Vois1 == Vois1_Q && Vois2 == Vois2_Q) && Vois3 == Vois3_Q))
  //   {
  //     printf(" ERROR NEIGHBOR AT %d \n", tri_k);
  //     printf(" Neighbors Q : %d, %d, %d \n", Vois1_Q, Vois2_Q, Vois3_Q);
  //     printf(" Neighbors : %d, %d, %d \n", Vois1, Vois2, Vois3);
  //   } 

  // }


  // printf("Quality evaluation :\n");
  // to = clock();
  // double* Qal = (double*)malloc(sizeof(double) * (Msh->NbrTri + 1));

  // int i1,i2,i3;

  // for (iTri = 1; iTri <= Msh->NbrTri; iTri++) {
  // // printf("calc %d \n",iTri);
  // i1 = Msh->Tri[iTri][0];        i2 = Msh->Tri[iTri][1];        i3 = Msh->Tri[iTri][2];
  // double2d P1 = {Msh->Crd[i1][0], Msh->Crd[i1][1] };
  // double2d P2 = {Msh->Crd[i2][0], Msh->Crd[i2][1] };
  // double2d P3 = {Msh->Crd[i3][0], Msh->Crd[i3][1] };
  // Qal[iTri] = quality(P1,P2,P3);
  // }

  // msh_write2dfield_Triangles("quality.sol", Msh->NbrTri, Qal);
  // ti = clock();
  // printf("quality evaluated in  %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  
  //================================================================================================
  //=================================         TP2        ===========================================
  //================================================================================================
  Maillage_Delauney(4);


  //================================================================================================
  //=============================         Exercice 3        ========================================
  //================================================================================================

  // double* color = calloc((Msh->NbrTri+1),sizeof(double));

  // color = (connex_comp(Msh));
  // // for(int i=0; i<Msh->NbrTri+1; i++){printf("color of element %d : %f \n",i,color[i]);}
  // msh_write2dfield_Triangles("color.sol", Msh->NbrTri, color);

  //--- TODO: compute metric field
  // double3d* Met = (double3d*)malloc(sizeof(double3d) * (Msh->NbrVer + 1));

  // for (iVer = 1; iVer <= Msh->NbrVer; iVer++) {
  //   Met[iVer][0] = 1.;
  //   Met[iVer][1] = 0.;
  //   Met[iVer][2] = 1.;
  // }

  // msh_write2dmetric("metric.sol", Msh->NbrVer, Met);

  // //--- Free memory
  // if (Qal != NULL) {
  //   free(Qal);
  //   Qal = NULL;
  // }
  // if (Met != NULL) {
  //   free(Met);
  //   Met = NULL;
  // }
  // if (color !=NULL){
  //   free(color);
  //   color = NULL;
  // }

  return 0;
}
