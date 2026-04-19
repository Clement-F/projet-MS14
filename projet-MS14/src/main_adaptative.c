#include "mesh.h"

int main(int argc, char* argv[])
{
  int    iTri, iVer;
  double to, ti;

  if (argc < 2) {
    printf(" usage : mesh file \n");
    return 0;
  }
  
  char* file = argv[1];
  char* file_mesh = malloc(1+4+strlen(file)); strcpy(file_mesh,file); strcat(file_mesh,".mesh"); 
  char* file_sol = malloc(1+4+strlen(file)); strcpy(file_sol,file); strcat(file_sol,".sol"); 
  // --- read a mesh
  to        = clock();
  Mesh* Msh = msh_read(file_mesh, 0);
  ti        = clock();
  
  printf("  Vertices   %10d \n", Msh->NbrVer);
  printf("  Triangles  %10d \n", Msh->NbrTri);
  printf("  time to read the mesh %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  
  printf("filename : %s \n",file_sol);
  double* sol =sol_read(file_sol,2,Msh->NbrVer);
  printf("sol recovered \n");

  double** gradu =grad_mesh(Msh,sol);

  printf(" gradient calculated \n");

  double* gradu_x = calloc(Msh->NbrVer,sizeof(double));

  for(int i=1;i<Msh->NbrVer;i++){ gradu_x[i] =gradu[i][0]; }

  
  msh_write2dfield_Vertices("Gradient.sol", Msh->NbrVer, gradu_x);

  return 0;
}