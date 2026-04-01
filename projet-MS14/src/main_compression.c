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

  Image* Raw_image = Imag_init();
  Raw_image->Msh = Msh;
  Raw_image->Sol = sol;

  if(sol==NULL){printf(" sol void \n");};

  to        = clock();
  Image* Comp_image = Compression(Raw_image);
  ti        = clock();

  printf("  Vertices   %10d \n", Raw_image->Msh->NbrVer);
  printf("  Triangles  %10d \n", Raw_image->Msh->NbrTri);
  printf("  time to create the compression %lg (s) \n", (ti - to) / CLOCKS_PER_SEC);

  
  // output mesh compressed and rough sol
  msh_write2dfield_Vertices("Compression.sol", Comp_image->Msh->NbrVer, Comp_image->Sol);
  msh_write(Comp_image->Msh,"Compression.mesh");

  return 0;
}