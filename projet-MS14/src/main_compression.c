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

  
  printf("filename : %s \n",file_sol);
  double* sol =sol_read(file_sol,2,Msh->NbrVer);

  if(sol==NULL){printf(" sol void \n");};
  Compression(Msh,sol,0.5);

  
  return 0;
}