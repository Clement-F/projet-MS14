#include "mesh.h"

int main(int argc, char* argv[])
{
  if (argc < 2) {
    printf(" usage : mesh file \n");
    return 0;
  }
  
  char* file = argv[1];
  char* file_mesh = malloc(1+4+strlen(file)); strcpy(file_mesh,file); strcat(file_mesh,".mesh"); 
  char* file_sol = malloc(1+4+strlen(file)); strcpy(file_sol,file); strcat(file_sol,".sol"); 

  // --- read a mesh
  Mesh* Msh = msh_read(file_mesh, 0);
    
  printf("filename : %s \n",file_sol);
  double* sol =sol_read(file_sol,2,Msh->NbrVer);

  Image* Raw_image = Imag_init();
  Raw_image->Msh = Msh;  Raw_image->Sol = sol;

  if(sol==NULL) {
    printf(" sol void \n"); 
    return 0;
  };

  Image* Comp_image  = Compression_alea(Raw_image,0.8);
  // Image* Comp_image  = Compression(Raw_image,100);

  double* sol_interpol        = Interpol_sol(Raw_image,Comp_image);

  // output mesh compressed and rough sol
  msh_write2dfield_Vertices("Compression.sol", Comp_image->Msh->NbrVer, Comp_image->Sol);
  msh_write2dfield_Vertices("Compression_interpolate.sol", Raw_image->Msh->NbrVer, sol_interpol);
    
  msh_write(Comp_image->Msh,"Compression.mesh");

  return 0;
}