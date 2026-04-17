  #include "mesh.h"
  extern int tri2edg[3][2];
int main(int argc, char* argv[])
{

    // carre_debug
    // test neighbours
    Mesh* Msh       = msh_read("../data/naca0012.mesh", 0);
    Mesh* Msh_ref   = msh_read("../data/naca0012.mesh", 0);
    msh_neighbors(Msh);     // crée la liste de voisins par table de hash
    msh_neighborsQ2(Msh_ref);   // crée la liste de voisins par brute force

    for(int tri_k=0;tri_k<Msh->NbrTri+1; tri_k++)
    {
        // printf(" i ");
        int Vois1 = Msh->TriVoi[tri_k][0]; int Vois2 = Msh->TriVoi[tri_k][1]; int Vois3 = Msh->TriVoi[tri_k][2];
        int Vois1_ref = Msh_ref->TriVoi[tri_k][0]; int Vois2_ref = Msh_ref->TriVoi[tri_k][1]; int Vois3_ref = Msh_ref->TriVoi[tri_k][2];

        if(!((Vois1 == Vois1_ref && Vois2 == Vois2_ref) && Vois3 == Vois3_ref))
        {
        printf(" ERROR NEIGHBOR AT %d \n", tri_k);
        printf(" Neighbors Q : %d, %d, %d \n", Vois1_ref, Vois2_ref, Vois3_ref);
        printf(" Neighbors : %d, %d, %d \n", Vois1, Vois2, Vois3);
        return 0;
        } 
    }
    // test Hash class and functions

    // hash_find
    int i_test=0;

    // limit cases
    i_test = hash_find(Msh->Hsh,0,1);              printf(" for (0,1)      :  i_test = %d \n",i_test);
    i_test = hash_find(Msh->Hsh,1,0);              printf(" for (1,0)      :  i_test = %d \n",i_test);
    i_test = hash_find(Msh->Hsh,Msh->NbrVer,1);    printf(" for (nmax,1)   :  i_test = %d \n",i_test);
    i_test = hash_find(Msh->Hsh,1,Msh->NbrVer);    printf(" for (1,nmax)   :  i_test = %d \n",i_test);

    // all edges test :
    int i1,i2, i_t1,i_t2;
    int A=0,B=0;
    for(int i_tri=1;i_tri<Msh->NbrTri;i_tri++)
    {
        for(int i_edg=0;i_edg<2;i_edg++)
        {
            i1 = Msh->Tri[i_tri][tri2edg[i_edg][0]];
            i2 = Msh->Tri[i_tri][tri2edg[i_edg][1]];
            i_test = hash_find(Msh->Hsh,i1,i2);     

            // vertex test
            i_t1 = Msh->Hsh->LstObj[i_test][0];  
            i_t2 = Msh->Hsh->LstObj[i_test][1];
            if((i1 != i_t1 && i2 != i_t2) && (i1 != i_t2 && i2 != i_t2))
            {
                printf("\n the object %d, linked to %d has wrong points associated with it \n should be : (%d,%d) instead of (%d,%d) \n", i_test,i_tri, i1,i2,i_t1,i_t2);
                hash_out_index(Msh->Hsh,i_test);  
            }

            // triangle test
            i_t1 = Msh->Hsh->LstObj[i_test][2];
            i_t2 = Msh->Hsh->LstObj[i_test][3];
            if(i_tri != i_t1 && i_tri != i_t2) 
            {
                printf("\n the object %d, linked to (%d,%d) has wrong triangle associated with it \n", i_test,i1,i2);
                hash_out_index(Msh->Hsh,i_test);  
            }

            //link two object disconnected
            i_t2 = 0; i_t1=0; A=0; B=0;
            A = Msh->Tri[i_tri][i_edg]; 
            
            if(Msh->Hsh->LstObj[i_test][2] == i_tri) i_t2 = Msh->Hsh->LstObj[i_test][3];
            else{                                    i_t2 = Msh->Hsh->LstObj[i_test][2];}

            if(i_t2 != 0){
            for(int j=0;j<=2;j++){ i_t1 = Msh->Tri[i_t2][j]; if( i_t1 != i1 && i_t1 != i2) B = i_t1;}
            i_test = hash_find(Msh->Hsh,A,B);
            if(i_test != 0)
            {
                printf("\n found a null object or an object that should not exist in position %d \n",i_test);
                printf("T1 : %d (%d,%d,%d), T2: %d (%d,%d,%d) \n",i_tri, Msh->Tri[i_tri][0],Msh->Tri[i_tri][1],Msh->Tri[i_tri][2], 
                        i_t2,Msh->Tri[i_t2][0],Msh->Tri[i_t2][1],Msh->Tri[i_t2][2] );
                printf("linking A %d, to B %d \n", A,B);
            }}

        }
    }
    
    // hash_add has no testing
    // hash_add needs no testing

    // hash_suppr
    int size_hash = Msh->Hsh->NbrObj-1;
    int i_h,iTri, n=0;
    // hash_out(Msh->Hsh);
    while(size_hash>1 && n<1000)
    {
        printf(" \n ===================== \n");
        i_h = rand()%size_hash;
        i1= Msh->Hsh->LstObj[i_h][0]; i2= Msh->Hsh->LstObj[i_h][1]; iTri= Msh->Hsh->LstObj[i_h][2];
        while(i1 == 0 || i2 ==0)
        {
        i_h = rand()%size_hash;
        i1= Msh->Hsh->LstObj[i_h][0]; i2= Msh->Hsh->LstObj[i_h][1]; iTri= Msh->Hsh->LstObj[i_h][2];
        }
        printf("deleting element (%d,%d) linked to %d \n", i1,i2,iTri);
        hash_suppr(Msh->Hsh,i1,i2,iTri);
        // size_hash -=1; 
        n+=1;
        printf(" \n  ==================== \n");
    }



    return 0;
}