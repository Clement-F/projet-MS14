#include "mesh.h"



int tri2edg[3][2] = { { 1, 2 }, { 2, 0 }, { 0, 1 } };

Mesh* msh_init()
{
  Mesh* Msh = malloc(sizeof(Mesh));
  if (!Msh) return NULL;

  Msh->Dim    = 0;
  Msh->NbrVer = 0;
  Msh->NbrTri = 0;
  Msh->NbrEfr = 0;
  Msh->NbrEdg = 0;

  Msh->NbrVerMax = 0;
  Msh->NbrTriMax = 0;
  Msh->NbrEfrMax = 0;
  Msh->NbrEdgMax = 0;

  Msh->Box[0] = 1.e30; // xmin
  Msh->Box[1] = -1.e30; // xmax
  Msh->Box[2] = 1.e30; // ymin
  Msh->Box[3] = -1.e30; // ymax

  //--- Data for the list of vertices
  Msh->Crd = NULL;

  //--- Data for the list of triangles
  Msh->Tri    = NULL;
  Msh->TriVoi = NULL;
  Msh->TriRef = NULL;
  Msh->TriMrk = NULL;

  //--- Data for the list of boundary edges
  Msh->Efr    = NULL;
  Msh->EfrVoi = NULL;
  Msh->EfrRef = NULL;

  //--- Data for the list of edges
  Msh->Edg = NULL;                    // <- serves nothing

  return Msh;
}

Mesh* msh_read(char* file, int readEfr)
{
  char   InpFil[1024];
  float  bufFlt[2];
  double bufDbl[2];
  int    i, bufTri[4], bufEfr[3];
  int    FilVer, ref;

  int fmsh = 0;

  if (!file) return NULL;

  Mesh* Msh = msh_init();

  //--- set file name
  strcpy(InpFil, file);
  if (strstr(InpFil, ".mesh")) {
    if (!(fmsh = GmfOpenMesh(InpFil, GmfRead, &FilVer, &Msh->Dim))) {
      return NULL;
    }
  }
  else {
    strcat(InpFil, ".meshb");
    if (!(fmsh = GmfOpenMesh(InpFil, GmfRead, &FilVer, &Msh->Dim))) {
      strcpy(InpFil, file);
      strcat(InpFil, ".mesh");
      if (!(fmsh = GmfOpenMesh(InpFil, GmfRead, &FilVer, &Msh->Dim))) {
        return NULL;
      }
    }
  }

  printf(" File %s opened Dimension %d Version %d \n", InpFil, Msh->Dim, FilVer);

  Msh->NbrVer = GmfStatKwd(fmsh, GmfVertices);
  Msh->NbrTri = GmfStatKwd(fmsh, GmfTriangles);

  Msh->NbrVerMax = Msh->NbrVer;
  Msh->NbrTriMax = Msh->NbrTri;

  //--- allocate arrays
  Msh->Crd    = calloc((Msh->NbrVerMax + 1), sizeof(double3d));
  Msh->Tri    = calloc((Msh->NbrTriMax + 1), sizeof(int3d));
  Msh->TriRef = calloc((Msh->NbrTriMax + 1), sizeof(int1d));
  Msh->TriMrk = calloc((Msh->NbrTriMax + 1), sizeof(int1d));

  //--- read vertices
  GmfGotoKwd(fmsh, GmfVertices);
  if (Msh->Dim == 2) {
    if (FilVer == GmfFloat) { // read 32 bits float
      for (i = 1; i <= Msh->NbrVer; ++i) {
        GmfGetLin(fmsh, GmfVertices, &bufFlt[0], &bufFlt[1], &ref);
        Msh->Crd[i][0] = (double)bufFlt[0];
        Msh->Crd[i][1] = (double)bufFlt[1];
      }
    }
    else { // read 64 bits float
      for (i = 1; i <= Msh->NbrVer; ++i) {
        GmfGetLin(fmsh, GmfVertices, &bufDbl[0], &bufDbl[1], &ref);
        Msh->Crd[i][0] = bufDbl[0];
        Msh->Crd[i][1] = bufDbl[1];
      }
    }
  }
  else {
    fprintf(stderr, "  ## ERROR: 3D is not implemented\n");
    exit(1);
  }

  //--- read triangles
  GmfGotoKwd(fmsh, GmfTriangles);
  for (i = 1; i <= Msh->NbrTri; ++i) {
    GmfGetLin(fmsh, GmfTriangles, &bufTri[0], &bufTri[1], &bufTri[2], &bufTri[3]);
    Msh->Tri[i][0] = bufTri[0];
    Msh->Tri[i][1] = bufTri[1];
    Msh->Tri[i][2] = bufTri[2];
    Msh->TriRef[i] = bufTri[3];
  }

  //--- read boundary edges
  if (readEfr == 1) {
    Msh->NbrEfr    = GmfStatKwd(fmsh, GmfEdges);
    Msh->NbrEfrMax = Msh->NbrEfr;

    Msh->Efr    = calloc((Msh->NbrEfrMax + 1), sizeof(int2d));
    Msh->EfrRef = calloc((Msh->NbrEfrMax + 1), sizeof(int1d));

    GmfGotoKwd(fmsh, GmfEdges);
    for (i = 1; i <= Msh->NbrEfr; ++i) {
      GmfGetLin(fmsh, GmfEdges, &bufEfr[0], &bufEfr[1], &bufEfr[2]);
      Msh->Efr[i][0] = bufEfr[0];
      Msh->Efr[i][1] = bufEfr[1];
      Msh->EfrRef[i] = bufEfr[2];
    }
  }

  GmfCloseMesh(fmsh);

  return Msh;
}

double* sol_read(char* file, int mshDim, int mshNbrSol)
{
  char   InpFil[1024];
  int    FilVer, SolTyp, NbrTyp, SolSiz, TypTab[GmfMaxTyp];
  float  bufFlt;
  double bufDbl;
  int    i, dim, nbrSol;

  int fsol = 0;

  if (!file) return NULL;

  double* sol = NULL;

  //--- set file name
  strcpy(InpFil, file);
  if (strstr(InpFil, ".sol")) {
    if (!(fsol = GmfOpenMesh(InpFil, GmfRead, &FilVer, &dim))) {
      return NULL;
    }
  }
  else {
    strcat(InpFil, ".solb");
    if (!(fsol = GmfOpenMesh(InpFil, GmfRead, &FilVer, &dim))) {
      strcpy(InpFil, file);
      strcat(InpFil, ".sol");
      if (!(fsol = GmfOpenMesh(InpFil, GmfRead, &FilVer, &dim))) {
        return NULL;
      }
    }
  }

  printf(" File %s opened Dimension %d Version %d \n", InpFil, dim, FilVer);

  SolTyp = GmfSolAtVertices; // read only sol at vertices
  nbrSol = GmfStatKwd(fsol, SolTyp, &NbrTyp, &SolSiz, TypTab);

  if (nbrSol == 0) {
    printf("  ## WARNING: No SolAtVertices in the solution file !\n");
    return NULL;
  }
  if (dim != mshDim) {
    printf("  ## WARNING: WRONG DIMENSION NUMBER. IGNORED\n");
    return NULL;
  }
  if (nbrSol != mshNbrSol) {
    printf("  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    return NULL;
  }
  if (NbrTyp != 1) {
    printf("  ## WARNING: WRONG FIELD NUMBER. IGNORED\n");
    return NULL;
  }
  if (TypTab[0] != GmfSca) {
    printf("  ## WARNING: WRONG FIELD TYPE. IGNORED\n");
    return NULL;
  }

  sol = (double*)calloc(nbrSol + 1, sizeof(double));

  GmfGotoKwd(fsol, SolTyp);

  for (i = 1; i <= nbrSol; ++i) {
    if (FilVer == GmfFloat) {
      GmfGetLin(fsol, SolTyp, &bufFlt);
      sol[i] = (double)bufFlt;
    }
    else {
      GmfGetLin(fsol, SolTyp, &bufDbl);
      sol[i] = bufDbl;
    }
  }

  if (!GmfCloseMesh(fsol)) {
    fprintf(stderr, "  ## ERROR: Cannot close solution file %s ! \n", InpFil);
    // myexit(1);
  }

  return sol;
}

int msh_boundingbox(Mesh* Msh)
{
  int1d iVer;
  double x_min =Msh->Crd[1][0], x_max =Msh->Crd[1][0];
  double y_min =Msh->Crd[1][1], y_max =Msh->Crd[1][1];
  //--- compute bounding box
  for (iVer = 2; iVer <= Msh->NbrVer; iVer++) {
    if(Msh->Crd[iVer][0]>x_max) x_max = Msh->Crd[iVer][0];
    if(Msh->Crd[iVer][0]<x_min) x_min = Msh->Crd[iVer][0];
    if(Msh->Crd[iVer][1]>y_max) y_max = Msh->Crd[iVer][1];
    if(Msh->Crd[iVer][1]<y_min) y_min = Msh->Crd[iVer][1];
  }
  Msh->Box[0]=x_min; Msh->Box[1]=x_max; Msh->Box[2]=y_min; Msh->Box[3]=y_max;
  return 1;
}

int msh_write(Mesh* Msh, char* file)
{
  int iVer, iTri, iEfr;
  int FilVer = 2;

  if (!Msh) return 0;
  if (!file) return 0;

  printf(" e ");
  int fmsh = GmfOpenMesh(file, GmfWrite, FilVer, Msh->Dim);
  printf(" p ");
  if (fmsh <= 0) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }

  printf(" a ");
  GmfSetKwd(fmsh, GmfVertices, Msh->NbrVer);
  for (iVer = 1; iVer <= Msh->NbrVer; iVer++)
    GmfSetLin(fmsh, GmfVertices, Msh->Crd[iVer][0], Msh->Crd[iVer][1], 0);

  printf(" b ");
  GmfSetKwd(fmsh, GmfTriangles, Msh->NbrTri);
  for (iTri = 1; iTri <= Msh->NbrTri; iTri++)
    GmfSetLin(fmsh, GmfTriangles, Msh->Tri[iTri][0], Msh->Tri[iTri][1], Msh->Tri[iTri][2], Msh->TriRef[iTri]);

  printf(" c ");
  if (Msh->NbrEfr > 0) {
    GmfSetKwd(fmsh, GmfEdges, Msh->NbrEfr);
    for (iEfr = 1; iEfr <= Msh->NbrEfr; iEfr++)
      GmfSetLin(fmsh, GmfEdges, Msh->Efr[iEfr][0], Msh->Efr[iEfr][1], Msh->EfrRef[iEfr]);
  }
  printf(" d ");
  GmfCloseMesh(fmsh);

  return 1;
}

int msh_neighborsQ2(Mesh* Msh)
{
  int iTri, iEdg, jTri, jEdg, iVer1, iVer2, jVer1, jVer2;
  double to, ti;

  if (!Msh) return 0;

  if (Msh->TriVoi == NULL)
    Msh->TriVoi = calloc((Msh->NbrTri + 1), sizeof(int3d));

  to = clock();  
  //--- Compute the neighbors using a quadratic-complexity algorithm
  for (iTri = 1; iTri <= Msh->NbrTri; iTri++) {
    if(iTri%5000 == 0){ ti = clock();   printf("--- task %d / %d full --- %lg (s) passed \n",iTri,Msh->NbrTri,(ti-to)/ CLOCKS_PER_SEC );}
    for (iEdg = 0; iEdg < 3; iEdg++) {
      if(Msh->TriVoi[iTri][iEdg] !=0){ continue;}
      iVer1 = Msh->Tri[iTri][tri2edg[iEdg][0]];
      iVer2 = Msh->Tri[iTri][tri2edg[iEdg][1]];

      //--- find the Tri different from iTri that has iVer1, iVer2 as vertices
      for (jTri = 1; jTri <= Msh->NbrTri; jTri++) {
        if(iTri == jTri){continue;}

        for (jEdg = 0; jEdg < 3; jEdg++) {
          jVer1 = Msh->Tri[jTri][tri2edg[jEdg][0]];
          jVer2 = Msh->Tri[jTri][tri2edg[jEdg][1]];


          if(((iVer1== jVer1) && (iVer2 ==jVer2)) || ((iVer1 == jVer2) && (iVer2 == jVer1)))
          {
            Msh->TriVoi[jTri][jEdg]=iTri;
            Msh->TriVoi[iTri][iEdg]=jTri;
          }

        }
      }
    }
  }
  
  // for (iTri = 1; iTri <= Msh->NbrTri; iTri++) {
  // printf("Trivoi of %d : %d, %d, %d \n", iTri, Msh->TriVoi[iTri][0], Msh->TriVoi[iTri][1], Msh->TriVoi[iTri][2]);
  // }

  return 1;
}

// ============================================================================
// ============================================================================

int msh_neighbors(Mesh* Msh)
{
  printf(" init neighbors\n");
  int iTri, iEdg, iVer1, iVer2, jVer1, jVer2, jTri, jEdg;

  if (!Msh) return 0;

  if (Msh->TriVoi == NULL)
    Msh->TriVoi = calloc((Msh->NbrTri + 1), sizeof(int3d));

  //--- initialize HashTable and set the hash table
  printf(" init hash\n");
  int SizHead = 2*(Msh->NbrVerMax);
  int NbrMaxObj = Msh->NbrVerMax + Msh->NbrTriMax ; // Euler caracteristique with a bit of security


  HashTable* hsh = hash_init(SizHead, NbrMaxObj); 
  for (iTri = 1; iTri <= Msh->NbrTri; iTri++) {
    for (iEdg = 0; iEdg < 3; iEdg++) {
      iVer1 = Msh->Tri[iTri][tri2edg[iEdg][0]];
      iVer2 = Msh->Tri[iTri][tri2edg[iEdg][1]];

      int j_hsh = hash_find(hsh,iVer1,iVer2);
      if(j_hsh !=0)
      {
        hash_add(hsh,iVer1,iVer2,iTri,j_hsh);
        jTri = hsh->LstObj[j_hsh][2];
        Msh->TriVoi[iTri][iEdg] = jTri;
        for(jEdg=0;jEdg<3;jEdg++)
        {
          jVer1 = Msh->Tri[jTri][tri2edg[jEdg][0]];
          jVer2 = Msh->Tri[jTri][tri2edg[jEdg][1]];
           if((jVer1 == iVer1 && jVer2 == iVer2) || (jVer1 == iVer2 && jVer2 == iVer1))
           {Msh->TriVoi[jTri][jEdg] = iTri;}
        }
      }
      if(j_hsh ==0)
      {
        hash_add(hsh,iVer1,iVer2,iTri,0);
      }

    }
  }

  Msh->Hsh = hsh;
  return 1;
}

HashTable* hash_init(int SizHead, int NbrMaxObj)
{
  // HashTable* hsh = NULL;
  HashTable* hsh = malloc(sizeof(HashTable));
  printf("init hash table \n");

  hsh->SizHead = SizHead; 
  hsh->NbrMaxObj = NbrMaxObj;
  hsh->NbrObj =0;
  
  printf("alloc hash table \n");

  hsh->LstObj = calloc((hsh->NbrMaxObj+1),sizeof(int5d));
  hsh->Head   = calloc((hsh->SizHead  +1),sizeof(int));

  printf("size of head :%10d \n", SizHead);
  printf("size of obj  :%10d \n", NbrMaxObj);
  for(int j=1;j<hsh->SizHead;j++){hsh->Head[j]=0;}
  printf("init finalized \n");

  return hsh;
}

int hash_find(HashTable* hsh, int iVer1, int iVer2)
{
  int key = iVer1 + iVer2;
  int j_hsh = hsh->Head[key];
  int n = 0; // security
  while(j_hsh!=0 && n<hsh->NbrMaxObj)
  {
    // printf("searching element %d if it has Vert (%d,%d) \n",j_hsh,iVer1,iVer2);
    if(hsh->LstObj[j_hsh][0]==iVer1 && hsh->LstObj[j_hsh][1]==iVer2){return j_hsh;}
    else{ if(hsh->LstObj[j_hsh][1]==iVer1 && hsh->LstObj[j_hsh][0]==iVer2){return j_hsh;}
          else {j_hsh = hsh->LstObj[j_hsh][4];}
    }
  }
  // printf("j_hsh = 0 \n");
  return 0;

}

int hash_add(HashTable* hsh, int iVer1, int iVer2, int iTri, int i_hsh)
{

  // i_hsh is a initial guess. it should be 0 by default, 
  // but if you have already run hash_find, you can bypass the find to imput the element at i_hsh
  if(i_hsh==0)  i_hsh = hash_find(hsh,iVer1,iVer2); // check if the element isn't in the hash_list already

  if(i_hsh==0)
  {
    if(hsh->NbrObj> hsh->NbrMaxObj+1) printf("  ## WARNING: HSH ELEMENT ALREADY FULL. IGNORED\n");
    // printf("adding element %d of Vertex (%d,%d) and Tri (%d,%d) to have next %d \n", hsh->NbrObj, iVer1,iVer2, iTri, hsh->LstObj[hsh->NbrObj][3], hsh->LstObj[hsh->NbrObj][4]);

    hsh->LstObj[hsh->NbrObj +1][0] = iVer1;
    hsh->LstObj[hsh->NbrObj +1][1] = iVer2;
    hsh->LstObj[hsh->NbrObj +1][2] = iTri; 
    hsh->LstObj[hsh->NbrObj +1][3] = 0;
    hsh->LstObj[hsh->NbrObj +1][4] = hsh->Head[iVer1+iVer2];  
    hsh->NbrObj +=1;
    hsh->Head[iVer1+iVer2] = hsh->NbrObj;   
    
  }

  if(i_hsh !=0)
  { 
    if(hsh->LstObj[i_hsh][3] ==0){hsh->LstObj[i_hsh][3] = iTri;} 
    else printf(" ## WARNING: HSH ELEMENT %d ALREADY COMPLETE (%d,%d); (%d,%d) -> %d. \n IGNORED\n",i_hsh,hsh->LstObj[i_hsh][0],hsh->LstObj[i_hsh][1],hsh->LstObj[i_hsh][2],hsh->LstObj[i_hsh][3],hsh->LstObj[i_hsh][4] );
    // printf("updating element %d of Vertex (%d,%d) and Tri (%d,%d) to have next %d \n", i_hsh, iVer1,iVer2, hsh->LstObj[hsh->NbrObj][3], iTri, hsh->LstObj[hsh->NbrObj][4]);
  }

  return 0;
}

// NEED TESTING
int hash_suppr(HashTable* hsh, int iVer1, int iVer2, int iTri)  // deletes an element of the hash table
{
  // printf(" \n ---------------------------------- \n ");
  int i_hsh;
  i_hsh = hash_find(hsh,iVer1,iVer2); // check if the element is in the hash_list

  if(i_hsh ==0){ printf("\n DELETING A NON EXISTING ELEMENT, IGNORED \n ");}
  if(i_hsh !=0)
  {
    // printf("deleting the edge %d :(%d,%d) of %d \n",i_hsh,iVer1,iVer2,iTri);
    int ToDelete=0;
    if((hsh->LstObj[i_hsh][2]!=0 && hsh->LstObj[i_hsh][3]==0) || (hsh->LstObj[i_hsh][3]!=0 && hsh->LstObj[i_hsh][2]==0)) ToDelete = 1;
    if(ToDelete==1)
    {
      // printf("deleting the object %d \n",i_hsh);
    // we redo the chain
    int key = iVer1 + iVer2; 
    int j_hsh = hsh->Head[key];
    int i_bef = 0;
    while(j_hsh!=0)
    {
      // printf("object %d in the chain of key : %d \n",j_hsh,key);
      // hash_out_index(hsh, j_hsh);
      if((hsh->LstObj[j_hsh][0]==iVer1 && hsh->LstObj[j_hsh][1]==iVer2) || (hsh->LstObj[j_hsh][1]==iVer1 && hsh->LstObj[j_hsh][0]==iVer2) ) 
      {
        // sewing it
        // printf("found the object %d in the chain \n",i_hsh);
        // printf("sewing the chain \n");
        if(i_bef==0){hsh->Head[key]       =hsh->LstObj[j_hsh][4];} 
        if(i_bef!=0){hsh->LstObj[i_bef][4]=hsh->LstObj[j_hsh][4];}

        if(j_hsh<hsh->NbrObj)
        {
          // permuting it with the last element
          int last_index = hsh->NbrObj;        
          // printf("permuting the position of the deleted object %d with the last element %d \n", j_hsh,last_index);
          for(int i=0;i<5;i++) hsh->LstObj[j_hsh][i] = hsh->LstObj[last_index][i];

          // printf("sewing the chain of the last element \n");
          // sewing the last element's chain 
          key = hsh->LstObj[last_index][0] + hsh->LstObj[last_index][1];
          int k_hsh= hsh->Head[key];
          
          int j_bef = 0;
          while(k_hsh!=0)
          {
            if(hsh->LstObj[k_hsh][0]==hsh->LstObj[j_hsh][0] && hsh->LstObj[k_hsh][1]==hsh->LstObj[j_hsh][1])
            {
              if(j_bef==0){hsh->Head[key]       =hsh->LstObj[k_hsh][4];} 
              if(j_bef!=0){hsh->LstObj[j_bef][4]=hsh->LstObj[k_hsh][4];}
              for(int i=0;i<5;i++) hsh->LstObj[last_index][i]=0;

            }
            else k_hsh = hsh->LstObj[j_bef][4];
          }
          hsh->NbrObj -=1;
        }
        
      }
      else j_hsh = hsh->LstObj[j_hsh][4];      
    } 
    }

    // we put an element to 0 and potentially put the element in a "good form"
    if(ToDelete==0)
    {
      // printf(" the element %d is decomposed as : \n",i_hsh);
      // for(int i=0;i<5;i++) printf(" %d ",hsh->LstObj[i_hsh][i]);
      // printf("\n");
      if(hsh->LstObj[i_hsh][3]==iTri) hsh->LstObj[i_hsh][3]=0;
      if(hsh->LstObj[i_hsh][2]==iTri) 
      {
        hsh->LstObj[i_hsh][2]= hsh->LstObj[i_hsh][3];
        hsh->LstObj[i_hsh][3]= 0;
      }  
      // printf(" the element %d is decomposed as : \n",i_hsh);
      // for(int i=0;i<5;i++) printf(" %d ",hsh->LstObj[i_hsh][i]);
      // printf("\n");
    }
  }
  
  // printf(" \n ---------------------------------- \n ");
  return 0;
}

int hash_out(HashTable* hsh)  // print the hash table
{
  for(int j=1;j<hsh->NbrMaxObj;j++){
    printf(" indices = %d ", j);
    printf(" Vertexes  : %d, %d ", hsh->LstObj[j][0], hsh->LstObj[j][1]);
    printf(" Triangles : %d, %d ", hsh->LstObj[j][2], hsh->LstObj[j][3]);
    printf(" next : %d\n", hsh->LstObj[j][4]);
    printf("----------- \n");
  }
  printf(" number of edge : %d / %d", hsh->NbrObj, hsh->NbrMaxObj);
  return 0;
}


int hash_out_index(HashTable* hsh, int index)  // print the element of the table
{

  printf(" indices = %d ", index);
  printf(" Vertexes  : %d, %d ", hsh->LstObj[index][0], hsh->LstObj[index][1]);
  printf(" Triangles : %d, %d ", hsh->LstObj[index][2], hsh->LstObj[index][3]);
  printf(" next : %d\n", hsh->LstObj[index][4]);
  printf("----------- \n");
  
  return 0;
}

int hash_bound(HashTable* hsh)
{
  int Nb_bound =0;

  for(int i_hsh=1; i_hsh <= hsh->NbrObj; i_hsh++ )
  {  
    // printf("Triangles : %d, %d ", hsh->LstObj[i_hsh][2], hsh->LstObj[i_hsh][3]);
    if(hsh->LstObj[i_hsh][3] == 0)
    {      
      // printf(" : edge ");
      Nb_bound +=1;
    }
    // printf(" \n ======= \n");
  }
  printf("Number of boundary edges : %d \n", Nb_bound);
  return 0;
}

int hash_collision(HashTable* hsh)
{
  int Max_col = 0,  col=0, length=0;
  double Mean_col = 0;
  int i_hsh;
  for(int i_key = 1; i_key<=hsh->SizHead; i_key++)
  {
    length =0;
    i_hsh = hsh->Head[i_key];
    if(i_hsh!=0){col+=1;}

    while(i_hsh!=0)
    {
      length+=1;
      i_hsh = hsh->LstObj[i_hsh][4];
    }

    Mean_col += length;
    if(length>Max_col) Max_col = length;

  }
  printf("number of collision : %d \n", col);
  printf("Maximum collision : %d \n", Max_col);
  printf("Mean collision : %lg \n",  Mean_col/col );
  return 0;
}

// ============================================================================
// ============================================================================

int msh_write2dfield_Vertices(char* file, int nfield, double* field)
{
  int iVer;

  int fmsh = GmfOpenMesh(file, GmfWrite, GmfDouble, 2);
  if (fmsh <= 0) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }

  int sizfld[1];
  sizfld[0] = GmfSca;

  GmfSetKwd(fmsh, GmfSolAtVertices, nfield, 1, sizfld);

  for (iVer = 1; iVer <= nfield; iVer++)
    GmfSetLin(fmsh, GmfSolAtVertices, &field[iVer]);

  GmfCloseMesh(fmsh);

  return 1;
}

int msh_write2dfield_Triangles(char* file, int nfield, double* field)
{
  int iTri;

  int fmsh = GmfOpenMesh(file, GmfWrite, GmfDouble, 2);
  if (fmsh <= 0) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }

  int sizfld[1];
  sizfld[0] = GmfSca;

  GmfSetKwd(fmsh, GmfSolAtTriangles, nfield, 1, sizfld);

  for (iTri = 1; iTri <= nfield; iTri++)
    GmfSetLin(fmsh, GmfSolAtTriangles, &field[iTri]);

  GmfCloseMesh(fmsh);

  return 1;
}

int msh_write2dmetric(char* file, int nmetric, double3d* metric)
{
  int iVer;

  int fmsh = GmfOpenMesh(file, GmfWrite, GmfDouble, 2);
  if (fmsh <= 0) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }

  int sizfld[1];
  sizfld[0] = GmfSymMat;

  GmfSetKwd(fmsh, GmfSolAtVertices, nmetric, 1, sizfld);

  for (iVer = 1; iVer <= nmetric; iVer++)
    GmfSetLin(fmsh, GmfSolAtVertices, &metric[iVer][0], &metric[iVer][1], &metric[iVer][2]);

  GmfCloseMesh(fmsh);

  return 1;
}

// ============================================================================
// ============================================================================

int valid_edge(Mesh* Msh, int iTri, int iEdg)
{
  int is_valid =1;
  int jEdg,jVer1,jVer2;
  int iVer1, iVer2;

  iVer1 = Msh->Tri[iTri][tri2edg[iEdg][0]];
  iVer2 = Msh->Tri[iTri][tri2edg[iEdg][1]];

  for (jEdg = 0; jEdg < Msh->NbrEfrMax; jEdg++) {
    jVer1 = Msh->Efr[jEdg][0];
    jVer2 = Msh->Efr[jEdg][1];

    if(((iVer1== jVer1) && (iVer2 ==jVer2)) || ((iVer1 == jVer2) && (iVer2 == jVer1))) is_valid = 0;
  }

  return is_valid;
}

double* connex_comp(Mesh* Msh)
{
  printf("creating neighbors list of the mesh \n");
  msh_neighbors(Msh) ;
  printf("neighbors list of the mesh done \n");
  printf("creating connex composante of the mesh \n");

  // init
  double* color_trig = calloc(Msh->NbrTri +1, sizeof(double)); // stored color of the triangles
  int* influence_ring = calloc(Msh->NbrTri +1 , sizeof(int));  // neighbohood of the seed, the pile   /!\ this might overload.
  int  nbr_influence = 1;                                      // top of the pile
  int  nbr_colored = 0;                                        // number of colored triangles
  int  id_last_seed;                                           // last seed
  int  color=1;                                                // color of the CC

  for(int i_trig=1;i_trig<Msh->NbrTri; i_trig++) color_trig[i_trig]=0;

  id_last_seed = 1;

  int trig_ring;

  printf("starting the loop on sub domains \n");
  while(id_last_seed<Msh->NbrTri && nbr_colored<Msh->NbrTri  ) // loop on subdomains
  {
    // init the seed
    influence_ring[1]= id_last_seed;  nbr_influence =1;           
    trig_ring = id_last_seed;

    // loop on the rings of neighbors  
    while (trig_ring !=0 && (nbr_influence<Msh->NbrTri+1 && nbr_influence>0))         
    { 
      trig_ring = influence_ring[nbr_influence];

      // if the element isn't colored, color it and add it's valid edges
      if(color_trig[trig_ring]==0)  
      {

        color_trig[trig_ring] = color; nbr_colored +=1;nbr_influence -=1;
        for(int i_edg=0; i_edg<3; i_edg++)
        {
          // an edge is valid if the edge isn't a border and the neighboring triangle isn't colored
          if(valid_edge(Msh,trig_ring,i_edg) && (color_trig[Msh->TriVoi[trig_ring][i_edg]] == 0 && Msh->TriVoi[trig_ring][i_edg] !=0))
          {
            // if the edge is valid, add on top of the pile
            nbr_influence +=1;
            influence_ring[nbr_influence] = Msh->TriVoi[trig_ring][i_edg];  
          }
        }
      }
      // if the element is colored, skip to the next element
      // the pile is one less high.
      else{ nbr_influence -=1;}
      

    }
    // either all neighbors have been explored or all triangles have a color
    // if all triangles have a color, the main loop will stop
    // if else create a new subdomain and restart the process
  
    color ++;
    for(int trig_ = 1; trig_<Msh->NbrTri; trig_ ++){ if(color_trig[trig_]==0){id_last_seed = trig_; break;}}
  }
  printf("sub domains created. there is %d sub domaines \n",color-1);
  return color_trig;
}

double surf(double2d P1, double2d P2, double2d P3)
{
  double K_surf = (0.5*((P2[0]-P1[0])*(P3[1]-P1[1]) - (P2[1]-P1[1])*(P3[0]-P1[0])));
  return K_surf;
}

double quality(double2d P1, double2d P2, double2d P3)
{
  
  double alpha = sqrt(3)/12;
  double K_surf;

  double a,b,c; 
  
  a = (P1[0]-P2[0])*(P1[0]-P2[0]) + (P1[1]-P2[1])*(P1[1]-P2[1]); 
  b = (P2[0]-P3[0])*(P2[0]-P3[0]) + (P2[1]-P3[1])*(P2[1]-P3[1]); 
  c = (P3[0]-P1[0])*(P3[0]-P1[0]) + (P3[1]-P1[1])*(P3[1]-P1[1]); 
  K_surf = fabs(surf(P1,P2,P3));
  
  if(K_surf < 1e-30){printf("\n ---- \n ERROR SURFACE NULLE (%f,%f),(%f,%f) (%f,%f) \n ---- \n",
                                                    P1[0],P1[1], P2[0],P2[1], P3[0],P3[1]);}
  
  double Qal = alpha *(a+b+c)/K_surf; 
  return Qal;
}

// ============================================================================
// ============================================================================

int Is_Inside_Circle_2(double2d Point, double2d P1,double2d P2, double2d P3)
{
  double R, dist;
  double a,b,c,K_surf; 
  a = (P1[0]-P2[0])*(P1[0]-P2[0]) + (P1[1]-P2[1])*(P1[1]-P2[1]); 
  b = (P2[0]-P3[0])*(P2[0]-P3[0]) + (P2[1]-P3[1])*(P2[1]-P3[1]); 
  c = (P3[0]-P1[0])*(P3[0]-P1[0]) + (P3[1]-P1[1])*(P3[1]-P1[1]); 
  K_surf = surf(P1,P2,P3);

  R = sqrt(a*b*c)/(4* (K_surf));
  double coord1,coord2,coord3;
  printf(" distance du triangle : (%f,%f,%f) \n", a,b,c);
  coord3 = b*(c+a -b); coord2 = c*(b+a -c); coord1 = a*(b+c -a);
  printf(" coords of the center in barycenter coord : (%f,%f,%f) \n", coord1, coord2, coord3);
  double2d P =  {coord1*P1[0]+ coord2*P2[0]+ coord3*P3[0],coord1*P1[1]+ coord2*P2[1]+ coord3*P3[1]};

  printf(" Center : (%f,%f) \n", P[0],P[1]);
  R = fmax( (P[0]-P1[0])*(P[0]-P1[0]) +  (P[1]-P1[1])*(P[1]-P1[1]), fmax( (P[0]-P2[0])*(P[0]-P2[0]) +  (P[1]-P2[1])*(P[1]-P2[1]), (P[0]-P3[0])*(P[0]-P3[0]) +  (P[1]-P3[1])*(P[1]-P3[1])));
  dist = (P[0]-Point[0])*(P[0]-Point[0]) +  (P[1]-Point[1])*(P[1]-Point[1]);
  
  printf(" Point : (%f,%f) \n", Point[0], Point[1]);
  printf(" P1 : (%f,%f) ", P1[0], P1[1]);
  printf(" P2 : (%f,%f) ", P2[0], P2[1]);
  printf(" P3 : (%f,%f) ", P3[0], P3[1]);
  printf("\n distance : %f vs R : %f \n", dist,R);
  if(dist>R*R) return 0;
  if(dist<=R*R) return 1;

  return -1;
}

int Is_Inside_Circle(double2d Point, double2d P1,double2d P2, double2d P3)
{
    // Let d be a determinant
    // | ax-x   ay-y  (ax-x)²+(ay-y)² |
    // | bx-x   by-y  (bx-x)²+(by-y)² |
    // | cx-x   cy-y  (cx-x)²+(cy-y)² |
    // if d = 0, then (x,y) is on the circle,
    // if d > 0, then (x,y) is in the circle,
    // if d < 0, then (x,y) is outside the circle

    double d11 = P1[0]-Point[0];
    double d12 = P1[1]-Point[1];
    double d13 = (P1[0]-Point[0])*(P1[0]-Point[0]) + (P1[1]-Point[1])* (P1[1]-Point[1]);

    double d21 = P2[0]-Point[0];
    double d22 = P2[1]-Point[1];
    double d23 = (P2[0]-Point[0])*(P2[0]-Point[0]) + (P2[1]-Point[1])*(P2[1]-Point[1]);

    double d31 = P3[0]-Point[0];
    double d32 = P3[1]-Point[1];
    double d33 = (P3[0]-Point[0])*(P3[0]-Point[0])+ ( P3[1]-Point[1])*( P3[1]-Point[1]) ;

    return d11*d22*d33 + d12*d23*d31 + d13*d21*d32 - d13*d22*d31 - d23*d32*d11 - d33*d12*d21 > 0.0;
}

int ajout_point(Mesh* Msh, double2d Point)
{
  //TODO
  // add the case -> on vertex

  // localisation of the point in the Mesh
  printf(" \n ================================ \n \n");
  printf(" INIT the method \n");
  printf(" ADDING POINT P : (%f, %f) to the mesh \n", Point[0], Point[1]);
  int in_trig = 0; // 0 if in the triangle, 1 if it is. 
  int on_edge = 0; // 0 if not on edge, 1 if it is.
  int compass; // check if the north( direction ) has been chosen
  int iTri =1, edge_Tri=0;
  int i1,i2,i3;
  double ax1,ax2,ax3, K;
  while(in_trig ==0 && on_edge ==0)
  {
    compass =0;
    printf(" In triangle : %d \n",iTri);
    i1 = Msh->Tri[iTri][0];        i2 = Msh->Tri[iTri][1];        i3 = Msh->Tri[iTri][2];
    double2d P1 = {Msh->Crd[i1][0], Msh->Crd[i1][1] };
    double2d P2 = {Msh->Crd[i2][0], Msh->Crd[i2][1] };
    double2d P3 = {Msh->Crd[i3][0], Msh->Crd[i3][1] };
    K =surf(P1,P2,P3);
    ax1 = surf(Point,P2,P3)/K; ax2 = surf(P1,Point, P3)/K; ax3 = surf(P1,P2,Point)/K;

    //check
    if(ax1<0 && (ax2<0 && ax3<0) ) printf("\n ERROR POINT OR TRIANGLE WRONGLY DEFINED, IGNORED \n");

    // 3 direct neighbours
    printf(" from the directions available : (%d,%d,%d) \n", Msh->TriVoi[iTri][0],Msh->TriVoi[iTri][1],Msh->TriVoi[iTri][2]);
    printf(" going in direction : (%d,%d,%d) towards : %d  \n", ax1<0, ax2<0, ax3<0, iTri);

    // we could change the choice of the direction
    // here we have a bias against the direction 1 

    if(fabs(ax1)<1e-30 && (ax2>0 && ax3>0)){ on_edge =1; edge_Tri = Msh->TriVoi[iTri][0]; compass =1;}
    if(fabs(ax2)<1e-30 && (ax1>0 && ax3>0)){ on_edge =1; edge_Tri = Msh->TriVoi[iTri][1]; compass =1;}
    if(fabs(ax3)<1e-30 && (ax1>0 && ax2>0)){ on_edge =1; edge_Tri = Msh->TriVoi[iTri][2]; compass =1;}

    if(ax1<0 && compass ==0){ iTri = Msh->TriVoi[iTri][0]; compass =1;}
    if(ax2<0 && compass ==0){ iTri = Msh->TriVoi[iTri][1]; compass =1;}
    if(ax3<0 && compass ==0){ iTri = Msh->TriVoi[iTri][2]; compass =1;}


    if((ax1>0 && (ax2>0 && ax3>0)) && on_edge ==0 ){ in_trig=1; printf(" Triangle Found ! \n");}
    // printf(" ---------------- \n");
  }

    
  printf(" \n -------------------------- \n \n");
  if(in_trig)  printf(" Point localised in Triangle : %d \n", iTri);
  if(on_edge)  printf(" Point localised on edge between %d and %d \n", iTri, edge_Tri);
  printf(" starting the cavity search \n");

  int id_tri = iTri;

  // calculting the elements to add to the cavity
  int* pile   = calloc(Msh->NbrTri +1, sizeof(int));  // pile of element to add to the mavity field 
  int* cavity = calloc(Msh->NbrTri +1, sizeof(int));  // elements to remove
  int  sizeof_pile = 1;                               // size of the pile / index of the top of the pile
  int  sizeof_cavity = 0;                             // size of the cavity/index of the last element added to the cavity 
  int jTri;
  int Is_In_Pile=0, Is_In_Cavity=0;

  if(in_trig ==1)  pile[1] = id_tri;
  if(on_edge ==1){ pile[1] = id_tri; pile[2] = edge_Tri; sizeof_pile +=1;}

  while(sizeof_pile>0)
  {
    printf(" \n -------------------------- \n");
    // add the element to the cavity and remove it from the pile
    iTri = pile[sizeof_pile]; sizeof_pile -=1;
    cavity[sizeof_cavity] = iTri; sizeof_cavity +=1;

    printf(" adding %d to the cavity \n",iTri);
    // printf(" checking it's neighbours \n \n");

    // check if the neighbours are in the cavity
    for(int jEdg=0;jEdg<3; jEdg++) 
    {
      Is_In_Cavity =0; Is_In_Pile =0;
      jTri = Msh->TriVoi[iTri][jEdg];
      // printf(" neighbour : %d \n", jTri);

      i1 = Msh->Tri[jTri][0];        i2 = Msh->Tri[jTri][1];        i3 = Msh->Tri[jTri][2];
      double2d P1 = {Msh->Crd[i1][0], Msh->Crd[i1][1]};
      double2d P2 = {Msh->Crd[i2][0], Msh->Crd[i2][1]};
      double2d P3 = {Msh->Crd[i3][0], Msh->Crd[i3][1]};

      // printf(" test \n");
      // check if the element hasn't already been added
      for(int i=0;i<=sizeof_cavity; i++){ if(cavity[i]==jTri ){ Is_In_Cavity = 1; break;} }               
      for(int i=0;i<=sizeof_pile  ; i++){ if(pile[i]  ==jTri ){ Is_In_Pile   = 1; break;} }

      // printf(" %d, %d, %d", Is_In_Cavity )
      if(Is_In_Cavity==0 && (Is_In_Pile==0 && jTri !=0))
      {
        // printf("the neighbour wasn't in the pile \n");
        // if it's inside the circle      
        if(Is_Inside_Circle(Point, P1,P2,P3)==1)
        {
          // add to the pile    
          // printf(" adding to the pile \n");
          sizeof_pile+=1;
          pile[sizeof_pile] = jTri;
        }
      }
      
    }
  }


  printf(" \n -------------------------- \n \n");
  // hash_out(Msh->Hsh);

  printf(" cavity search done \n");
  // printf(" we will be deleting %d elements \n", sizeof_cavity);
  // printf(" we will delete the elements : \n");
  for(int i=0; i<sizeof_cavity;i++){printf(" a "); printf(" Tri %d : (%d,%d,%d) \n", cavity[i], Msh->Tri[cavity[i]][0], Msh->Tri[cavity[i]][1], Msh->Tri[cavity[i]][2] );}

  printf(" b");
  // deleting the elements of the cavity
  int Tri_neigh,i_obj,iVer1,iVer2;

  int sizeof_boundary=0;
  int2d* boundary_cavity;
  
  printf(" %d", Msh->NbrVer);
  
  boundary_cavity = (int2d*)calloc(2*(Msh->NbrVer) +1 , sizeof(int2d));

  printf(" c ");
  for(int i_cavity=0; i_cavity<sizeof_cavity; i_cavity ++)
  {
    //----------------------------------------------------------------------------------
    // deleting it from it's neighbours's list of neighbours
    // printf(" deleting from the neighbours' list \n");
    for(int i=0; i<3; i++)
    {
      Tri_neigh =Msh->TriVoi[cavity[i_cavity]][i]; 
      // printf(" %d ",Tri_neigh ); // check all neighbours 
      for(int j=0; j<3; j++)
      {
        // printf(" checking neighbours %d",Tri_neigh);
        if(Msh->TriVoi[Tri_neigh][j]==cavity[i_cavity]) Msh->TriVoi[Tri_neigh][j] = 0; // removing it from their neighbours
      }
    }
    printf("\n");

    //----------------------------------------------------------------------------------
    // deleting the edge from the Hashtable and calculating the boundary of the cavity
    for(int i_Edg=0;i_Edg<3;i_Edg++)
    {
      printf(" deleting from the hashtable element %d's edge %d \n", cavity[i_cavity],i_Edg);
      iVer1 = Msh->Tri[cavity[i_cavity]][tri2edg[i_Edg][0]];
      iVer2 = Msh->Tri[cavity[i_cavity]][tri2edg[i_Edg][1]];
      printf("deleting (%d,%d) \n",iVer1,iVer2); 
      i_obj=hash_find(Msh->Hsh,iVer1,iVer2);   
      printf(" i.e object : %d \n",i_obj);        
      if(i_obj!=0) hash_suppr(Msh->Hsh,iVer1,iVer2,cavity[i_cavity]);
    }
    //----------------------------------------------------------------------------------   
  } 
  printf(" \n -------------------------- \n");
  printf(" creating the boundary of the cavity \n");
  for(int i_cavity=0; i_cavity<sizeof_cavity; i_cavity ++)
  {
    printf("    checking element : %d \n", cavity[i_cavity]);
    for(int i_Edg=0;i_Edg<3;i_Edg++)
    {
      iVer1 = Msh->Tri[cavity[i_cavity]][tri2edg[i_Edg][0]];
      iVer2 = Msh->Tri[cavity[i_cavity]][tri2edg[i_Edg][1]];
      i_obj=hash_find(Msh->Hsh,iVer1,iVer2);
      Is_In_Cavity =0;

      for(int i=0;i<sizeof_boundary;i++) 
      {
        if((boundary_cavity[i][0] ==iVer1 && boundary_cavity[i][1]==iVer2 ) || (boundary_cavity[i][1] ==iVer1 && boundary_cavity[i][0]==iVer2) )
        { 
          // if it is in the boundary -> remove it  
          Is_In_Cavity= 1;
          boundary_cavity[i][0]=boundary_cavity[sizeof_boundary][0];
          boundary_cavity[i][1]=boundary_cavity[sizeof_boundary][1];
          // boundary_cavity[sizeof_boundary]=0;
          printf(" removing element (%d,%d) \n", iVer1, iVer2);
          // sizeof_boundary -=1;
          break; 
        }
      }     
        
      // if it is not in the boundary -> adding it
      if(Is_In_Cavity==0)
      {
      printf(" adding element (%d,%d) \n", iVer1, iVer2);
      boundary_cavity[sizeof_boundary][0]=iVer1;
      boundary_cavity[sizeof_boundary][1]=iVer2;
      sizeof_boundary +=1;
      }
    } 
  } 
  
  printf(" \n -------------------------- \n");
  printf(" the cavity is of size %d \n", sizeof_cavity);
  for(int i=0; i<sizeof_cavity; i++){printf(" %d ", cavity[i]);}
  printf("\n");

  
  printf(" sorting the boundary \n");
  int NumOfZero=0;
  for(int i =0;i<sizeof_boundary-NumOfZero;i++)
  {
    if(boundary_cavity[i][0]==0)
    {
      // printf(" (%d,%d) ", boundary_cavity[i][0], boundary_cavity[i][1]);
      // printf(" i: %d, Nz :%d, Sb: %d \n", i,NumOfZero, sizeof_boundary);
      boundary_cavity[i][0] = boundary_cavity[sizeof_boundary-NumOfZero-1][0];  boundary_cavity[sizeof_boundary-NumOfZero-1][0]=0;
      boundary_cavity[i][1] = boundary_cavity[sizeof_boundary-NumOfZero-1][1];  boundary_cavity[sizeof_boundary-NumOfZero-1][1]=0;
      NumOfZero +=1;
    }
  }
  sizeof_boundary -= NumOfZero;
  

  printf(" the boundary is of size %d \n", sizeof_boundary);
  for(int i=0; i<=sizeof_boundary; i++){printf(" (%d,%d) ", boundary_cavity[i][0], boundary_cavity[i][1]);}
  printf("\n");
  
  printf(" \n -------------------------- \n \n");
  printf(" cavity made \n");


  printf(" \n -------------------------- \n \n");

  // Filling the cavity with the star centered on P
  
  printf(" reallocating the memory \n");
  // ===========================================================================
  // realloc the memory
  // NbrTri & Tri & TriVoi & TriRef & Crd
  Msh->NbrTriMax = Msh->NbrTriMax + sizeof_boundary - sizeof_cavity;
  Msh->NbrTri = Msh->NbrTri - sizeof_cavity;

  printf(" there will be at most %d triangles at the end \n", Msh->NbrTriMax);

  printf("allocating memory for Tri \n");
  int3d* temp_Tri = realloc(Msh->Tri, (Msh->NbrTriMax+1  )*sizeof(int3d));
  if (temp_Tri == NULL) {
    // If reallocation fails
    printf("ERROR. Unable to resize memory");
  } else {
    // If reallocation is successful
    Msh->Tri = temp_Tri;  // Update ptr1 to point to the newly allocated memory
    // free(temp_Tri); 
    // temp_Tri= NULL;
  } 
  
  
  printf("allocating memory for TriRef \n");
  int1d* temp_Ref = realloc(Msh->TriRef, (Msh->NbrTriMax+1 )*sizeof(int1d));
  if (temp_Ref == NULL) {
    // If reallocation fails
    printf("ERROR. Unable to resize memory");
  } else {
    // If reallocation is successful
    Msh->TriRef = temp_Ref;  // Update ptr1 to point to the newly allocated memory
    // free(temp_Ref); 
    // temp_Ref= NULL;
  } 

  printf("allocating memory for TriVoi \n");
  int3d* temp_Voi = realloc(Msh->TriVoi, (Msh->NbrTriMax+1 )*sizeof(int3d));
  if (temp_Voi == NULL) {
    // If reallocation fails
    printf("ERROR. Unable to resize memory");
  } else {
    // If reallocation is successful
    Msh->TriVoi = temp_Voi;  // Update ptr1 to point to the newly allocated memory
    
    // free(temp_Voi); 
    // temp_Voi= NULL;
  } 

  Msh->NbrVer +=1; Msh->NbrVerMax +=1;
  
  printf("allocating memory for Crd \n");
  double2d* temp_crd = realloc(Msh->Crd , (Msh->NbrVerMax+1 )*sizeof(double3d));
  if (temp_crd == NULL) {
    // If reallocation fails
    printf("ERROR. Unable to resize memory");
  } else {
    // If reallocation is successful
    Msh->Crd = temp_crd;  // Update ptr1 to point to the newly allocated memory
    Msh->Crd[Msh->NbrVer][0] = Point[0];
    Msh->Crd[Msh->NbrVer][1] = Point[1];
    
    // free(temp_crd); 
    // temp_crd= NULL;
    } 

  // ---------------------------------------- Hash Table -----------------------------------------
  Msh->Hsh->NbrMaxObj += (1+(sizeof_boundary-sizeof_cavity)) ; 

  printf("allocating memory for LstObj \n");
  int5d* temp_Lst =  realloc(Msh->Hsh->LstObj , (Msh->Hsh->NbrMaxObj )*sizeof(int5d)); 
  if (temp_Lst == NULL) {
    // If reallocation fails
    printf("ERROR. Unable to resize memory");
  } else {
    // If reallocation is successful
    Msh->Hsh->LstObj = temp_Lst;  // Update ptr1 to point to the newly allocated memory
    for(int j=Msh->Hsh->SizHead + 1;    j<2*(Msh->Hsh->NbrMaxObj);   j++)
    {for(int i=0;i<5;i++) Msh->Hsh->LstObj[j][i]=0;}
    
    // free(temp_Lst); 
    // temp_Lst= NULL;
  } 
  
  printf("allocating memory for Head \n");
  int* temp_Hd =  realloc(Msh->Hsh->Head ,  2*(Msh->NbrVerMax )*sizeof(int)); 
  if (temp_Hd == NULL) {
    // If reallocation fails
    printf("ERROR. Unable to resize memory");
  } else {
    // If reallocation is successful
    Msh->Hsh->Head = temp_Hd;  // Update ptr1 to point to the newly allocated memory
    for(int j=Msh->Hsh->SizHead + 1;    j<2*(Msh->NbrVerMax );   j++){Msh->Hsh->Head[j]=0;}

    Msh->Hsh->SizHead = 2*(Msh->NbrVerMax);
    // free(temp_Hd); 
    // temp_Hd= NULL;
  } 
  // ===========================================================================

  printf(" allocated the memory \n");
  int Tri_added;
  int j_hsh,jVer1,jVer2;
  for(int i_boundary=0; i_boundary<sizeof_boundary; i_boundary++)
  { 
    while( (boundary_cavity[i_boundary][0]==0 && boundary_cavity[i_boundary][1]==0) && i_boundary<sizeof_boundary){i_boundary +=1;}
    // if(i_boundary == sizeof_boundary) break;

    // printf(" test : %d \n", cavity[sizeof_cavity-1]);
    if(sizeof_cavity>0){Tri_added = cavity[sizeof_cavity-1];}
    else{Tri_added = Msh->NbrTri+1;}

    // printf(" adding tri : %d \n",Tri_added);
    
    // add to Tri
    iVer1 = boundary_cavity[i_boundary][0]; 
    iVer2 = boundary_cavity[i_boundary][1]; 
    // printf(" triangle linked to the boundary : (%d,%d) \n", iVer1,iVer2 );
    
    if(surf(Msh->Crd[iVer1],Msh->Crd[iVer2],Point)<0)
    {
      // printf(" NEGATIVE SURFACE \n");
    Msh->Tri[Tri_added][0]= iVer2;
    Msh->Tri[Tri_added][1]= iVer1;
    Msh->Tri[Tri_added][2]= Msh->NbrVer;
    Msh->TriRef[Tri_added] =1;
    }

    else
    {      
      // printf(" POSITIVE SURFACE \n");
      Msh->Tri[Tri_added][0]= iVer1;
      Msh->Tri[Tri_added][1]= iVer2;
      Msh->Tri[Tri_added][2]= Msh->NbrVer;
    Msh->TriRef[Tri_added] =1;
    }
    
    // printf(" added to the Mesh to add to the HashTable \n");
    // update the hashtable and TriVoi
    for(int iEdg=0;iEdg<3;iEdg++)
    {
      iVer1 = Msh->Tri[Tri_added][tri2edg[iEdg][0]];
      iVer2 = Msh->Tri[Tri_added][tri2edg[iEdg][1]];

      // printf("test \n");
      j_hsh = hash_find(Msh->Hsh,iVer1,iVer2);
      // printf(" %d, (%d,%d) (x,%d) \n",j_hsh, iVer1,iVer2,Tri_added);
      hash_add(Msh->Hsh,iVer1,iVer2,Tri_added,j_hsh);
      jTri = Msh->Hsh->LstObj[j_hsh][2];

      Msh->TriVoi[Tri_added][iEdg] = jTri;
      for(int jEdg=0;jEdg<3;jEdg++)
        {
          jVer1 = Msh->Tri[jTri][tri2edg[jEdg][0]];
          jVer2 = Msh->Tri[jTri][tri2edg[jEdg][1]];
            if((jVer1 == iVer1 && jVer2 == iVer2) || (jVer1 == iVer2 && jVer2 == iVer1))
            {Msh->TriVoi[jTri][jEdg] = Tri_added;}
        }
      
    }

    sizeof_cavity -=1; Msh->NbrTri +=1;
  }
  printf(" \n ====================================  \n");

  printf(" POINT ADDED !! \n Point (%f,%f) will be named as %d thereafter. \n", Point[0],Point[1],Msh->NbrVer);

  // printf(" releasing memory \n");

  // if(boundary_cavity !=NULL)
  // {
  //   free(boundary_cavity); 
  //   boundary_cavity = NULL;
  // }
  // if(cavity !=NULL)
  // {
  //   free(cavity); 
  //   cavity= NULL;
  // }
  // if(pile !=NULL)
  // {
  //   free(pile); 
  //   pile= NULL;
  // }

  // if(temp_crd !=NULL)
  // {
  //   free(temp_crd); 
  //   temp_crd= NULL;
  // }

  // if(temp_Hd !=NULL)
  // {
  //   free(temp_Hd); 
  //   temp_Hd= NULL;
  // }
  // printf(" bl ");
  // if(temp_Ref !=NULL)
  // {
  //   free(temp_Ref); 
  //   temp_Ref= NULL;
  // }
  
  // printf(" bl ");
  // if(temp_Voi !=NULL)
  // {
  //   free(temp_Voi); 
  //   temp_Voi= NULL;
  // }
  // printf(" bl ");
  // if(temp_Lst !=NULL)
  // {
  //   free(temp_Lst); 
  //   temp_Lst= NULL;
  // }
  
  // printf(" memory released \n");




  return 0;
}

int Maillage_Delauney(int Nb_Point)
{
  srand(10);
  //--- read a mesh
  Mesh* Msh = msh_read("../data/carre_base.mesh", 0);

  double Crd_x,Crd_y;

  msh_neighbors(Msh);
  hash_out(Msh->Hsh);

  for(int it=0;it<Nb_Point;it++)
  {
    Crd_x=(double)(rand())/RAND_MAX;
    while(Crd_x<1e-2 || Crd_x>1-1e-2) Crd_x=(double)(rand())/RAND_MAX;
    Crd_y=(double)(rand())/RAND_MAX;
    while(Crd_y<1e-2 || Crd_y>1-1e-2) Crd_y=(double)(rand())/RAND_MAX;

    ajout_point(Msh, (double2d){Crd_x,Crd_y} );
    hash_out(Msh->Hsh);
  }

  printf("\n data : tri %d/%d, Ver %d/%d \n", Msh->NbrTri, Msh->NbrTriMax, Msh->NbrVer, Msh->NbrVerMax);


  printf(" \n outing the mesh \n");
  msh_write(Msh,"TEST.mesh");
  return 0;
}
// ============================================================================
// ============================================================================