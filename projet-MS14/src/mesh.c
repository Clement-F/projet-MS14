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

  int fmsh = GmfOpenMesh(file, GmfWrite, FilVer, Msh->Dim);
  if (fmsh <= 0) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }

  GmfSetKwd(fmsh, GmfVertices, Msh->NbrVer);
  for (iVer = 1; iVer <= Msh->NbrVer; iVer++)
    GmfSetLin(fmsh, GmfVertices, Msh->Crd[iVer][0], Msh->Crd[iVer][1], 0);

  GmfSetKwd(fmsh, GmfTriangles, Msh->NbrTri);
  for (iTri = 1; iTri <= Msh->NbrTri; iTri++)
    GmfSetLin(fmsh, GmfTriangles, Msh->Tri[iTri][0], Msh->Tri[iTri][1], Msh->Tri[iTri][2], Msh->TriRef[iTri]);

  if (Msh->NbrEfr > 0) {
    GmfSetKwd(fmsh, GmfEdges, Msh->NbrEfr);
    for (iEfr = 1; iEfr <= Msh->NbrEfr; iEfr++)
      GmfSetLin(fmsh, GmfEdges, Msh->Efr[iEfr][0], Msh->Efr[iEfr][1], Msh->EfrRef[iEfr]);
  }

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

  // hash_out(hsh);
  hash_bound(hsh);
  hash_collision(hsh);

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
    else printf("  ## WARNING: HSH ELEMENT ALREADY COMPLETE. IGNORED\n");
    // printf("updating element %d of Vertex (%d,%d) and Tri (%d,%d) to have next %d \n", i_hsh, iVer1,iVer2, hsh->LstObj[hsh->NbrObj][3], iTri, hsh->LstObj[hsh->NbrObj][4]);
  }

  return 0;
}

// NEED TEST
int hash_suppr(HashTable* hsh, int iVer1, int iVer2, int iTri)  // deletes an element of the hash table
{
  
  int i_hsh;
  i_hsh = hash_find(hsh,iVer1,iVer2); // check if the element is in the hash_list


  if(i_hsh ==0){ printf("\n DELETING A NON EXISTING ELEMENT, IGNORED \n ");}
  if(i_hsh !=0)
  {
    int ToDelete=0;
    if((hsh->LstObj[i_hsh][2]!=0 && hsh->LstObj[i_hsh][3]==0) || (hsh->LstObj[i_hsh][3]!=0 && hsh->LstObj[i_hsh][2]==0)) ToDelete = 1;
    if(ToDelete==1)
    {
    // we redo the chain
    int key = iVer1 + iVer2; 
    int j_hsh = hsh->Head[key];
    int i_bef = 0;
    while(j_hsh!=0)
    {
      if(hsh->LstObj[j_hsh][0]==iVer1 && hsh->LstObj[j_hsh][1]==iVer2)
      {
        // sewing it
        if(i_bef==0){hsh->Head[key]       =hsh->LstObj[j_hsh][4];} 
        if(i_bef!=0){hsh->LstObj[i_bef][4]=hsh->LstObj[j_hsh][4];}

        // permuting it with the last element
        int last_index = hsh->NbrObj;
        hsh->LstObj[j_hsh][0] = hsh->LstObj[last_index][0]; hsh->LstObj[j_hsh][1] =hsh->LstObj[last_index][1];
        hsh->LstObj[j_hsh][2] = hsh->LstObj[last_index][2]; hsh->LstObj[j_hsh][3] =hsh->LstObj[last_index][3];
        hsh->LstObj[j_hsh][4] = hsh->LstObj[last_index][4];

        // sewing the last element's chain 
        key = hsh->LstObj[last_index][0] + hsh->LstObj[last_index][1];
        int k_hsh= hsh->Head[key];
        
        int j_bef = 0;
        while(k_hsh!=0)
        {
          if(hsh->LstObj[k_hsh][0]==iVer1 && hsh->LstObj[k_hsh][1]==iVer2)
          {
            if(j_bef==0){hsh->Head[key]       =hsh->LstObj[k_hsh][4];} 
            if(j_bef!=0){hsh->LstObj[j_bef][4]=hsh->LstObj[k_hsh][4];}
          }
          else k_hsh = hsh->LstObj[j_bef][4];
        }
        
      }
      else j_hsh = hsh->LstObj[j_hsh][4];      
    } 
    }

    // we put an element to 0 and potentially put the element in a "good form"
    if(ToDelete==0)
    {
      if(hsh->LstObj[i_hsh][3]==iTri) hsh->LstObj[i_hsh][3]=0;
      if(hsh->LstObj[i_hsh][2]==iTri) 
      {
        hsh->LstObj[i_hsh][2]= hsh->LstObj[i_hsh][3];
        hsh->LstObj[i_hsh][3]= 0;
      }  
    }
  }
  return 0;
}

int hash_out(HashTable* hsh)  // print the hash table
{
  for(int j=0;j<hsh->NbrMaxObj;j++){
    printf("indices = %d \n", j);
    printf("Vertexes  : %d, %d \n", hsh->LstObj[j][0], hsh->LstObj[j][1]);
    printf("Triangles : %d, %d \n", hsh->LstObj[j][2], hsh->LstObj[j][3]);
    printf("next : %d\n", hsh->LstObj[j][4]);
    printf("----------- \n");
  }
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
  double K_surf = fabs(0.5*((P2[0]-P1[0])*(P3[1]-P1[1]) - (P2[1]-P1[1])*(P3[0]-P1[0])));
  if(K_surf < 1e-30){printf("ERROR SURFACE NULLE");}
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
  K_surf = surf(P1,P2,P3);
  
  double Qal = alpha *(a+b+c)/K_surf; 
  return Qal;
}

// ============================================================================
// ============================================================================

int Is_Inside_Circle(double2d Point, double2d P1,double2d P2, double2d P3)
{
  double R, dist;
  double a,b,c,K_surf; 
  a = (P1[0]-P2[0])*(P1[0]-P2[0]) + (P1[1]-P2[1])*(P1[1]-P2[1]); 
  b = (P2[0]-P3[0])*(P2[0]-P3[0]) + (P2[1]-P3[1])*(P2[1]-P3[1]); 
  c = (P3[0]-P1[0])*(P3[0]-P1[0]) + (P3[1]-P1[1])*(P3[1]-P1[1]); 
  K_surf = surf(P1,P2,P3);

  R = sqrt(a*b*c)/(4*K_surf);
  double coord1,coord2,coord3;
  coord1 = a*(b+c-a); coord2 = b*(c+a-b); coord3 = c*(b+a-c);
  double2d P =  {coord1*P1[0]+ coord2*P2[0]+ coord3*P3[0],coord1*P1[1]+ coord2*P2[1]+ coord3*P3[1]};
  dist = (P[0]-Point[0])* (P[0]-Point[0]) +  (P[1]-Point[1])*(P[1]-Point[1]);
  
  if(dist>R*R) return 0;
  if(dist<=R*R) return 1;

  return -1;
}

int ajout_point(Mesh* Msh, double2d Point)
{
  // localisation of the point in the Mesh
  int in_trig = 0; // 0 if in the triangle, 1 if it is. 
  int iTri =1;
  int i1,i2,i3;
  double ax1,ax2,ax3;
  while(in_trig ==0)
  {
    i1 = Msh->Tri[iTri][0];        i2 = Msh->Tri[iTri][1];        i3 = Msh->Tri[iTri][2];
    double2d P1 = {Msh->Crd[i1][0], Msh->Crd[i1][1] };
    double2d P2 = {Msh->Crd[i2][0], Msh->Crd[i2][1] };
    double2d P3 = {Msh->Crd[i3][0], Msh->Crd[i3][1] };

    ax1 = surf(Point,P2,P3); ax2 = surf(P1,Point, P3); ax3 = surf(P1,P2,Point);

    //check
    if(ax1<0 && (ax2<0 && ax3<0) ) printf("\n ERROR POINT OR TRIANGLE WRONGLY DEFINED, IGNORED \n");

    // 3 direct neighbours
    if(ax1<0) iTri = Msh->TriVoi[iTri][0];
    if(ax2<0) iTri = Msh->TriVoi[iTri][1];
    if(ax3<0) iTri = Msh->TriVoi[iTri][2];

    // 3 non direct neighbours
    // if(ax1<0 && ax2<0) iTri = Msh->TriVoi[Msh->TriVoi[iTri][0]][0];
    // if(ax2<0 && ax3<0) iTri = Msh->TriVoi[Msh->TriVoi[iTri][0]][0];
    // if(ax3<0 && ax1<0) iTri = Msh->TriVoi[Msh->TriVoi[iTri][0]][0];

    if(ax1>0 && (ax2>0 && ax3>0)) in_trig=1;
  }

  int id_tri = iTri;

  // calculting the elements to add to the cavity
  int* pile   = calloc(Msh->NbrTri +1, sizeof(int));  // pile of element to add to the mavity field 
  int* cavity = calloc(Msh->NbrTri +1, sizeof(int));  // elements to remove
  int  sizeof_pile = 1;                               // size of the pile / index of the top of the pile
  int  sizeof_cavity = 0;                             // size of the cavity/index of the last element added to the cavity 
  
  pile[0] = id_tri;
  int jTri;
  int Is_In_Pile, Is_In_Cavity;

  while(sizeof_pile>0)
  {
    // add the element to the cavity and remove it from the pile
    iTri = pile[sizeof_pile]; sizeof_pile -=1;
    cavity[sizeof_cavity] = iTri; sizeof_cavity +=1;

    // check if the neighbours are in the cavity
    for(int jEdg=0;jEdg<3; jEdg++) 
    {
      jTri = Msh->TriVoi[iTri][jEdg];
      i1 = Msh->Tri[iTri][0];        i2 = Msh->Tri[iTri][1];        i3 = Msh->Tri[iTri][2];
      double2d P1 = {Msh->Crd[i1][0], Msh->Crd[i1][1]};
      double2d P2 = {Msh->Crd[i2][0], Msh->Crd[i2][1]};
      double2d P3 = {Msh->Crd[i3][0], Msh->Crd[i3][1]};


      // check if the element hasn't already been added
      for(int i=0;i<sizeof_cavity; i++) if(cavity[i]==jTri ){ Is_In_Cavity = 1; break;}
      for(int i=0;i<sizeof_pile  ; i++) if(pile[i]  ==jTri ){ Is_In_Pile   = 1; break;}

      if(Is_In_Cavity==0 && Is_In_Pile==0)
      {
        // if it's inside the circle      
        if(Is_Inside_Circle(Point, P1,P2,P3)==1)
        {
          // add to the pile    
          sizeof_pile+=1;
          pile[sizeof_pile] = jTri;
        }
      }
    }
  }

  // deleting the elements of the cavity
  int Tri_neigh,i_obj,iVer1,iVer2;

  int sizeof_boundary=0;
  int* boundary_cavity; 
  boundary_cavity = calloc(2*(Msh->NbrVer) +1 , sizeof(int));

  for(int i_cavity=0; i_cavity<sizeof_cavity; i_cavity ++)
  {
    //----------------------------------------------------------------------------------
    // deleting it from it's neighbours's list of neighbours
    for(int i=0; i<3; i++)
    {
    Tri_neigh =Msh->TriVoi[i_cavity][i]; // check all neighbours 
    for(int j=0; j<3; j++) if(Msh->TriVoi[Tri_neigh][j]==i_cavity) Msh->TriVoi[Tri_neigh][j] = 0; // removing it from their neighbours
    }
    //----------------------------------------------------------------------------------
    // deleting the edge from the Hashtable and calculating the boundary of the cavity
    for(int i_Edg=0;i_Edg<3;i_Edg++)
    {
      iVer1 = Msh->Tri[i_cavity][tri2edg[i_Edg][0]];
      iVer2 = Msh->Tri[i_cavity][tri2edg[i_Edg][1]];
      i_obj=hash_find(Msh->Hsh,iVer1,iVer2);             
      hash_suppr(Msh->Hsh,iVer1,iVer2,i_cavity);
    }
    //----------------------------------------------------------------------------------   
  } 
  
  for(int i_cavity=0; i_cavity<sizeof_cavity; i_cavity ++)
  {
    for(int i_Edg=0;i_Edg<3;i_Edg++)
    {
      iVer1 = Msh->Tri[i_cavity][tri2edg[i_Edg][0]];
      iVer2 = Msh->Tri[i_cavity][tri2edg[i_Edg][1]];
      i_obj=hash_find(Msh->Hsh,iVer1,iVer2);
      Is_In_Cavity =0;
      for(int i=0;i<sizeof_boundary;i++) 
      {
        if((boundary_cavity[i]==i_obj))
        { 
          // if it is in the boundary -> remove it  
          Is_In_Cavity= 1;
          boundary_cavity[i]=boundary_cavity[sizeof_boundary];
          boundary_cavity[sizeof_boundary]=0;
          sizeof_boundary -=1;
          break;
        }
      }     
        
      // if it is not in the boundary -> adding it
      if(Is_In_Cavity==0)
      {
      boundary_cavity[sizeof_boundary]=i_obj;
      sizeof_boundary +=1;
      }
    } 
  } 

  // Filling the cavity with the star centered on P
  
  // ===========================================================================
  // realloc the memory
  // NbrTri & Tri & TriVoi & TriRef 
  Msh->NbrTriMax = Msh->NbrTriMax + sizeof_boundary - sizeof_cavity;
  Msh->NbrTri = Msh->NbrTri + sizeof_boundary - sizeof_cavity;

  int3d* temp_Tri = realloc(Msh->Tri, Msh->NbrTriMax*sizeof(int3d));
  if (temp_Tri == NULL) {
    // If reallocation fails
    printf("ERROR. Unable to resize memory");
  } else {
    // If reallocation is successful
    Msh->Tri = temp_Tri;  // Update ptr1 to point to the newly allocated memory
    free(temp_Tri); 
    temp_Tri= NULL;
  } 
  
  
  int1d* temp_Ref = realloc(Msh->TriRef, Msh->NbrTriMax*sizeof(int1d));
  if (temp_Ref == NULL) {
    // If reallocation fails
    printf("ERROR. Unable to resize memory");
  } else {
    // If reallocation is successful
    Msh->TriRef = temp_Ref;  // Update ptr1 to point to the newly allocated memory
    
    free(temp_Ref); 
    temp_Ref= NULL;
  } 

  int3d* temp_Voi = realloc(Msh->TriVoi, Msh->NbrTriMax*sizeof(int3d));
  if (temp_Voi == NULL) {
    // If reallocation fails
    printf("ERROR. Unable to resize memory");
  } else {
    // If reallocation is successful
    Msh->TriVoi = temp_Voi;  // Update ptr1 to point to the newly allocated memory
    
    free(temp_Voi); 
    temp_Voi= NULL;
  } 

  Msh->NbrVer +=1; Msh->NbrVerMax +=1;
  
  double2d* temp_crd = realloc(Msh->Crd, Msh->NbrVerMax*sizeof(double2d));
  if (temp_crd == NULL) {
    // If reallocation fails
    printf("ERROR. Unable to resize memory");
  } else {
    // If reallocation is successful
    Msh->Crd = temp_crd;  // Update ptr1 to point to the newly allocated memory
    Msh->Crd[Msh->NbrVer][0] = Point[0];
    Msh->Crd[Msh->NbrVer][1] = Point[1];
    
    free(temp_crd); 
    temp_crd= NULL;
    } 


  // ===========================================================================

  for(int i_boundary=0; i_boundary<sizeof_boundary; i_boundary++)
  {  
    // add to Tri
    iVer1 = Msh->Hsh->LstObj[boundary_cavity[i_boundary]][0]; 
    iVer2 = Msh->Hsh->LstObj[boundary_cavity[i_boundary]][1]; 
    
    if(surf(Msh->Crd[iVer1],Msh->Crd[iVer2],Point)<0)
    {
    Msh->Tri[cavity[sizeof_cavity]][0]= iVer2;
    Msh->Tri[cavity[sizeof_cavity]][1]= iVer1;
    Msh->Tri[cavity[sizeof_cavity]][2]= Msh->NbrVer;
    }
    else
    {      
      Msh->Tri[cavity[sizeof_cavity]][0]= iVer1;
      Msh->Tri[cavity[sizeof_cavity]][1]= iVer2;
      Msh->Tri[cavity[sizeof_cavity]][2]= Msh->NbrVer;
    }
    
    // update the hashtable
    hash_add(Msh->Hsh,iVer1,Msh->NbrVer,cavity[sizeof_cavity],boundary_cavity[i_boundary]);
    hash_add(Msh->Hsh,iVer2,Msh->NbrVer,cavity[sizeof_cavity],boundary_cavity[i_boundary]);
    
  }
  return 0;
}

// ============================================================================
// ============================================================================