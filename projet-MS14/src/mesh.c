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
  Msh->Edg = NULL;

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
  double x_min,x_max, y_min, y_max;
  //--- compute bounding box
  for (iVer = 1; iVer <= Msh->NbrVer; iVer++) {
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

  if (!Msh) return 0;

  if (Msh->TriVoi == NULL)
    Msh->TriVoi = calloc((Msh->NbrTri + 1), sizeof(int3d));

  //--- Compute the neighbors using a quadratic-complexity algorithm
  for (iTri = 1; iTri <= Msh->NbrTri; iTri++) {
    for (iEdg = 0; iEdg < 3; iEdg++) {
      iVer1 = Msh->Tri[iTri][tri2edg[iEdg][0]];
      iVer2 = Msh->Tri[iTri][tri2edg[iEdg][1]];

      //--- find the Tri different from iTri that has iVer1, iVer2 as vertices
      for (jTri = 1; jTri <= Msh->NbrTri; jTri++) {
        if (iTri == jTri)
          continue;

        for (jEdg = 0; jEdg < 3; jEdg++) {
          jVer1 = Msh->Tri[jTri][tri2edg[jEdg][0]];
          jVer2 = Msh->Tri[jTri][tri2edg[jEdg][1]];

          if((iVer1== jVer1 && iVer2 ==jVer2 )|| (iVer1 == jVer2 && iVer1 == jVer1)){
            Msh->TriVoi[jTri][jEdg]=iTri;
            Msh->TriVoi[iTri][iEdg]=jTri;
          }

        }
      }
    }
  }

  return 1;
}

int msh_neighbors(Mesh* Msh)
{
  printf(" init neighbors\n");
  int iTri, iEdg, iVer1, iVer2;

  if (!Msh) return 0;

  if (Msh->TriVoi == NULL)
    Msh->TriVoi = calloc((Msh->NbrTri + 1), sizeof(int3d));

  //--- initialize HashTable and set the hash table
  printf(" init hash\n");
  int SizHead = 2*(Msh->NbrVerMax);
  int NbrMaxObj = Msh->NbrVer + Msh->NbrTri -2 ;


  HashTable* hsh = hash_init(SizHead, NbrMaxObj); 

  //--- Compute the neighbors using the hash table
  for (iTri = 1; iTri <= Msh->NbrTri; iTri++) {
    for (iEdg = 0; iEdg < 3; iEdg++) {
      iVer1 = Msh->Tri[iTri][tri2edg[iEdg][0]];
      iVer2 = Msh->Tri[iTri][tri2edg[iEdg][1]];

      // TODO:
      // compute the key : iVer1+iVer2
      int key = iVer1 + iVer2;
      // do we have objects as that key hash_find () */
      int j = hsh->Head[key];
      //  if yes ===> look among objects and potentially update TriVoi */
      int i_buf,n;
      
      n=0;
      while(j !=0 && n<10)  // look through the chain
      {
        i_buf=j;
        // if the object is the one we want
        if((hsh->LstObj[j][0]==iVer1 && hsh->LstObj[j][1]==iVer2) || (hsh->LstObj[j][1]==iVer1 && hsh->LstObj[j][0]==iVer2))
        {
          printf(" ------------ (%d,%d) --------- \n", iVer1, iVer2);
          printf("Checking the Neighbor place j= %d \n",j);
          // and the triangle is the one we are on
          if(hsh->LstObj[j][2]==iTri){ Msh->TriVoi[iTri][iEdg] = hsh->LstObj[j][3]; j=0;}     // we add it to the neighbor list
          // else if(hsh->LstObj[j][3]==iTri){ Msh->TriVoi[iTri][iEdg] = hsh->LstObj[j][2]; j=0;}// we add it to the neighbor list
          else printf("  ## WARNING: HSH TABLE WRONG. IGNORED\n");
        }
        // if it isn't the one we want
        else
        {
          printf("elsevier \n");
          j = hsh->LstObj[j][4]; // we look for the next element in the list, and start again
        }
        
        n+=1;
      }

      printf(" ============= %d ; %d =============== \n", iTri, iEdg);

      //  if no  ===> add to hash table   hash_add()   */

      hash_add(hsh,iVer1,iVer2,iTri);
      hsh->LstObj[hsh->NbrObj][4] = hsh->Head[key];   //  add the element as the beginning of the key
      printf("changed the element at position %d to %d \n",hsh->NbrObj, hsh->Head[key]);
      hsh->Head[key]= hsh->NbrObj;                    //  
      printf("changed the Head at position %d to %d \n",key, hsh->Head[key]);
      printf("---------------- \n"); 
    }
  printf("================== \n");
  }

  hash_out(hsh);

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

  hsh->LstObj = calloc(hsh->NbrMaxObj+1,sizeof(int5d));
  hsh->Head   = calloc(hsh->SizHead  +1,1);

  printf("size of head :%10d \n", SizHead);
  printf("size of obj  :%10d \n", NbrMaxObj);
  for(int j=1;j<SizHead+1;j++){hsh->Head[j]=0;}

  return hsh;
}

int hash_find(HashTable* hsh, int iVer1, int iVer2)
{

  int j_hsh;
  for(j_hsh=1; j_hsh <hsh->NbrMaxObj; j_hsh++){
    if(hsh->LstObj[j_hsh][0]==iVer1 && hsh->LstObj[j_hsh][1]==iVer2){return j_hsh;}
    if(hsh->LstObj[j_hsh][1]==iVer1 && hsh->LstObj[j_hsh][0]==iVer2){return j_hsh;}
  }
  return 0;

}

int hash_add(HashTable* hsh, int iVer1, int iVer2, int iTri)
{
  int i_hsh;
  i_hsh = hash_find(hsh,iVer1,iVer2); // check if the element isn't in the hash_list already

  if(i_hsh==0)
  {
    hsh->LstObj[hsh->NbrObj +1][0] = iVer1;
    hsh->LstObj[hsh->NbrObj +1][1] = iVer2;
    hsh->LstObj[hsh->NbrObj +1][2] = iTri; 
    hsh->LstObj[hsh->NbrObj +1][3] = 0;
    hsh->LstObj[hsh->NbrObj +1][4] = 0;     
    hsh->NbrObj +=1;
    if(hsh->NbrObj> hsh->NbrMaxObj) printf("  ## WARNING: HSH ELEMENT ALREADY FULL. IGNORED\n");
  }
  if(i_hsh !=0)
  { 
    if(hsh->LstObj[i_hsh][3] ==0){hsh->LstObj[i_hsh][3] = iTri;} 
    else printf("  ## WARNING: HSH ELEMENT ALREADY COMPLETE. IGNORED\n");
  }

  return 0;
}

int hash_suppr(HashTable* hsh, int iVer1, int iVer2, int iTri)
{

  int i_hsh;
  i_hsh = hash_find(hsh,iVer1,iVer2); // check if the element is in the hash_list

  if(i_hsh ==0) return 0;
  if(hsh->LstObj[hsh->NbrObj +1][2] ==iTri) hsh->LstObj[hsh->NbrObj +1][2] =0;
  else if(hsh->LstObj[hsh->NbrObj +1][3] ==iTri) hsh->LstObj[hsh->NbrObj +1][3] =0;
  else printf("  ## WARNING: HSH ELEMENT DOESNT HAVE TRIANGLE. IGNORED\n");

  return 0;
}

int hash_out(HashTable* hsh)
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

// double Q1_el(Mesh* Msh, int3d* Tri)
// {
//   double alpha = sqrt(3)/12;
//   double2d P1;
//   double2d P2;
//   double2d P3; 
//   int i1,i2,i3; double a,b,c;
//   i1 = Tri[0];        i2 = Tri[1];        i3 = Tri[2];
//   P1 = Msh->Crd[i1];  P2 = Msh->Crd[i2];  P3 = Msh->Crd[i3];
//   a = (P1[0]-P2[0])*(P1[0]-P2[0]) + (P1[1]-P2[1])*(P1[1]-P2[1]); 
//   b = (P2[0]-P3[0])*(P2[0]-P3[0]) + (P2[1]-P3[1])*(P2[1]-P3[1]); 
//   c = (P3[0]-P1[0])*(P3[0]-P1[0]) + (P3[1]-P1[1])*(P3[1]-P1[1]); 
//   double K_surf = 0.5*((P2[0]-P1[0])*(P3[1]-P1[1]) - (P1[1]-P2[1])*(P3[0]-P1[0]));
//   return (a+b+c)/K_surf;
// }

// double Q2_el(Mesh* Msh, int3d* Tri)
// {
//   double alpha = sqrt(3)/12;
//   double2d P1;
//   double2d P2;
//   double2d P3; 
//   int i1,i2,i3; double a,b,c;
//   i1 = Tri[0]; i2= Tri[1]; i3 = Tri[2];
//   P1 = Msh->Crd[i1];  P2 = Msh->Crd[i2];  P3 = Msh->Crd[i3];
//   a = sqrt((P1[0]-P2[0])*(P1[0]-P2[0]) + (P1[1]-P2[1])*(P1[1]-P2[1])); 
//   b = sqrt((P2[0]-P3[0])*(P2[0]-P3[0]) + (P2[1]-P3[1])*(P2[1]-P3[1])); 
//   c = sqrt((P3[0]-P1[0])*(P3[0]-P1[0]) + (P3[1]-P1[1])*(P3[1]-P1[1])); 
//   double h = max(a,b,c); 
//   double r = a*b*c/sqrt((a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c));
// }
