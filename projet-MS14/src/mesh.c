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

  //printf(" File %s opened Dimension %d Version %d \n", InpFil, Msh->Dim, FilVer);

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
    printf(" %s ",InpFil);
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

  //printf(" File %s opened Dimension %d Version %d \n", InpFil, dim, FilVer);

  SolTyp = GmfSolAtVertices; // read only sol at vertices
  nbrSol = GmfStatKwd(fsol, SolTyp, &NbrTyp, &SolSiz, TypTab);

  if (nbrSol == 0) {
    //printf("  ## WARNING: No SolAtVertices in the solution file !\n");
    return NULL;
  }
  if (dim != mshDim) {
    //printf("  ## WARNING: WRONG DIMENSION NUMBER. IGNORED\n");
    return NULL;
  }
  if (nbrSol != mshNbrSol) {
    //printf("  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    return NULL;
  }
  if (NbrTyp != 1) {
    //printf("  ## WARNING: WRONG FIELD NUMBER. IGNORED\n");
    return NULL;
  }
  if (TypTab[0] != GmfSca) {
    //printf("  ## WARNING: WRONG FIELD TYPE. IGNORED\n");
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
    //printf("  ## ERROR: CANNOT CREATE FILE \n");
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
    if(iTri%5000 == 0){ ti = clock();   //printf("--- task %d / %d full --- %lg (s) passed \n",iTri,Msh->NbrTri,(ti-to)/ CLOCKS_PER_SEC );
      }
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
  // //printf("Trivoi of %d : %d, %d, %d \n", iTri, Msh->TriVoi[iTri][0], Msh->TriVoi[iTri][1], Msh->TriVoi[iTri][2]);
  // }

  return 1;
}

// ============================================================================
// ============================================================================

int msh_neighbors(Mesh* Msh)
{
  int iTri, iEdg, iVer1, iVer2, jVer1, jVer2, jTri, jEdg;

  if (!Msh) return 0;

  if (Msh->TriVoi == NULL)
    Msh->TriVoi = calloc((Msh->NbrTri + 1), sizeof(int3d));

  printf(" creating a neighbour list \n");

  //--- initialize HashTable and set the hash table
  int SizHead = 2*(Msh->NbrVerMax);
  int NbrMaxObj = Msh->NbrVerMax + Msh->NbrTriMax ; // Euler caracteristique with a bit of security

  printf(" max : %d", Msh->NbrTri);

  HashTable* hsh = hash_init(SizHead, NbrMaxObj); 
  for (iTri = 1; iTri <= Msh->NbrTri; iTri++) {
    printf(" \n a %d ",iTri);
    for (iEdg = 0; iEdg < 3; iEdg++) {
      iVer1 = Msh->Tri[iTri][tri2edg[iEdg][0]];
      iVer2 = Msh->Tri[iTri][tri2edg[iEdg][1]];

      int j_hsh = hash_find(hsh,iVer1,iVer2);
      if(j_hsh !=0)
      {
        printf("complete ");
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
        printf("add ");
        hash_add(hsh,iVer1,iVer2,iTri,0);
      }

    }
  }
  printf(" end ");

  free(hsh);
  hsh = NULL;
  // Msh->Hsh = hsh;
  return 1;
}

HashTable* hash_init(int SizHead, int NbrMaxObj)
{
  // HashTable* hsh = NULL;
  HashTable* hsh = malloc(sizeof(HashTable));

  hsh->SizHead = SizHead; 
  hsh->NbrMaxObj = NbrMaxObj;
  hsh->NbrObj =0;

  hsh->LstObj = calloc((hsh->NbrMaxObj+1),sizeof(int5d));
  hsh->Head   = calloc((hsh->SizHead  +1),sizeof(int));

  for(int j=1;j<hsh->SizHead;j++){hsh->Head[j]=0;}

  return hsh;
}

int hash_find(HashTable* hsh, int iVer1, int iVer2)
{
  int key = iVer1 + iVer2;
  int j_hsh = hsh->Head[key];
  int n = 0; // security
  while(j_hsh!=0 && n<hsh->NbrMaxObj)
  {
    // //printf("searching element %d if it has Vert (%d,%d) \n",j_hsh,iVer1,iVer2);
    if(hsh->LstObj[j_hsh][0]==iVer1 && hsh->LstObj[j_hsh][1]==iVer2){return j_hsh;}
    else{ if(hsh->LstObj[j_hsh][1]==iVer1 && hsh->LstObj[j_hsh][0]==iVer2){return j_hsh;}
          else {j_hsh = hsh->LstObj[j_hsh][4];}
    }
  }
  // //printf("j_hsh = 0 \n");
  return 0;

}

int hash_add(HashTable* hsh, int iVer1, int iVer2, int iTri, int i_hsh)
{

  // i_hsh is a initial guess. it should be 0 by default, 
  // but if you have already run hash_find, you can bypass the find to imput the element at i_hsh
  if(i_hsh==0)  i_hsh = hash_find(hsh,iVer1,iVer2); // check if the element isn't in the hash_list already

  if(i_hsh==0)
  {
    if(hsh->NbrObj> hsh->NbrMaxObj+1) //printf("  ## WARNING: HSH ELEMENT ALREADY FULL. IGNORED\n");
    // //printf("adding element %d of Vertex (%d,%d) and Tri (%d,%d) to have next %d \n", hsh->NbrObj, iVer1,iVer2, iTri, hsh->LstObj[hsh->NbrObj][3], hsh->LstObj[hsh->NbrObj][4]);

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
    else printf(" ## WARNING: HSH ELEMENT %d / %d ALREADY COMPLETE (%d,%d); (%d,%d) -> %d. \n IGNORED\n",i_hsh, hsh->NbrMaxObj,hsh->LstObj[i_hsh][0],hsh->LstObj[i_hsh][1],hsh->LstObj[i_hsh][2],hsh->LstObj[i_hsh][3],hsh->LstObj[i_hsh][4] );
    // printf("updating element %d of Vertex (%d,%d) and Tri (%d,%d) to have next %d \n", i_hsh, iVer1,iVer2, hsh->LstObj[hsh->NbrObj][3], iTri, hsh->LstObj[hsh->NbrObj][4]);
  }

  return 0;
}

// NEED TESTING
int hash_suppr(HashTable* hsh, int iVer1, int iVer2, int iTri)  // deletes an element of the hash table
{
  int i_hsh;
  i_hsh = hash_find(hsh,iVer1,iVer2); // check if the element is in the hash_list

  if(i_hsh ==0){ printf("\n DELETING A NON EXISTING ELEMENT, IGNORED \n ");
    }
  if(i_hsh !=0)
  {
    printf(" Deleting element %d : (%d,%d) linked to %d",i_hsh, iVer1,iVer2, iTri);
    printf(" Deleting element %d : (%d,%d) linked to %d",i_hsh, hsh->LstObj[i_hsh][0],hsh->LstObj[i_hsh][1], hsh->LstObj[i_hsh][2]);
    int ToDelete=0;
    if((hsh->LstObj[i_hsh][2]!=0 && hsh->LstObj[i_hsh][3]==0) || (hsh->LstObj[i_hsh][3]!=0 && hsh->LstObj[i_hsh][2]==0)) ToDelete = 1;
    
    if(ToDelete==1)
    {
    printf("deleting the element %d",i_hsh);
    hash_out_index(hsh,i_hsh);

    // we redo the chain
    int key = iVer1 + iVer2; 
    int j_hsh = hsh->Head[key];
    int i_bef = 0;
    while(j_hsh!=0)
    {
      if((hsh->LstObj[j_hsh][0]==iVer1 && hsh->LstObj[j_hsh][1]==iVer2) || (hsh->LstObj[j_hsh][1]==iVer1 && hsh->LstObj[j_hsh][0]==iVer2) ) 
      {
        // sewing it
        if(i_bef==0){hsh->Head[key]       =hsh->LstObj[j_hsh][4]; printf("sewing the element to the head \n");} 
        if(i_bef!=0){hsh->LstObj[i_bef][4]=hsh->LstObj[j_hsh][4]; printf("sewing the element to element %d \n",i_bef);}

        // if(j_hsh<hsh->NbrObj)
        {
          // permuting it with the last element
          int last_index = hsh->NbrObj;        
          for(int i=0;i<5;i++) hsh->LstObj[j_hsh][i] = hsh->LstObj[last_index][i];

          // sewing the last element's chain 
          key = hsh->LstObj[last_index][0] + hsh->LstObj[last_index][1];
          int k_hsh= hsh->Head[key];
          
          int j_bef = 0;
          while(k_hsh!=0)
          {
            printf(" test %d vs %d : (%d,%d) vs (%d,%d)    \n",k_hsh,j_hsh,hsh->LstObj[k_hsh][0],hsh->LstObj[k_hsh][1], hsh->LstObj[j_hsh][0],hsh->LstObj[j_hsh][1]);
            if(hsh->LstObj[k_hsh][0]==hsh->LstObj[j_hsh][0] && hsh->LstObj[k_hsh][1]==hsh->LstObj[j_hsh][1])
            {
              if(j_bef==0){hsh->Head[key]       =hsh->LstObj[k_hsh][4];} 
              if(j_bef!=0){hsh->LstObj[j_bef][4]=hsh->LstObj[k_hsh][4];}
              for(int i=0;i<5;i++) hsh->LstObj[last_index][i]=0;
            }
            else{j_bef = k_hsh; k_hsh = hsh->LstObj[j_bef][4]; printf("chain secondaire : %d -> %d \n",j_bef,k_hsh);}
          }
          hsh->NbrObj -=1;
        }
        
      }
      else{i_bef = j_hsh; j_hsh = hsh->LstObj[j_hsh][4]; printf("chain primaire : %d -> %d \n",i_bef,j_hsh);}      
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
  for(int j=1;j<=hsh->NbrMaxObj;j++){
    hash_out_index(hsh,j);
  }
  printf(" \n ------------------ \n");

  for(int i=1;i<=hsh->SizHead;i++)
  {
    printf("Head %d : %d",i,hsh->Head[i]);
    if(i%5 ==0)printf(" \n");
  }

  return 0;
}

int hash_out_index(HashTable* hsh, int index)  // print the element of the table
{

  printf(" indices = %d ", index);
  printf(" Vertexes  : (%d,%d) ", hsh->LstObj[index][0], hsh->LstObj[index][1]);
  printf(" Triangles : (%d,%d) ", hsh->LstObj[index][2], hsh->LstObj[index][3]);
  printf(" key : %d ",hsh->LstObj[index][0]+ hsh->LstObj[index][1]);
  printf(" next : %d\n", hsh->LstObj[index][4]);
  // printf("----------- \n");
  
  return 0;
}

int hash_bound(HashTable* hsh)
{
  int Nb_bound =0;

  for(int i_hsh=1; i_hsh <= hsh->NbrObj; i_hsh++ )
  {  
    // //printf("Triangles : %d, %d ", hsh->LstObj[i_hsh][2], hsh->LstObj[i_hsh][3]);
    if(hsh->LstObj[i_hsh][3] == 0)
    {      
      // //printf(" : edge ");
      Nb_bound +=1;
    }
    // //printf(" \n ======= \n");
  }
  //printf("Number of boundary edges : %d \n", Nb_bound);
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

int Check_hash(Mesh* Msh) // check if the hashtable is correct for the mesh.
{
  HashTable* hsh = Msh->Hsh;
  // int Headed=0;     // the collision links are on the same key
  // int Neighbours=0; // the edge are connected to the right triangle
  // int ancestry=0;   // there are no branch
  
  // check the keys
  int i_obj=0;
  int iVer1,iVer2,iTri1,iTri2;
  for(int key=1;key<hsh->SizHead;key++)
  {
    printf("check key %d \n",key);
    i_obj = hsh->Head[key];
    while(i_obj !=0)
    {
      hash_out_index(hsh,i_obj);
      iVer1 = hsh->LstObj[i_obj][0];
      iVer2 = hsh->LstObj[i_obj][1];
      if(iVer1 + iVer2 != key){printf("The element %d is associated to the wrong key \n",i_obj);}

      i_obj = hsh->LstObj[i_obj][4];
    }
  }
  
  // check the Neighbours
  int loc_neigh1=0, loc_neigh2 =0;
  for(int key=0;key<hsh->SizHead;key++)
  {
    i_obj = hsh->Head[key];
    while(i_obj !=0)
    {
      iVer1 = hsh->LstObj[i_obj][0];
      iVer2 = hsh->LstObj[i_obj][1];
      iTri1 = hsh->LstObj[i_obj][2];
      iTri2 = hsh->LstObj[i_obj][3];

      // hash_out_index(hsh,i_obj);

      if(Msh->Tri[iTri1][0]==iVer1 && Msh->Tri[iTri1][1]==iVer2) loc_neigh1 =1;
      if(Msh->Tri[iTri1][1]==iVer1 && Msh->Tri[iTri1][0]==iVer2) loc_neigh1 =1;
      if(Msh->Tri[iTri1][2]==iVer1 && Msh->Tri[iTri1][1]==iVer2) loc_neigh1 =1;
      if(Msh->Tri[iTri1][1]==iVer1 && Msh->Tri[iTri1][2]==iVer2) loc_neigh1 =1;
      if(Msh->Tri[iTri1][0]==iVer1 && Msh->Tri[iTri1][2]==iVer2) loc_neigh1 =1;
      if(Msh->Tri[iTri1][2]==iVer1 && Msh->Tri[iTri1][0]==iVer2) loc_neigh1 =1;
      
      // printf("Tri : %d,  (%d,%d,%d), ref : %d \n", iTri1, Msh->Tri[iTri1][0], Msh->Tri[iTri1][1], Msh->Tri[iTri1][2],Msh->TriRef[iTri1]);
      
      if(Msh->Tri[iTri2][0]==iVer1 && Msh->Tri[iTri2][1]==iVer2) loc_neigh2 =1;
      if(Msh->Tri[iTri2][1]==iVer1 && Msh->Tri[iTri2][0]==iVer2) loc_neigh2 =1;
      if(Msh->Tri[iTri2][2]==iVer1 && Msh->Tri[iTri2][1]==iVer2) loc_neigh2 =1;
      if(Msh->Tri[iTri2][1]==iVer1 && Msh->Tri[iTri2][2]==iVer2) loc_neigh2 =1;
      if(Msh->Tri[iTri2][0]==iVer1 && Msh->Tri[iTri2][2]==iVer2) loc_neigh2 =1;
      if(Msh->Tri[iTri2][2]==iVer1 && Msh->Tri[iTri2][0]==iVer2) loc_neigh2 =1;  
      // printf("Tri : %d,  (%d,%d,%d), ref : %d \n", iTri2, Msh->Tri[iTri2][0], Msh->Tri[iTri2][1], Msh->Tri[iTri2][2],Msh->TriRef[iTri2]);

      if((loc_neigh1!=1 && iTri1!=0 )|| (loc_neigh2 !=1 && iTri2!=0) ){printf("The element %d is associated to the wrong triangles \n",i_obj);}
      
      i_obj = hsh->LstObj[i_obj][4];
    }
  }
  return 0;
}

// ============================================================================
// ============================================================================

int msh_write2dfield_Vertices(char* file, int nfield, double* field)
{
  int iVer;

  int fmsh = GmfOpenMesh(file, GmfWrite, GmfDouble, 2);
  if (fmsh <= 0) {
    //printf("  ## ERROR: CANNOT CREATE FILE \n");
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
    //printf("  ## ERROR: CANNOT CREATE FILE \n");
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
    //printf("  ## ERROR: CANNOT CREATE FILE \n");
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

  for (jEdg = 0; jEdg < Msh->NbrEfrMax; jEdg++)  // NbrEfrMax <- non def
  {
    jVer1 = Msh->Efr[jEdg][0];
    jVer2 = Msh->Efr[jEdg][1];

    if(((iVer1== jVer1) && (iVer2 ==jVer2)) || ((iVer1 == jVer2) && (iVer2 == jVer1))) is_valid = 0;
  }

  return is_valid;
}

double* connex_comp(Mesh* Msh)
{
  //printf("creating neighbors list of the mesh \n");
  msh_neighbors(Msh) ;
  //printf("neighbors list of the mesh done \n");
  //printf("creating connex composante of the mesh \n");

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
  //printf("sub domains created. there is %d sub domaines \n",color-1);
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
                                                    P1[0],P1[1], P2[0],P2[1], P3[0],P3[1]); return 0;}
  
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

    return (d11*d22*d33 + d12*d23*d31 + d13*d21*d32 - d13*d22*d31 - d23*d32*d11 - d33*d12*d21 > 1e-10);
}

double* Coord_bary(double2d Point, Mesh* Msh, int iTri)
{
  
  int i1,i2,i3;
  double K;
  
  double* crd = calloc(3,sizeof(double)); crd[0]= 0; crd[1]=0; crd[2]=0;

  double* P1; double* P2; double* P3;
  i1 = Msh->Tri[iTri][0];        i2 = Msh->Tri[iTri][1];        i3 = Msh->Tri[iTri][2];
  P1 = (double2d){Msh->Crd[i1][0], Msh->Crd[i1][1] };
  P2 = (double2d){Msh->Crd[i2][0], Msh->Crd[i2][1] };
  P3 = (double2d){Msh->Crd[i3][0], Msh->Crd[i3][1] };
  K =(surf(P1,P2,P3));
  crd[0] = surf(Point,P2,P3)/K; crd[1] = surf(P1,Point, P3)/K; crd[2] = surf(P1,P2,Point)/K;

  
  // free(P1); P1=NULL;
  // free(P2); P2=NULL;
  // free(P3); P3=NULL;

  return crd;

}

// localise_Tri 
// NULL                         if on point <-- CHANGED
// iTri, edge_Tri, 0, on_edge;  if on edge 
// iTri, 0, in_Trig, 0;         if in Tri

int* Localise_Tri(Mesh* Msh, double2d Point)
{

  // printf(" --------------------------------- \n");
  int in_trig = 0; // 0 if in the triangle, 1 if it is. 
  int on_edge = 0; // 0 if not on edge, 1 if it is.
  int compass; // check if the north( direction ) has been chosen
  int iTri =1, edge_Tri=0;
  double ax1,ax2,ax3;

  // edge cases 

  // naive non convex mesh
  int edge_case =1;           

  // Point is on Point;
  int on_Point = 0; int Ver=0, jTri=0;
  
  while(in_trig ==0 && on_edge ==0)
  {
    compass =0;
    double* coord = Coord_bary(Point,Msh,iTri);
    ax1 = coord[0]; ax2 = coord[1]; ax3= coord[2];

    // printf("In the triangle %d, with neighbours (%d,%d,%d) the values are (%10f,%10f,%10f) \n",iTri,Msh->TriVoi[iTri][0],Msh->TriVoi[iTri][1],Msh->TriVoi[iTri][2],ax1,ax2,ax3);
    //check
    if(ax1<0 && (ax2<0 && ax3<0) ) //printf("\n ERROR POINT OR TRIANGLE WRONGLY DEFINED, IGNORED \n");

    // 3 direct neighbours
    // printf(" from the directions available : (%d,%d,%d) \n", Msh->TriVoi[iTri][0],Msh->TriVoi[iTri][1],Msh->TriVoi[iTri][2]);
    

    // we could change the choice of the direction
    // here we have a bias against the direction 1 

    if(fabs(ax1-1)<1e-30 && (fabs(ax2)<1e-30 && fabs(ax3)<1e-30) ){ on_Point=1; Ver=Msh->Tri[iTri][0]; break;} 
    if(fabs(ax2-1)<1e-30 && (fabs(ax1)<1e-30 && fabs(ax3)<1e-30) ){ on_Point=1; Ver=Msh->Tri[iTri][1]; break;} 
    if(fabs(ax3-1)<1e-30 && (fabs(ax1)<1e-30 && fabs(ax2)<1e-30) ){ on_Point=1; Ver=Msh->Tri[iTri][2]; break;}

    if(fabs(ax1)<1e-30 && (ax2>0 && ax3>0)){ on_edge =1; edge_Tri = Msh->TriVoi[iTri][0]; compass =1; break;}
    if(fabs(ax2)<1e-30 && (ax1>0 && ax3>0)){ on_edge =1; edge_Tri = Msh->TriVoi[iTri][1]; compass =1; break;}
    if(fabs(ax3)<1e-30 && (ax1>0 && ax2>0)){ on_edge =1; edge_Tri = Msh->TriVoi[iTri][2]; compass =1; break;}

    if(ax1<0 && compass ==0){ iTri = Msh->TriVoi[iTri][0]; compass =1;}
    if(ax2<0 && compass ==0){ iTri = Msh->TriVoi[iTri][1]; compass =1;}
    if(ax3<0 && compass ==0){ iTri = Msh->TriVoi[iTri][2]; compass =1;}

    // printf(" going in direction : (%d,%d,%d) towards : %d  \n", ax1<0, ax2<0, ax3<0, iTri);

    if((ax1>0 && (ax2>0 && ax3>0)) && on_edge ==0 ){ in_trig=1;}
    
    if(iTri == 0){edge_case +=1; iTri = edge_case;}

  }

  // returns an array of the triangle linked to the point
  if(on_Point)
  {
    // printf(" ------ Ver : %d ----------- \n", Ver);
    // return NULL;
    int* domain = calloc(2+20,sizeof(int)); // how many triangle could a point possibly have ?
    domain[0] =0; domain[1]=0;
    for(int iEdg=0;iEdg<3;iEdg++) 
    {
      if(Msh->Tri[iTri][tri2edg[iEdg][0]] == Ver)       jTri = Msh->TriVoi[iTri][tri2edg[iEdg][1]]; 
      else if (Msh->Tri[iTri][tri2edg[iEdg][1]] == Ver) jTri = Msh->TriVoi[iTri][tri2edg[iEdg][0]];
    }
    
    // printf(" iTri : %d (%d,%d,%d) \n",iTri, Msh->Tri[iTri][0],  Msh->Tri[iTri][1],  Msh->Tri[iTri][2]);
      
    int sizeof_domain =2;
    while((jTri != iTri && jTri !=0) && sizeof_domain<22 )
    {
      // printf(" jTri : %d (%d,%d,%d) \n",jTri, Msh->Tri[jTri][0],  Msh->Tri[jTri][1],  Msh->Tri[jTri][2]);
      
      for(int iEdg=0;iEdg<3;iEdg++) 
      {
        if(Msh->Tri[jTri][tri2edg[iEdg][0]] == Ver)       jTri = Msh->TriVoi[iTri][tri2edg[iEdg][1]]; 
        else if (Msh->Tri[jTri][tri2edg[iEdg][1]] == Ver) jTri = Msh->TriVoi[iTri][tri2edg[iEdg][0]];
      }
      // printf(" to %d \n",jTri);
      domain[sizeof_domain] = jTri; sizeof_domain++;

    

    }
    return domain;
  }

  // border case 
  if(on_edge && iTri == 0){on_edge=0; iTri = edge_Tri; in_trig =1;}
  if(on_edge && edge_Tri == 0){on_edge=0; in_trig=1;}

  // if(on_edge) printf("The Point has been located on an edge between %d and %d \n", iTri,edge_Tri);
  // if(in_trig) printf("The Point has been located in triangle %d \n",iTri);

  int* Result = (int*)calloc(4,sizeof(int));
  Result[0]= in_trig;   Result[1]= on_edge;
  Result[2]= iTri;    Result[3] = edge_Tri;

  return Result;

} 

int sizeof_cavity = 0; // size of the cavity/index of the last element added to the cavity 
int sizeof_boundary =0; // var global

int* digging_a_hole(Mesh* Msh, double2d Point) //diggy diggy hole
{

  int id_tri,edge_Tri,in_trig,on_edge;
  int* Result= Localise_Tri(Msh, Point);

  if(Result == NULL || (Result[0]==0 && Result[1]==0)) return NULL;
  in_trig = Result[0]; on_edge = Result[1];
  id_tri = Result[2]; edge_Tri = Result[3]; 


  free(Result); Result = NULL;

  // calculting the elements to add to the cavity
  int i1,i2,i3;
  int* pile   = calloc(Msh->NbrTri +1, sizeof(int));  // pile of element to add to the mavity field 
  int* cavity = calloc(Msh->NbrTri +1, sizeof(int));  // elements to remove
  int  sizeof_pile = 1;                               // size of the pile / index of the top of the pile                        
 
  int Is_In_Pile=0, Is_In_Cavity=0;

  int iTri = id_tri; int jTri;

  // initializing the pile
  if(in_trig ==1)  pile[1] = id_tri;
  if(on_edge ==1){ pile[1] = id_tri; pile[2] = edge_Tri; sizeof_pile +=1;}

  while(sizeof_pile>0)
  {
    // add the element to the cavity and remove it from the pile
    iTri = pile[sizeof_pile]; sizeof_pile -=1;
    cavity[sizeof_cavity] = iTri; sizeof_cavity +=1;

    // printf("adding tri %d to the pile \n",iTri);

    // check if the neighbours are in the cavity
    for(int jEdg=0;jEdg<3; jEdg++) 
    {
      Is_In_Cavity =0; Is_In_Pile =0;
      jTri = Msh->TriVoi[iTri][jEdg];

      i1 = Msh->Tri[jTri][0];        i2 = Msh->Tri[jTri][1];        i3 = Msh->Tri[jTri][2];
      double2d P1 = {Msh->Crd[i1][0], Msh->Crd[i1][1]};
      double2d P2 = {Msh->Crd[i2][0], Msh->Crd[i2][1]};
      double2d P3 = {Msh->Crd[i3][0], Msh->Crd[i3][1]};

      // check if the element hasn't already been added
      for(int i=0;i<=sizeof_cavity; i++){ if(cavity[i]==jTri ){ Is_In_Cavity = 1; break;} }               
      for(int i=0;i<=sizeof_pile  ; i++){ if(pile[i]  ==jTri ){ Is_In_Pile   = 1; break;} }

      if(Is_In_Cavity==0 && (Is_In_Pile==0 && jTri !=0))
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
  
  free(pile); pile =NULL;
  return cavity;
}

int deleting_Cavity(Mesh* Msh, int* cavity)  // deleting the elements of the cavity
{
  int Tri_neigh;
  int Tri_cavity;
  
  for(int i_cavity=0; i_cavity<sizeof_cavity; i_cavity ++)
  {
    Tri_cavity = cavity[i_cavity];
    // deleting it from it's neighbours's list of neighbours
    for(int i=0; i<3; i++)
    {
      Tri_neigh =Msh->TriVoi[Tri_cavity][i];// check all neighbours 
      for(int j=0; j<3; j++)      
        if(Msh->TriVoi[Tri_neigh][j]==Tri_cavity) Msh->TriVoi[Tri_neigh][j] = 0; // removing it from their neighbours
    } 
  } 

  return 0;
}

int3d* boundary_hole(Mesh* Msh, int* cavity)
{
  int3d* boundary_cavity; 
  boundary_cavity = calloc((Msh->NbrTri+ Msh->NbrVer) +1 , sizeof(int3d));

  int Is_On_Boundary=0;
  int iVer1,iVer2;

  for(int i_cavity=0; i_cavity<sizeof_cavity; i_cavity ++)  // on check tout les triangles de la cavité
  {
    for(int i_Edg=0;i_Edg<3;i_Edg++)                        // toutes les arrêtes
    {
      iVer1 = Msh->Tri[cavity[i_cavity]][tri2edg[i_Edg][0]];
      iVer2 = Msh->Tri[cavity[i_cavity]][tri2edg[i_Edg][1]];
      // i_obj=hash_find(Msh->Hsh,iVer1,iVer2);
      Is_On_Boundary =0;

      for(int i=0;i<sizeof_boundary;i++)                    // On regarde  
      {
        if((boundary_cavity[i][0] ==iVer1 && boundary_cavity[i][1]==iVer2 ) || (boundary_cavity[i][1] ==iVer1 && boundary_cavity[i][0]==iVer2) )
        { 
          // if it is in the boundary -> remove it  
          Is_On_Boundary= 1;
          boundary_cavity[i][0]=boundary_cavity[sizeof_boundary][0];
          boundary_cavity[i][1]=boundary_cavity[sizeof_boundary][1];
          boundary_cavity[i][2]=boundary_cavity[sizeof_boundary][2];
          // printf(" removing boundary : (%d,%d) linked to %d \n", boundary_cavity[sizeof_boundary][0],boundary_cavity[sizeof_boundary][1],boundary_cavity[sizeof_boundary][2]);
          break; 
        }
      }     
        
      // if it is not in the boundary -> adding it
      if(Is_On_Boundary==0)
      {
      boundary_cavity[sizeof_boundary][0]=iVer1;
      boundary_cavity[sizeof_boundary][1]=iVer2;
      boundary_cavity[sizeof_boundary][2]=Msh->TriVoi[cavity[i_cavity]][i_Edg];
      // printf(" adding boundary : (%d,%d) linked to %d \n", iVer1,iVer2,Msh->TriVoi[cavity[i_cavity]][i_Edg]);
      sizeof_boundary +=1;
      }
    } 
  } 
  
  int NumOfZero=0;
  for(int i =0;i<sizeof_boundary-NumOfZero;i++)
  {
    if(boundary_cavity[i][0]==0)
    {
      boundary_cavity[i][0] = boundary_cavity[sizeof_boundary-NumOfZero-1][0];  boundary_cavity[sizeof_boundary-NumOfZero-1][0]=0;
      boundary_cavity[i][1] = boundary_cavity[sizeof_boundary-NumOfZero-1][1];  boundary_cavity[sizeof_boundary-NumOfZero-1][1]=0;
      boundary_cavity[i][2] = boundary_cavity[sizeof_boundary-NumOfZero-1][2];  boundary_cavity[sizeof_boundary-NumOfZero-1][2]=0;
      NumOfZero +=1;
    }
  }
  sizeof_boundary -= NumOfZero;
  
  return boundary_cavity;
}

int memory_allocation_start(Mesh* Msh, int esti_Ver, int esti_Tri)
{
  int security = 1;
  Msh->NbrTriMax = Msh->NbrTri + esti_Tri;
  Msh->NbrVerMax = Msh->NbrVer + esti_Ver;

  Msh->Tri = realloc(Msh->Tri, (Msh->NbrTriMax+1  +security)*sizeof(int3d));
  if (Msh->Tri == NULL) {
    // If reallocation fails
    printf("ERROR. Unable to resize memory of Tri \n");
    exit(0);
  }
  
  Msh->TriRef = realloc(Msh->TriRef, (Msh->NbrTriMax+1 +security)*sizeof(int1d));
  if (Msh->TriRef == NULL) {
    // If reallocation fails
    printf("ERROR. Unable to resize memory of Tri_Ref \n");
  }

  Msh->TriVoi = realloc(Msh->TriVoi, (Msh->NbrTriMax+1 +security)*sizeof(int3d));
  if (Msh->TriVoi == NULL) {
    // If reallocation fails
    printf("ERROR. Unable to resize memory of Tri_Voi \n");
  } 
  
  Msh->Crd = realloc(Msh->Crd , (Msh->NbrVerMax+1 +security)*sizeof(double2d));
  if (Msh->Crd == NULL) {
    // If reallocation fails
    printf("ERROR. Unable to resize memory of Crd \n");
  } 

  return 0;
}

int memory_allocation_point(Mesh* Msh)
{

  int security =0;
  Msh->NbrTriMax = Msh->NbrTriMax + sizeof_boundary - sizeof_cavity ;
  Msh->NbrTri = Msh->NbrTri - sizeof_cavity;


  Msh->Tri = realloc(Msh->Tri, (Msh->NbrTriMax+1  +security)*sizeof(int3d));
  if (Msh->Tri == NULL) {
    // If reallocation fails
    printf("ERROR. Unable to resize memory of Tri \n");
    exit(0);
  }

  Msh->TriRef = realloc(Msh->TriRef, (Msh->NbrTriMax+1 +security)*sizeof(int1d));
  if (Msh->TriRef == NULL) {
    // If reallocation fails
    printf("ERROR. Unable to resize memory of Tri_Ref \n");
  }

  Msh->TriVoi = realloc(Msh->TriVoi, (Msh->NbrTriMax+1 +security)*sizeof(int3d));
  if (Msh->TriVoi == NULL) {
    // If reallocation fails
    printf("ERROR. Unable to resize memory of Tri_Voi \n");
  } 

  Msh->NbrVer +=1; Msh->NbrVerMax +=1 ;
  
  Msh->Crd = realloc(Msh->Crd , (Msh->NbrVerMax+1 +security)*sizeof(double2d));
  if (Msh->Crd == NULL) {
    // If reallocation fails
    printf("ERROR. Unable to resize memory of Crd \n");
  } 

  return 0;
}

int Starring_Point(Mesh* Msh, int* cavity, int3d* boundary, double2d Point)
{
  int Tri_added, Tri_bound;
  int iVer1,iVer2,jVer1,jVer2;
  double surface;
  double* P1 =Msh->Crd[1]; double* P2=Msh->Crd[1];

  // starring the point to the cavity's boundary
  for(int i_boundary=0; i_boundary<sizeof_boundary; i_boundary++)
  { 
    iVer1 = boundary[i_boundary][0]; 
    iVer2 = boundary[i_boundary][1]; 
    P1 = Msh->Crd[iVer1]; P2 = Msh->Crd[iVer2];

    while( ((boundary[i_boundary][0]==0 || boundary[i_boundary][1]==0) && i_boundary<sizeof_boundary) || fabs(surf(P1,P2,Point)) <1e-30)
    {
    i_boundary +=1;
    iVer1 = boundary[i_boundary][0]; 
    iVer2 = boundary[i_boundary][1]; 
    P1 = Msh->Crd[iVer1]; P2 = Msh->Crd[iVer2];
    
    if(i_boundary > sizeof_boundary) break;
    }


    if(i_boundary > sizeof_boundary) break;

    if(sizeof_cavity>0){Tri_added = cavity[sizeof_cavity-1];}
    else{Tri_added = Msh->NbrTri+1; cavity[i_boundary]= Tri_added;}
    
    // add to Tri
    iVer1 = boundary[i_boundary][0]; 
    iVer2 = boundary[i_boundary][1]; 
    surface = surf(Msh->Crd[iVer1],Msh->Crd[iVer2],Point);
    if(surface<0)
    {
    Msh->Tri[Tri_added][0]= iVer2;
    Msh->Tri[Tri_added][1]= iVer1;
    Msh->Tri[Tri_added][2]= Msh->NbrVer;
    Msh->TriRef[Tri_added] =1;
    Msh->TriVoi[Tri_added][0] =0; Msh->TriVoi[Tri_added][1] =0;
    }
    else
    {      
      Msh->Tri[Tri_added][0]= iVer1;
      Msh->Tri[Tri_added][1]= iVer2;
      Msh->Tri[Tri_added][2]= Msh->NbrVer;
    Msh->TriRef[Tri_added] =1;
    Msh->TriVoi[Tri_added][0] =0; Msh->TriVoi[Tri_added][1] =0;
    }

    // linking the edge of the cavity to the outside
    Tri_bound = boundary[i_boundary][2];
    Msh->TriVoi[Tri_added][2] =Tri_bound;
    for(int iEdg=0;iEdg<3;iEdg++)
    {
    jVer1 = Msh->Tri[Tri_bound][tri2edg[iEdg][0]];
    jVer2 = Msh->Tri[Tri_bound][tri2edg[iEdg][1]];
    if((iVer1==jVer1 && iVer2 == jVer2) || (iVer2 == jVer1 && iVer1 == jVer2)){Msh->TriVoi[Tri_bound][iEdg]=Tri_added;}
    }

    sizeof_cavity -=1; Msh->NbrTri +=1;
  }
  
  // free(P1); P1 = NULL;
  // free(P2); P2 = NULL;

  // updating the neighbours of the star using a quadratic algorithm.
  int iTri,jTri;
  for(int i_boundary=0;i_boundary<sizeof_boundary;i_boundary++)
  {
    // or just stock it inside an array
    iTri = cavity[i_boundary];
    iVer1 = Msh->Tri[iTri][0];
    iVer2 = Msh->Tri[iTri][1];
    // check which Tri has a point in common with iTri
    for(int j_bound=0;j_bound<sizeof_boundary;j_bound++)
    {
      
      jTri = cavity[j_bound];
      
      jVer1 = Msh->Tri[jTri][tri2edg[2][0]];
      jVer2 = Msh->Tri[jTri][tri2edg[2][1]];
      if(iTri != jTri )
      {
      if(jVer2 == iVer1){ Msh->TriVoi[jTri][0] = iTri; Msh->TriVoi[iTri][1] = jTri; }
      if(iVer2 == jVer1){ Msh->TriVoi[iTri][0] = jTri; Msh->TriVoi[jTri][1] = iTri; }
      if(jVer1 == iVer1 || jVer2 == iVer2) printf("ERROR TRIANGLE NOT DEF CORRECTLY \n");
      }
    }
  }
  return 0;
}

int ajout_point(Mesh* Msh, double2d Point)
{
  sizeof_cavity=0; sizeof_boundary =0;

  int* cavity = digging_a_hole(Msh,Point);
  if(cavity == NULL){ printf("point already there \n"); return 0;}

  deleting_Cavity(Msh,cavity);

  int3d* boundary = boundary_hole(Msh, cavity);
  

  // if((Msh->NbrTri + sizeof_boundary - sizeof_cavity +1> Msh->NbrTriMax ) || (Msh->NbrVer +1 > Msh->NbrVerMax))
  memory_allocation_point(Msh);
  
  Msh->Crd[Msh->NbrVer][0] = Point[0];
  Msh->Crd[Msh->NbrVer][1] = Point[1];

  Starring_Point(Msh,cavity,boundary,Point);

  free(cavity);
  free(boundary);
  cavity =NULL; boundary =NULL;
  return 1;
}

/// ------------------------------------------------------------

Mesh* Maillage_Delaunay(int Nb_Point, Mesh* Msh)
{
  printf("Creating a Mesh based on Delaunay method \n");
  srand(10);

  double Crd_x,Crd_y;
  msh_neighbors(Msh);

  for(int it=0;it<=Nb_Point;it++)
  {
    printf("------------ adding point %d to the mesh ----------- \n",it+1);
    Crd_x=(double)(rand())/RAND_MAX;
    while(Crd_x<1e-2 || Crd_x>1-1e-2) Crd_x=(double)(rand())/RAND_MAX;
    Crd_y=(double)(rand())/RAND_MAX;
    while(Crd_y<1e-2 || Crd_y>1-1e-2) Crd_y=(double)(rand())/RAND_MAX;

    ajout_point(Msh, (double2d){Crd_x,Crd_y} );

    printf(" point added \n ");

   
  }

  double regul = (double)(Nb_Point)/(double)64;
  // if we have enough points, add an edge point
  for(int jt=0;jt+1<regul; jt++)
  {
    printf("------------ adding point %d to the mesh's edges  ----------- \n",4*jt+1);
    ajout_point(Msh, (double2d){(double)(jt+1)/regul ,0});
    printf("------------ adding point %d to the mesh's edges  ----------- \n",4*jt+2);
    ajout_point(Msh, (double2d){(double)(jt+1)/regul,1});
    printf("------------ adding point %d to the mesh's edges  ----------- \n",4*jt+3);
    ajout_point(Msh, (double2d){0,(double)(jt+1)/regul}); 
    printf("------------ adding point %d to the mesh's edges  ----------- \n",4*jt+4);
    ajout_point(Msh, (double2d){1,(double)(jt+1)/regul});
  }


  // printf("------------ adding point %d to the mesh's edges  ----------- \n",1);
  // ajout_point(Msh, (double2d){0.5 ,0});
  printf(" mesh completed \n");
  

  msh_write(Msh,"Delaunay.mesh");
  return Msh;
  }

int Mesh_out(Mesh* Msh)
{
  printf("The mesh is of dimension %d \n", Msh->Dim);

  // printf("The mesh is in a box of x: (%f,%f) y:(%f,%f) \n", Msh->Box[0],Msh->Box[1],Msh->Box[2],Msh->Box[3]);



  printf("The mesh is made out of \n");
  printf("Vertexes : \n");
  for(int i=1;i<=Msh->NbrVerMax;i++) printf("Vert : %d, (%f,%f) \n", i,Msh->Crd[i][0], Msh->Crd[i][1]);
  printf(" \n ============== \n");
  printf("Triangles : \n");
  for(int i=1;i<=Msh->NbrTriMax;i++) printf("Tri : %d,  (%d,%d,%d), ref : %d, Voi :(%d,%d,%d) \n", i, Msh->Tri[i][0], Msh->Tri[i][1], Msh->Tri[i][2],Msh->TriRef[i], Msh->TriVoi[i][0],Msh->TriVoi[i][1],Msh->TriVoi[i][2]);
  printf(" \n ============== \n");

  printf(" there are in total : %d Triangle and a maximum of %d \n", Msh->NbrTri, Msh->NbrTriMax);
  printf(" there are in total : %d Vertex   and a maximum of %d \n", Msh->NbrVer  , Msh->NbrVerMax);
  // printf(" there are in total : %d Edges    and a maximum of %d \n", Msh->Hsh->NbrObj, Msh->Hsh->NbrMaxObj);

  return 0;

}

int Sort_Mesh(Mesh* Msh)
{
  for(int iTri=0; iTri<Msh->NbrTriMax;iTri ++)
  {
    if(Msh->Tri[iTri][0] ==0 || (Msh->Tri[iTri][1] ==0 || Msh->Tri[iTri][2] ==0)) 
    {

    }
  }

  return 0;
}

int Compression(Mesh* Msh, double* sol, double factor)
{
  srand(8);
  msh_boundingbox(Msh);

  Mesh* Msh_red = msh_read("../data/carre_base.mesh", 0);
  double* sol_red = calloc(Msh->NbrVerMax+1,sizeof(double));

  printf("bounding box : x (%10f,%10f), y (%10f,%10f) \n", Msh->Box[0],Msh->Box[1],Msh->Box[2],Msh->Box[3] );
  
  Msh_red->Crd[1][0] = Msh->Box[0]; Msh_red->Crd[1][1] = Msh->Box[2];
  Msh_red->Crd[2][0] = Msh->Box[0]; Msh_red->Crd[2][1] = Msh->Box[3];
  Msh_red->Crd[3][0] = Msh->Box[1]; Msh_red->Crd[3][1] = Msh->Box[2];
  Msh_red->Crd[4][0] = Msh->Box[1]; Msh_red->Crd[4][1] = Msh->Box[3];

  msh_neighbors(Msh_red);

  int jVer=1, newpoint=0;
  for(int iVer=1;iVer<Msh->NbrVer; iVer++) 
  {  if(iVer%1000 ==0) printf("-------- work done %d/%d : %10f ------------ \n", iVer,Msh->NbrVer, (double)iVer/(double)Msh->NbrVer);
    if((double)(rand())/RAND_MAX<factor)
    {
      // printf(" adding a point \n");
      newpoint =ajout_point(Msh_red, Msh->Crd[iVer]);

      if(newpoint ==1 && sol !=NULL)
      {
      sol_red[jVer] = sol[iVer];
      jVer+=1;
      }
    }
  }
  printf("Reduction done \n");
  // output mesh compressed and rough sol
  msh_write2dfield_Vertices("Compression.sol", Msh_red->NbrVer, sol_red);
  msh_write(Msh_red,"Compression.mesh");

  double* sol_interpol = Interpol_sol(Msh,Msh_red,sol,sol_red);

  // calcul qc : 
  double qc=0;
  for(int i=0;i<Msh->NbrVer;i++) qc+= fabs(sol_interpol[i]- sol[i]);
  printf(" qc = %10f \n", qc/Msh->NbrVer);

  // output mesh base et interpole
  msh_write(Msh,"Base_mesh.mesh");
  msh_write2dfield_Vertices("Compressed_sol.sol", Msh->NbrVer, sol_interpol);

  printf("Mesh written \n");
  return 0;
}

double* Interpol_sol(Mesh* Msh, Mesh* Msh_red, double* sol, double* sol_red)
{
  int iTri;
  double* Point; double* Crd;
  double* sol_interpol = calloc(Msh->NbrVer , sizeof(double));
  for(int iVer=1;iVer<Msh->NbrVer;iVer++)
  {

    Point = Msh->Crd[iVer];                   // take a point of the base mesh
    int* Result =Localise_Tri(Msh_red,Point); // locate it in the compressed mesh

    if((Result[0]!=0 && Result[1] !=0)&& Result != NULL) 
    {
      printf(" sol in triangle or edge \n");
      iTri = Result[2];             
      Crd = Coord_bary(Point,Msh_red,iTri);       
      for(int i=0;i<3;i++)  sol_interpol[iVer] += Crd[i]*sol_red[Msh_red->Tri[iTri][i]];
      
      if(Result[1]!=0)
      {
        iTri = Result[3];             
        Crd = Coord_bary(Point,Msh_red,iTri);   
        for(int i=0;i<3;i++)  sol_interpol[iVer] += Crd[i]*sol_red[Msh_red->Tri[iTri][i]];
      }

    }

    if((Result[0]==0 && Result[1] ==0)&& Result != NULL)
    {
      // printf(" sol on vertex ");
      int Ver=0;
      iTri= Result[2]; int jTri = Result[3];
      // printf(" shared with %d and %d ",iTri,jTri);
      for(int i=0;i<3;i++) for(int j=0;j<3;j++) if(Msh_red->Tri[iTri][i] == Msh_red->Tri[jTri][j]) Ver = Msh_red->Tri[iTri][i];
      
      sol_interpol[iVer] += sol_red[Ver];
    }
    
    
  }
  return sol_interpol;
}

// ============================================================================
// ============================================================================