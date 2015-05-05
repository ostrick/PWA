#include "Parse_HELI.h"

//-----------------------------------------------------------------------------

void Parse_HELI(Char_t* Path)
{
  Char_t Buffer[1024];
  Double_t Energy, Re, Im,nz,Egamma;
  Double_t t, z, ReH1, ImH1, ReH2, ImH2, ReH3, ImH3, ReH4, ImH4;
  FILE* HELI_H1[LBINS];
  FILE* HELI_H2[LBINS];
  FILE* HELI_H3[LBINS];
  FILE* HELI_H4[LBINS];
  FILE* HELI;
  
  EPSILON = 1e38; //Initialize for search

  /* 
  printf("------------------------------------------------------------------------------------\n");
  printf("Loading Hedim's helicity aplitudes... ");

  //Open file
  sprintf(Buffer, "%s/Helcon_SE.dat", Path); HELI = fopen(Buffer, "r");
  
  //Load Helicity amplitudes
  heli_bin = 0;
   
  while(!feof(HELI))
  {
    if(fscanf(HELI, "%lf %lf %lf\n", &Energy, &Egamma, &nz)!=3) break;
    heli_nz[heli_bin] = nz;
    for(Int_t i=0; i<nz; i++)
    {
     fscanf(HELI, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &t, &z, &ReH1, &ImH1, &ReH2, &ImH2, &ReH3, &ImH3, &ReH4, &ImH4);
      
     heli_H1[i][heli_bin] = TComplex(ReH1/1.41382, ImH1/1.41382);
     heli_H2[i][heli_bin] = TComplex(ReH2/1.41382, ImH2/1.41382);
     heli_H3[i][heli_bin] = TComplex(ReH3/1.41382, ImH3/1.41382);
     heli_H4[i][heli_bin] = TComplex(ReH4/1.41382, ImH4/1.41382);
    
     heli_z[i][heli_bin] = z;
     heli_en[heli_bin] = Egamma*1000.;
     
    }
    printf("%d, %d, %lf \n", heli_bin, heli_nz[heli_bin],  heli_en[heli_bin]);
    heli_bin++;
  }
  
  */ 

  printf("------------------------------------------------------------------------------------\n");
  printf("Loading Sven's helicity amplitudes... ");

  //Open file
  sprintf(Buffer, "%s/Heli_Sven.dat", Path); HELI = fopen(Buffer, "r");
  
  //Load Helicity amplitudes
  heli_bin = 0;
  nz = 40; 
  while(!feof(HELI))
  {
    heli_nz[heli_bin] = nz;
    for(Int_t i=0; i<nz; i++)
    {
      if(fscanf(HELI, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &Egamma, &z, &ReH1, &ImH1, &ReH2, &ImH2, &ReH3, &ImH3, &ReH4, &ImH4)!=10) break;;
      
     heli_H1[i][heli_bin] = TComplex(ReH1, ImH1);
     heli_H2[i][heli_bin] = TComplex(ReH2, ImH2);
     heli_H3[i][heli_bin] = TComplex(ReH3, ImH3);
     heli_H4[i][heli_bin] = TComplex(ReH4, ImH4);
    
     heli_z[i][heli_bin] = z;
     heli_en[heli_bin] = Egamma;
     
    }
    printf("%d, %d, %lf \n", heli_bin, heli_nz[heli_bin], Egamma);
    heli_bin++;
  }

  

  //Close file
  fclose(HELI);

}
//-----------------------------------------------------------------------------

Int_t GetEnergyBin_heli()
{
  //Get energy bin for sigma0 for given global energy
  Double_t Min = 1e38;
  Int_t eM = 0;

  for(Int_t e=0; e<heli_bin; e++)
    if(fabs(heli_en[e] - gEnergy) < Min)
    {
      Min = fabs(heli_en[e] - gEnergy);
      eM = e;
    }
  return eM;
}

//-----------------------------------------------------------------------------

