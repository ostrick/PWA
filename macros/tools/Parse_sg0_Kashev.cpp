static const Int_t EBINS  = 2048;
static const Int_t THBINS = 256;

//-----------------------------------------------------------------------------

Double_t Obs_val[EBINS][THBINS];
Double_t Obs_err[EBINS][THBINS];
Double_t Obs_unc[EBINS][THBINS];
Double_t Obs_th[EBINS][THBINS];
Double_t Obs_lo[EBINS];
Double_t Obs_en[EBINS];
Double_t Obs_hi[EBINS];
Double_t Obs_sy[EBINS];
Int_t Obs_pre[EBINS];
Int_t Obs_pts[EBINS];
Int_t Obs_bin;

//-----------------------------------------------------------------------------

void Parse_sg0(Char_t* VikFile, Double_t Sys, Char_t* ID, Int_t CthBins=30)
{
  Char_t Buffer[1024];
  Int_t ThetaBin;
  Int_t Dummy1, Dummy2;
  Double_t Energy, EnergyLo, EnergyHi, W, CosTheta, Theta, sigma0, Dsigma0, Esigma0;
  FILE* File_sg0_in;
  FILE* File_sg0_out;

  File_sg0_in = fopen(VikFile, "r");

  Obs_bin = 0;
  //Skip 9 uninteresting lines
  for(Int_t n=0; n<9; n++)
    fgets(Buffer, sizeof(Buffer), File_sg0_in);

  while(!feof(File_sg0_in))
  {
    for(Int_t t=0; t<CthBins; t++)
    {
      if(fscanf(File_sg0_in, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
         &Energy, &EnergyLo, &EnergyHi, &W, &CosTheta, &Theta, &sigma0, &Dsigma0, &Esigma0)!=9) break;
      Obs_th[Obs_bin][t] = Theta;
      //Obs_th[Obs_bin][t] = TMath::ACos(CosTheta)*TMath::RadToDeg();
      Obs_val[Obs_bin][t] = sigma0;
      Obs_err[Obs_bin][t] = Dsigma0;
      Obs_unc[Obs_bin][t] = TMath::Sqrt(Esigma0*Esigma0 - Dsigma0*Dsigma0);
      Obs_en[Obs_bin] = Energy;
      Obs_lo[Obs_bin] = EnergyLo;
      Obs_hi[Obs_bin] = EnergyHi;
      Obs_sy[Obs_bin] = Sys;
      Obs_pre[Obs_bin] = 0;
      Obs_pts[Obs_bin] = CthBins;
    }
    Obs_bin++;
  }
  fclose(File_sg0_in);

  sprintf(Buffer, "out_%s", VikFile);
  File_sg0_out = fopen(Buffer, "w");
  for(Int_t e=0; e<Obs_bin; e++)
  {
    fprintf(File_sg0_out, "E = %8.3f MeV, E_lo = %8.3f MeV, E_hi = %8.3f MeV\n", Obs_en[e], Obs_lo[e], Obs_hi[e]);
    fprintf(File_sg0_out, "Systematic = %5.3f, Preliminary = 0, %s\n", Obs_sy[e], ID);
    for(Int_t t=0; t<Obs_pts[e]; t++)
      fprintf(File_sg0_out, "%7.3f   %10.6f   %10.6f   %10.6f\n", Obs_th[e][t], Obs_val[e][t], Obs_err[e][t], Obs_unc[e][t]);
      //fprintf(File_sg0_out, "%7.3f   %10.6f   %10.6f\n", Obs_th[e][t], Obs_val[e][t], Obs_err[e][t]);
    fprintf(File_sg0_out, "-----------------------------------------------------------\n");
  }
  fclose(File_sg0_out);
}

//-----------------------------------------------------------------------------
