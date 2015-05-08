static const Int_t EBINS  = 2048;
static const Int_t THBINS = 256;

//-----------------------------------------------------------------------------

Double_t O1_val[EBINS][THBINS];
Double_t O1_err[EBINS][THBINS];
Double_t O1_unc[EBINS][THBINS];
Double_t O1_th[EBINS][THBINS];
Double_t O1_lo[EBINS];
Double_t O1_en[EBINS];
Double_t O1_hi[EBINS];
Double_t O1_sy[EBINS];
Int_t O1_pre[EBINS];
Int_t O1_pts[EBINS];
Int_t O1_bin;

Double_t O2_[EBINS][THBINS];
Double_t O2_err[EBINS][THBINS];
Double_t O2_unc[EBINS][THBINS];
Double_t O2_th[EBINS][THBINS];
Double_t O2_lo[EBINS];
Double_t O2_en[EBINS];
Double_t O2_hi[EBINS];
Double_t O2_sy[EBINS];
Int_t O2_pre[EBINS];
Int_t O2_pts[EBINS];
Int_t O2_bin;

//-----------------------------------------------------------------------------

void Parse_MAID(Char_t* VikFile, Double_t Sys, Char_t* ID, Int_t ThBins=18)
{
  Char_t Buffer[1024];
  Int_t ThetaBin;
  Int_t Dummy1, Dummy2;
  Double_t Energy, EnergyWd, Theta;
  Double_t Cx, DCx, Cz, DCz;
  FILE* File_in;
  FILE* File_out;

  File_in = fopen(VikFile, "r");

  O1_bin = 0;
  O2_bin = 0;
  //Skip 1 uninteresting line
  fgets(Buffer, sizeof(Buffer), File_in);

  while(!feof(File_in))
  {
    for(Int_t t=0; t<ThBins; t++)
    {
      if(fscanf(File_in, "%lf %lf %lf %lf %lf %lf %lf\n",
         &Energy, &EnergyWd, &Theta, &Cx, &DCx, &Cz, &DCz)!=7) break;
      O1_val[O1_bin][t] = Cx;
      O1_err[O1_bin][t] = DCx;
      O1_th[O1_bin][t] = Theta;
      O1_en[O1_bin] = Energy;
      O1_lo[O1_bin] = Energy-EnergyWd;
      O1_hi[O1_bin] = Energy+EnergyWd;
      O1_pre[O1_bin] = 0;
      O1_pts[O1_bin] = ThBins;
      O1_sy[O1_bin] = Sys;

      O2_[O2_bin][t] = Cz;
      O2_err[O2_bin][t] = DCz;
      O2_th[O2_bin][t] = Theta;
      O2_en[O2_bin] = Energy;
      O2_lo[O2_bin] = Energy-EnergyWd;
      O2_hi[O2_bin] = Energy+EnergyWd;
      O2_pre[O2_bin] = 0;
      O2_pts[O2_bin] = ThBins;
      O2_sy[O2_bin] = Sys;
    }
    O1_bin++;
    O2_bin++;
  }
  fclose(File_in);

  File_out = fopen("O1", "w");
  for(Int_t e=0; e<O1_bin; e++)
  {
    fprintf(File_out, "E = %8.3f MeV, E_lo = %8.3f MeV, E_hi = %8.3f MeV\n", O1_en[e], O1_lo[e], O1_hi[e]);
    fprintf(File_out, "Systematic = %5.3f, Preliminary = 0, %s\n", O1_sy[e], ID);
    for(Int_t t=0; t<O1_pts[e]; t++)
      //fprintf(File_out, "%7.3f   %10.6f   %10.6f   %10.6f\n", O1_th[e][t], O1_val[e][t], O1_err[e][t], O1_unc[e][t]);
      fprintf(File_out, "%7.3f   %10.6f   %10.6f\n", O1_th[e][t], O1_val[e][t], O1_err[e][t]);
    fprintf(File_out, "-----------------------------------------------------------\n");
  }
  fclose(File_out);


  File_out = fopen("O2", "w");
  for(Int_t e=0; e<O2_bin; e++)
  {
    fprintf(File_out, "E = %8.3f MeV, E_lo = %8.3f MeV, E_hi = %8.3f MeV\n", O2_en[e], O2_lo[e], O2_hi[e]);
    fprintf(File_out, "Systematic = %5.3f, Preliminary = 0, %s\n", O2_sy[e], ID);
    for(Int_t t=0; t<O2_pts[e]; t++)
      //fprintf(File_out, "%7.3f   %10.6f   %10.6f   %10.6f\n", O2_th[e][t], O2_[e][t], O2_err[e][t], O2_unc[e][t]);
      fprintf(File_out, "%7.3f   %10.6f   %10.6f\n", O2_th[e][t], O2_[e][t], O2_err[e][t]);
    fprintf(File_out, "-----------------------------------------------------------\n");
  }
  fclose(File_out);
}

//-----------------------------------------------------------------------------
