#include "Fitter.h"

Int_t ierflg = 0;
Double_t arglist[10];
static Double_t vstart[8*(LBINS-1) + OBS];
static Double_t step[8*(LBINS-1) + OBS];

//-----------------------------------------------------------------------------

void fcn_main(Double_t *par)
{
  Int_t Index;

  //Force some parameters to be positive
  //A parameter range definition in TMinuit::mnparm() leads to some strange fitting problems
  if(FIX_EP_PHASE[0]) par[0] = fabs(par[0]); //E0+ magnitude
  if(FIX_EP_PHASE[1]) par[2] = fabs(par[2]); //E1+ magnitude
  if(FIX_MP_PHASE[1]) par[4] = fabs(par[4]); //M1+ magnitude
  if(FIX_MM_PHASE[1]) par[6] = fabs(par[6]); //M1- magnitude

  //Set fittable mulipoles to given parameters
  //s,p waves
  Ep[0] = TComplex(par[0],  par[1], FIX_EP_PHASE[0]);
  Ep[1] = TComplex(par[2],  par[3], FIX_EP_PHASE[1]);
  Mp[1] = TComplex(par[4],  par[5], FIX_MP_PHASE[1]);
  Mm[1] = TComplex(par[6],  par[7], FIX_MM_PHASE[1]);

  //Higher waves
  for(Int_t l=2; l<L_MAX+1; l++)
  {
    Index = 8*(l-1);
    //Force some parameters to be positive
    //A parameter range definition in TMinuit::mnparm() leads to some strange fitting problems
    if(FIX_EP_PHASE[l]) par[Index+0] = fabs(par[Index+0]); //El+ magnitude
    if(FIX_MP_PHASE[l]) par[Index+2] = fabs(par[Index+2]); //Ml+ magnitude
    if(FIX_EM_PHASE[l]) par[Index+4] = fabs(par[Index+4]); //El- magnitude
    if(FIX_MM_PHASE[l]) par[Index+6] = fabs(par[Index+6]); //Ml- magnitude

    //Set fittable mulipoles to given parameters
    //d,f,... waves
    Ep[l] = TComplex(par[Index+0], par[Index+1], FIX_EP_PHASE[l]);
    Mp[l] = TComplex(par[Index+2], par[Index+3], FIX_MP_PHASE[l]);
    Em[l] = TComplex(par[Index+4], par[Index+5], FIX_EM_PHASE[l]);
    Mm[l] = TComplex(par[Index+6], par[Index+7], FIX_MM_PHASE[l]);
  }

  //Set observable scale factors to given parameters
  for(Int_t o=0; o<OBS+1; o++)
    f_obs[o] = par[8*L_MAX + o];
}

//-----------------------------------------------------------------------------

void fcn_chi_pen(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  //Call main part of fcn routine
  fcn_main(par);

  //Return sum of chi^2 and penalties from multipoles and observable scaling to fitter
  f = ChiSq() + Penalty() + Scale();
}

//-----------------------------------------------------------------------------

void fcn_chi(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  //Call main part of fcn routine
  fcn_main(par);

  //Return only chi^2 to fitter
  f = ChiSq();
}

//-----------------------------------------------------------------------------

Double_t ChiSq()
{
  Double_t ChiSq_0,  ChiSq_S,  ChiSq_T,  ChiSq_P;
  Double_t ChiSq_E,  ChiSq_F,  ChiSq_G,  ChiSq_H;
  Double_t ChiSq_Cx, ChiSq_Cz, ChiSq_Ox, ChiSq_Oz;
  Double_t ChiSq_Lx, ChiSq_Lz, ChiSq_Tx, ChiSq_Tz;

  //Get chi^2 for all observables for given global energy
  ChiSq_0  = GetChiSq_sg0();
  ChiSq_S  = GetChiSq_sgS()  + GetChiSq_S();
  ChiSq_T  = GetChiSq_sgT()  + GetChiSq_T();
  ChiSq_P  = GetChiSq_sgP()  + GetChiSq_P();
  ChiSq_E  = GetChiSq_sgE()  + GetChiSq_E();
  ChiSq_F  = GetChiSq_sgF()  + GetChiSq_F();
  ChiSq_G  = GetChiSq_sgG()  + GetChiSq_G();
  ChiSq_H  = GetChiSq_sgH()  + GetChiSq_H();
  ChiSq_Cx = GetChiSq_sgCx() + GetChiSq_Cx();
  ChiSq_Cz = GetChiSq_sgCz() + GetChiSq_Cz();
  ChiSq_Ox = GetChiSq_sgOx() + GetChiSq_Ox();
  ChiSq_Oz = GetChiSq_sgOz() + GetChiSq_Oz();
  ChiSq_Lx = GetChiSq_sgLx() + GetChiSq_Lx();
  ChiSq_Lz = GetChiSq_sgLz() + GetChiSq_Lz();
  ChiSq_Tx = GetChiSq_sgTx() + GetChiSq_Tx();
  ChiSq_Tz = GetChiSq_sgTz() + GetChiSq_Tz();

  if(ONLY_CROSS_S) return (ChiSq_0 + ChiSq_S);
  if(ONLY_CROSS_F) return (ChiSq_0 + ChiSq_F);
  return (ChiSq_0  + ChiSq_S  + ChiSq_T  + ChiSq_P  +
          ChiSq_E  + ChiSq_F  + ChiSq_G  + ChiSq_H  +
          ChiSq_Cx + ChiSq_Cz + ChiSq_Ox + ChiSq_Oz +
          ChiSq_Lx + ChiSq_Lz + ChiSq_Tx + ChiSq_Tz);
}

//-----------------------------------------------------------------------------

Double_t Scale()
{
  if(FIX_SCALES) return 0.0;

  Double_t Scale_0,  Scale_S,  Scale_T,  Scale_P;
  Double_t Scale_E,  Scale_F,  Scale_G,  Scale_H;
  Double_t Scale_Cx, Scale_Cz, Scale_Ox, Scale_Oz;
  Double_t Scale_Lx, Scale_Lz, Scale_Tx, Scale_Tz;

  Scale_0  = GetScale_sg0();
  Scale_S  = GetScale_sgS()  + GetScale_S();
  Scale_T  = GetScale_sgT()  + GetScale_T();
  Scale_P  = GetScale_sgP()  + GetScale_P();
  Scale_E  = GetScale_sgE()  + GetScale_E();
  Scale_F  = GetScale_sgF()  + GetScale_F();
  Scale_G  = GetScale_sgG()  + GetScale_G();
  Scale_H  = GetScale_sgH()  + GetScale_H();
  Scale_Cx = GetScale_sgCx() + GetScale_Cx();
  Scale_Cz = GetScale_sgCz() + GetScale_Cz();
  Scale_Ox = GetScale_sgOx() + GetScale_Ox();
  Scale_Oz = GetScale_sgOz() + GetScale_Oz();
  Scale_Lx = GetScale_sgLx() + GetScale_Lx();
  Scale_Lz = GetScale_sgLz() + GetScale_Lz();
  Scale_Tx = GetScale_sgTx() + GetScale_Tx();
  Scale_Tz = GetScale_sgTz() + GetScale_Tz();

  if(ONLY_CROSS_S) return SCALING*(Scale_0 + Scale_S)/NSca();
  if(ONLY_CROSS_F) return SCALING*(Scale_0 + Scale_F)/NSca();

  //Without additional weighting (SCALING), the penalty should be comparable to chi^2.
  //So, it is divided by the number of scaling parameters (NSca) that are taken into
  //account (i.e. normalise the penalty to a single parameter case). The number of data
  //points is already considered when calculating the individual contributions.
  return SCALING*(Scale_0  + Scale_S  + Scale_T  + Scale_P  +
                  Scale_E  + Scale_F  + Scale_G  + Scale_H  +
                  Scale_Cx + Scale_Cz + Scale_Ox + Scale_Oz +
                  Scale_Lx + Scale_Lz + Scale_Tx + Scale_Tz)/NSca();
}

//-----------------------------------------------------------------------------

Double_t Penalty()
{
  switch(PENALTY_MODE)
  {
   case MLP1:
    return PenaltyMLP1();
   case MLP2:
    return PenaltyMLP2();
   case MLP3:
    return PenaltyMLP1() + PenaltyMLP2();
   case CGLN:
    return PenaltyCGLN();
   case HELI:
    return PenaltyHELI();
   case HELI2:
    return PenaltyHELI2();
   case INVA:
    return PenaltyINVA();
   case INVB:
    return PenaltyINVB();
   case CONT:
    return PenaltyCONT();
   case NONE:
   default:
    return 0.0;
  }
}

//-----------------------------------------------------------------------------

Double_t PenaltyMLP1()
{
  Int_t eM = GetEnergyBin_maid();
  Double_t SumSq = 0.0;
  Double_t MagSq = 0.0;

  //Calculate deviation from model for s,p wave real parts
  SumSq+=(Ep[0].Re() - maid_Ep[0][eM].Re())*(Ep[0].Re() - maid_Ep[0][eM].Re());
  SumSq+=(Ep[1].Re() - maid_Ep[1][eM].Re())*(Ep[1].Re() - maid_Ep[1][eM].Re());
  SumSq+=(Mp[1].Re() - maid_Mp[1][eM].Re())*(Mp[1].Re() - maid_Mp[1][eM].Re());
  SumSq+=(Mm[1].Re() - maid_Mm[1][eM].Re())*(Mm[1].Re() - maid_Mm[1][eM].Re());
  if(ONLY_CROSS_S || ONLY_CROSS_F) //When p waves are kept real...
  {
    //Calculate deviation from unitary value for s wave imaginary parts
    SumSq+=(Ep[0].Im() - ImE0p())*(Ep[0].Im() - ImE0p());
    //Calculate deviation from 0 for p wave imaginary parts
    SumSq+=Ep[1].Im()*Ep[1].Im();
    SumSq+=Mp[1].Im()*Mp[1].Im();
    SumSq+=Mm[1].Im()*Mm[1].Im();
  }
  else
  {
    //Calculate deviation from model for s,p wave imaginary parts
    SumSq+=(Ep[1].Im() - maid_Ep[1][eM].Im())*(Ep[1].Im() - maid_Ep[1][eM].Im());
    SumSq+=(Mp[1].Im() - maid_Mp[1][eM].Im())*(Mp[1].Im() - maid_Mp[1][eM].Im());
    SumSq+=(Mm[1].Im() - maid_Mm[1][eM].Im())*(Mm[1].Im() - maid_Mm[1][eM].Im());
    SumSq+=(Ep[0].Im() - maid_Ep[0][eM].Im())*(Ep[0].Im() - maid_Ep[0][eM].Im());
  }

  //Calculate deviations for d,f,... waves
  for(Int_t l=2; l<L_MAX+1; l++)
  {
    //Calculate deviation from model for d,f,... waves
    SumSq+=(Ep[l].Re() - maid_Ep[l][eM].Re())*(Ep[l].Re() - maid_Ep[l][eM].Re());
    SumSq+=(Mp[l].Re() - maid_Mp[l][eM].Re())*(Mp[l].Re() - maid_Mp[l][eM].Re());
    SumSq+=(Em[l].Re() - maid_Em[l][eM].Re())*(Em[l].Re() - maid_Em[l][eM].Re());
    SumSq+=(Mm[l].Re() - maid_Mm[l][eM].Re())*(Mm[l].Re() - maid_Mm[l][eM].Re());
    SumSq+=(Ep[l].Im() - maid_Ep[l][eM].Im())*(Ep[l].Im() - maid_Ep[l][eM].Im());
    SumSq+=(Mp[l].Im() - maid_Mp[l][eM].Im())*(Mp[l].Im() - maid_Mp[l][eM].Im());
    SumSq+=(Em[l].Im() - maid_Em[l][eM].Im())*(Em[l].Im() - maid_Em[l][eM].Im());
    SumSq+=(Mm[l].Im() - maid_Mm[l][eM].Im())*(Mm[l].Im() - maid_Mm[l][eM].Im());
  }

  if(ONLY_CROSS_S || ONLY_CROSS_F) //When p waves are kept real...
  {
    //Sum up magnitude^2 of s wave model multipole with unitary s wave
    MagSq+=(maid_Ep[0][eM].Re()*maid_Ep[0][eM].Re() + ImE0p()*ImE0p());
    //Sum up magnitudes^2 of p wave model multipole with real p waves
    MagSq+=maid_Ep[1][eM].Re()*maid_Ep[1][eM].Re();
    MagSq+=maid_Mp[1][eM].Re()*maid_Mp[1][eM].Re();
    MagSq+=maid_Mm[1][eM].Re()*maid_Mm[1][eM].Re();
  }
  else
  {
    //Sum up magnitudes^2 of s,p wave model multipoles
    MagSq+=maid_Ep[0][eM].Rho2();
    MagSq+=maid_Ep[1][eM].Rho2();
    MagSq+=maid_Mp[1][eM].Rho2();
    MagSq+=maid_Mm[1][eM].Rho2();
  }
  for(Int_t l=2; l<L_MAX+1; l++)
  {
    //Sum up magnitudes^2 of d,f,... wave model multipoles
    MagSq+=maid_Ep[l][eM].Rho2();
    MagSq+=maid_Mp[l][eM].Rho2();
    MagSq+=maid_Em[l][eM].Rho2();
    MagSq+=maid_Mm[l][eM].Rho2();
  }
  //Calculate average magnitude^2 of a single real or imaginary part of model multipoles
  MagSq/=(8.0*L_MAX);

  //Without additional weighting (PENALTY), the penalty should be comparable to chi^2.
  //So, it is divided by the number of multipole parameters (NMlp) that are taken into
  //account (i.e. normalise the penalty to a single parameter case) and multiplied by
  //the number of data points (NPts) that are used for chi^2. Also, it is normalised
  //to an average magnitude^2 of a single real (or imaginary) part of model multipoles.
  return PENALTY[MLP1]*SumSq*NPts()/(MagSq*NMlp());
}

//-----------------------------------------------------------------------------

Double_t PenaltyMLP2()
{
  Int_t eM = GetEnergyBin_maid();
  Double_t SumSq = 0.0;

  //Calculate deviation from model for s,p wave real parts
  SumSq+=(Ep[0].Re() - maid_Ep[0][eM].Re())*(Ep[0].Re() - maid_Ep[0][eM].Re())/((maid_Ep[0][eM].Re()*maid_Ep[0][eM].Re() + EPSILON2));
  SumSq+=(Ep[1].Re() - maid_Ep[1][eM].Re())*(Ep[1].Re() - maid_Ep[1][eM].Re())/((maid_Ep[1][eM].Re()*maid_Ep[1][eM].Re() + EPSILON2));
  SumSq+=(Mp[1].Re() - maid_Mp[1][eM].Re())*(Mp[1].Re() - maid_Mp[1][eM].Re())/((maid_Mp[1][eM].Re()*maid_Mp[1][eM].Re() + EPSILON2));
  SumSq+=(Mm[1].Re() - maid_Mm[1][eM].Re())*(Mm[1].Re() - maid_Mm[1][eM].Re())/((maid_Mm[1][eM].Re()*maid_Mm[1][eM].Re() + EPSILON2));
  if(ONLY_CROSS_S || ONLY_CROSS_F) //When p waves are kept real...
  {
    //Calculate deviation from unitary value for s wave imaginary parts
    SumSq+=(Ep[0].Im() - ImE0p())*(Ep[0].Im() - ImE0p())/((ImE0p()*ImE0p()+EPSILON2));
    //Calculate deviation from 0 for p wave imaginary parts
    SumSq+=Ep[1].Im()*Ep[1].Im()/EPSILON2;
    SumSq+=Mp[1].Im()*Mp[1].Im()/EPSILON2;
    SumSq+=Mm[1].Im()*Mm[1].Im()/EPSILON2;
  }
  else
  {
    //Calculate deviation from model for s,p wave imaginary parts
    SumSq+=(Ep[0].Im() - maid_Ep[0][eM].Im())*(Ep[0].Im() - maid_Ep[0][eM].Im())/((maid_Ep[0][eM].Im()*maid_Ep[0][eM].Im() + EPSILON2));
    SumSq+=(Ep[1].Im() - maid_Ep[1][eM].Im())*(Ep[1].Im() - maid_Ep[1][eM].Im())/((maid_Ep[1][eM].Im()*maid_Ep[1][eM].Im() + EPSILON2));
    SumSq+=(Mp[1].Im() - maid_Mp[1][eM].Im())*(Mp[1].Im() - maid_Mp[1][eM].Im())/((maid_Mp[1][eM].Im()*maid_Mp[1][eM].Im() + EPSILON2));
    SumSq+=(Mm[1].Im() - maid_Mm[1][eM].Im())*(Mm[1].Im() - maid_Mm[1][eM].Im())/((maid_Mm[1][eM].Im()*maid_Mm[1][eM].Im() + EPSILON2));
  }

  //Calculate deviations for d,f,... waves
  for(Int_t l=2; l<L_MAX+1; l++)
  {
    //Calculate deviation from model for d,f,... waves
    SumSq+=(Ep[l].Re() - maid_Ep[l][eM].Re())*(Ep[l].Re() - maid_Ep[l][eM].Re())/((maid_Ep[l][eM].Re()*maid_Ep[l][eM].Re() + EPSILON2));
    SumSq+=(Mp[l].Re() - maid_Mp[l][eM].Re())*(Mp[l].Re() - maid_Mp[l][eM].Re())/((maid_Mp[l][eM].Re()*maid_Mp[l][eM].Re() + EPSILON2));
    SumSq+=(Em[l].Re() - maid_Em[l][eM].Re())*(Em[l].Re() - maid_Em[l][eM].Re())/((maid_Em[l][eM].Re()*maid_Em[l][eM].Re() + EPSILON2));
    SumSq+=(Mm[l].Re() - maid_Mm[l][eM].Re())*(Mm[l].Re() - maid_Mm[l][eM].Re())/((maid_Mm[l][eM].Re()*maid_Mm[l][eM].Re() + EPSILON2));
    SumSq+=(Ep[l].Im() - maid_Ep[l][eM].Im())*(Ep[l].Im() - maid_Ep[l][eM].Im())/((maid_Ep[l][eM].Im()*maid_Ep[l][eM].Im() + EPSILON2));
    SumSq+=(Mp[l].Im() - maid_Mp[l][eM].Im())*(Mp[l].Im() - maid_Mp[l][eM].Im())/((maid_Mp[l][eM].Im()*maid_Mp[l][eM].Im() + EPSILON2));
    SumSq+=(Em[l].Im() - maid_Em[l][eM].Im())*(Em[l].Im() - maid_Em[l][eM].Im())/((maid_Em[l][eM].Im()*maid_Em[l][eM].Im() + EPSILON2));
    SumSq+=(Mm[l].Im() - maid_Mm[l][eM].Im())*(Mm[l].Im() - maid_Mm[l][eM].Im())/((maid_Mm[l][eM].Im()*maid_Mm[l][eM].Im() + EPSILON2));
  }

  //Without additional weighting (PENALTY), the penalty should be comparable to chi^2.
  //So, it is divided by the number of multipole parameters (NMlp) that are taken into
  //account (i.e. normalise the penalty to a single parameter case) and multiplied by
  //the number of data points (NPts) that are used for chi^2.
  return PENALTY[MLP2]*SumSq*NPts()/NMlp();
}

//-----------------------------------------------------------------------------

Double_t PenaltyCGLN()
{
  Int_t eM = GetEnergyBin_maid();
  Int_t Npts = NPts();
  Int_t Nsteps;
  Double_t CosTheta, DCosTheta;
  Double_t SumSq = 0.0;

  Nsteps = Npts;
  if(Nsteps > NSTEPS) Nsteps = NSTEPS; //Limit number of step for reasonable performance

  DCosTheta = 2.0/Nsteps; //Calculate costheta step size
  for(Int_t n=0; n<Nsteps; n++)
  {
    CosTheta = DCosTheta*n + DCosTheta/2.0 - 1.0; //Calculate central position of current costheta step
    SumSq+=EvaluateCGLN(CosTheta, eM);
  }

  //Without additional weighting, the penalty should be comparable to chi^2.
  //Here it is summed up to the number of theta points from all observables
  //(i.e. NPts costheta positions are used) and EvaluateCGLN() normalises to
  //the sum of F1...F4 magnitudes. Hence, no further normalisation is required.
  //Note: The number of steps is limited, hence for large number of points,
  //a correction is required.
  return SumSq*Npts/Nsteps;
}

//-----------------------------------------------------------------------------

Double_t PenaltyHELI()
{
  Int_t eM = GetEnergyBin_maid();
  Int_t Npts = NPts();
  Int_t Nsteps;
  Double_t CosTheta, DCosTheta;
  Double_t SumSq = 0.0;

  Nsteps = Npts;
  if(Nsteps > NSTEPS) Nsteps = NSTEPS; //Limit number of step for reasonable performance

  DCosTheta = 2.0/Nsteps; //Calculate costheta step size
  for(Int_t n=0; n<Nsteps; n++)
  {
    CosTheta = DCosTheta*n + DCosTheta/2.0 - 1.0; //Calculate central position of current costheta step
    SumSq+=EvaluateHELI(CosTheta, eM);
  }

  //Without additional weighting, the penalty should be comparable to chi^2.
  //Here it is summed up to the number of theta points from all observables
  //(i.e. NPts costheta positions are used) and EvaluateHELI() normalises to
  //the sum of H1...H4 magnitudes. Hence, no further normalisation is required.
  //Note: The number of steps is limited, hence for large number of points,
  //a correction is required.
  return SumSq*Npts/Nsteps;
}

//-----------------------------------------------------------------------------


Double_t PenaltyHELI2()
{
  Int_t eM = GetEnergyBin_heli();
  Double_t CosTheta, DCosTheta;
  Double_t SumSq = 0.0;

  for(Int_t n=0; n < heli_nz[eM]; n++)
  {
    SumSq+=(heli_H1[n][eM].Re()-H1(heli_z[n][eM]).Re())*(heli_H1[n][eM].Re()-H1(heli_z[n][eM]).Re());
    SumSq+=(heli_H2[n][eM].Re()-H2(heli_z[n][eM]).Re())*(heli_H2[n][eM].Re()-H2(heli_z[n][eM]).Re());
    SumSq+=(heli_H3[n][eM].Re()-H3(heli_z[n][eM]).Re())*(heli_H3[n][eM].Re()-H3(heli_z[n][eM]).Re());
    SumSq+=(heli_H4[n][eM].Re()-H4(heli_z[n][eM]).Re())*(heli_H4[n][eM].Re()-H4(heli_z[n][eM]).Re());
    
    SumSq+=(heli_H1[n][eM].Im()-H1(heli_z[n][eM]).Im())*(heli_H1[n][eM].Im()-H1(heli_z[n][eM]).Im());
    SumSq+=(heli_H2[n][eM].Im()-H2(heli_z[n][eM]).Im())*(heli_H2[n][eM].Im()-H2(heli_z[n][eM]).Im());
    SumSq+=(heli_H3[n][eM].Im()-H3(heli_z[n][eM]).Im())*(heli_H3[n][eM].Im()-H3(heli_z[n][eM]).Im());
    SumSq+=(heli_H4[n][eM].Im()-H4(heli_z[n][eM]).Im())*(heli_H4[n][eM].Im()-H4(heli_z[n][eM]).Im());
  
  }

  return SumSq;
}




//-----------------------------------------------------------------------------

Double_t PenaltyINVA()
{
  Int_t eM = GetEnergyBin_maid();
  Int_t Npts = NPts();
  Int_t Nsteps;
  Double_t CosTheta, DCosTheta;
  Double_t SumSq = 0.0;

  Nsteps = Npts;
  if(Nsteps > NSTEPS) Nsteps = NSTEPS; //Limit number of step for reasonable performance

  DCosTheta = 2.0/Nsteps; //Calculate costheta step size
  for(Int_t n=0; n<Nsteps; n++)
  {
    CosTheta = DCosTheta*n + DCosTheta/2.0 - 1.0; //Calculate central position of current costheta step
    SumSq+=EvaluateINVA(CosTheta, eM);
  }

  //Without additional weighting, the penalty should be comparable to chi^2.
  //Here it is summed up to the number of theta points from all observables
  //(i.e. NPts costheta positions are used) and EvaluateINVA() normalises to
  //the sum of A1...A4 magnitudes. Hence, no further normalisation is required.
  //Note: The number of steps is limited, hence for large number of points,
  //a correction is required.
  return SumSq*Npts/Nsteps;
}

//-----------------------------------------------------------------------------

Double_t PenaltyINVB()
{
  Int_t eM = GetEnergyBin_maid();
  Int_t Npts = NPts();
  Int_t Nsteps;
  Double_t CosTheta, DCosTheta;
  Double_t SumSq = 0.0;

  Nsteps = Npts;
  if(Nsteps > NSTEPS) Nsteps = NSTEPS; //Limit number of step for reasonable performance

  DCosTheta = 2.0/Nsteps; //Calculate costheta step size
  for(Int_t n=0; n<Nsteps; n++)
  {
    CosTheta = DCosTheta*n + DCosTheta/2.0 - 1.0; //Calculate central position of current costheta step
    SumSq+=EvaluateINVB(CosTheta, eM);
  }

  //Without additional weighting, the penalty should be comparable to chi^2.
  //Here it is summed up to the number of theta points from all observables
  //(i.e. NPts costheta positions are used) and EvaluateINVB() normalises to
  //the sum of B1...B8 magnitudes. Hence, no further normalisation is required.
  //Note: The number of steps is limited, hence for large number of points,
  //a correction is required.
  return SumSq*Npts/Nsteps;
}



//-----------------------------------------------------------------------------

Double_t PenaltyCONT()
{
  Int_t eM = GetEnergyBin_maid();
  Int_t Npts = NPts();
  Int_t Nsteps;
  Double_t CosTheta, DCosTheta;
  Double_t SumSq = 0.0;

  Nsteps = Npts;
  if(Nsteps > NSTEPS) Nsteps = NSTEPS; //Limit number of step for reasonable performance

  DCosTheta = 2.0/Nsteps; //Calculate costheta step size
  for(Int_t n=0; n<Nsteps; n++)
  {
    CosTheta = DCosTheta*n + DCosTheta/2.0 - 1.0; //Calculate central position of current costheta step
    SumSq+=EvaluateCONT(CosTheta, eM);
  }

  //Without additional weighting, the penalty should be comparable to chi^2.
  //Here it is summed up to the number of theta points from all observables
  //(i.e. NPts costheta positions are used) and EvaluateCGLN() normalises to
  //the sum of F1...F4 magnitudes. Hence, no further normalisation is required.
  //Note: The number of steps is limited, hence for large number of points,
  //a correction is required.
  return SumSq*Npts/Nsteps;
}

//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------

Int_t NPts()
{
  Int_t NPts_0,  NPts_S,  NPts_T,  NPts_P;
  Int_t NPts_E,  NPts_F,  NPts_G,  NPts_H;
  Int_t NPts_Cx, NPts_Cz, NPts_Ox, NPts_Oz;
  Int_t NPts_Lx, NPts_Lz, NPts_Tx, NPts_Tz;

  NPts_0  = GetNPts_sg0();
  NPts_S  = GetNPts_sgS()  + GetNPts_S();
  NPts_T  = GetNPts_sgT()  + GetNPts_T();
  NPts_P  = GetNPts_sgP()  + GetNPts_P();
  NPts_E  = GetNPts_sgE()  + GetNPts_E();
  NPts_F  = GetNPts_sgF()  + GetNPts_F();
  NPts_G  = GetNPts_sgG()  + GetNPts_G();
  NPts_H  = GetNPts_sgH()  + GetNPts_H();
  NPts_Cx = GetNPts_sgCx() + GetNPts_Cx();
  NPts_Cz = GetNPts_sgCz() + GetNPts_Cz();
  NPts_Ox = GetNPts_sgOx() + GetNPts_Ox();
  NPts_Oz = GetNPts_sgOz() + GetNPts_Oz();
  NPts_Lx = GetNPts_sgLx() + GetNPts_Lx();
  NPts_Lz = GetNPts_sgLz() + GetNPts_Lz();
  NPts_Tx = GetNPts_sgTx() + GetNPts_Tx();
  NPts_Tz = GetNPts_sgTz() + GetNPts_Tz();

  if(ONLY_CROSS_S) return (NPts_0 + NPts_S);
  if(ONLY_CROSS_F) return (NPts_0 + NPts_F);

  return (NPts_0  + NPts_S  + NPts_T  + NPts_P  +
          NPts_E  + NPts_F  + NPts_G  + NPts_H  +
          NPts_Cx + NPts_Cz + NPts_Ox + NPts_Oz +
          NPts_Lx + NPts_Lz + NPts_Tx + NPts_Tz);
}

//-----------------------------------------------------------------------------

Int_t NMlp()
{
  Int_t NMlp = L_MAX*8; //We start with 4 multipoles (8 real numbers) per angular momentum (L=0,1 combined), ...

  if(FIX_EP[0]) NMlp-=2; //... the E0+ multipole can be fixed, ...
  if(FIX_EP[1]) NMlp-=2; //... the E1+ multipole can be fixed, ...
  if(FIX_MP[1]) NMlp-=2; //... the M1+ multipole can be fixed, ...
  if(FIX_MM[1]) NMlp-=2; //... the M1- multipole can be fixed, ...
  //Fixed phases may decrease parameter count only, if complete wave not already fixed
  if(FIX_EP_PHASE[0] && !FIX_EP[0]) NMlp-=1; //... the phase of the E0+ multipole can be fixed, ...
  if(FIX_EP_PHASE[1] && !FIX_EP[1]) NMlp-=1; //... the phase of the E1+ multipole can be fixed, ...
  if(FIX_MP_PHASE[1] && !FIX_MP[1]) NMlp-=1; //... the phase of the M1+ multipole can be fixed, ...
  if(FIX_MM_PHASE[1] && !FIX_MM[1]) NMlp-=1; //... the phase of the M1- multipole can be fixed, ...
  for(Int_t l=2; l<L_MAX+1; l++)
  {
    if(FIX_EP[l]) NMlp-=2; //... the El+ multipole can be fixed, ...
    if(FIX_MP[l]) NMlp-=2; //... the Ml+ multipole can be fixed, ...
    if(FIX_EM[l]) NMlp-=2; //... the El- multipole can be fixed, ...
    if(FIX_MM[l]) NMlp-=2; //... the Ml- multipole can be fixed, ...
    //Fixed phases may decrease parameter count only, if complete wave not already fixed
    if(FIX_EP_PHASE[l] && !FIX_EP[l]) NMlp-=1; //... the El+ phase can be fixed, ...
    if(FIX_MP_PHASE[l] && !FIX_MP[l]) NMlp-=1; //... the Ml+ phase can be fixed, ...
    if(FIX_EM_PHASE[l] && !FIX_EM[l]) NMlp-=1; //... the El- phase can be fixed, ...
    if(FIX_MM_PHASE[l] && !FIX_MM[l]) NMlp-=1; //... the Ml- phase can be fixed, ...
  }
  if(ONLY_CROSS_S || ONLY_CROSS_F) NMlp-=3; //... the imaginary parts of three p waves can be fixed, ...
  if(FIX_RE_E0P) NMlp-=1; //... the real part of one s wave can be fixed.
  if(FIX_IM_E0P) NMlp-=1; //... the imaginary part of one s wave can be fixed.

  return NMlp;
}

//-----------------------------------------------------------------------------

Int_t NSca()
{
  if(FIX_SCALES) return 0;

  Int_t NSca = 0;

  //Handle special cases first
  if(ONLY_CROSS_S)
  {
    if(GetNPts_sg0()) NSca+=1;
    if(GetNPts_sgS()) NSca+=1;
    if(GetNPts_S())   NSca+=1;
    return NSca;
  }
  if(ONLY_CROSS_F)
  {
    if(GetNPts_sg0()) NSca+=1;
    if(GetNPts_sgF()) NSca+=1;
    if(GetNPts_F())   NSca+=1;
    return NSca;
  }

  //If observable scalings are not fixed, we get one parameter for each scalable observable (i.e. observable with data)
  if(GetNPts_sg0())  NSca+=1;
  if(GetNPts_sgS())  NSca+=1;
  if(GetNPts_sgT())  NSca+=1;
  if(GetNPts_sgP())  NSca+=1;
  if(GetNPts_sgE())  NSca+=1;
  if(GetNPts_sgF())  NSca+=1;
  if(GetNPts_sgG())  NSca+=1;
  if(GetNPts_sgH())  NSca+=1;
  if(GetNPts_sgCx()) NSca+=1;
  if(GetNPts_sgCz()) NSca+=1;
  if(GetNPts_sgOx()) NSca+=1;
  if(GetNPts_sgOz()) NSca+=1;
  if(GetNPts_sgLx()) NSca+=1;
  if(GetNPts_sgLz()) NSca+=1;
  if(GetNPts_sgTx()) NSca+=1;
  if(GetNPts_sgTz()) NSca+=1;
  if(GetNPts_S())    NSca+=1;
  if(GetNPts_T())    NSca+=1;
  if(GetNPts_P())    NSca+=1;
  if(GetNPts_E())    NSca+=1;
  if(GetNPts_F())    NSca+=1;
  if(GetNPts_G())    NSca+=1;
  if(GetNPts_H())    NSca+=1;
  if(GetNPts_Cx())   NSca+=1;
  if(GetNPts_Cz())   NSca+=1;
  if(GetNPts_Ox())   NSca+=1;
  if(GetNPts_Oz())   NSca+=1;
  if(GetNPts_Lx())   NSca+=1;
  if(GetNPts_Lz())   NSca+=1;
  if(GetNPts_Tx())   NSca+=1;
  if(GetNPts_Tz())   NSca+=1;

  return NSca;
}

//-----------------------------------------------------------------------------

Int_t NPar()
{
  return NMlp()+NSca();
}

//-----------------------------------------------------------------------------

Int_t NDF()
{
  return NPts()-NPar();
}

//-----------------------------------------------------------------------------

Double_t VariateRel()
{
  return (1.0 - VARIATION[REL] + gRandom->Rndm()*2.0*VARIATION[REL]);
}

//-----------------------------------------------------------------------------

Double_t VariateAbs()
{
  return (gRandom->Rndm()*2.0*VARIATION[ABS] - VARIATION[ABS]);
}

//-----------------------------------------------------------------------------

Bool_t SingleFit()
{
  Double_t Param[L_MAX*8 + OBS];
  Double_t Error[L_MAX*8 + OBS];

  //Initialise Minuit fitter and use fcn with penalty
  InitMinuit();
  gMinuit->SetFCN(fcn_chi_pen);

  //Set starting values for parameters (start with MAID solution)
  SetSPWaves(); //s,p waves
  for(Int_t l=2; l<L_MAX+1; l++) //Higher waves
    SetWaves(l);

  //When fitting threshold data, keep p waves real...
  if(ONLY_CROSS_S || ONLY_CROSS_F) SetRealPWaves();
  //...and parametrise E0+ imaginary part using unitarity
  if(ONLY_CROSS_S || ONLY_CROSS_F) SetUnitarity();

  //Set fitting start parameters in Minuit
  SetParameters();

  //Fix s,p waves or phases during the fitting process
  FixSPPhases();
  FixSPWaves();
  //Fix d,f,... waves or phases during the fitting process when requested (fit contains not enough data to constrain all waves simultaneously)
  for(Int_t l=2; l<L_MAX+1; l++)
  {
    FixPhases(l);
    FixWaves(l);
  }
  //Fix observable scalings
  FixScalings();

  //When fitting threshold data, fix remaining p waves to be real,...
  if(ONLY_CROSS_S || ONLY_CROSS_F) FixRealPWaves();
  //...,fix imaginary part of E0+ to its unitarity value,...
  if(ONLY_CROSS_S || ONLY_CROSS_F) FixUnitarity();
  //...fix real part of E0+ to model value
  FixThreshold();

  //Perform minimisation with MIGRAD technique
  arglist[0] = 1e+06;
  arglist[1] = 1e-03;
  gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  if(ierflg) return false; //MIGRAD failed, fit was unsuccessful

  //Calculate covariance matrix
  arglist[0] = 1e+06;
  gMinuit->mnexcm("HESSE", arglist, 1, ierflg);

  //Try to improve minimisation result
  arglist[0] = 1e+06;
  gMinuit->mnexcm("IMPROVE", arglist, 1, ierflg);

  //Get fit parameters and use them for multipoles
  GetParameters(Param, Error);
  UseParameters(Param, Error);

  /*Double_t CosTheta;
  Double_t E_i, E_f;
  Double_t m, m2, M_i, M2_i, M_f, M2_f;
  Double_t W, q, k, w;
  Double_t s, t;

  W = W_cm(gEnergy)/1e3;
  q = q_mes(gEnergy)/1e3;
  k = omega_cm(gEnergy)/1e3;

  m   = MASS_MESON/1e3;   m2   = MASS2_MESON/1e6;
  M_i = MASS_INITIAL/1e3; M2_i = MASS2_INITIAL/1e6;
  M_f = MASS_INITIAL/1e3; M2_f = MASS2_INITIAL/1e6;

  w = Sqrt(m2 + q*q);

  s = W*W;

  for(CosTheta = -1.0; CosTheta < +1.0; CosTheta+=0.01)
  {
    t = m2 - 2.0*k*(w - q*CosTheta);
    if((t < -0.149) && (t > -0.151))
    {
      printf("IA Ax: %f  %f     %f  %f     %f  %f     %f  %f     %f  %f\n", W, t,
             A1(CosTheta).Re(), A1(CosTheta).Im(), A2(CosTheta).Re(), A2(CosTheta).Im(),
             A3(CosTheta).Re(), A3(CosTheta).Im(), A4(CosTheta).Re(), A4(CosTheta).Im());
      printf("IA Bx: %f  %f     %f  %f     %f  %f     %f  %f     %f  %f\n", W, t,
             B1(CosTheta).Re(), B1(CosTheta).Im(), B2(CosTheta).Re(), B2(CosTheta).Im(),
             B6(CosTheta).Re(), B6(CosTheta).Im(), B8(CosTheta).Re(), B8(CosTheta).Im());
    }
  }*/

  return true; //Fit was successful
}

//-----------------------------------------------------------------------------

Double_t GetErrors(Double_t* Par, Double_t* Err)
{
  Double_t Mean[2];

  switch(ERROR_MODE)
  {
   case CHI2ONLY:
    return ErrorChi2(Par, Err);
   case CHI2PENALTY:
    return ErrorChi2Penalty(Par, Err);
   case ADAPTIVE:
   default:
    Mean[0] = ErrorChi2(Par, Err);
    Mean[1] = ErrorChi2Penalty(Par, Err);
    if(Mean[1] < Mean[0])
      return ErrorChi2Penalty(Par, Err);
    else
      return ErrorChi2(Par, Err);
  }
}

//-----------------------------------------------------------------------------

Double_t ErrorBase(Double_t* Par, Double_t* Err)
{
  Char_t Buffer[16];
  Double_t Mean = 0.0;

  //Use previously found multipole parameters from fit with penalty
  for(Int_t l=0; l<L_MAX*8; l++)
  {
    sprintf(Buffer, "Param%d", l);
    gMinuit->mnparm(l, Buffer, Par[l], step[l], 0, 0, ierflg);
  }
  //Use previously found observable scale parameters from fit with penalty
  for(Int_t o=0; o<OBS; o++)
  {
    sprintf(Buffer, "Scale%d", o);
    gMinuit->mnparm(L_MAX*8 + o, Buffer, Par[L_MAX*8 + o], step[L_MAX*8 + o], 0, 0, ierflg);
  }

  //Fix s,p waves ore phases during the fitting process
  FixSPPhases();
  FixSPWaves();
  //Fix d,f,... waves or phases during the fitting process when requested (fit contains not enough data to constrain all waves simultaneously)
  for(Int_t l=2; l<L_MAX+1; l++)
  {
    FixPhases(l);
    FixWaves(l);
  }
  //Fix observabe scalings
  FixScalings();

  //When fitting threshold data, fix remaining p waves to be real,...
  if(ONLY_CROSS_S || ONLY_CROSS_F) FixRealPWaves();
  //...,fix imaginary part of E0+ to its unitarity value,...
  if(ONLY_CROSS_S || ONLY_CROSS_F) FixUnitarity();
  //...fix real part of E0+ to model value
  FixThreshold();

  //Calculate covariance matrix
  arglist[0] = 1e+06;
  gMinuit->mnexcm("HESSE", arglist, 1, ierflg);

  //Get fit parameters and use them for multipoles
  GetParameters(Par, Err);
  UseParameters(Par, Err);

  //Calculate mean error over all multipole parameters
  for(Int_t l=0; l<L_MAX*8; l++)
    Mean+=Err[l];

  return Mean/NMlp();
}

//-----------------------------------------------------------------------------

void Fit()
{
  Double_t Quality[ITERATIONS];
  Double_t MeanErr[ITERATIONS];
  Double_t Param[ITERATIONS][L_MAX*8 + OBS];
  Double_t Error[ITERATIONS][L_MAX*8 + OBS];
  Int_t Attempt[ITERATIONS];
  Int_t Sol[SOLUTIONS+1];
  Int_t Index = 0; //Indicates which iterations have already been processed
  Int_t Unique = 1; //Counts number of found unique solutions

  //Perform multiple fits for given energy point

  for(Int_t n=0; n<ITERATIONS; n++)
  {
    //Try to perform successful single fit with penalty (stop after finite number of attempts)
    for(Int_t m=0; m<FINITE; m++)
      if(SingleFit()) break;

    //Get parameters, true errors (only chi^2), and fit quality (chi^2+penalty) for single fit
    GetParameters(Param[n], Error[n]);
    MeanErr[n] = GetErrors(Param[n], Error[n]);
    Quality[n] = ChiSq() + Penalty();
    Attempt[n] = n;
  }
  //Sort fits w.r.t. fit quality (i.e. chi^2+penalty)
  SortFits(Quality, MeanErr, Attempt, 0, ITERATIONS-1);

  //Pick up requested number of unique solutions, each one with lowest mean error
  for(Int_t t=0; t<SOLUTIONS; t++) Sol[t] = 0; //Start with absolutely lowest chi^2+penalty iteration
  for(Int_t s=0; s<SOLUTIONS; s++) //Collect requested number of unique solutions
  {
    //Check all iterations that may belong to current solution or upcoming next-best solutions
    for(Int_t n=Index; n<ITERATIONS; n++)
    {
      //Consider only iterations close to best chi^2+penalty of current solution
      if(Quality[n]-Quality[Sol[s]] < 1e-07)
      {
        if(MeanErr[n] < MeanErr[Sol[s]]) //An iteration with lower mean error has been found
          for(Int_t t=s; t<SOLUTIONS; t++) Sol[t] = n; //Select this iteration for current and all upcoming solutions (in case of no other solutions to be found)
      }
      //Fit quality changes significantly in a new group of iterations
      if(Quality[n]-Quality[Sol[s]] > 1e-02)
      {
        for(Int_t t=s+1; t<SOLUTIONS; t++) Sol[t] = n; //Use this iteration as preliminary result for the next and all upcoming solutions (in case of no other solutions to be found)
        Index = n; //In next loop, begin at first iteration with significantly changed quality...
        break; //...and start over with processsing iterations for the next-best solution
      }
    }
  }
  //Determine number of found unique solutions
  for(Int_t n=1; n<SOLUTIONS; n++)
    if(Sol[n]!=Sol[n-1]) Unique++;

  //Apply and store unique solutions (in reverse order, so that best solution is applied at the end).
  //This will make the best solution in Ep[], Em[], Mp[], Mm[] available for further use after fitting is finished.
  for(Int_t s=Unique-1; s>-1; s--)
  {
    //Use fit parameters and errors from current solution for multipoles
    UseParameters(Param[Attempt[Sol[s]]], Error[Attempt[Sol[s]]]);
    //Store fit parameters and errors from current solution for graph plots
    StoreFit(s);
  }
  //Store model parameters for graph plots
  StoreModel();
  for(Int_t l=0; l<L_MAX+1; l++)
    {
      Ep_prev[l]=Ep[l];
      Mp_prev[l]=Mp[l];
      Em_prev[l]=Em[l];      
      Mm_prev[l]=Mm[l];
    }


  return;

  //Debug output
  for(Int_t n=0; n<ITERATIONS; n++)
    printf("%4d %4d: chi^2+penalty = %19.15f, distance = %e, error = %11.8f\n", n, Attempt[n], Quality[n], Quality[n]-Quality[0], MeanErr[n]);
  for(Int_t s=0; s<SOLUTIONS; s++)
    printf("Solution %d is %d\n", s, Sol[s]);
}

//-----------------------------------------------------------------------------

void SetSPWaves()
{
  Int_t eM = GetEnergyBin_maid();

  //Set starting values for parameters (start with MAID solution)
  if(FIX_EP_PHASE[0])
  {
    vstart[0] = maid_Ep[0][eM].Rho();   //s wave, expressed in magnitude/phase
    vstart[1] = maid_Ep[0][eM].Theta(); //s wave, expressed in magnitude/phase
  }
  else
  {
    vstart[0] = maid_Ep[0][eM].Re(); //s wave
    vstart[1] = maid_Ep[0][eM].Im(); //s wave
  }

  if(FIX_EP_PHASE[1])
  {
    vstart[2] = maid_Ep[1][eM].Rho();   //p wave, expressed in magnitude/phase
    vstart[3] = maid_Ep[1][eM].Theta(); //p wave, expressed in magnitude/phase
  }
  else
  {
    vstart[2] = maid_Ep[1][eM].Re(); //p wave
    vstart[3] = maid_Ep[1][eM].Im(); //p wave
  }

  if(FIX_MP_PHASE[1])
  {
    vstart[4] = maid_Mp[1][eM].Rho();   //p wave, expressed in magnitude/phase
    vstart[5] = maid_Mp[1][eM].Theta(); //p wave, expressed in magnitude/phase
  }
  else
  {
    vstart[4] = maid_Mp[1][eM].Re(); //p wave
    vstart[5] = maid_Mp[1][eM].Im(); //p wave
  }

  if(FIX_MM_PHASE[1])
  {
    vstart[6] = maid_Mm[1][eM].Rho();   //p wave, expressed in magnitude/phase
    vstart[7] = maid_Mm[1][eM].Theta(); //p wave, expressed in magnitude/phase
  }
  else
  {
    vstart[6] = maid_Mm[1][eM].Re(); //p wave
    vstart[7] = maid_Mm[1][eM].Im(); //p wave
  }
}

//-----------------------------------------------------------------------------

void SetWaves(Int_t l)
{
  Int_t eM = GetEnergyBin_maid();
  Int_t Index = 8*(l-1);

  if(FIX_EP_PHASE[l])
  {
    vstart[Index+0] = maid_Ep[l][eM].Rho();
    vstart[Index+1] = maid_Ep[l][eM].Theta();
  }
  else
  {
    vstart[Index+0] = maid_Ep[l][eM].Re();
    vstart[Index+1] = maid_Ep[l][eM].Im();
  }

  if(FIX_MP_PHASE[l])
  {
    vstart[Index+2] = maid_Mp[l][eM].Rho();
    vstart[Index+3] = maid_Mp[l][eM].Theta();
  }
  else
  {
    vstart[Index+2] = maid_Mp[l][eM].Re();
    vstart[Index+3] = maid_Mp[l][eM].Im();
  }

  if(FIX_EM_PHASE[l])
  {
    vstart[Index+4] = maid_Em[l][eM].Rho();
    vstart[Index+5] = maid_Em[l][eM].Theta();
  }
  else
  {
    vstart[Index+4] = maid_Em[l][eM].Re();
    vstart[Index+5] = maid_Em[l][eM].Im();
  }

  if(FIX_MM_PHASE[l])
  {
    vstart[Index+6] = maid_Mm[l][eM].Rho();
    vstart[Index+7] = maid_Mm[l][eM].Theta();
  }
  else
  {
    vstart[Index+6] = maid_Mm[l][eM].Re();
    vstart[Index+7] = maid_Mm[l][eM].Im();
  }
}

//-----------------------------------------------------------------------------

void SetRealPWaves()
{
  //Keep p waves real (E1+, M1+, M1-) for threshold data fits
  vstart[3] = 0.0; //ImE1+
  vstart[5] = 0.0; //ImM1+
  vstart[7] = 0.0; //ImM1-
}

//-----------------------------------------------------------------------------

void SetUnitarity()
{
  //Parametrise imaginary part of E0+ using unitarity for threshold data fits
  vstart[1] = ImE0p();
}

//-----------------------------------------------------------------------------

void FixSPWaves()
{
  //Fix phases (global phase needs to be fixed) during the fitting process to MAID values
  if(FIX_EP[0]) { gMinuit->FixParameter(0); gMinuit->FixParameter(1); } //E0+
  if(FIX_EP[1]) { gMinuit->FixParameter(2); gMinuit->FixParameter(3); } //E1+
  if(FIX_MP[1]) { gMinuit->FixParameter(4); gMinuit->FixParameter(5); } //M1+
  if(FIX_MM[1]) { gMinuit->FixParameter(6); gMinuit->FixParameter(7); } //M1-
}

//-----------------------------------------------------------------------------

void FixSPPhases()
{
  //Fix phases (global phase needs to be fixed) during the fitting process to MAID values
  if(FIX_EP_PHASE[0]) gMinuit->FixParameter(1); //E0+ phase
  if(FIX_EP_PHASE[1]) gMinuit->FixParameter(3); //E1+ phase
  if(FIX_MP_PHASE[1]) gMinuit->FixParameter(5); //M1+ phase
  if(FIX_MM_PHASE[1]) gMinuit->FixParameter(7); //M1- phase
}

//-----------------------------------------------------------------------------

void FixWaves(Int_t l)
{
  Int_t Index = 8*(l-1);
  //Fix L=l waves during the fitting process
  if(FIX_EP[l]) { gMinuit->FixParameter(Index+0); gMinuit->FixParameter(Index+1); } //El+
  if(FIX_MP[l]) { gMinuit->FixParameter(Index+2); gMinuit->FixParameter(Index+3); } //Ml+
  if(FIX_EM[l]) { gMinuit->FixParameter(Index+4); gMinuit->FixParameter(Index+5); } //El-
  if(FIX_MM[l]) { gMinuit->FixParameter(Index+6); gMinuit->FixParameter(Index+7); } //Ml-
}

//-----------------------------------------------------------------------------

void FixPhases(Int_t l)
{
  Int_t Index = 8*(l-1);
  //Fix L=l phases during the fitting process
  if(FIX_EP_PHASE[l]) gMinuit->FixParameter(Index+1); //El+ phase
  if(FIX_MP_PHASE[l]) gMinuit->FixParameter(Index+3); //Ml+ phase
  if(FIX_EM_PHASE[l]) gMinuit->FixParameter(Index+5); //El- phase
  if(FIX_MM_PHASE[l]) gMinuit->FixParameter(Index+7); //Ml- phase
}

//-----------------------------------------------------------------------------

void FixRealPWaves()
{
  //Fix imaginary parts of E1+, M1+, and M1- for threshold data fits
  gMinuit->FixParameter(3); //ImE1+
  gMinuit->FixParameter(5); //ImM1+
  gMinuit->FixParameter(7); //ImM1-
}

//-----------------------------------------------------------------------------

void FixUnitarity()
{
  //Fix imaginary part of E0+ for threshold data fits
  if(FIX_IM_E0P) gMinuit->FixParameter(1);
}

//-----------------------------------------------------------------------------

void FixThreshold()
{
  //Fix real part of E0+ for threshold data fits
  if(FIX_RE_E0P) gMinuit->FixParameter(0);
}

//-----------------------------------------------------------------------------

void SetParameters()
{
  Double_t VariatRel[2];
  Double_t VariatAbs[2];
  Char_t Buffer[64];
  Int_t Index;

  //Set start values for s,p waves (depending on Re/Im or magnitude/phase parametrisation)
  if(FIX_EP[0]) //Random variations if this multipole is fitted
  { VariatRel[0] = 1.0;          VariatRel[1] = 1.0;          VariatAbs[0] = 0.0;          VariatAbs[1] = 0.0;          }
  else
  { VariatRel[0] = VariateRel(); VariatRel[1] = VariateRel(); VariatAbs[0] = VariateAbs(); VariatAbs[1] = VariateAbs(); }
  if(FIX_EP_PHASE[0])
  {
    gMinuit->mnparm(0, "MagE0p", vstart[0]*VariatRel[0]+VariatAbs[0], step[0], 0, 0, ierflg); //s wave
    gMinuit->mnparm(1, "PhsE0p", vstart[1],                           step[1], 0, 0, ierflg);
  }
  else
  {
    gMinuit->mnparm(0, "ReE0p",  vstart[0]*VariatRel[0]+VariatAbs[0], step[0], 0, 0, ierflg); //s wave
    gMinuit->mnparm(1, "ImE0p",  vstart[1]*VariatRel[1]+VariatAbs[1], step[1], 0, 0, ierflg);
  }

  if(FIX_EP[1]) //Random variations if this multipole is fitted
  { VariatRel[0] = 1.0;          VariatRel[1] = 1.0;          VariatAbs[0] = 0.0;          VariatAbs[1] = 0.0;          }
  else
  { VariatRel[0] = VariateRel(); VariatRel[1] = VariateRel(); VariatAbs[0] = VariateAbs(); VariatAbs[1] = VariateAbs(); }
  if(FIX_EP_PHASE[1])
  {
    gMinuit->mnparm(2, "MagE1p", vstart[2]*VariatRel[0]+VariatAbs[0], step[2], 0, 0, ierflg); //p wave
    gMinuit->mnparm(3, "PhsE1p", vstart[3],                           step[3], 0, 0, ierflg);
  }
  else
  {
    gMinuit->mnparm(2, "ReE1p",  vstart[2]*VariatRel[0]+VariatAbs[0], step[2], 0, 0, ierflg); //p wave
    gMinuit->mnparm(3, "ImE1p",  vstart[3]*VariatRel[1]+VariatAbs[1], step[3], 0, 0, ierflg);
  }

  if(FIX_MP[1]) //Random variations if this multipole is fitted
  { VariatRel[0] = 1.0;          VariatRel[1] = 1.0;          VariatAbs[0] = 0.0;          VariatAbs[1] = 0.0;          }
  else
  { VariatRel[0] = VariateRel(); VariatRel[1] = VariateRel(); VariatAbs[0] = VariateAbs(); VariatAbs[1] = VariateAbs(); }
  if(FIX_MP_PHASE[1])
  {
    gMinuit->mnparm(4, "MagM1p", vstart[4]*VariatRel[0]+VariatAbs[0], step[4], 0, 0, ierflg); //p wave
    gMinuit->mnparm(5, "PhsM1p", vstart[5],                           step[5], 0, 0, ierflg);
  }
  else
  {
    gMinuit->mnparm(4, "ReM1m",  vstart[4]*VariatRel[0]+VariatAbs[0], step[4], 0, 0, ierflg); //p wave
    gMinuit->mnparm(5, "ImM1m",  vstart[5]*VariatRel[1]+VariatAbs[1], step[5], 0, 0, ierflg);
  }

  if(FIX_MM[1]) //Random variations if this multipole is fitted
  { VariatRel[0] = 1.0;          VariatRel[1] = 1.0;          VariatAbs[0] = 0.0;          VariatAbs[1] = 0.0;          }
  else
  { VariatRel[0] = VariateRel(); VariatRel[1] = VariateRel(); VariatAbs[0] = VariateAbs(); VariatAbs[1] = VariateAbs(); }
  if(FIX_MM_PHASE[1])
  {
    gMinuit->mnparm(6, "MagM1m", vstart[6]*VariatRel[0]+VariatAbs[0], step[6], 0, 0, ierflg); //p wave
    gMinuit->mnparm(7, "PhsM1m", vstart[7],                           step[7], 0, 0, ierflg);
  }
  else
  {
    gMinuit->mnparm(6, "ReM1m",  vstart[6]*VariatRel[0]+VariatAbs[0], step[6], 0, 0, ierflg); //p wave
    gMinuit->mnparm(7, "ImM1m",  vstart[7]*VariatRel[1]+VariatAbs[1], step[7], 0, 0, ierflg);
  }

  //No variation on imaginary parts of s,p waves for threshold fits
  if(FIX_RE_E0P) //ReE0+ fixed
    gMinuit->mnparm(0, "ReE0p",  vstart[0], step[0], 0, 0, ierflg);
  if(FIX_IM_E0P) //ImE0+ fixed to unitary parametrisation
    gMinuit->mnparm(1, "ImE0p",  vstart[1], step[1], 0, 0, ierflg);
  if(ONLY_CROSS_S || ONLY_CROSS_F) //Real p waves
  {
    gMinuit->mnparm(3, "ImE1p",  vstart[3], step[3], 0, 0, ierflg);
    gMinuit->mnparm(5, "ImM1p",  vstart[5], step[5], 0, 0, ierflg);
    gMinuit->mnparm(7, "ImM1m",  vstart[7], step[7], 0, 0, ierflg);
  }

  //Set start values for d,f,... waves (depending on Re/Im or magnitude/phase parametrisation)
  for(Int_t l=2; l<L_MAX+1; l++)
  {
    Index = 8*(l-1);

    if(FIX_EP[l]) //Random variations if this multipole is fitted
    { VariatRel[0] = 1.0;          VariatRel[1] = 1.0;          VariatAbs[0] = 0.0;          VariatAbs[1] = 0.0;          }
    else
    { VariatRel[0] = VariateRel(); VariatRel[1] = VariateRel(); VariatAbs[0] = VariateAbs(); VariatAbs[1] = VariateAbs(); }
    if(FIX_EP_PHASE[l])
    {
      sprintf(Buffer, "MagE%dp", l); gMinuit->mnparm(0+Index, Buffer, vstart[0+Index]*VariatRel[0]+VariatAbs[0], step[0+Index], 0, 0, ierflg);
      sprintf(Buffer, "PhsE%dp", l); gMinuit->mnparm(1+Index, Buffer, vstart[1+Index],                           step[1+Index], 0, 0, ierflg);
    }
    else
    {
      sprintf(Buffer, "ReE%dp",  l); gMinuit->mnparm(0+Index, Buffer, vstart[0+Index]*VariatRel[0]+VariatAbs[0], step[0+Index], 0, 0, ierflg);
      sprintf(Buffer, "ImE%dp",  l); gMinuit->mnparm(1+Index, Buffer, vstart[1+Index]*VariatRel[1]+VariatAbs[1], step[1+Index], 0, 0, ierflg);
    }

    if(FIX_MP[l]) //Random variations if this multipole is fitted
    { VariatRel[0] = 1.0;          VariatRel[1] = 1.0;          VariatAbs[0] = 0.0;          VariatAbs[1] = 0.0;          }
    else
    { VariatRel[0] = VariateRel(); VariatRel[1] = VariateRel(); VariatAbs[0] = VariateAbs(); VariatAbs[1] = VariateAbs(); }
    if(FIX_MP_PHASE[l])
    {
      sprintf(Buffer, "MagM%dp", l); gMinuit->mnparm(2+Index, Buffer, vstart[2+Index]*VariatRel[0]+VariatAbs[0], step[2+Index], 0, 0, ierflg);
      sprintf(Buffer, "PhsM%dp", l); gMinuit->mnparm(3+Index, Buffer, vstart[3+Index],                           step[3+Index], 0, 0, ierflg);
    }
    else
    {
      sprintf(Buffer, "ReM%dp",  l); gMinuit->mnparm(2+Index, Buffer, vstart[2+Index]*VariatRel[0]+VariatAbs[0], step[2+Index], 0, 0, ierflg);
      sprintf(Buffer, "ImM%dp",  l); gMinuit->mnparm(3+Index, Buffer, vstart[3+Index]*VariatRel[1]+VariatAbs[1], step[3+Index], 0, 0, ierflg);
    }

    if(FIX_EM[l]) //Random variations if this multipole is fitted
    { VariatRel[0] = 1.0;          VariatRel[1] = 1.0;          VariatAbs[0] = 0.0;          VariatAbs[1] = 0.0;          }
    else
    { VariatRel[0] = VariateRel(); VariatRel[1] = VariateRel(); VariatAbs[0] = VariateAbs(); VariatAbs[1] = VariateAbs(); }
    if(FIX_EM_PHASE[l])
    {
      sprintf(Buffer, "MagE%dm", l); gMinuit->mnparm(4+Index, Buffer, vstart[4+Index]*VariatRel[0]+VariatAbs[0], step[4+Index], 0, 0, ierflg);
      sprintf(Buffer, "PhsE%dm", l); gMinuit->mnparm(5+Index, Buffer, vstart[5+Index],                           step[5+Index], 0, 0, ierflg);
    }
    else
    {
      sprintf(Buffer, "ReE%dm",  l); gMinuit->mnparm(4+Index, Buffer, vstart[4+Index]*VariatRel[0]+VariatAbs[0], step[4+Index], 0, 0, ierflg);
      sprintf(Buffer, "ImE%dm",  l); gMinuit->mnparm(5+Index, Buffer, vstart[5+Index]*VariatRel[1]+VariatAbs[1], step[5+Index], 0, 0, ierflg);
    }

    if(FIX_MM[l]) //Random variations if this multipole is fitted
    { VariatRel[0] = 1.0;          VariatRel[1] = 1.0;          VariatAbs[0] = 0.0;          VariatAbs[1] = 0.0;          }
    else
    { VariatRel[0] = VariateRel(); VariatRel[1] = VariateRel(); VariatAbs[0] = VariateAbs(); VariatAbs[1] = VariateAbs(); }
    if(FIX_MM_PHASE[l])
    {
      sprintf(Buffer, "MagM%dm", l); gMinuit->mnparm(6+Index, Buffer, vstart[6+Index]*VariatRel[0]+VariatAbs[0], step[6+Index], 0, 0, ierflg);
      sprintf(Buffer, "PhsM%dm", l); gMinuit->mnparm(7+Index, Buffer, vstart[7+Index],                           step[7+Index], 0, 0, ierflg);
    }
    else
    {
      sprintf(Buffer, "ReM%dm",  l); gMinuit->mnparm(6+Index, Buffer, vstart[6+Index]*VariatRel[0]+VariatAbs[0], step[6+Index], 0, 0, ierflg);
      sprintf(Buffer, "ImM%dm",  l); gMinuit->mnparm(7+Index, Buffer, vstart[7+Index]*VariatRel[1]+VariatAbs[1], step[7+Index], 0, 0, ierflg);
    }
  }

  //Set start values for observable scale parameters
  for(Int_t o=0; o<OBS; o++)
  {
    vstart[8*L_MAX + o] = 1.0;
    sprintf(Buffer, "f_%d", o); gMinuit->mnparm(8*L_MAX + o, Buffer, vstart[8*L_MAX + o], step[8*L_MAX + o], 0, 0, ierflg);
  }
}

//-----------------------------------------------------------------------------

void GetParameters(Double_t* Param, Double_t* Error)
{
  for(Int_t p=0; p<L_MAX*8; p++) gMinuit->GetParameter(p, Param[p], Error[p]);
  for(Int_t o=0; o<OBS; o++) gMinuit->GetParameter(8*L_MAX + o, Param[8*L_MAX + o], Error[8*L_MAX + o]);
}

//-----------------------------------------------------------------------------

void UseParameters(Double_t* Param, Double_t* Error)
{
  Int_t Index;

  //Get s,p wave parameters (depending on Re/Im or magnitude/phase parametrisation)
  if(FIX_EP_PHASE[0])
  {
    Ep[0]  = TComplex(fabs(Param[0]), Param[1], true);
    DEp[0] = TComplex(fabs(Error[0]), Param[1], true); //Only valid if phase has no error!
  }
  else
  {
    Ep[0]  = TComplex(Param[0], Param[1]);
    DEp[0] = TComplex(fabs(Error[0]), fabs(Error[1]));
  }

  if(FIX_EP_PHASE[1])
  {
    Ep[1]  = TComplex(fabs(Param[2]), Param[3], true);
    DEp[1] = TComplex(fabs(Error[2]), Param[3], true); //Only valid if phase has no error!
  }
  else
  {
    Ep[1]  = TComplex(Param[2], Param[3]);
    DEp[1] = TComplex(fabs(Error[2]), fabs(Error[3]));
  }

  if(FIX_MP_PHASE[1])
  {
    Mp[1]  = TComplex(fabs(Param[4]), Param[5], true);
    DMp[1] = TComplex(fabs(Error[4]), Param[5], true); //Only valid if phase has no error!
  }
  else
  {
    Mp[1]  = TComplex(Param[4], Param[5]);
    DMp[1] = TComplex(fabs(Error[4]), fabs(Error[5]));
  }

  if(FIX_MM_PHASE[1])
  {
    Mm[1]  = TComplex(fabs(Param[6]), Param[7], true);
    DMm[1] = TComplex(fabs(Error[6]), Param[7], true); //Only valid if phase has no error!
  }
  else
  {
    Mm[1]  = TComplex(Param[6], Param[7]);
    DMm[1] = TComplex(fabs(Error[6]), fabs(Error[7]));
  }
  //Make errors on multipoles positive
  DEp[0] = TComplex(fabs(DEp[0].Re()), fabs(DEp[0].Im()));
  DEp[1] = TComplex(fabs(DEp[1].Re()), fabs(DEp[1].Im()));
  DMp[1] = TComplex(fabs(DMp[1].Re()), fabs(DMp[1].Im()));
  DMm[1] = TComplex(fabs(DMm[1].Re()), fabs(DMm[1].Im()));

  //Get d,f,... wave parameters
  for(Int_t l=2; l<L_MAX+1; l++)
  {
    Index = 8*(l-1);
    if(FIX_EP_PHASE[l])
    {
      Ep[l]  = TComplex(fabs(Param[0+Index]), Param[1+Index], true);
      DEp[l] = TComplex(fabs(Error[0+Index]), Param[1+Index], true); //Only valid if phase has no error!
    }
    else
    {
      Ep[l]  = TComplex(Param[0+Index], Param[1+Index]);
      DEp[l] = TComplex(fabs(Error[0+Index]), fabs(Error[1+Index]));
    }
    if(FIX_MP_PHASE[l])
    {
      Mp[l]  = TComplex(fabs(Param[2+Index]), Param[3+Index], true);
      DMp[l] = TComplex(fabs(Error[2+Index]), Param[3+Index], true); //Only valid if phase has no error!
    }
    else
    {
      Mp[l]  = TComplex(Param[2+Index], Param[3+Index]);
      DMp[l] = TComplex(fabs(Error[2+Index]), fabs(Error[3+Index]));
    }
    if(FIX_EM_PHASE[l])
    {
      Em[l]  = TComplex(fabs(Param[4+Index]), Param[5+Index], true);
      DEm[l] = TComplex(fabs(Error[4+Index]), Param[5+Index], true); //Only valid if phase has no error!
    }
    else
    {
      Em[l]  = TComplex(Param[4+Index], Param[5+Index]);
      DEm[l] = TComplex(fabs(Error[4+Index]), fabs(Error[5+Index]));
    }
    if(FIX_MM_PHASE[l])
    {
      Mm[l]  = TComplex(fabs(Param[6+Index]), Param[7+Index], true);
      DMm[l] = TComplex(fabs(Error[6+Index]), Param[7+Index], true); //Only valid if phase has no error!
    }
    else
    {
      Mm[l]  = TComplex(Param[6+Index], Param[7+Index]);
      DMm[l] = TComplex(fabs(Error[6+Index]), fabs(Error[7+Index]));
    }
    //Make errors on multipoles positive
    DEp[l] = TComplex(fabs(DEp[l].Re()), fabs(DEp[l].Im()));
    DMp[l] = TComplex(fabs(DMp[l].Re()), fabs(DMp[l].Im()));
    DEm[l] = TComplex(fabs(DEm[l].Re()), fabs(DEm[l].Im()));
    DMm[l] = TComplex(fabs(DMm[l].Re()), fabs(DMm[l].Im()));
  }

  //Get observable scale parameters
  for(Int_t o=0; o<OBS; o++)
  {
    f_obs[o]  = Param[8*L_MAX + o];
    Df_obs[o] = Error[8*L_MAX + o];
  }
}

//-----------------------------------------------------------------------------

void InitMinuit()
{
  //Delete existing instance and create new Minuit fitter
  if(gMinuit) gMinuit->Delete();
  gMinuit = new TMinuit(L_MAX*8 + OBS);

  //Set to quiet mode
  gMinuit->SetPrintLevel(-1);

  //Set step sizes up to maximum number of multipole parameters
  for(Int_t l=0; l<L_MAX*8; l++)
    step[l] = 1e-06;
  //Set step sizes up to maximum number of observable scale parameters
  for(Int_t o=0; o<OBS; o++)
    step[L_MAX*8 + o] = 1e-06;

  //Set error definition (1.0 is recommended in MINUIT manual for chi^2 fits)
  arglist[0] = 1.0;
  gMinuit->mnexcm("SET ERRORDEF", arglist, 1, ierflg);

  //Set minimisation strategy (2 should provide more reliable results)
  arglist[0] = 2;
  gMinuit->mnexcm("SET STRATEGY", arglist, 1, ierflg);
}

//-----------------------------------------------------------------------------

void SortFits(Double_t* Quality, Double_t* MeanErr, Int_t* Attempt, Int_t l, Int_t r) //Quicksort implementation
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

    for(;;)
    {
      while(Quality[++i] < Quality[r]);
      while((Quality[--j] > Quality[r]) && (j>i));
      if(i>=j) break;
      Swap(&Attempt[i], &Attempt[j]);
      Swap(&Quality[i], &Quality[j]);
      Swap(&MeanErr[i], &MeanErr[j]);
    }
    Swap(&Attempt[i], &Attempt[r]);
    Swap(&Quality[i], &Quality[r]);
    Swap(&MeanErr[i], &MeanErr[r]);

    SortFits(Quality, MeanErr, Attempt, l, i-1);
    SortFits(Quality, MeanErr, Attempt, i+1, r);
  }
}

//-----------------------------------------------------------------------------

void StoreFit(Int_t Sol)
{
  Int_t Index;

  Fit_en[Sol][Fit_pts[Sol]] = gEnergy;

  //Store fit data for graph plots
  Fit_val[Sol][0][Fit_pts[Sol]] = Ep[0].Re(); Fit_val[Sol][1][Fit_pts[Sol]] = Ep[0].Im(); Fit_err[Sol][0][Fit_pts[Sol]] = DEp[0].Re(); Fit_err[Sol][1][Fit_pts[Sol]] = DEp[0].Im();
  Fit_val[Sol][2][Fit_pts[Sol]] = Ep[1].Re(); Fit_val[Sol][3][Fit_pts[Sol]] = Ep[1].Im(); Fit_err[Sol][2][Fit_pts[Sol]] = DEp[1].Re(); Fit_err[Sol][3][Fit_pts[Sol]] = DEp[1].Im();
  Fit_val[Sol][4][Fit_pts[Sol]] = Mp[1].Re(); Fit_val[Sol][5][Fit_pts[Sol]] = Mp[1].Im(); Fit_err[Sol][4][Fit_pts[Sol]] = DMp[1].Re(); Fit_err[Sol][5][Fit_pts[Sol]] = DMp[1].Im();
  Fit_val[Sol][6][Fit_pts[Sol]] = Mm[1].Re(); Fit_val[Sol][7][Fit_pts[Sol]] = Mm[1].Im(); Fit_err[Sol][6][Fit_pts[Sol]] = DMm[1].Re(); Fit_err[Sol][7][Fit_pts[Sol]] = DMm[1].Im();
  for(Int_t l=2; l<L_MAX+1; l++)
  {
    Index = 8*(l-1);
    Fit_val[Sol][Index+0][Fit_pts[Sol]] = Ep[l].Re(); Fit_val[Sol][Index+1][Fit_pts[Sol]] = Ep[l].Im(); Fit_err[Sol][Index+0][Fit_pts[Sol]] = DEp[l].Re(); Fit_err[Sol][Index+1][Fit_pts[Sol]] = DEp[l].Im();
    Fit_val[Sol][Index+2][Fit_pts[Sol]] = Mp[l].Re(); Fit_val[Sol][Index+3][Fit_pts[Sol]] = Mp[l].Im(); Fit_err[Sol][Index+2][Fit_pts[Sol]] = DMp[l].Re(); Fit_err[Sol][Index+3][Fit_pts[Sol]] = DMp[l].Im();
    Fit_val[Sol][Index+4][Fit_pts[Sol]] = Em[l].Re(); Fit_val[Sol][Index+5][Fit_pts[Sol]] = Em[l].Im(); Fit_err[Sol][Index+4][Fit_pts[Sol]] = DEm[l].Re(); Fit_err[Sol][Index+5][Fit_pts[Sol]] = DEm[l].Im();
    Fit_val[Sol][Index+6][Fit_pts[Sol]] = Mm[l].Re(); Fit_val[Sol][Index+7][Fit_pts[Sol]] = Mm[l].Im(); Fit_err[Sol][Index+6][Fit_pts[Sol]] = DMm[l].Re(); Fit_err[Sol][Index+7][Fit_pts[Sol]] = DMm[l].Im();
  }
  Fit_chi[Sol][Fit_pts[Sol]] = ChiSq()/NDF();
  Fit_pen[Sol][Fit_pts[Sol]] = Penalty()/NDF();
  Fit_pts[Sol]++;
}

//-----------------------------------------------------------------------------

void StoreModel()
{
  Int_t Index;
  Int_t eM = GetEnergyBin_maid();

  //Store model data for graph plots
  Model_val[0][Model_pts] = maid_Ep[0][eM].Re(); Model_val[1][Model_pts] = maid_Ep[0][eM].Im();
  Model_val[2][Model_pts] = maid_Ep[1][eM].Re(); Model_val[3][Model_pts] = maid_Ep[1][eM].Im();
  Model_val[4][Model_pts] = maid_Mp[1][eM].Re(); Model_val[5][Model_pts] = maid_Mp[1][eM].Im();
  Model_val[6][Model_pts] = maid_Mm[1][eM].Re(); Model_val[7][Model_pts] = maid_Mm[1][eM].Im();
  for(Int_t l=2; l<L_MAX+1; l++)
  {
    Index = 8*(l-1);
    Model_val[Index+0][Model_pts] = maid_Ep[l][eM].Re(); Model_val[Index+1][Model_pts] = maid_Ep[l][eM].Im();
    Model_val[Index+2][Model_pts] = maid_Mp[l][eM].Re(); Model_val[Index+3][Model_pts] = maid_Mp[l][eM].Im();
    Model_val[Index+4][Model_pts] = maid_Em[l][eM].Re(); Model_val[Index+5][Model_pts] = maid_Em[l][eM].Im();
    Model_val[Index+6][Model_pts] = maid_Mm[l][eM].Re(); Model_val[Index+7][Model_pts] = maid_Mm[l][eM].Im();
  }
  Model_en[Model_pts] = maid_en[eM];
  Model_pts++;
}

//-----------------------------------------------------------------------------

void FixScalings()
{
  //Fix either all observable scales (when requested)...
  if(FIX_SCALES)
    for(Int_t o=0; o<OBS; o++) gMinuit->FixParameter(8*L_MAX + o);
  //...or fix observable scales for observables without data
  else
  {
    if(!sg0_pts[GetEnergyBins_sg0()])   gMinuit->FixParameter(8*L_MAX + SIG_0);
    if(!sgS_pts[GetEnergyBins_sgS()])   gMinuit->FixParameter(8*L_MAX + SIG_S);
    if(!sgT_pts[GetEnergyBins_sgT()])   gMinuit->FixParameter(8*L_MAX + SIG_T);
    if(!sgP_pts[GetEnergyBins_sgP()])   gMinuit->FixParameter(8*L_MAX + SIG_P);
    if(!sgE_pts[GetEnergyBins_sgE()])   gMinuit->FixParameter(8*L_MAX + SIG_E);
    if(!sgF_pts[GetEnergyBins_sgF()])   gMinuit->FixParameter(8*L_MAX + SIG_F);
    if(!sgG_pts[GetEnergyBins_sgG()])   gMinuit->FixParameter(8*L_MAX + SIG_G);
    if(!sgH_pts[GetEnergyBins_sgH()])   gMinuit->FixParameter(8*L_MAX + SIG_H);
    if(!sgCx_pts[GetEnergyBins_sgCx()]) gMinuit->FixParameter(8*L_MAX + SIG_CX);
    if(!sgCz_pts[GetEnergyBins_sgCz()]) gMinuit->FixParameter(8*L_MAX + SIG_CZ);
    if(!sgOx_pts[GetEnergyBins_sgOx()]) gMinuit->FixParameter(8*L_MAX + SIG_OX);
    if(!sgOz_pts[GetEnergyBins_sgOz()]) gMinuit->FixParameter(8*L_MAX + SIG_OZ);
    if(!sgCx_pts[GetEnergyBins_sgLx()]) gMinuit->FixParameter(8*L_MAX + SIG_LX);
    if(!sgCz_pts[GetEnergyBins_sgLz()]) gMinuit->FixParameter(8*L_MAX + SIG_LZ);
    if(!sgOx_pts[GetEnergyBins_sgTx()]) gMinuit->FixParameter(8*L_MAX + SIG_TX);
    if(!sgOz_pts[GetEnergyBins_sgTz()]) gMinuit->FixParameter(8*L_MAX + SIG_TZ);
    if(!S_pts[GetEnergyBins_S()])       gMinuit->FixParameter(8*L_MAX + ASY_S);
    if(!T_pts[GetEnergyBins_T()])       gMinuit->FixParameter(8*L_MAX + ASY_T);
    if(!P_pts[GetEnergyBins_P()])       gMinuit->FixParameter(8*L_MAX + ASY_P);
    if(!E_pts[GetEnergyBins_E()])       gMinuit->FixParameter(8*L_MAX + ASY_E);
    if(!F_pts[GetEnergyBins_F()])       gMinuit->FixParameter(8*L_MAX + ASY_F);
    if(!G_pts[GetEnergyBins_G()])       gMinuit->FixParameter(8*L_MAX + ASY_G);
    if(!H_pts[GetEnergyBins_H()])       gMinuit->FixParameter(8*L_MAX + ASY_H);
    if(!Cx_pts[GetEnergyBins_Cx()])     gMinuit->FixParameter(8*L_MAX + ASY_CX);
    if(!Cz_pts[GetEnergyBins_Cz()])     gMinuit->FixParameter(8*L_MAX + ASY_CZ);
    if(!Ox_pts[GetEnergyBins_Ox()])     gMinuit->FixParameter(8*L_MAX + ASY_OX);
    if(!Oz_pts[GetEnergyBins_Oz()])     gMinuit->FixParameter(8*L_MAX + ASY_OZ);
    if(!Cx_pts[GetEnergyBins_Lx()])     gMinuit->FixParameter(8*L_MAX + ASY_LX);
    if(!Cz_pts[GetEnergyBins_Lz()])     gMinuit->FixParameter(8*L_MAX + ASY_LZ);
    if(!Ox_pts[GetEnergyBins_Tx()])     gMinuit->FixParameter(8*L_MAX + ASY_TX);
    if(!Oz_pts[GetEnergyBins_Tz()])     gMinuit->FixParameter(8*L_MAX + ASY_TZ);
  }
}

//-----------------------------------------------------------------------------

