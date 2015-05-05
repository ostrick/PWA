#ifndef __Invariant__
#define __Invariant__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"
#include "Parse_MAID.h"

//-----------------------------------------------------------------------------

inline TComplex X(Int_t i, Int_t j, Double_t PlMi, Double_t CosTheta)
{
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
  t = m2 - 2.0*k*(w - q*CosTheta);

  E_i = W - k;
  E_f = W - w;

  return (1.0/(E_i-M_i))*F(i, CosTheta) + PlMi*(E_f-M_f)*F(j, CosTheta)/(k*q);
}

//-----------------------------------------------------------------------------

inline TComplex Y(Int_t i, Int_t j, Double_t PlMi, Double_t CosTheta)
{
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
  t = m2 - 2.0*k*(w - q*CosTheta);

  E_i = W - k;
  E_f = W - w;

  return ((W-M_i)/(E_i-M_i))*F(i, CosTheta) + PlMi*(E_f-M_f)*(W+M_f)*F(j, CosTheta)/(k*q);
}

//-----------------------------------------------------------------------------

inline TComplex A1(Double_t CosTheta)
{
  Double_t E_i, E_f;
  Double_t m, m2, M_i, M2_i, M_f, M2_f;
  Double_t W, q, k, w;
  Double_t s, t;
  Double_t N;

  W = W_cm(gEnergy)/1e3;
  q = q_mes(gEnergy)/1e3;
  k = omega_cm(gEnergy)/1e3;

  m   = MASS_MESON/1e3;   m2   = MASS2_MESON/1e6;
  M_i = MASS_INITIAL/1e3; M2_i = MASS2_INITIAL/1e6;
  M_f = MASS_INITIAL/1e3; M2_f = MASS2_INITIAL/1e6;

  w = Sqrt(m2 + q*q);

  s = W*W;
  t = m2 - 2.0*k*(w - q*CosTheta);

  E_i = W - k;
  E_f = W - w;

  N = TwoPi()/(W*Sqrt((E_i+M_i)*(E_f+M_f)));

  return N*((s-M2_i)*X(1,2,-1.0,CosTheta) + M_i*(t-m2)*X(3,4,+1.0,CosTheta)/q);
}

//-----------------------------------------------------------------------------

inline TComplex A2(Double_t CosTheta)
{
  Double_t E_i, E_f;
  Double_t m, m2, M_i, M2_i, M_f, M2_f;
  Double_t W, q, k, w;
  Double_t s, t;
  Double_t N;

  W = W_cm(gEnergy)/1e3;
  q = q_mes(gEnergy)/1e3;
  k = omega_cm(gEnergy)/1e3;

  m   = MASS_MESON/1e3;   m2   = MASS2_MESON/1e6;
  M_i = MASS_INITIAL/1e3; M2_i = MASS2_INITIAL/1e6;
  M_f = MASS_INITIAL/1e3; M2_f = MASS2_INITIAL/1e6;

  w = Sqrt(m2 + q*q);

  s = W*W;
  t = m2 - 2.0*k*(w - q*CosTheta);

  E_i = W - k;
  E_f = W - w;

  N = TwoPi()/(W*Sqrt((E_i+M_i)*(E_f+M_f)));

  return N*Y(3,4,-1.0,CosTheta)/q;
}

//-----------------------------------------------------------------------------

inline TComplex A3(Double_t CosTheta)
{
  Double_t E_i, E_f;
  Double_t m, m2, M_i, M2_i, M_f, M2_f;
  Double_t W, q, k, w;
  Double_t s, t;
  Double_t N;

  W = W_cm(gEnergy)/1e3;
  q = q_mes(gEnergy)/1e3;
  k = omega_cm(gEnergy)/1e3;

  m   = MASS_MESON/1e3;   m2   = MASS2_MESON/1e6;
  M_i = MASS_INITIAL/1e3; M2_i = MASS2_INITIAL/1e6;
  M_f = MASS_INITIAL/1e3; M2_f = MASS2_INITIAL/1e6;

  w = Sqrt(m2 + q*q);

  s = W*W;
  t = m2 - 2.0*k*(w - q*CosTheta);

  E_i = W - k;
  E_f = W - w;

  N = TwoPi()/(W*Sqrt((E_i+M_i)*(E_f+M_f)));

  return N*(Y(1,2,+1.0,CosTheta) + (t-m2-4.0*k*W)*X(3,4,+1.0,CosTheta)/(2.0*q));
}

//-----------------------------------------------------------------------------

inline TComplex A4(Double_t CosTheta)
{
  Double_t E_i, E_f;
  Double_t m, m2, M_i, M2_i, M_f, M2_f;
  Double_t W, q, k, w;
  Double_t s, t;
  Double_t N;

  W = W_cm(gEnergy)/1e3;
  q = q_mes(gEnergy)/1e3;
  k = omega_cm(gEnergy)/1e3;

  m   = MASS_MESON/1e3;   m2   = MASS2_MESON/1e6;
  M_i = MASS_INITIAL/1e3; M2_i = MASS2_INITIAL/1e6;
  M_f = MASS_INITIAL/1e3; M2_f = MASS2_INITIAL/1e6;

  w = Sqrt(m2 + q*q);

  s = W*W;
  t = m2 - 2.0*k*(w - q*CosTheta);

  E_i = W - k;
  E_f = W - w;

  N = TwoPi()/(W*Sqrt((E_i+M_i)*(E_f+M_f)));

  return N*(Y(1,2,+1.0,CosTheta) + (t-m2)*X(3,4,+1.0,CosTheta)/(2.0*q));
}

//-----------------------------------------------------------------------------

//Wrapper for explicit invariant amplitudes A_i
inline TComplex A(Int_t i, Double_t CosTheta)
{
  switch(i)
  {
    case 1: return A1(CosTheta);
    case 2: return A2(CosTheta);
    case 3: return A3(CosTheta);
    case 4: return A4(CosTheta);
    default: return TComplex(0.0, 0.0);
  }
}

//-----------------------------------------------------------------------------

inline TComplex A1cc(Double_t CosTheta){ return TComplex::Conjugate(A1(CosTheta)); }
inline TComplex A2cc(Double_t CosTheta){ return TComplex::Conjugate(A2(CosTheta)); }
inline TComplex A3cc(Double_t CosTheta){ return TComplex::Conjugate(A3(CosTheta)); }
inline TComplex A4cc(Double_t CosTheta){ return TComplex::Conjugate(A4(CosTheta)); }

//-----------------------------------------------------------------------------

//Wrapper for explicit invariant amplitudes A_i*
inline TComplex Acc(Int_t i, Double_t CosTheta)
{
  switch(i)
  {
    case 1: return A1cc(CosTheta);
    case 2: return A2cc(CosTheta);
    case 3: return A3cc(CosTheta);
    case 4: return A4cc(CosTheta);
    default: return TComplex(0.0, 0.0);
  }
}

//-----------------------------------------------------------------------------

inline TComplex B1(Double_t CosTheta)
{
  Double_t M_i = MASS_INITIAL/1e3;

  return A1(CosTheta) - 2.0*M_i*A4(CosTheta);
}

//-----------------------------------------------------------------------------

inline TComplex B2(Double_t CosTheta)
{
  Double_t m2 = MASS2_MESON/1e6;
  Double_t q = q_mes(gEnergy)/1e3;
  Double_t k = omega_cm(gEnergy)/1e3;
  Double_t w = Sqrt(m2 + q*q);
  Double_t t = m2 - 2.0*k*(w - q*CosTheta);

  return k*(w-q*CosTheta)*A2(CosTheta);
}

//-----------------------------------------------------------------------------

inline TComplex B6(Double_t CosTheta)
{
  Double_t M_i = MASS_INITIAL/1e3;

  return -2.0*A4(CosTheta);
}

//-----------------------------------------------------------------------------

inline TComplex B8(Double_t CosTheta)
{
  Double_t M_i = MASS_INITIAL/1e3;

  return -A3(CosTheta);
}

//-----------------------------------------------------------------------------

//Wrapper for explicit invariant amplitudes B_i
inline TComplex B(Int_t i, Double_t CosTheta)
{
  switch(i)
  {
    case 1: return B1(CosTheta);
    case 2: return B2(CosTheta);
    case 6: return B6(CosTheta);
    case 8: return B8(CosTheta);
    default: return TComplex(0.0, 0.0);
  }
}

//-----------------------------------------------------------------------------

inline TComplex B1cc(Double_t CosTheta){ return TComplex::Conjugate(B1(CosTheta)); }
inline TComplex B2cc(Double_t CosTheta){ return TComplex::Conjugate(B2(CosTheta)); }
inline TComplex B6cc(Double_t CosTheta){ return TComplex::Conjugate(B6(CosTheta)); }
inline TComplex B8cc(Double_t CosTheta){ return TComplex::Conjugate(B8(CosTheta)); }

//-----------------------------------------------------------------------------

//Wrapper for explicit invariant amplitudes A_i*
inline TComplex Bcc(Int_t i, Double_t CosTheta)
{
  switch(i)
  {
    case 1: return B1cc(CosTheta);
    case 2: return B2cc(CosTheta);
    case 6: return B6cc(CosTheta);
    case 8: return B8cc(CosTheta);
    default: return TComplex(0.0, 0.0);
  }
}

//-----------------------------------------------------------------------------

inline TComplex maid_X(Int_t i, Int_t j, Double_t PlMi, Double_t CosTheta, Int_t e)
{
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
  t = m2 - 2.0*k*(w - q*CosTheta);

  E_i = W - k;
  E_f = W - w;

  return (1.0/(E_i-M_i))*maid_F(i, CosTheta, e) + PlMi*(E_f-M_f)*maid_F(j, CosTheta, e)/(k*q);
}

//-----------------------------------------------------------------------------

inline TComplex maid_Y(Int_t i, Int_t j, Double_t PlMi, Double_t CosTheta, Int_t e)
{
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
  t = m2 - 2.0*k*(w - q*CosTheta);

  E_i = W - k;
  E_f = W - w;

  return ((W-M_i)/(E_i-M_i))*maid_F(i, CosTheta, e) + PlMi*(E_f-M_f)*(W+M_f)*maid_F(j, CosTheta, e)/(k*q);
}

//-----------------------------------------------------------------------------

inline TComplex maid_A1(Double_t CosTheta, Int_t e)
{
  Double_t E_i, E_f;
  Double_t m, m2, M_i, M2_i, M_f, M2_f;
  Double_t W, q, k, w;
  Double_t s, t;
  Double_t N;

  W = W_cm(gEnergy)/1e3;
  q = q_mes(gEnergy)/1e3;
  k = omega_cm(gEnergy)/1e3;

  m   = MASS_MESON/1e3;   m2   = MASS2_MESON/1e6;
  M_i = MASS_INITIAL/1e3; M2_i = MASS2_INITIAL/1e6;
  M_f = MASS_INITIAL/1e3; M2_f = MASS2_INITIAL/1e6;

  w = Sqrt(m2 + q*q);

  s = W*W;
  t = m2 - 2.0*k*(w - q*CosTheta);

  E_i = W - k;
  E_f = W - w;

  N = TwoPi()/(W*Sqrt((E_i+M_i)*(E_f+M_f)));

  return N*((s-M2_i)*maid_X(1,2,-1.0,CosTheta,e) + M_i*(t-m2)*maid_X(3,4,+1.0,CosTheta,e)/q);
}

//-----------------------------------------------------------------------------

inline TComplex maid_A2(Double_t CosTheta, Int_t e)
{
  Double_t E_i, E_f;
  Double_t m, m2, M_i, M2_i, M_f, M2_f;
  Double_t W, q, k, w;
  Double_t s, t;
  Double_t N;

  W = W_cm(gEnergy)/1e3;
  q = q_mes(gEnergy)/1e3;
  k = omega_cm(gEnergy)/1e3;

  m   = MASS_MESON/1e3;   m2   = MASS2_MESON/1e6;
  M_i = MASS_INITIAL/1e3; M2_i = MASS2_INITIAL/1e6;
  M_f = MASS_INITIAL/1e3; M2_f = MASS2_INITIAL/1e6;

  w = Sqrt(m2 + q*q);

  s = W*W;
  t = m2 - 2.0*k*(w - q*CosTheta);

  E_i = W - k;
  E_f = W - w;

  N = TwoPi()/(W*Sqrt((E_i+M_i)*(E_f+M_f)));

  return N*maid_Y(3,4,-1.0,CosTheta,e)/q;
}

//-----------------------------------------------------------------------------

inline TComplex maid_A3(Double_t CosTheta, Int_t e)
{
  Double_t E_i, E_f;
  Double_t m, m2, M_i, M2_i, M_f, M2_f;
  Double_t W, q, k, w;
  Double_t s, t;
  Double_t N;

  W = W_cm(gEnergy)/1e3;
  q = q_mes(gEnergy)/1e3;
  k = omega_cm(gEnergy)/1e3;

  m   = MASS_MESON/1e3;   m2   = MASS2_MESON/1e6;
  M_i = MASS_INITIAL/1e3; M2_i = MASS2_INITIAL/1e6;
  M_f = MASS_INITIAL/1e3; M2_f = MASS2_INITIAL/1e6;

  w = Sqrt(m2 + q*q);

  s = W*W;
  t = m2 - 2.0*k*(w - q*CosTheta);

  E_i = W - k;
  E_f = W - w;

  N = TwoPi()/(W*Sqrt((E_i+M_i)*(E_f+M_f)));

  return N*(maid_Y(1,2,+1.0,CosTheta,e) + (t-m2-4.0*k*W)*maid_X(3,4,+1.0,CosTheta,e)/(2.0*q));
}

//-----------------------------------------------------------------------------

inline TComplex maid_A4(Double_t CosTheta, Int_t e)
{
  Double_t E_i, E_f;
  Double_t m, m2, M_i, M2_i, M_f, M2_f;
  Double_t W, q, k, w;
  Double_t s, t;
  Double_t N;

  W = W_cm(gEnergy)/1e3;
  q = q_mes(gEnergy)/1e3;
  k = omega_cm(gEnergy)/1e3;

  m   = MASS_MESON/1e3;   m2   = MASS2_MESON/1e6;
  M_i = MASS_INITIAL/1e3; M2_i = MASS2_INITIAL/1e6;
  M_f = MASS_INITIAL/1e3; M2_f = MASS2_INITIAL/1e6;

  w = Sqrt(m2 + q*q);

  s = W*W;
  t = m2 - 2.0*k*(w - q*CosTheta);

  E_i = W - k;
  E_f = W - w;

  N = TwoPi()/(W*Sqrt((E_i+M_i)*(E_f+M_f)));

  return N*(maid_Y(1,2,+1.0,CosTheta,e) + (t-m2)*maid_X(3,4,+1.0,CosTheta,e)/(2.0*q));
}

//-----------------------------------------------------------------------------

//Wrapper for explicit invariant MAID amplitudes B_i
inline TComplex maid_A(Int_t i, Double_t CosTheta, Int_t e)
{
  switch(i)
  {
    case 1: return maid_A1(CosTheta, e);
    case 2: return maid_A2(CosTheta, e);
    case 3: return maid_A3(CosTheta, e);
    case 4: return maid_A4(CosTheta, e);
    default: return TComplex(0.0, 0.0);
  }
}

//-----------------------------------------------------------------------------

inline TComplex maid_A1cc(Double_t CosTheta, Int_t e){ return TComplex::Conjugate(maid_A1(CosTheta, e)); }
inline TComplex maid_A2cc(Double_t CosTheta, Int_t e){ return TComplex::Conjugate(maid_A2(CosTheta, e)); }
inline TComplex maid_A3cc(Double_t CosTheta, Int_t e){ return TComplex::Conjugate(maid_A3(CosTheta, e)); }
inline TComplex maid_A4cc(Double_t CosTheta, Int_t e){ return TComplex::Conjugate(maid_A4(CosTheta, e)); }

//-----------------------------------------------------------------------------

//Wrapper for explicit invariant MAID amplitudes A_i*
inline TComplex maid_Acc(Int_t i, Double_t CosTheta, Int_t e)
{
  switch(i)
  {
    case 1: return maid_A1cc(CosTheta, e);
    case 2: return maid_A2cc(CosTheta, e);
    case 3: return maid_A3cc(CosTheta, e);
    case 4: return maid_A4cc(CosTheta, e);
    default: return TComplex(0.0, 0.0);
  }
}

//-----------------------------------------------------------------------------

inline TComplex maid_B1(Double_t CosTheta, Int_t e)
{
  Double_t M_i = MASS_INITIAL/1e3;

  return maid_A1(CosTheta, e) - 2.0*M_i*maid_A4(CosTheta, e);
}

//-----------------------------------------------------------------------------

inline TComplex maid_B2(Double_t CosTheta, Int_t e)
{
  Double_t m2 = MASS2_MESON/1e6;
  Double_t q = q_mes(gEnergy)/1e3;
  Double_t k = omega_cm(gEnergy)/1e3;
  Double_t w = Sqrt(m2 + q*q);
  Double_t t = m2 - 2.0*k*(w - q*CosTheta);

  return k*(w-q*CosTheta)*maid_A2(CosTheta, e);
}

//-----------------------------------------------------------------------------

inline TComplex maid_B6(Double_t CosTheta, Int_t e)
{
  Double_t M_i = MASS_INITIAL/1e3;

  return -2.0*maid_A4(CosTheta, e);
}

//-----------------------------------------------------------------------------

inline TComplex maid_B8(Double_t CosTheta, Int_t e)
{
  Double_t M_i = MASS_INITIAL/1e3;

  return -maid_A3(CosTheta, e);
}

//-----------------------------------------------------------------------------

//Wrapper for explicit invariant MAID amplitudes B_i
inline TComplex maid_B(Int_t i, Double_t CosTheta, Int_t e)
{
  switch(i)
  {
    case 1: return maid_B1(CosTheta, e);
    case 2: return maid_B2(CosTheta, e);
    case 6: return maid_B6(CosTheta, e);
    case 8: return maid_B8(CosTheta, e);
    default: return TComplex(0.0, 0.0);
  }
}

//-----------------------------------------------------------------------------

inline TComplex maid_B1cc(Double_t CosTheta, Int_t e){ return TComplex::Conjugate(maid_B1(CosTheta, e)); }
inline TComplex maid_B2cc(Double_t CosTheta, Int_t e){ return TComplex::Conjugate(maid_B2(CosTheta, e)); }
inline TComplex maid_B6cc(Double_t CosTheta, Int_t e){ return TComplex::Conjugate(maid_B6(CosTheta, e)); }
inline TComplex maid_B8cc(Double_t CosTheta, Int_t e){ return TComplex::Conjugate(maid_B8(CosTheta, e)); }

//-----------------------------------------------------------------------------

//Wrapper for explicit invariant MAID amplitudes A_i*
inline TComplex maid_Bcc(Int_t i, Double_t CosTheta, Int_t e)
{
  switch(i)
  {
    case 1: return maid_B1cc(CosTheta, e);
    case 2: return maid_B2cc(CosTheta, e);
    case 6: return maid_B6cc(CosTheta, e);
    case 8: return maid_B8cc(CosTheta, e);
    default: return TComplex(0.0, 0.0);
  }
}

//-----------------------------------------------------------------------------

#endif
