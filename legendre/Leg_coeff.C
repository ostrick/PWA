#include "TMinuit.h"

#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"

#include <Riostream.h>
#include "TLegend.h"
#include "TLegendEntry.h"

#include "Math/IFunction.h"
#include <cmath>
#include "TSystem.h"


//program to calculate the legendre coefficients for a specific order
  const char *GetID(Int_t type)  //changes Int to String
{
	return Form("%d", type);
}

//Definition of legendre polynoms
//----------------------------------------------------------------
inline Double_t L(Int_t l, Double_t x)
{
  if(l==0) return 1.0; //Legendre polynom P_0
  if(l==1) return x;   //Legendre polynom P_1
  if(l > 1) return (1.0/l) *((2.0*l - 1.0)*x*L(l-1, x) - (l - 1.0)*L(l-2, x));
  return 0.0; //Any other case, e.g. for l < 0
}

//-----------------------------------------------------------------------------

inline Double_t DL(Int_t l, Double_t x)
{
  if(l==0) return 0.0; //Legendre polynom 1st derivative P'_0
  if(l==1) return 1.0; //Legendre polynom 1st derivative P'_1
  if(l > 1) return (1.0/(l - 1.0)) * ((2.0*l - 1.0)*x*DL(l-1, x) - l*DL(l-2, x));
  return 0.0; //Any other case, e.g. for l < 0
}

//-----------------------------------------------------------------------------

inline Double_t D2L(Int_t l, Double_t x)
{
  if(l==0) return 0.0; //Legendre polynom 2nd derivative P"_0
  if(l==1) return 0.0; //Legendre polynom 2nd derivative P"_1
  if(l==2) return 3.0; //Legendre polynom 2nd derivative P"_2
  if(l > 2) return (1.0/(l - 2.0)) * ((2.0*l - 1.0)*x*D2L(l-1, x) - (l + 1.0)*D2L(l-2, x));
  return 0.0; //Any other case, e.g. for l < 0
}

//-----------------------------------------------------------------------------





  Double_t sigma1[30]={}, cross1[30]={}, err1[30]={};  //arrays for data
  const Int_t polynumber=0;  //order of legendre polynom
//----------------------------------------------------------------

 //Definition of fit function, the legendre polynomial
 
Double_t legendre(Double_t* var,Double_t* par)
{

	Double_t PI(3.1415926535);
    Double_t value=0;
    
    for(Int_t k=0; k<=polynumber; k++)
      {value = value + par[k]*L(k,cos(var[0]*PI/180));}
      
      return value;
                  
}

Double_t func(float x,Double_t* par)
{
	Double_t PI(3.1415926535);
	 
    Double_t value=0;
    
    for(Int_t k=0; k<=polynumber;k++)
      {value = value + par[k]*L(k,cos(x*PI/180));}

	 
 
 return value;
}


 //Calculation of chi square
//--------------------------------------------------------------------------

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   Double_t a=0;
   const Int_t nbins = 30;  //amount of datapoints per data set

   

   Double_t chisq = 0;
   Double_t delta;
   for (Int_t i=0;i<nbins; i++) {
     delta  = (cross1[i]-
     func(sigma1[i],par))/err1[i];
     chisq += delta*delta;
   }
   f = chisq;
}

// reduced chi square
Double_t chisquare(Double_t *x, Double_t *mess, Double_t *err, Double_t *outpar)

{ Double_t chi=0, delta=0;
	for(Int_t i=0; i<30; i++)
	{  delta = (mess[i]-func(x[i],outpar))/err[i];
		chi += delta*delta;
	}
	  chi = chi/(30-(polynumber+1));
	return chi;
}

  //main function
//-------------------------------------------------------
void Leg_coeff()
{
  cout<<"Order of Legendre polynom"<<endl;
  cin>>polynumber;
  

	
 FILE* file_prakhov;   
 
 float Energy, Dummy1, Dummy2, Dummy3;
 float sigma, crossc, err;
 float ENERGY[246]={}, SIGMA[7380]={}, CROSS[7380]={}, ERROR[7380]={};
 char Buffer[1024];

 
 int eng(0), zeilen(0);  //variables to count the amount of data (energy and cross section)
 
 file_prakhov=fopen("sg0_CBMAMI_2014_Prakhov.txt","r");  //open data file
 
 int evalue=0;   
 int value=0;    
 
 
 while(!feof(file_prakhov)) //read in energy
 
   {    
	   if(fscanf(file_prakhov,"E =  %f MeV, E_lo =  %f MeV, E_hi =  %f MeV\n",&Energy, &Dummy1, &Dummy2)==3)
	   { 
		   ENERGY[evalue]=Energy;  
		   eng++;
		   evalue++;
	   }
	   
	    //skip one uninteresting line
	   fgets(Buffer, sizeof(Buffer), file_prakhov);  
	 
	   //read in cross section, angle and error of cross section
	   
	  while(fscanf(file_prakhov,"%f     %f     %f     %f\n", &sigma, &crossc, &err, &Dummy3)==4)
	  {
		  SIGMA[value]=sigma;
		  CROSS[value]=crossc;    
		  ERROR[value]=err;
		  value++;
		  zeilen++;
	  }
	    //Skip one line again
	  fgets(Buffer, sizeof(Buffer), file_prakhov);
	  
  }
   //amount of energy values and lines of data
  cout<<"Anzahl Energiewerte: "<<eng<<" Zeilenzahl: "<<zeilen<<endl; 

  fclose(file_prakhov);
  
  const Int_t Anz=246;   //amount of evaluated data sets
  Double_t coeff[(polynumber+1)*Anz]={}, coefferr[(polynumber+1)*Anz]={};
  Double_t cs[Anz]={};

for (Int_t i=1; i<=Anz; i++)

{
  //plot of first measurement (lowest energy)
//----------------------------------------------------------------------
	 
	 for (Int_t k=0;k<30;k++)
	 {
	  sigma1[k]= SIGMA[30*(i-1)+k];
	  cross1[k]=CROSS[30*(i-1)+k];
	  err1[k]=ERROR[30*(i-1)+k];
     }
  
  



 //definition of the fit with Minuit
//---------------------------------------------------------------------- 
  TMinuit *gMinuit = new TMinuit(polynumber+1);
  gMinuit->SetFCN(fcn);  //opens the funtion fcn for chi square calculation
  
     Double_t arglist[10];
     arglist[0]=1;
    
     Int_t iflag=1;
   
   
   
   gMinuit->SetPrintLevel(-1);  //-1 no output, 1 output
   
   gMinuit->mnexcm("SET ERR", arglist,1,iflag);  //chi square 1, neg log likelihood 0.5
   
   // Set starting values and step sizes for parameters
   //the last two numbers define the limits, 0 0 -> no limits
   
    TString str = "a";
    for(Int_t I=0; I<=polynumber;I++)
    {   str=str+GetID(I);
		gMinuit->mnparm(I, str, 1, 1, 0,0,iflag);
		str="a";}
   
   gMinuit->mnexcm("CALL FCN", arglist,1,iflag);  //chi square calculation
   
   arglist[0]=1;  //2 better error estimation, 1 faster
   gMinuit->mnexcm("SET STR",arglist,1,iflag);   //minimization strategy
   
   arglist[0]=750;  //500 iterations of migrad
   gMinuit->mnexcm("MIGRAD",arglist,1,iflag);    //mimization process

   //gMinuit->mnexcm("HESSE",arglist,1,iflag);     //error matrix 
   
     //Plot the fit function
   Double_t outpar[polynumber+1], error[polynumber+1];
   for (Int_t j=0; j<=polynumber; j++)
   {
	   gMinuit->GetParameter(j,outpar[j],error[j]);

	   
	   coeff[(polynumber+1)*(i-1)+j]=outpar[j];
		coefferr[(polynumber+1)*(i-1)+j]=error[j];
   }
  

  cout<<i<<endl;
   
   
   
   cs[i-1]=chisquare(sigma1,cross1,err1,outpar);
   
}


//output of coefficients	 

     TCanvas* c = new TCanvas("c","L-Coeff",1900,1200);
     c->Divide(3,4);
     Double_t x[Anz]={};
     
		for(Int_t i=0;i<Anz;i++)
		{ x[i]=ENERGY[i];}
		
     for(Int_t j=0;j<=polynumber;j++)
     {		
		 Double_t y[Anz]={}, erry[Anz]={};
		 for(Int_t n=0; n<Anz;n++)
		 {
			 y[n]=coeff[(polynumber+1)*n+j];
			 erry[n]=coeff[(polynumber+1)*n+j];
			 

		 }
		 	 TGraphErrors *f1 = new TGraphErrors(Anz,x,y,NULL,NULL);

		 	 f1->SetTitle(Form("A%d",j));       //"A"Form("%d",i)
			 f1->SetMarkerStyle(6);
			 f1->GetXaxis()->SetTitle("Energy [MeV]");
			 f1->GetYaxis()->SetTitle("Leg Coeff");
			 c->cd(j+1);
			 f1->Draw();
	 }
			 c->cd(polynumber+2);
		  
             TGraph *chi = new TGraph(Anz,x,cs);
             chi->SetTitle("Red. chi square");
			 chi->GetXaxis()->SetTitle("Energy [MeV]");
			 chi->GetYaxis()->SetTitle("chi square");             
             chi->SetMarkerStyle(6);
             chi->SetLineColor(2);
             chi->Draw();
             
             for(Int_t t=0; t<Anz; t++)
             {  Double_t chisq[Anz]={};
				   chisq[t]=cs[t];
				   cout<<chisq[t]<<endl;
				 }	
     
     c->Print("legendre_poly.pdf");
}
