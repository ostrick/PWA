 
 
 const Int_t EBINS(300),THBINS(40);
  Float_t sgE_th2[THBINS], sgE_val2[THBINS], sgE_err2[THBINS];
  Int_t sgE_pts[EBINS];
  Int_t l=0;
  Int_t dataset=0;  //gives the right bin number for the fcn function



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
//----------------------------------------------------------------------------

 //Definition of fit function, the legendre polynomial
 
Double_t legendre(Double_t* var,Double_t* par)
{

	Double_t PI(3.1415926535);
    Double_t value=0;
    Int_t polynumber=2*l;
    
    for(Int_t k=0; k<=polynumber; k++)
      {value = value + par[k]*L(k,cos(var[0]*PI/180));}
      
      return value;
                  
}



//-----------------------------------------------------------------------------
Double_t func(float x,Double_t* par)
{
	Double_t PI(3.1415926535);
	Int_t polynumber=2*l;
    Double_t value=0;
    
    for(Int_t k=0; k<=polynumber;k++)
      {value = value + par[k]*L(k,cos(x*PI/180));}

	 
 
 return value;
}

//chi square reduction 
//-----------------------------------------------------------------------------

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

   Int_t nbins =sgE_pts[dataset]; //gE_pts[0];  //amount of datapoints

   Double_t chisq = 0;
   Double_t delta;
   for (Int_t i=0;i<nbins; i++) {
     delta  = (sgE_val2[i]-
     func(sgE_th2[i],par))/sgE_err2[i];
     chisq += delta*delta;
   }
   f = chisq;
}

// reduced chi square
Double_t chisq(Float_t *x, Float_t *mess, Float_t *err, Double_t *para)

{   Int_t polynumber=2*l;
	Int_t nbins = sgE_pts[dataset];
	Double_t chi=0, delta=0;
	for(Int_t i=0; i<nbins; i++)
	{  delta = (mess[i]-func(x[i],para))/err[i];
		chi += delta*delta;
	}
	 // chi = chi/(nbins-(polynumber+1));
	return chi;
}


//fit definition with minuit
//-------------------------------------------------------------------------
void minuit(Float_t* sgE_en, Float_t* sgE_th2, Float_t* sgE_val2, Float_t* sgE_err2, Int_t* sgE_pts,Int_t number,Double_t* parameter, Double_t* err_parameter, Int_t l)

{
	
	const Int_t bin=number;
	
	/*TCanvas *c = new TCanvas("c","c",1200,1000);
	TGraphErrors *graph = new TGraphErrors(sgE_pts[bin],sgE_th2,sgE_val2,NULL,sgE_err2);
	graph->SetTitle(Form("Energy: %.0f MeV",sgE_en[number]));
	graph->GetXaxis()->SetTitle("Angle [#circ]");
	graph->GetYaxis()->SetTitle("sgE");
	graph->Draw("AP");*/
	
	
	const Int_t polynumber=2*l;
	 
	 

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
   
   arglist[0]=500;  //500 iterations of migrad
   gMinuit->mnexcm("MIGRAD",arglist,1,iflag);    //mimization process


     //Plot the fit function
//----------------------------------------------------------------------


   
   for (Int_t i=0; i<=polynumber; i++)
   {
	   gMinuit->GetParameter(i,parameter[i],err_parameter[i]);
   }
   
    
    
   // TF1* fitfunction = new TF1("fitfunction",legendre,0,180,polynumber+1);   
	//for(Int_t i=0; i<=polynumber; i++)
	//{fitfunction->SetParameter(i,outpar[i]);}

							   
    
	//fitfunction->Draw("SAME");
	//c->Print(Form("sgE_%.0f.png",sgE_en[number]));
}


//main file for data import
//-----------------------------------------------------------------------------
void obs_read_in()

{  
	Int_t l_low=0, l_hi=0;
	
	cout<<"value for l: "<<endl;
  cin>>l;
  
//	cout<<"max value for l: "<<endl;
//  cin>>l_hi;

  
 FILE* file_sgE;   
 
 Float_t Energy, Dummy1, Dummy2;
 Float_t Theta, sigmaE, err_sigmaE;
   Int_t ThetaBin, sgE_bin;



Float_t sgE_val[EBINS][THBINS];          //EBINS  THBINS
Float_t sgE_err[EBINS][THBINS];
Float_t sgE_th[EBINS][THBINS];


Float_t sgE_en[EBINS];
 
 file_sgE=fopen("sgE_DAPHNE_2001_Preobrajenski.txt","r");  //open data file
 
 
 
 while(!feof(file_sgE)) //read in energy
 
   {    
	   if(fscanf(file_sgE,"E =  %f MeV, E_lo =  %f MeV, E_hi =  %f MeV\n",&Energy, &Dummy1, &Dummy2)!=3) break;

    
    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
      ReadLine_sgE(file_sgE, &Theta, &sigmaE, &err_sigmaE);
      while(ReadLine_sgE(file_sgE, &Theta, &sigmaE, &err_sigmaE)>=2)
    {
      sgE_val[sgE_bin][ThetaBin] = sigmaE;
      sgE_err[sgE_bin][ThetaBin] = err_sigmaE;
      sgE_th[sgE_bin][ThetaBin]  = Theta;
      if(err_sigmaE!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
  }

    //Store data for this energy bin
    sgE_pts[sgE_bin] = ThetaBin;
    sgE_en[sgE_bin] = Energy;
   // strcpy(sgE_id[sgE_bin], Ident);

    //Increase energy bin counter
    sgE_bin++;
  }

  fclose(file_sgE);
    //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgE_bin; t++) n+=sgE_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgE_bin; t++) if(sgE_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  
 /* //Debug output
  printf("EBins: %d\n", sgE_bin);
  for(Int_t e=0; e<sgE_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgE_en[e], sgE_pts[e]);
    for(Int_t th=0; th<sgE_pts[e]; th++)
      printf("%f %f %f\n", sgE_th[e][th], sgE_val[e][th], sgE_err[e][th]);
  }*/
  
 // for(Int_t l_val=l_low; l_val<=l_hi; l_val++) {  
  const Int_t polynumber=2*l, a = sgE_bin;
  Double_t parameter[polynumber+1], err_parameter[polynumber+1];
  Double_t para_arr[a][polynumber+1], err_para_arr[a][polynumber+1];
  Double_t chisquare[a];
  

 for(Int_t k=0; k<sgE_bin; k++){
  for(Int_t i=0; i<sgE_pts[k];i++){
	  sgE_th2[i]=sgE_th[k]i];
	  sgE_val2[i]=sgE_val[k][i];
	  sgE_err2[i]=sgE_err[k][i];}

     minuit(sgE_en, sgE_th2, sgE_val2, sgE_err2, sgE_pts, k, parameter, err_parameter,l);
     for(Int_t j=0; j<=polynumber; j++){
		 para_arr[k][j]=parameter[j];
		 err_para_arr[k][j]=err_parameter[j];}
		 
		 chisquare[k]=chisq(sgE_th2, sgE_val2, sgE_err2, parameter);
     
     dataset++;
     cout<<k<<endl;}
     
 
     
     TCanvas *c = new TCanvas("c","c",(polynumber+2)*600,600);
     c->Divide(polynumber+1,1);
     

     for(Int_t i=0; i<polynumber+1;i++){
		 c->cd(i+1);
     Double_t para[a],err_para[a], en[a];
     for(Int_t j=0; j<sgE_bin; j++){
		 para[j]=para_arr[j][i];
		 err_para[j]=err_para_arr[j][i];
		 en[j]=sgE_en[j];
		 cout<<para[j]<<"	"<<err_para[j]<<"	"<<sgE_en[j]<<endl;}
		 cout<<sgE_bin<<endl;
		 
		 
		 TGraphErrors *f1 = new TGraphErrors(sgE_bin,en,para,NULL,err_para);
		 f1->SetTitle(Form("A_{%d}",i));
		 f1->Draw();    }
		 
		/* c->cd(1);
     TGraphErrors *chi = new TGraphErrors(sgE_bin,en,chisquare,NULL,NULL);
     chi->SetTitle("#chi^{2}");
     chi->Draw();*/
     
     c->Print(Form("Coefficients_l%d.png",l));
     

    
}

/*void obs_read_in()
{
 Int_t l_low=0, l_hi=0;
 
/* cout<<"Min value for l: "<<endl;
 cin>>l_low;
 
 cout<<"Max value for l: "<<endl;
 cin>>l_hi;
 
 for(Int_t i=l_low; i<=l_hi; i++){
	 main_func(i); }
	 main_func;
	 
}*/
  
Int_t ReadLine_sgE(FILE* file_sgE, Float_t* Theta, Float_t* sigmaE, Float_t* err_sigmaE)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), file_sgE);
  return sscanf(Buffer, "%f  %f    %f", Theta, sigmaE, err_sigmaE);
}


  



