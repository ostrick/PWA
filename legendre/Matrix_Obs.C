 
 
 const Int_t EBINS(300),THBINS(40);
  Float_t sgE_th2[THBINS], sgE_val2[THBINS], sgE_err2[THBINS];
  Int_t sgE_pts[EBINS];
  Int_t l_angular=0;
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


//-----------------------------------------------------------------------------
Double_t func(float x,Double_t* par)
{
	Double_t PI(3.1415926535);
	Int_t polynumber=l_angular;
    Double_t value=0;
    
    for(Int_t k=0; k<=polynumber;k++)
      {value = value + par[k]*L(k,cos(x*PI/180));}

	 
 
 return value;
}
 //Definition of fit function, the legendre polynomial
 
Double_t legendre(Double_t* var,Double_t* par)
{

	Double_t PI(3.1415926535);
    Double_t value=0;
    Int_t polynumber=l_angular;
    
    for(Int_t k=0; k<=polynumber; k++)
      {value = value + par[k]*L(k,cos(var[0]*PI/180));}
      
      return value;
                  
}



// reduced chi square
Double_t chisq(Float_t *x, Float_t *mess, Float_t *err, Double_t *para)

{   Int_t polynumber=l_angular;
	Int_t nbins = sgE_pts[dataset];
	Double_t chi=0, delta=0;
	for(Int_t i=0; i<nbins; i++)
	{  delta = (mess[i]-func(x[i],para))/err[i];
		chi += delta*delta;
	}
	  chi = chi/(nbins-(polynumber+1));
	return chi;
}





//main file for data import
//-----------------------------------------------------------------------------
void Matrix_Obs()

{ 
	
	cout<<"value for polynom order: "<<endl;
  cin>>l_angular;
  
	Double_t PI(3.1415926535);
  
 FILE* file_sgE;   
 
 Float_t Energy, Dummy1, Dummy2;
 Float_t Theta, sigmaE, err_sigmaE;
   Int_t ThetaBin, sgE_bin;



Float_t sgE_val[EBINS][THBINS];          //EBINS  THBINS
Float_t sgE_err[EBINS][THBINS];
Float_t sgE_th[EBINS][THBINS];


Float_t sgE_en[EBINS];
 
 file_sgE=fopen("sg0_CBMAMI_2014_Prakhov.txt","r");  //open data file
 
 
 
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
  
  const Int_t polynumber=l_angular, a = sgE_bin;
  
  Double_t parameter[polynumber+1], err_parameter[polynumber+1];
  
  Double_t para_arr[a][polynumber+1], err_para_arr[a][polynumber+1];
  Double_t chisquare[a];
  
  
     TF1* fitfunction = new TF1("fitfunction",legendre,0,180,polynumber+1);
     
 for(Int_t k=0; k<sgE_bin; k++){
  for(Int_t i=0; i<sgE_pts[k];i++){
	  sgE_th2[i]=sgE_th[k]i];
	  sgE_val2[i]=sgE_val[k][i];
	  sgE_err2[i]=sgE_err[k][i];}

      const Int_t bin = sgE_pts[k];
 
     TMatrixD theta(bin,1);
     TMatrixD value(bin,1);

	   for(Int_t j=0;j<sgE_pts[k]; j++) {
		   theta(j,0)=sgE_th2[j];
		   value(j,0)=sgE_val2[j]; }
		 //  value.Print();
	   TMatrixD A(bin,polynumber+1);
	  for(Int_t m=0;m<=polynumber;m++) {
	   for(Int_t j=0;j<bin;j++){
		   A[j][m]=L(m,cos(sgE_th2[j]*PI/180)); }          //L(m,cos(sgE_th2[j]*PI/180))
	   }		    
        TMatrixD W(bin,bin);
        for(Int_t j=0;j<sgE_pts[k];j++) {
			W(j,j)=1/sgE_err2[j]; }
			
			
			
        TMatrixD WA(W,TMatrixD::kMult,A);
        TMatrixD AWA(A,TMatrixD::kTransposeMult,WA);
      //  AWA.Print();
        
       // Double_t det = AWA.Determinant();
       // cout<<det<<endl;
        
        
          TDecompSVD svd(AWA);
          TMatrixD AWAinv = svd.Invert();
          
           TMatrixD WY(W,TMatrixD::kMult,value);
           TMatrixD AWY(A,TMatrixD::kTransposeMult,WY);
           TMatrixD para(AWAinv,TMatrixD::kMult,AWY);
        

            for(Int_t j=0; j<=polynumber; j++){
			   parameter[j]=para(j,0);
			   err_parameter[j]=sqrt(AWA(j,j)); }    
     
     
     for(Int_t j=0; j<=polynumber; j++){
		 para_arr[k][j]=parameter[j];
		 err_para_arr[k][j]=err_parameter[j];}
		 
		 
		 
		 chisquare[k]=chisq(sgE_th2, sgE_val2, sgE_err2, parameter);
     
     
     
     
     
     dataset++;
     cout<<k<<endl;}
     
    /* 
     const Int_t bin = sgE_pts[5];
     
     for(Int_t i=0;i<=polynumber;i++){
		 parameter[i]=para_arr[5][i];
		 cout<<parameter[i]<<endl;}
		 
		 for(Int_t i=0; i<sgE_pts[5];i++){
			 sgE_val2[i]=sgE_val[5]i];
			 sgE_err2[i]=sgE_err[5][i];
			 sgE_th2[i]=sgE_th[5][i]; }
			 
			 TGraphErrors *f1 = new TGraphErrors(sgE_pts[5],sgE_th2,sgE_val2,NULL,sgE_err2);
			 f1->Draw("AP");
			 
			 
			     TF1* fitfunction = new TF1("fitfunction",legendre,sgE_th2[0],sgE_th2[bin-1],polynumber+1);   
				for(Int_t i=0; i<=polynumber; i++)
				{fitfunction->SetParameter(i,parameter[i]);}
				fitfunction->Draw("SAME");   */
     
     
     
     if(polynumber<4){
     TCanvas *c = new TCanvas("c","c",(polynumber+1)*600,600);
     c->Divide(polynumber+1,1);}
     else{ if(polynumber<8){
		          TCanvas *c = new TCanvas("c","c",4*600,1200);
					c->Divide(4,2);}
					else{
						     TCanvas *c = new TCanvas("c","c",4*600,1800);
								c->Divide(4,3);}
							}
     

     for(Int_t i=0; i<polynumber+1;i++){
		 c->cd(i+1);
     Double_t para_fit[a],err_para[a], en[a];
     for(Int_t j=0; j<sgE_bin; j++){
		 para_fit[j]=para_arr[j][i];
		 err_para[j]=err_para_arr[j][i];
		 en[j]=sgE_en[j];}

		 
		 
		 TGraphErrors *f1 = new TGraphErrors(sgE_bin,en,para_fit,NULL,NULL);
		 f1->SetTitle(Form("A_{%d}",i));
		 f1->Draw();    }
		 
		 c->cd(polynumber+2);
     TGraphErrors *chi = new TGraphErrors(sgE_bin,en,chisquare,NULL,NULL);
     chi->SetTitle("#chi^{2}");
     chi->Draw();   
     
     
     c->Print(Form("Coefficients_l%d_sg0.pdf",l_angular));
     
     
     

}


  
Int_t ReadLine_sgE(FILE* file_sgE, Float_t* Theta, Float_t* sigmaE, Float_t* err_sigmaE)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), file_sgE);
  return sscanf(Buffer, "%f  %f    %f", Theta, sigmaE, err_sigmaE);
}


  



