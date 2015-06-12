void ExtractPlot(Int_t SOLUTION=0)
{
Int_t  nsol = SOLUTION+1;

  Extract("peta",5,nsol);

  Multipole("E0p",SOLUTION,SOLUTION,true,true);
 c1->Close();
Multipole("E1p",SOLUTION,SOLUTION,true,true);
 c1->Close();
Multipole("M1p",SOLUTION,SOLUTION,true,true);
 c1->Close();
Multipole("M1m",SOLUTION,SOLUTION,true,true);
 c1->Close();

Multipole("E2p",SOLUTION,SOLUTION,true,true);
 c1->Close();
Multipole("E2m",SOLUTION,SOLUTION,true,true);
 c1->Close();
Multipole("M2p",SOLUTION,SOLUTION,true,true);
 c1->Close();
Multipole("M2m",SOLUTION,SOLUTION,true,true);
 c1->Close();


Multipole("E3p",SOLUTION,SOLUTION,true,true);
 c1->Close();
Multipole("E3m",SOLUTION,SOLUTION,true,true);
 c1->Close();
Multipole("M3p",SOLUTION,SOLUTION,true,true);
 c1->Close();
Multipole("M3m",SOLUTION,SOLUTION,true,true);
 c1->Close();

Chi2(SOLUTION);
Penalty(SOLUTION);
 c1->Close();

}
