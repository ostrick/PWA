void ExtractPlot()
{
Extract("peta",5);

Multipole("E0p",0,0,true,true);
 c1->Close();
Multipole("E1p",0,0,true,true);
 c1->Close();
Multipole("M1p",0,0,true,true);
 c1->Close();
Multipole("M1m",0,0,true,true);
 c1->Close();

Multipole("E2p",0,0,true,true);
 c1->Close();
Multipole("E2m",0,0,true,true);
 c1->Close();
Multipole("M2p",0,0,true,true);
 c1->Close();
Multipole("M2m",0,0,true,true);
 c1->Close();


Multipole("E3p",0,0,true,true);
 c1->Close();
Multipole("E3m",0,0,true,true);
 c1->Close();
Multipole("M3p",0,0,true,true);
 c1->Close();
Multipole("M3m",0,0,true,true);
 c1->Close();

Chi2();
 c1->Close();

}
