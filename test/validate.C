
//make plots, compare numbers
void validate(char *fname) {
  
  TFile *fin=new TFile(fname);

  //first deal with the numbers:
  TH1F *theNumbers=(TH1F*)fin->Get("theNumbers_");
  int nEventsWithDuplicates=theNumbers->GetBinContent(1);
  int nEleMisMatchId=theNumbers->GetBinContent(2);
  int nTopElectrons=theNumbers->GetBinContent(3);
  int nTopElectronsPassEleId=theNumbers->GetBinContent(4);
  int nTopElectronsPassEleIdPassTqaf_tkIso=theNumbers->GetBinContent(5);
  int nTopElectronsPassEleIdPassTqaf_caloIso=theNumbers->GetBinContent(6);
  int nPMGsfElectrons=theNumbers->GetBinContent(7);
  int nPMGsfElectronsPassEleId=theNumbers->GetBinContent(8);


  cout <<"Number of events with duplicate TopElectrons (should be zero): "<<nEventsWithDuplicates<<endl;
  cout <<"Number of TopElectrons that can't be matched to PMGsfElectrons with same electronId result (should be zero): "
       <<nEleMisMatchId<<endl;

  cout <<"Fraction of electrons passing electron id cuts: "<<((float) nTopElectronsPassEleId)/nTopElectrons<<endl;
  cout <<"Fraction of electrons that, having passed electron id cuts, pass tqaf tracker iso < 3 GeV cut: "<<
    ((float) nTopElectronsPassEleIdPassTqaf_tkIso)/nTopElectronsPassEleId<<endl;
  cout <<"Fraction of electrons that, having passed electron id cuts, pass tqaf calo iso < 6 GeV cut: "<<
    ((float) nTopElectronsPassEleIdPassTqaf_caloIso)/nTopElectronsPassEleId<<endl;
  cout <<"total number of PMGsfElectrons (after removing duplicates): "<<nPMGsfElectrons<<endl;
  cout <<"total number of TopElectrons: "<<nTopElectrons<<endl;
  cout <<"number of PMGsfElectrons passing electronId (after removing duplicates): "<<nPMGsfElectronsPassEleId<<endl;
  cout <<"number of TopElectrons passing electronId: "<<nTopElectronsPassEleId<<endl;
  


  //then the histograms:
  TH1F *allGenEleVsEt=(TH1F*)fin->Get("allGenEleVsEt_");
  TH1F *allGenEleVsEta=(TH1F*)fin->Get("allGenEleVsEta_");
  TH1F *matchedGenEleVsEt=(TH1F*)fin->Get("matchedGenEleVsEt_");
  TH1F *matchedGenEleVsEta=(TH1F*)fin->Get("matchedGenEleVsEta_");
  
  TH1F *effVsEt=matchedGenEleVsEt->Clone();
  effVsEt->Divide(allGenEleVsEt);
  TH1F *effVsEta=matchedGenEleVsEta->Clone();
  effVsEta->Divide(allGenEleVsEta);

  TCanvas *c0=new TCanvas();
  effVsEt->Draw();
  
  TCanvas *c1=new TCanvas();
  effVsEta->Draw();
  
}
