

void validate_electrons(char *fname) {
  
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
  TH1F *recoMatchedGenEleVsEt=(TH1F*)fin->Get("recoMatchedGenEleVsEt_");
  TH1F *recoMatchedGenEleVsEta=(TH1F*)fin->Get("recoMatchedGenEleVsEta_");
  TH1F *selMatchedGenEleVsEt=(TH1F*)fin->Get("selMatchedGenEleVsEt_");
  TH1F *selMatchedGenEleVsEta=(TH1F*)fin->Get("selMatchedGenEleVsEta_");
  TH1F *isoMatchedGenEleVsEt=(TH1F*)fin->Get("isoMatchedGenEleVsEt_");
  TH1F *isoMatchedGenEleVsEta=(TH1F*)fin->Get("isoMatchedGenEleVsEta_");

  TH1F *deltaRReco=(TH1F*)fin->Get("deltaRReco_");
  TH1F *deltaRSel=(TH1F*)fin->Get("deltaRSel_");
  TH1F *deltaRIso=(TH1F*)fin->Get("deltaRIso_");


  
  TH1F *recoEffVsEt=recoMatchedGenEleVsEt->Clone();
  recoEffVsEt->Divide(allGenEleVsEt);
  TH1F *recoEffVsEta=recoMatchedGenEleVsEta->Clone();
  recoEffVsEta->Divide(allGenEleVsEta);
  
  TH1F *selEffVsEt=selMatchedGenEleVsEt->Clone();
  selEffVsEt->Divide(allGenEleVsEt);
  TH1F *selEffVsEta=selMatchedGenEleVsEta->Clone();
  selEffVsEta->Divide(allGenEleVsEta);

  TH1F *isoEffVsEt=isoMatchedGenEleVsEt->Clone();
  isoEffVsEt->Divide(allGenEleVsEt);
  TH1F *isoEffVsEta=isoMatchedGenEleVsEta->Clone();
  isoEffVsEta->Divide(allGenEleVsEta);


  TCanvas *c0=new TCanvas();
  recoEffVsEt->SetXTitle("gen. level ET");
  recoEffVsEt->SetYTitle("eff.");
  recoEffVsEt->SetLineWidth(3);
  selEffVsEt->SetLineWidth(3);
  isoEffVsEt->SetLineWidth(3);
  recoEffVsEt->Draw();
  selEffVsEt->SetLineColor(4);
  selEffVsEt->Draw("same");
  isoEffVsEt->SetLineColor(2);
  isoEffVsEt->Draw("same");
  TLegend *leg0=new TLegend(0.4,0.6,0.9,0.9);
  leg0->AddEntry(recoEffVsEt,"reco.");
  leg0->AddEntry(selEffVsEt,"selected, \"loose\"");
  leg0->AddEntry(isoEffVsEt,"selected, \"loose\" and isolated");
  leg0->Draw("same");

  TCanvas *c1=new TCanvas();
  recoEffVsEta->SetXTitle("gen. level eta");
  recoEffVsEta->SetYTitle("eff.");
  recoEffVsEta->SetLineWidth(3);
  selEffVsEta->SetLineWidth(3);
  isoEffVsEta->SetLineWidth(3);
  recoEffVsEta->Draw();
  selEffVsEta->SetLineColor(4);
  selEffVsEta->Draw("same");
  isoEffVsEta->SetLineColor(2);
  isoEffVsEta->Draw("same");
  TLegend *leg1=new TLegend(0.4,0.6,0.9,0.9);
  leg1->AddEntry(recoEffVsEta,"reco.");
  leg1->AddEntry(selEffVsEta,"selected, \"loose\"");
  leg1->AddEntry(isoEffVsEta,"selected, \"loose\" and isolated");
  leg1->Draw("same");
  
  TCanvas *c2=new TCanvas();
  deltaRReco_->SetLineWidth(3);
  deltaRReco_->DrawNormalized();
  deltaRSel_->SetLineWidth(3);
  deltaRSel_->SetLineColor(4);
  deltaRIso_->SetLineWidth(3);
  deltaRIso_->SetLineColor(2);
  deltaRSel->DrawNormalized("same");
  deltaRIso->DrawNormalized("same");
  TLegend *leg2=new TLegend(0.4,0.6,0.9,0.9);
  leg2->AddEntry(deltaRReco,"reco.");
  leg2->AddEntry( deltaRSel,"selected, \"loose\"");
  leg2->AddEntry(deltaRIso,"selected, \"loose\" and isolated");
  leg2->Draw("same");
  
  
}



void validate_muons(char *fname) {
  
  TFile *fin=new TFile(fname);

  TH1F *allGenMuonVsEt=(TH1F*)fin->Get("allGenMuonVsEt_");
  TH1F *allGenMuonVsEta=(TH1F*)fin->Get("allGenMuonVsEta_");
  TH1F *recoMatchedGenMuonVsEt=(TH1F*)fin->Get("recoMatchedGenMuonVsEt_");
  TH1F *recoMatchedGenMuonVsEta=(TH1F*)fin->Get("recoMatchedGenMuonVsEta_");
  TH1F *isoMatchedGenMuonVsEt=(TH1F*)fin->Get("isoMatchedGenMuonVsEt_");
  TH1F *isoMatchedGenMuonVsEta=(TH1F*)fin->Get("isoMatchedGenMuonVsEta_");

  TH1F *deltaRRecoMuon=(TH1F*)fin->Get("deltaRRecoMuon_");
  TH1F *deltaRIsoMuon=(TH1F*)fin->Get("deltaRIsoMuon_");

  TH1F *recoMuonEffVsEt=recoMatchedGenMuonVsEt->Clone();
  recoMuonEffVsEt->Divide(allGenMuonVsEt);
  TH1F *recoMuonEffVsEta=recoMatchedGenMuonVsEta->Clone();
  recoMuonEffVsEta->Divide(allGenMuonVsEta);
  
  TH1F *isoMuonEffVsEt=isoMatchedGenMuonVsEt->Clone();
  isoMuonEffVsEt->Divide(allGenMuonVsEt);
  TH1F *isoMuonEffVsEta=isoMatchedGenMuonVsEta->Clone();
  isoMuonEffVsEta->Divide(allGenMuonVsEta);
  
  TCanvas *c0=new TCanvas();
  recoMuonEffVsEt->SetXTitle("gen. level ET");
  recoMuonEffVsEt->SetYTitle("eff.");
  recoMuonEffVsEt->SetLineWidth(3);
  isoMuonEffVsEt->SetLineWidth(3);
  recoMuonEffVsEt->Draw();
  isoMuonEffVsEt->SetLineColor(2);
  isoMuonEffVsEt->Draw("same");
  TLegend *leg0=new TLegend(0.4,0.6,0.9,0.9);
  leg0->AddEntry(recoMuonEffVsEt,"reco.");
  leg0->AddEntry(isoMuonEffVsEt,"isolated");
  leg0->Draw("same");
  
}



void validate_met(char *fname) {

  TFile *fin=new TFile(fname);

  TH1F *genMET=(TH1F*)fin->Get("genMET_");
  TH1F *tqafGenMET=(TH1F*)fin->Get("tqafGenMET_");
  TH1F *corrCaloMET=(TH1F*)fin->Get("corrCaloMET_");
  TH1F *topMET=(TH1F*)fin->Get("topMET_");
  
  genMET->SetLineWidth(3);
  tqafGenMET->SetLineWidth(3);
  topMET->SetLineWidth(3);
  corrCaloMET->SetLineWidth(3);


  TCanvas *c0=new TCanvas();
  genMET->Draw();
  tqafGenMET->SetLineColor(4);
  tqafGenMET->Draw("same");
  topMET->SetLineColor(2);
  topMET->Draw("same");
  corrCaloMET->SetLineColor(6);
  corrCaloMET->Draw("same");
  TLegend *leg0=new TLegend(0.4,0.6,0.9,0.9);
  leg0->AddEntry(genMET,"reco::GenMET");
  leg0->AddEntry(tqafGenMET,"TQAF Gen-level MET");
  leg0->AddEntry(topMET,"TopMET");
  leg0->AddEntry(corrCaloMET,"Type1-Corrected reco::CaloMET");
  leg0->Draw("same");


}
