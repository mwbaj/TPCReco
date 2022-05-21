/*
Usage:

$ root -l
root [0] .L makePlots_report.cpp
root [1] makePlots_report("Histos.root", 11.5);

*/

#ifndef __ROOTLOGON__
R__ADD_INCLUDE_PATH(../../DataFormats/include)
R__ADD_LIBRARY_PATH(../lib)
#endif

#include <tuple>
#include "GeometryTPC.h"

void makePlots_report(std::string fileNameHistos, float energyMeV, std::string cutInfo=""){

  if (!gROOT->GetClass("GeometryTPC")){
    R__LOAD_LIBRARY(libDataFormats.so);
  }
  GeometryTPC *geo=new GeometryTPC("geometry_ELITPC.dat", false);
  if(!geo) exit(-1);
  double xmin, xmax, ymin, ymax;
  std::tie(xmin, xmax, ymin, ymax)=geo->rangeXY();
  TGraph *area1=new TGraph(geo->GetActiveAreaConvexHull(0));
  TGraph *area2=new TGraph(geo->GetActiveAreaConvexHull(5.0));
  area1->GetXaxis()->SetLimits(xmin-0.025*(xmax-xmin), xmax+0.025*(xmax-xmin));
  area1->GetYaxis()->SetLimits(ymin-0.100*(ymax-ymin), ymax+0.100*(ymax-ymin));
  area1->SetLineColor(kBlack);
  area1->SetLineWidth(2);
  area1->SetLineStyle(kSolid);
  area2->SetLineColor(kRed);
  area2->SetLineWidth(2);
  area2->SetLineStyle(kSolid);

  TFile *aFile = new TFile(fileNameHistos.c_str(), "OLD");
  if(!aFile) exit(-1);

  gStyle->SetOptFit(111); // PCEV : prob=disabled / chi2,errors,variables=enabled
  gStyle->SetOptStat(10); // KISOURMEN : entries=enabled, name=disabled
  gStyle->SetPadLeftMargin(0.125);
  gStyle->SetPadBottomMargin(0.125);
  gStyle->SetTitleOffset(1.6, "X");
  gStyle->SetTitleOffset(1.7, "Y");
  
  TCanvas *c=new TCanvas("c","c",800,600);

  std::string prefix=Form("figures_%.1fMeV", energyMeV);
  c->Print(((string)(prefix)+".pdf[").c_str());
  gPad->SetLogy(false);
  auto *f=TFile::Open("Histos.root");
  f->cd();

  TH1D *h1;
  TH2D *h2;
  TProfile *hp;
  TPaveStats *st;
  auto statX1=0.75;
  auto statY1=0.825;
  auto statX2=statX1-0.1;
  auto statY2=statY1-0.05;
  auto statX3=statX1-0.04;
  auto statY3=statY1;
  auto statX4=0.65;
  auto statY4=0.25;
  auto statX5=0.625;
  auto statY5=0.20;
  auto statX6=statX1-0.1;
  auto statY6=statY1;
  auto lineH=0.075;
  auto lineW=0.15;
  std::string energyInfo=Form("E_{#gamma}=%.1f MeV : ", energyMeV);
  
  ///// control plot: number of tracks
  h1=(TH1D*)f->Get("h_ntracks");
  h1->SetTitle((energyInfo+h1->GetTitle()).c_str());
  gStyle->SetOptStat(10); // KISOURMEN : show entries only
  h1->SetStats(true);
  h1->Draw("HIST");
  h1->Draw("TEXT00 SAME");
  h1->UseCurrentStyle();
  h1->SetLineWidth(2);
  h1->GetXaxis()->SetNdivisions(5);
  c->Update();
  st = (TPaveStats *)c->GetPrimitive("stats");
  if(st) {
    st->SetX1NDC(statX1); st->SetX2NDC(statX1+lineW); st->SetY1NDC(statY1); st->SetY2NDC(statY1+lineH);
  }
  c->Update();
  c->Modified();
  c->Print(((string)(prefix)+".pdf").c_str());

  //// ALL categories, any track endpoint X vs Y
  c->Clear();
  h2=(TH2D*)f->Get("h_all_endXY")->Clone(); // copy for modifications of the same histogram
  h2->SetTitle((energyInfo+"Fully contained tracks").c_str());
  h2->SetXTitle("Track endpoint X_{DET} [mm]");
  h2->SetYTitle("Track endpoint Y_{DET} [mm]");
  gStyle->SetOptStat(10); // KISOURMEN : show entries only
  h2->Rebin2D(5,5);
  //  h2->Rebin2D(10,10);
  h2->UseCurrentStyle();
  area1->GetXaxis()->SetTitle(h2->GetXaxis()->GetTitle());
  area1->GetYaxis()->SetTitle(h2->GetYaxis()->GetTitle());
  area1->SetTitle(h2->GetTitle());
  area1->Draw("AL"); // UVW perimeter
  h2->Draw("COLZ SAMES"); // draw stat box
  area2->Draw("L SAME"); // UVW perimeter with 5mm veto band
  h2->SetStats(true);
  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.14);
  area1->GetXaxis()->SetTitleOffset(1.6);
  area1->GetYaxis()->SetTitleOffset(1.4);
  h2->SetTitleOffset(1.2, "Z");
  gPad->Update();
  st = (TPaveStats *)c->GetPrimitive("stats");
  if(st) {
    st->SetFillColor(kWhite);
    st->SetFillStyle(1001); // solid
    st->SetX1NDC(statX3); st->SetX2NDC(statX3+lineW); st->SetY1NDC(statY3); st->SetY2NDC(statY3+lineH);
    st->Draw();
  }
  c->Update();
  c->Modified();
  c->Print(((string)(prefix)+".pdf").c_str());

  //// 2-prong, alpha endpoint X vs Y
  c->Clear();
  h2=(TH2D*)f->Get("h_2prong_alpha_endXY")->Clone(); // copy for modifications of the same histogram
  h2->SetTitle((energyInfo+h2->GetTitle()).c_str());
  gStyle->SetOptStat(10); // KISOURMEN : show entries only
  h2->Rebin2D(5,5);
  //  h2->Rebin2D(10,10);
  h2->UseCurrentStyle();
  area1->GetXaxis()->SetTitle(h2->GetXaxis()->GetTitle());
  area1->GetYaxis()->SetTitle(h2->GetYaxis()->GetTitle());
  area1->SetTitle(h2->GetTitle());
  area1->Draw("AL"); // UVW perimeter
  h2->Draw("COLZ SAMES"); // draw stat box
  area2->Draw("L SAME"); // UVW perimeter with 5mm veto band
  h2->SetStats(true);
  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.14);
  area1->GetXaxis()->SetTitleOffset(1.6);
  area1->GetYaxis()->SetTitleOffset(1.4);
  h2->SetTitleOffset(1.2, "Z");
  gPad->Update();
  st = (TPaveStats *)c->GetPrimitive("stats");
  if(st) {
    st->SetFillColor(kWhite);
    st->SetFillStyle(1001); // solid
    st->SetX1NDC(statX3); st->SetX2NDC(statX3+lineW); st->SetY1NDC(statY3); st->SetY2NDC(statY3+lineH);
    st->Draw();
  }
  c->Update();
  c->Modified();
  c->Print(((string)(prefix)+".pdf").c_str());

  //// 2-prong, carbon endpoint X vs Y
  c->Clear();
  h2=(TH2D*)f->Get("h_2prong_carbon_endXY")->Clone(); // copy for modifications of the same histogram
  h2->SetTitle((energyInfo+h2->GetTitle()).c_str());
  gStyle->SetOptStat(10); // KISOURMEN : show entries only
  h2->Rebin2D(5,5);
  //  h2->Rebin2D(10,10);
  h2->UseCurrentStyle();
  area1->GetXaxis()->SetTitle(h2->GetXaxis()->GetTitle());
  area1->GetYaxis()->SetTitle(h2->GetYaxis()->GetTitle());
  area1->SetTitle(h2->GetTitle());
  area1->Draw("AL"); // UVW perimeter
  h2->Draw("COLZ SAMES"); // draw stat box
  area2->Draw("L SAME"); // UVW perimeter with 5mm veto band
  h2->SetStats(true);
  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.14);
  area1->GetXaxis()->SetTitleOffset(1.6);
  area1->GetYaxis()->SetTitleOffset(1.4);
  h2->SetTitleOffset(1.2, "Z");
  gPad->Update();
  st = (TPaveStats *)c->GetPrimitive("stats");
  if(st) {
    st->SetFillColor(kWhite);
    st->SetFillStyle(1001); // solid
    st->SetX1NDC(statX3); st->SetX2NDC(statX3+lineW); st->SetY1NDC(statY3); st->SetY2NDC(statY3+lineH);
    st->Draw();
  }
  c->Update();
  c->Modified();
  c->Print(((string)(prefix)+".pdf").c_str());

  //// 2-prong, vertex X vs Y (ALL) 
  c->Clear();
  h2=(TH2D*)f->Get("h_2prong_vertexXY")->Clone(); // copy for modifications of the same histogram
  h2->SetTitle((energyInfo+h2->GetTitle()).c_str());
  gStyle->SetOptStat(10); // KISOURMEN : show entries only
  h2->Rebin2D(3,3);
  //  h2->Rebin2D(10,10);
  h2->UseCurrentStyle();
  area1->GetXaxis()->SetTitle(h2->GetXaxis()->GetTitle());
  area1->GetYaxis()->SetTitle(h2->GetYaxis()->GetTitle());
  area1->SetTitle(h2->GetTitle());
  area1->Draw("AL"); // UVW perimeter
  h2->Draw("COLZ SAMES"); // draw stat box
  area2->Draw("L SAME"); // UVW perimeter with 5mm veto band
  h2->SetStats(true);
  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.14);
  area1->GetXaxis()->SetTitleOffset(1.6);
  area1->GetYaxis()->SetTitleOffset(1.4);
  h2->SetTitleOffset(1.2, "Z");
  gPad->Update();
  st = (TPaveStats *)c->GetPrimitive("stats");
  if(st) {
    st->SetFillColor(kWhite);
    st->SetFillStyle(1001); // solid
    st->SetX1NDC(statX3); st->SetX2NDC(statX3+lineW); st->SetY1NDC(statY3); st->SetY2NDC(statY3+lineH);
    st->Draw();
  }
  c->Update();
  c->Modified();
  c->Print(((string)(prefix)+".pdf").c_str());

  //// 2-prong, vertex X vs Y, zoomed Y=[-15mm, 15mm]
  c->Clear();
  h2=(TH2D*)f->Get("h_2prong_vertexXY")->Clone(); // copy for modifications of the same histogram
  h2->SetTitle((energyInfo+h2->GetTitle()).c_str());
  gStyle->SetOptStat(1110); // KISOURMEN : show entries + mean + rms
  h2->Rebin2D(5,1);
  h2->UseCurrentStyle();
  h2->SetStats(true);
  h2->Draw("COLZ");
  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.14);
  h2->SetTitleOffset(1.6, "X");
  h2->SetTitleOffset(1.4, "Y");
  h2->SetTitleOffset(1.2, "Z");
  h2->GetYaxis()->SetRangeUser(-15.0, 15.0);
  gPad->Update();
  st = (TPaveStats *)c->GetPrimitive("stats");
  if(st) {
    st->SetX1NDC(statX3); st->SetX2NDC(statX3+lineW); st->SetY1NDC(statY3-lineH*2); st->SetY2NDC(statY3+lineH);
  }
  c->Update();
  c->Modified();
  c->Print(((string)(prefix)+".pdf").c_str());

  //// 2-prong, profile histogram: vertex X vs average vertex Y
  c->Clear();
  hp=(TProfile*)f->Get("h_2prong_vertexXY_prof");
  hp->SetTitle((energyInfo+hp->GetTitle()).c_str());
  gStyle->SetOptStat(10); // KISOURMEN : show entries only
  hp->Rebin(1);
  hp->Draw();
  hp->UseCurrentStyle();
  hp->SetStats(true);
  hp->SetLineWidth(2);
  TF1 *tiltFunc=new TF1("tiltFunc", "[0]+[1]*x", hp->GetXaxis()->GetXmin()+20, hp->GetXaxis()->GetXmax()-20);
  tiltFunc->SetParNames("Offset", "Slope");
  hp->Fit(tiltFunc);
  hp->Draw("E1 SAME");
  gPad->SetLeftMargin(0.125);
  gPad->SetRightMargin(0.1);
  gPad->Update();
  st = (TPaveStats *)c->GetPrimitive("stats");
  if(st) {
    st->SetX1NDC(statX5); st->SetX2NDC(statX5+lineW*1.5); st->SetY1NDC(statY5); st->SetY2NDC(statY5+lineH*2.5);
  }
  //  if(st) {
  //    st->SetX1NDC(statX4); st->SetX2NDC(statX4+lineW); st->SetY1NDC(statY4-lineH); st->SetY2NDC(statY4+lineH);
  //  }
  c->Update();
  c->Modified();
  c->Print(((string)(prefix)+".pdf").c_str());

  //// 2-prong, vertex X distribution (along beam)
  c->Clear();
  h1=(TH1D*)f->Get("h_2prong_vertexX");
  h1->SetTitle((energyInfo+h1->GetTitle()).c_str());
  gStyle->SetOptStat(1110); // KISOURMEN : show entries + mean + rms
  h1->Rebin(5);
  h1->Draw();
  h1->UseCurrentStyle();
  h1->SetStats(true);
  h1->SetLineWidth(2);
  gPad->SetLeftMargin(0.125);
  gPad->SetRightMargin(0.1);
  gPad->Update();
  st = (TPaveStats *)c->GetPrimitive("stats");
  if(st) {
    st->SetX1NDC(statX4); st->SetX2NDC(statX4+lineW); st->SetY1NDC(statY4-lineH); st->SetY2NDC(statY4+lineH);
  }
  c->Update();
  c->Modified();
  c->Print(((string)(prefix)+".pdf").c_str());

  //// 2-prong, vertex Y distribution, zoomed Y=[-15mm, 15mm], no correction for beam tilt
  c->Clear();
  h1=(TH1D*)f->Get("h_2prong_vertexY");
  h1->SetTitle((energyInfo+h1->GetTitle()).c_str());
  gStyle->SetOptStat(1110); // KISOURMEN : show entries + mean + rms
  h1->Draw();
  h1->UseCurrentStyle();
  h1->SetStats(true);
  h1->SetLineWidth(2);
  h1->GetXaxis()->SetRangeUser(-15.0, 15.0);
  gPad->SetLeftMargin(0.125);
  gPad->SetRightMargin(0.1);
  gPad->Update();
  st = (TPaveStats *)c->GetPrimitive("stats");
  if(st) {
    st->SetX1NDC(statX1); st->SetX2NDC(statX1+lineW); st->SetY1NDC(statY1-lineH); st->SetY2NDC(statY1+lineH);
  }
  c->Update();
  c->Modified();
  c->Print(((string)(prefix)+".pdf").c_str());

  //// 2-prong, Sum of alpha + C track lengths, zoomed X=[30mm, 130mm], log Y scale
  c->Clear();
  h1=(TH1D*)f->Get("h_2prong_alpha_len");
  h1->SetTitle((energyInfo+h1->GetTitle()).c_str());
  gStyle->SetOptStat(10); // KISOURMEN : show entries only
  h1->Draw();
  h1->UseCurrentStyle();
  h1->SetStats(true);
  h1->SetLineWidth(2);
  h1->GetXaxis()->SetRangeUser(30.0, 130.0);
  gPad->SetLeftMargin(0.125);
  gPad->SetRightMargin(0.1);
  gPad->SetLogy(true);
  gPad->Update();
  st = (TPaveStats *)c->GetPrimitive("stats");
  if(st) {
    st->SetX1NDC(statX1); st->SetX2NDC(statX1+lineW); st->SetY1NDC(statY1); st->SetY2NDC(statY1+lineH);
  }
  c->Update();
  c->Modified();
  c->Print(((string)(prefix)+".pdf").c_str());

  //// 2-prong, Sum of alpha + C track lengths, zoomed X=[60mm, 80mm], linear Y scale
  c->Clear();
  h1=(TH1D*)f->Get("h_2prong_alpha_len");
  h1->SetTitle((energyInfo+h1->GetTitle()).c_str());
  gStyle->SetOptStat(10); // KISOURMEN : show entries only
  h1->Draw();
  h1->UseCurrentStyle();
  h1->SetStats(true);
  h1->SetLineWidth(2);
  h1->GetXaxis()->SetRangeUser(60.0, 80.0);
  gPad->SetLeftMargin(0.125);
  gPad->SetRightMargin(0.1);
  gPad->SetLogy(false);
  gPad->Update();
  st = (TPaveStats *)c->GetPrimitive("stats");
  if(st) {
    st->SetX1NDC(statX1); st->SetX2NDC(statX1+lineW); st->SetY1NDC(statY1); st->SetY2NDC(statY1+lineH);
  }
  c->Update();
  c->Modified();
  c->Print(((string)(prefix)+".pdf").c_str());

  //// 2-prong, Alpha track length (X) vs Carbon track length (Y), zoomed X=[20mm, 120mm] x Y=[4mm, 18mm]
  c->Clear();
  h2=(TH2D*)f->Get("h_2prong_alpha_len_carbon_len")->Clone(); // copy for modifications of the same histogram
  h2->SetTitle((energyInfo+h2->GetTitle()).c_str());
  gStyle->SetOptStat(10); // KISOURMEN : show entries only
  h2->Rebin2D(2,1);
  h2->UseCurrentStyle();
  h2->SetStats(true);
  h2->Draw("COLZ");
  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.14);
  h2->SetTitleOffset(1.6, "X");
  h2->SetTitleOffset(1.4, "Y");
  h2->SetTitleOffset(1.2, "Z");
  h2->GetXaxis()->SetRangeUser(20, 120);
  h2->GetYaxis()->SetRangeUser(4, 18);
  gPad->Update();
  st = (TPaveStats *)c->GetPrimitive("stats");
  if(st) {
    st->SetX1NDC(statX3); st->SetX2NDC(statX3+lineW); st->SetY1NDC(statY3); st->SetY2NDC(statY3+lineH);
  }
  c->Update();
  c->Modified();
  c->Print(((string)(prefix)+".pdf").c_str());

  //// define 2D cut to select Oxygen-16 only for next 3 plots
  Double_t theta[5]  = {-1,    1,       1,    -1, -1};
  Double_t length[5] = {58, 58+5, 58+5+20, 58+20, 58};
  TCutG *mycut=new TCutG("select_O16", 5, theta, length);
  mycut->SetLineColor(kRed);
  mycut->SetLineWidth(3);
 
  //// 2-prong, Alpha track length (X) vs Alpha cos(theta_BEAM_LAB), zoomed Y=[20mm, 120mm]
  c->Clear();
  h2=(TH2D*)f->Get("h_2prong_alpha_cosThetaBEAM_len_LAB")->Clone(); // copy for modifications of the same histogram
  h2->SetTitle((energyInfo+h2->GetTitle()).c_str());
  gStyle->SetOptStat(10); // KISOURMEN : show entries only
  h2->Rebin2D(2,2);
  h2->UseCurrentStyle();
  h2->SetStats(true);
  h2->Draw("COLZ");
  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.14);
  h2->SetTitleOffset(1.6, "X");
  h2->SetTitleOffset(1.4, "Y");
  h2->SetTitleOffset(1.2, "Z");
  h2->GetYaxis()->SetRangeUser(20, 120);
  gPad->Update();
  st = (TPaveStats *)c->GetPrimitive("stats");
  if(st) {
    st->SetX1NDC(statX3); st->SetX2NDC(statX3+lineW); st->SetY1NDC(statY3); st->SetY2NDC(statY3+lineH);
  }
  mycut->Draw("SAME");
  c->Update();
  c->Modified();
  c->Print(((string)(prefix)+".pdf").c_str());

  //// 2-prong, Alpha track length (X) vs Alpha cos(theta_BEAM_LAB), zoomed Y=[50mm, 90mm], tagged O-16
  c->Clear();
  h2=(TH2D*)f->Get("h_2prong_alpha_cosThetaBEAM_len_LAB")->Clone(); // copy for modifications of the same histogram
  h2->SetTitle((energyInfo+h2->GetTitle()).c_str());
  //  gStyle->SetOptStat(10); // KISOURMEN : show entries only
  h2->Rebin2D(2,2);
  h2->UseCurrentStyle();
  h2->SetStats(false);
  h2->Draw("COLZ [select_O16]");
  h2->SetTitle(Form("%s, ^{16}O candidates", h2->GetTitle())); // modify the title
  mycut->Draw("SAME");
  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.14);
  h2->SetTitleOffset(1.6, "X");
  h2->SetTitleOffset(1.4, "Y");
  h2->SetTitleOffset(1.2, "Z");
  h2->GetYaxis()->SetRangeUser(50, 90);
  gPad->Update();
  //  st = (TPaveStats *)c->GetPrimitive("stats");
  //  if(st) {
  //    st->SetX1NDC(statX3); st->SetX2NDC(statX3+lineW); st->SetY1NDC(statY3); st->SetY2NDC(statY3+lineH);
  //  }
  c->Update();
  c->Modified();
  c->Print(((string)(prefix)+".pdf").c_str());

  //// 2-prong, Alpha cos(theta_BEAM_LAB), tagged Oxygen-16
  h2=(TH2D*)f->Get("h_2prong_alpha_cosThetaBEAM_len_LAB");
  h2->SetTitle((energyInfo+h2->GetTitle()).c_str());
  h1=h2->ProjectionX("_px",1, h2->GetNbinsY(),"[select_O16]");
  h1->Rebin(4);
  h1->SetYTitle(h2->GetZaxis()->GetTitle());
  h1->SetTitle(Form("%s, ^{16}O candidates", h2->GetTitle())); // modify the title
  h1->Draw("HIST E1");
  h1->UseCurrentStyle();
  h1->SetStats(true);
  h1->SetLineWidth(2);
  h1->GetXaxis()->SetRangeUser(60.0, 80.0);
  h1->GetMinimum(0);
  gPad->SetLeftMargin(0.125);
  gPad->SetRightMargin(0.1);
  gPad->SetLogy(false);
  gPad->Update();
  st = (TPaveStats *)c->GetPrimitive("stats");
  if(st) {
    st->SetX1NDC(statX1); st->SetX2NDC(statX1+lineW); st->SetY1NDC(statY1); st->SetY2NDC(statY1+lineH);
  }
  c->Update();
  c->Modified();
  c->Print(((string)(prefix)+".pdf").c_str());

  //// 2-prong, Alpha phi_BEAM_LAB, no correction for beam tilt
  c->Clear();
  h1=(TH1D*)f->Get("h_2prong_alpha_phiBEAM_LAB");
  h1->SetTitle((energyInfo+h1->GetTitle()).c_str());
  gStyle->SetOptStat(10); // KISOURMEN : show entries only
  h1->Rebin(5);
  h1->Draw();
  h1->UseCurrentStyle();
  h1->SetStats(true);
  h1->SetLineWidth(2);
  h1->SetMinimum(0);
  gPad->SetLeftMargin(0.125);
  gPad->SetRightMargin(0.1);
  gPad->SetLogy(false);
  TF1 *phiFunc=new TF1("phiFunc", "[0]*(1+[1]*cos(2*(x+[2])))", h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
  phiFunc->SetParNames("A", "f", "#delta");
  h1->Fit(phiFunc);
  h1->Draw("E1 SAME");
  gPad->Update();
  st = (TPaveStats *)c->GetPrimitive("stats");
  if(st) {
    st->SetX1NDC(statX5); st->SetX2NDC(statX5+lineW*1.5); st->SetY1NDC(statY5); st->SetY2NDC(statY5+lineH*2.5);
  }
  TLatex *tl=new TLatex;
  tl->SetTextSize(0.030);
  tl->SetTextAlign(11); // align at bottom, left
  tl->DrawLatexNDC(0.2, 0.275, "Fit :  #font[12]{#frac{dN}{d#phi} = A #[]{1 + f cos 2#(){#phi + #delta}}}");
  c->Update();
  c->Modified();
  c->Print(((string)(prefix)+".pdf").c_str());
  
  //// 2-prong, (Alpha track enpoint_Y vs vertex_Y) vs (Alpha track enpoint_Z vs vertex_Z)
  c->Clear();
  h2=(TH2D*)f->Get("h_2prong_alpha_deltaYZ")->Clone(); // copy for modifications of the same histogram
  h2->SetTitle((energyInfo+h2->GetTitle()).c_str());
  gStyle->SetOptStat(10); // KISOURMEN : show entries only
  h2->Rebin2D(10,10);
  h2->UseCurrentStyle();
  h2->SetStats(true);
  h2->Draw("COLZ");
  gPad->SetLeftMargin(0.200);
  gPad->SetRightMargin(0.200);
  h2->SetTitleOffset(1.6, "X");
  h2->SetTitleOffset(1.4, "Y");
  h2->SetTitleOffset(1.25, "Z");
  h2->SetXTitle("#(){#alpha track end wrt vertex}_{Y_{DET}} [mm]");
  h2->SetYTitle("#(){#alpha track end wrt vertex}_{Z_{DET}} [mm]");
  gPad->Update();
  st = (TPaveStats *)c->GetPrimitive("stats");
  if(st) {
    st->SetX1NDC(statX6); st->SetX2NDC(statX6+lineW); st->SetY1NDC(statY6); st->SetY2NDC(statY6+lineH);
  }
  c->Update();
  c->Modified();
  c->Print(((string)(prefix)+".pdf").c_str());
  
  //////////////////////////////////
  f->Close();
  c->Print(((string)(prefix)+".pdf]").c_str());

  return;
}
