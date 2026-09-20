void plot_comparison() {
    double E[] = {20, 50, 100, 300, 600, 1000, 5000, 10000};
    double M1_sim[]   = {0.770, 2.212, 4.317, 5.297, 2.850, 1.843, 8.833, 4.923};
    double M1_paper[] = {0.770, 2.210, 4.300, 5.350, 2.900, 1.830, 7.100, 4.020};

    TGraph* gSim = new TGraph(8, E, M1_sim);
    TGraph* gPaper = new TGraph(8, E, M1_paper);
    gSim->SetMarkerStyle(20); gSim->SetMarkerColor(kBlue); gSim->SetLineColor(kBlue);
    gPaper->SetMarkerStyle(21); gPaper->SetMarkerColor(kRed); gPaper->SetLineColor(kRed);

    TCanvas* c = new TCanvas("c2", "M1: simulato vs paper", 900, 600);
    c->SetLogx();
    gSim->SetTitle("M1: nanoICSD vs paper (Geant4-DNA opt.2);Energia (eV);M1");
    gSim->Draw("APL");
    gPaper->Draw("PL SAME");
    auto leg = new TLegend(0.15, 0.7, 0.4, 0.85);
    leg->AddEntry(gSim, "nanoICSD (mio, G4 11.4)", "lp");
    leg->AddEntry(gPaper, "Paper 2026 (Tab.3)", "lp");
    leg->Draw();
    c->SaveAs("m1_comparison.png");
}
