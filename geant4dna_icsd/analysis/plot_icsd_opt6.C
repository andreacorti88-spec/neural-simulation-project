void plot_icsd_opt6() {
    int energies[] = {20, 50, 100, 300, 600, 1000, 5000, 10000};
    TCanvas* c = new TCanvas("c6", "ICSD per energia - opzione 6", 1600, 900);
    c->Divide(4, 2);
    for (int i = 0; i < 8; i++) {
        int E = energies[i];
        c->cd(i+1);
        TString fname = Form("results_opt6/icsd_%deV.csv", E);
        TTree* t = new TTree("t6", "t6");
        t->ReadFile(fname, "n/I", ',');
        int maxN = t->GetMaximum("n") + 2;
        TH1D* h = new TH1D(Form("h6_%d", E), Form("%d eV (opt.6);ionizzazioni;eventi", E), maxN, 0, maxN);
        t->Draw(Form("n>>h6_%d", E), "", "goff");
        h->SetFillColor(kGreen+2);
        h->Draw();
        delete t;
    }
    c->SaveAs("icsd_histograms_opt6.png");
}
