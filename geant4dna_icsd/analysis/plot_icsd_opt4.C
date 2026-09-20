void plot_icsd_opt4() {
    int energies[] = {20, 50, 100, 300, 600, 1000, 5000, 10000};
    TCanvas* c = new TCanvas("c4", "ICSD per energia - opzione 4", 1600, 900);
    c->Divide(4, 2);
    for (int i = 0; i < 8; i++) {
        int E = energies[i];
        c->cd(i+1);
        TString fname = Form("results_opt4/icsd_%deV.csv", E);
        TTree* t = new TTree("t4", "t4");
        t->ReadFile(fname, "n/I", ',');
        int maxN = t->GetMaximum("n") + 2;
        TH1D* h = new TH1D(Form("h4_%d", E), Form("%d eV (opt.4);ionizzazioni;eventi", E), maxN, 0, maxN);
        t->Draw(Form("n>>h4_%d", E), "", "goff");
        h->SetFillColor(kOrange+1);
        h->Draw();
        delete t;
    }
    c->SaveAs("icsd_histograms_opt4.png");
}
