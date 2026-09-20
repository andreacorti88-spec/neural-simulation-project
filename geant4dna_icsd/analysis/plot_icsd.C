void plot_icsd() {
    int energies[] = {20, 50, 100, 300, 600, 1000, 5000, 10000};
    TCanvas* c = new TCanvas("c", "ICSD per energia", 1600, 900);
    c->Divide(4, 2);

    for (int i = 0; i < 8; i++) {
        int E = energies[i];
        c->cd(i+1);
        TString fname = Form("results_opt2/icsd_%deV.csv", E);
        TTree* t = new TTree("t", "t");
        t->ReadFile(fname, "n/I", ',');
        int maxN = t->GetMaximum("n") + 2;
        TH1D* h = new TH1D(Form("h%d", E), Form("%d eV;ionizzazioni;eventi", E), maxN, 0, maxN);
        t->Draw(Form("n>>h%d", E), "", "goff");
        h->SetFillColor(kAzure-4);
        h->Draw();
    }
    c->SaveAs("icsd_histograms.png");
}
