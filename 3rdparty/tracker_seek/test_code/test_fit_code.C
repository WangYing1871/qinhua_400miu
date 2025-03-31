void test_fit_code(){

    TCanvas *c = new TCanvas("c","c",0,0,800,800);
    TCanvas *c1 = new TCanvas("c1","c1",0,0,800,800);
    TCanvas *c2 = new TCanvas("c2","c2",500,0,800,800);
    TCanvas *c3 = new TCanvas("c3","c3",500,0,800,800);
    
    TH1D *h = new TH1D("h","h",100,-10,800);
    TH1D *rd = new TH1D("rd","rd",100,0,100);
    TH1D *sloph = new TH1D("sloph","sloph",100,-2,2);
    TH1D *interh = new TH1D("inter","inter",100,0,100);
    TH1D *res[6];
    for(int i=0; i<6; i++){
      res[i] = new TH1D(Form("res_%d",i),Form("res_%d",i),100,-20,20);
    }
    int point = 6;
    double zi[6]={0.,150.,300.,450.,600.,750.};
    double xi[6]={50.,80.,110.,140.,170.,200.};
    double slop;
    double intercept;

    TRandom3 rm;

    TF1 *f = new TF1("gaus","[0]*x+[1]",-10,800);

    TGraph *g;
    for(int n=0;n<1000;n++){
   
      xi[1] = rm.Gaus(80,8);
      xi[2] = rm.Gaus(110,10);
      xi[3] = rm.Gaus(140,10);
      xi[4] = rm.Gaus(170,10);
      g = new TGraph(point,zi,xi);
      
      g->Fit(f,"Q");
      slop = f->GetParameter(0);
      intercept = f->GetParameter(1);
      
      sloph->Fill(slop);
      interh->Fill(intercept);

      double mse = 0;
      double prox[6] = {0};
      for(int i=0;i<point;i++){

             mse += (xi[i]-slop*zi[i]-intercept)*(xi[i]-slop*zi[i]-intercept);
             prox[i] = slop*zi[i]+intercept-xi[i];
             res[i]->Fill(prox[i]);
      }
      double rmse = sqrt(mse/point);
      rd->Fill(rmse);
      cout<<"rmse:  "<<rmse<<endl;
    }
      c->cd();
      h->Draw();
      h->SetMaximum(300);
      g->Draw("samePL");
      g->SetMarkerStyle(20);
      g->SetMarkerSize(2);

      c1->cd();
      rd->Draw();
      c2->Divide(3,2);
      for(int i=0; i<6; i++){
      c2->cd(i+1);
      res[i]->Draw();
      }

      c3->Divide(2,1);
      c3->cd(1);
      sloph->Draw();
      c3->cd(2);
      interh->Draw();

}
