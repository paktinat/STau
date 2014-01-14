if (HBHENoiseFilterResult == 1) {
  int n[72]={0};
  for (int i = 0; i < 72; i++) {
    for (HBHERecHitCollection::const_iterator HBHEiter = hbhe_rechit->begin(); HBHEiter != hbhe_rechit->end(); ++ HBHEiter) {
      HcalDetId ID = HBHEiter->id();
      int RBXIndexHit = HcalHPDRBXMap::indexRBX(ID);
      if (i == RBXIndexHit) n[i]++;
    }
    nrRecHitsPerRBX->Fill(n[i]);
    nrRecHitsVsRBXIndex->Fill(i,n[i]);
  }
 }
