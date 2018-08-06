const char* genToken_ = "genParticles";
char delim = ',';

ostream ofs("protons.csv");

void getProton(edm::Event iEvent) {
    edm::Handle<edm::View<reco::GenParticle> > genPartColl;
    iEvent.getByToken( genToken_, genPartColl );

    for ( unsigned int i = 0; i < genPartColl->size(); ++i ) {
        const edm::Ptr<reco::GenParticle> genPart = genPartColl->ptrAt( i );

        if ( !genPart->isPromptFinalState() ) continue;

        // generated outgoing protons
        if ( genPart->pdgId() == 2212) {
            ofs << genPart->pt() << delim;
            ofs << genPart->eta() << delim;
            ofs << genPart->phi() << delim;
            ofs << genPart->energy() << delim;
            ofs << "\n";
        }
}


void extractProton(const char* proton_gensim = "../../../SDGEN.root") {
    TFile* file = new TFile(proton_gensim);
    TTree* tree = (TTree*)file->Get("Events");
    edm::Event iEvent;
    tree->SetBranchAddress("recoGenParticles_genParticles__GEN", iEvent);
    int N = tree->GetEntriesFast();

    for(int i = 0; i < N; i++) {
        cout << i << endl;
        getProton()
    }
};