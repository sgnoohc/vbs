// NanoCORE
#include "Nano.h"
#include "Config.h"
#include "ElectronSelections.h"
#include "MuonSelections.h"
#include "tqdm.h"

#include "core/collections.h"
#include "core/cuts.h"
#include "vvv/cuts.h"
// RAPIDO
#include "arbusto.h"
#include "looper.h"
#include "cutflow.h"
#include "utilities.h"
// ROOT
#include "TString.h"
#include "Math/VectorUtil.h"

typedef std::vector<LorentzVector> LorentzVectors;
typedef std::vector<double> Doubles;
typedef std::vector<int> Integers;
typedef std::vector<unsigned int> Indices;

int main(int argc, char** argv) 
{
    // CLI
    HEPCLI cli = HEPCLI(argc, argv);

    // Initialize Looper
    Looper looper = Looper(cli);
    // Initialize Arbusto
    Arbusto arbusto = Arbusto(
        cli,
        {
            "Electron*", "nElectron",
            "Muon*", "nMuon",
            "Tau*", "nTau", 
            "Jet*", "nJet", 
            "FatJet*", "nFatJet", 
            "GenPart*", "nGenPart",
            "GenJet*", "nGenJet", 
            "Generator*",
            "MET*",
            "event*",
            "run*",
            "luminosityBlock*",
            "genWeight*",
            "btagWeight*",
            "LHE*",
            "*Weight*",
            "Flag*",
            "SubJet*",
            "HLT_*",
            "Pileup*",
            "fixedGridRhoFastjetAll"
        }
    );

    // Initialize Cutflow
    Cutflow cutflow = Cutflow(cli.output_name+"_Cutflow");
    cutflow.globals.newVar<LorentzVectors>("veto_lep_p4s", {});
    cutflow.globals.newVar<LorentzVectors>("tight_lep_p4s", {});
    cutflow.globals.newVar<LorentzVectors>("jet_p4s", {});
    cutflow.globals.newVar<LorentzVectors>("fatjet_p4s", {});

    Core::Skimmer skimmer = Core::Skimmer(arbusto, nt, cli, cutflow);

    /* --- Assemble cutflow --- */

    Cut* base = new LambdaCut("Base", [&]() { return true; });
    cutflow.setRoot(base);

    Cut* find_leps = new VVV::FindLeptons("FindLeptons", skimmer);
    cutflow.insert(base, find_leps, Right);

    Cut* find_fatjets = new VVV::FindFatJets("FindFatJets", skimmer);
    cutflow.insert(find_leps, find_fatjets, Right);

    Cut* skim_cuts = new LambdaCut(
        "SkimSelection", 
        [&]() 
        { 
            float HT = 0;
            for (auto& jet_p4 : nt.Jet_p4()) HT += jet_p4.pt();
            for (auto& fatjet_p4 : nt.FatJet_p4()) HT += fatjet_p4.pt();
            if (HT < 1100)
                return false;
            if (cutflow.globals.getVal<LorentzVectors>("veto_lep_p4s").size() == 0)
            {
                return (cutflow.globals.getVal<LorentzVectors>("fatjet_p4s").size() >= 2 and cutflow.globals.getVal<LorentzVectors>("fatjet_p4s").at(0).pt() > 500.);
            } 
            else
            {
                return (cutflow.globals.getVal<LorentzVectors>("fatjet_p4s").size() >= 1);
            }
        }
    );
    cutflow.insert(find_fatjets, skim_cuts, Right);

    /* ------------------------ */

    // Run looper
    tqdm bar;
    looper.run(
        [&](TTree* ttree)
        {
            nt.Init(ttree);
            skimmer.init(ttree);
        },
        [&](int entry) 
        {
            if (cli.debug && looper.n_events_processed == 10000) { looper.stop(); }
            else
            {
                // Reset branches and globals
                arbusto.resetBranches();
                cutflow.globals.resetVars();
                // Run cutflow
                nt.GetEntry(entry);
                bool passed = cutflow.run("SkimSelection");
                if (passed) { arbusto.fill(entry); }
                bar.progress(looper.n_events_processed, looper.n_events_total);
            }
        }
    );

    // Wrap up
    if (!cli.is_data) { cutflow.print(); }

    skimmer.write();

    return 0;
}
