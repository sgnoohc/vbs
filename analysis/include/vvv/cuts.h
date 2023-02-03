#ifndef VVV_CUTS_H
#define VVV_CUTS_H

// RAPIDO
#include "arbol.h"
#include "cutflow.h"
#include "utilities.h"
// VBS
#include "core/collections.h"
#include "core/cuts.h"
#include "core/pku.h"
#include "core/vvv.h"
#include "corrections/all.h"

namespace VVV
{

class FindLeptons : public Core::SkimmerCut
{
public:
    FindLeptons(std::string name, Core::Skimmer& skimmer) : Core::SkimmerCut(name, skimmer) 
    {
        // Do nothing
    };

    virtual bool passesVetoElecID(int elec_i)
    {
        return VVV::passesElecID(elec_i, VVV::IDveto);
    };

    virtual bool passesVetoMuonID(int muon_i)
    {
        return VVV::passesMuonID(muon_i, VVV::IDveto);
    };

    virtual bool passesTightElecID(int elec_i)
    {
        return VVV::passesElecID(elec_i, VVV::IDtight);
    };

    virtual bool passesTightMuonID(int muon_i)
    {
        return VVV::passesMuonID(muon_i, VVV::IDtight);
    };

    bool evaluate()
    {
        LorentzVectors veto_lep_p4s;
        LorentzVectors tight_lep_p4s;
        for (unsigned int elec_i = 0; elec_i < nt.nElectron(); elec_i++)
        {
            LorentzVector lep_p4 = nt.Electron_p4().at(elec_i);
            if (passesVetoElecID(elec_i)) { veto_lep_p4s.push_back(lep_p4); }
            if (passesTightElecID(elec_i)) { tight_lep_p4s.push_back(lep_p4); }
        }
        for (unsigned int muon_i = 0; muon_i < nt.nMuon(); muon_i++)
        {
            LorentzVector lep_p4 = nt.Muon_p4().at(muon_i);
            if (passesVetoMuonID(muon_i)) { veto_lep_p4s.push_back(lep_p4); }
            if (passesTightMuonID(muon_i)) { tight_lep_p4s.push_back(lep_p4); }
        }
        globals.setVal<LorentzVectors>("veto_lep_p4s", veto_lep_p4s);
        globals.setVal<LorentzVectors>("tight_lep_p4s", tight_lep_p4s);
        return true;
    };
};

class FindFatJets : public Core::SkimmerCut
{
public:
    FindFatJets(std::string name, Core::Skimmer& skimmer) : Core::SkimmerCut(name, skimmer) 
    {
        // Do nothing
    };

    virtual bool passesFatJetID(int fatjet_i)
    {
        return VVV::passesFatJetID(fatjet_i);
    };

    bool evaluate()
    {
        LorentzVectors lep_p4s = globals.getVal<LorentzVectors>("veto_lep_p4s");
        LorentzVectors fatjet_p4s;
        for (unsigned int fatjet_i = 0; fatjet_i < nt.nFatJet(); fatjet_i++)
        {
            LorentzVector fatjet_p4 = nt.FatJet_p4().at(fatjet_i);
            if (not (passesFatJetID(fatjet_i)))
                continue;
            bool is_overlap = false;
            for (auto lep_p4 : lep_p4s)
            {
                if (ROOT::Math::VectorUtil::DeltaR(lep_p4, fatjet_p4) < 0.8)
                {
                    is_overlap = true;
                    break;
                }
            }
            if (not is_overlap)
            {
                fatjet_p4s.push_back(fatjet_p4);
            }
        }
        globals.setVal<LorentzVectors>("fatjet_p4s", fatjet_p4s);
        return true;
    };
};

} // End namespace VVV;

#endif
