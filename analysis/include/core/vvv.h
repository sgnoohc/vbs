#ifndef VVV_H
#define VVV_H

// NanoCORE
#include "Nano.h"

namespace VVV
{

enum IDLevel 
{
    IDveto = 0,
    IDtight = 1
};

bool passesElecID(unsigned int elec_i, IDLevel id_level)
{
    // Selections
    if (not (nt.Electron_mvaFall17V2Iso_WP90()[elec_i])) return false;
    if (not (nt.Electron_p4()[elec_i].pt()       > 10.)) return false;
    if (not (abs(nt.Electron_p4()[elec_i].eta()) < 2.5)) return false;
    if (abs(nt.Electron_p4()[elec_i].eta()) < 1.566 && abs(nt.Electron_p4()[elec_i].eta()) > 1.444) return false; 
    return true;
};

bool passesMuonID(unsigned int muon_i, IDLevel id_level)
{
    // Selections
    if (not (nt.Muon_mediumId()[muon_i]             )) return false; // TODO: What is Muon_mediumPromptId in NanoAOD?
    if (not (nt.Muon_p4()[muon_i].pt()        > 10. )) return false;
    if (not (nt.Muon_pfRelIso04_all()[muon_i] < 0.25)) return false; // i.e. Loose from https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonSelection#Particle_Flow_isolation
    if (not (abs(nt.Muon_p4()[muon_i].eta())  < 2.4 )) return false;
    return true;
};

bool passesFatJetID(unsigned int fatjet_i)
{
    if (not (nt.FatJet_p4()[fatjet_i].pt()         > 200)) return false;
    if (not (fabs(nt.FatJet_p4()[fatjet_i].eta())  < 2.5)) return false;
    if (not (nt.FatJet_msoftdrop()[fatjet_i]       > 40 )) return false;
    if (not (nt.FatJet_msoftdrop()[fatjet_i]       < 150)) return false;
    return true;
};

};

#endif
