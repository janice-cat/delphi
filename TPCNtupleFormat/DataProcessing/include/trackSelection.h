#ifndef TRACKSELECTION_H
#define TRACKSELECTION_H

#include <vector>

#include "TMath.h"
#include "TLorentzVector.h"
#include "DataProcessing/include/particleData.h"
// #include "DataProcessing/include/alephTrkEfficiency.h"

/*example do getting highPurity requirement for ith particle from pData:
TrackSelection trkSel = TrackSelection();
particleData pData;
int i;
bool isGoodTrk = TrackSelection.highPurity(&pData,i);
*/

class TrackSelection{
 public:
   TrackSelection();
   ~TrackSelection();

   void setnTPCCut(Short_t cut);
   void setThetaCutLow(float cut);
   void setThetaCutHigh(float cut);
   void setAbsCosThCut(float cut);
   void setPCut(float cut);
   void setPtCut(float cut);
   void setD0Cut(float cut);
   void setZ0Cut(float cut);

   inline bool passesNTPC(Short_t ntpc);
   inline bool passesTheta(float theta);
   inline bool passesAbsCosThCut(float theta);
   inline bool passesP(float p);
   inline bool passesPt(float pt);
   inline bool passesD0(float d0);
   inline bool passesZ0(float z0);
   inline bool passesPWFlag(Short_t pwflag);
   // inline bool nonzeroTrkEff(float theta, float phi, float pt, int nTrk);

   bool isConversionElectron(particleData * p, int indx);

   bool highPurity(particleData * p, int indx);
   bool highPurityBit(particleData * p, int indx);
   void fillHighPurity(particleData * p); 
 private:
   Short_t nTPCcut = 4;
   float thetaCutLow = 20.*TMath::Pi()/180.;       //currently not used in highPurity
   float thetaCutHigh = 160.*TMath::Pi()/180.;     //currently not used in highPurity
   float pCut = 0.2;                               //currently not used in highPurity
   float absCosThCut = 0.94;                       //maximum abs(cos(th)) of charged tracks
   float ptCut = 0.2;
   float d0Cut = 2;
   float z0Cut = 10;
   // alephTrkEfficiency effCorrector;

   float conversionDPhi = 0.05;
   float conversionDTheta = 0.05;
};

bool TrackSelection::highPurity(particleData * p, int indx){
  if(!passesPWFlag(p->pwflag[indx])) return false;
  if(!passesAbsCosThCut(p->theta[indx]))   return false;
  if(!passesPt(p->pt[indx]))         return false;
  if(!passesD0(p->d0[indx]))         return false;
  if(!passesZ0(p->z0[indx]))         return false;
  if(!passesNTPC(p->ntpc[indx]))     return false;

  return true;
}

bool TrackSelection::highPurityBit(particleData * p, int indx){
  if(!(p->highPurity[indx])) return false;
  return true;
}


void TrackSelection::fillHighPurity(particleData * p){
  for(int i = 0; i<p->nParticle; i++){
    if(highPurity(p,i)) p->highPurity[i] = true;
    else p->highPurity[i] = false;
  }
}

inline bool TrackSelection::passesNTPC(Short_t ntpc){ return ntpc>=nTPCcut; }
inline bool TrackSelection::passesTheta(float theta){ return (theta>=thetaCutLow) && (theta<=thetaCutHigh); }
inline bool TrackSelection::passesAbsCosThCut(float theta){ return TMath::Abs(cos(theta))<=absCosThCut; }
inline bool TrackSelection::passesP(float p){ return p>=pCut;}
inline bool TrackSelection::passesPt(float pt){ return pt>=ptCut;}
inline bool TrackSelection::passesD0(float d0){ return TMath::Abs(d0)<=d0Cut;}
inline bool TrackSelection::passesZ0(float z0){ return TMath::Abs(z0)<=z0Cut;}
inline bool TrackSelection::passesPWFlag(Short_t pwflag){ return (pwflag>=0 && pwflag<=2);}
// inline bool TrackSelection::nonzeroTrkEff(float theta,float phi, float pt, int nTrk){ return effCorrector.efficiency(theta,phi,pt,nTrk) > 0; }



void TrackSelection::setnTPCCut(Short_t cut){ nTPCcut = cut;}
void TrackSelection::setThetaCutLow(float cut){ thetaCutLow = cut;}
void TrackSelection::setThetaCutHigh(float cut){ thetaCutHigh = cut;}
void TrackSelection::setAbsCosThCut(float cut){ absCosThCut = cut;}
void TrackSelection::setPCut(float cut){ pCut = cut;}
void TrackSelection::setPtCut(float cut){ ptCut = cut;}
void TrackSelection::setD0Cut(float cut){ d0Cut = cut;}
void TrackSelection::setZ0Cut(float cut){ z0Cut = cut;}

//checks the particle against the previous particle in the list of particles and sees if the pair looks like a conversion to e+e-
bool TrackSelection::isConversionElectron(particleData * p, int indx){ 
  //deal with the fact that there is nothing before 0 by comparing it to index 1 for the 'relative' comparisons
  if(indx==0) indx++;

  //both electrons
  if( p->pwflag[indx] != 2 ) return false;
  if( p->pwflag[indx-1] != 2 ) return false;
 
  //opposite charge required 
  if(p->charge[indx] != -(p->charge[indx-1])) return false;

  //dtheta and dphi matching
  if( TMath::Abs(p->theta[indx] - p->theta[indx-1]) > conversionDTheta) return false;
  if( TMath::ACos(TMath::Cos(p->phi[indx] - p->phi[indx-1])) > conversionDPhi) return false;

  return true;
}


TrackSelection::TrackSelection(){}

TrackSelection::~TrackSelection(){}

#endif
