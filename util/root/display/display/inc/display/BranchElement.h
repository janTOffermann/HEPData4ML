/*
 *  This code is a modification of code from the Delphes fast detector
 *  simulation framework:
 *
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef BranchElement_h
#define BranchElement_h

#include "CaloData.h" // if "display/CaloData.h", works during build but breaks when loading header with gInterpeter in PyROOT... -Jan
#include "JetElementList.h"

//standard library includes
#include <exception>
#include <iostream>

//ROOT includes
#include "TClass.h"
#include "TClonesArray.h"
#include "TColor.h"
#include "TEveElement.h"
#include "TEveTrack.h"
#include "TString.h"

// Delphes includes
// #include "display/CaloData.h"

namespace Display{

  // virtual class to represent objects from a Delphes-tree branch
  class BranchBase
  {
  public:
    BranchBase(const char *name = "", const enum EColor color = kBlack, Float_t maxPt = 50.) :
      name_(name), maxPt_(maxPt), color_(color) {} // branch_(branch) ; TClonesArray *branch = NULL
    virtual ~BranchBase() {}
    const char *GetName() const { return (const char *)name_; }
    // const char *GetType() const { return branch_ ? branch_->GetClass()->GetName() : "None"; }
    virtual const char *GetClassName() = 0;
    enum EColor GetColor() const { return color_; }
    void SetColor(const enum EColor color){color_ = color;}
    virtual void Reset() = 0;
    virtual void SetTrackingVolume(Float_t r, Float_t l, Float_t Bz = 0.)
    {
      tkRadius_ = r;
      tkHalfLength_ = l;
      tk_Bz_ = Bz;
    }

    std::vector<Float_t> GetTrackingVolume(){
      return std::vector<Float_t>{tkRadius_,tkHalfLength_,tk_Bz_};
    }

    // virtual void ReadBranch() = 0;
    // virtual std::vector<TLorentzVector> GetVectors() = 0;
    // TClonesArray *GetBranch() const {return branch_;}

  protected:
    TString name_;
    Float_t maxPt_;
    // TClonesArray *branch_;
    enum EColor color_;
    Float_t tkRadius_, tkHalfLength_, tk_Bz_;
  };

  // concrete implementations. EveContainer can be a TrackList, ElementList, JetElementList or CaloData.
  template <typename EveContainer>
  class BranchElement: public BranchBase
  {
  public:
    // constructor
    BranchElement(const char *name = "", const enum EColor color = kBlack, Float_t maxPt = 50.) :
      BranchBase(name, color, maxPt)
    {
      throw std::exception();
    }

    BranchElement(EveContainer* container, const char *name = "", const enum EColor color = kBlack, Float_t maxPt = 50.) :
      BranchBase(name, color, maxPt)
    {
      data_ = container;
    }


    // destructor
    virtual ~BranchElement() { delete data_; }

    // get the container (ElementList, TrackList, or CaloData)
    EveContainer *GetContainer() { return data_; }

    // externally set hte container
    void SetContainer(EveContainer* container){data_ = container;};

    // tracking volume
    virtual void SetTrackingVolume(Float_t r, Float_t l, Float_t Bz = 0.)
    {
      tkRadius_ = r;
      tkHalfLength_ = l;
      tk_Bz_ = Bz;
    }

    // resets the collection (before moving to the next event)
    virtual void Reset(){};

    // template class name
    virtual const char *GetClassName() { return data_->ClassName(); }

    // // read the branch and fill elements for display
    // virtual void ReadBranch() {}

    // // return the vector for all elements
    // virtual std::vector<TLorentzVector> GetVectors()
    // {
    //   std::vector<TLorentzVector> v;
    //   return v;
    // }

  private:
    EveContainer *data_;
  };

  #if !defined(__CINT__) && !defined(__CLING__)

  // special case for calo towers
  template <>
  BranchElement<CaloData>::BranchElement(const char *name, const enum EColor color, Float_t maxPt);
  template <>
  void BranchElement<CaloData>::Reset();
  // template <>
  // void BranchElement<CaloData>::ReadBranch();
  // template <>
  // std::vector<TLorentzVector> BranchElement<CaloData>::GetVectors();

  // special case for element lists
  template <>
  BranchElement<TEveElementList>::BranchElement(const char *name, const enum EColor color, Float_t maxPt);
  template <>
  void BranchElement<TEveElementList>::Reset();
  // template <>
  // void BranchElement<TEveElementList>::ReadBranch();
  // template <>
  // std::vector<TLorentzVector> BranchElement<TEveElementList>::GetVectors();

  // special case for jet element lists
  template <>
  BranchElement<JetElementList>::BranchElement(const char *name, const enum EColor color, Float_t maxPt);
  template <>
  void BranchElement<JetElementList>::Reset();
  // template <>
  // void BranchElement<TEveElementList>::ReadBranch();
  // template <>
  // std::vector<TLorentzVector> BranchElement<TEveElementList>::GetVectors();


  // special case for track lists
  template <>
  BranchElement<TEveTrackList>::BranchElement(const char *name, const enum EColor color, Float_t maxPt);
  template <>
  void BranchElement<TEveTrackList>::SetTrackingVolume(Float_t r, Float_t l, Float_t Bz);
  template <>
  void BranchElement<TEveTrackList>::Reset();
  // template <>
  // void BranchElement<TEveTrackList>::ReadBranch();
  // template <>
  // std::vector<TLorentzVector> BranchElement<TEveTrackList>::GetVectors();

  #endif // CINT, CLING

}


#endif //BranchElement_h
