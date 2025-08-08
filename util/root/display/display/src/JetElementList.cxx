#include <iostream>

#include "display/JetElementList.h"
#include "display/JetCone.h"

namespace Display{

  Bool_t JetElementList::SetRnrChildren(Bool_t rnr){
    Bool_t result = TEveElementList::SetRnrChildren(rnr);

    // Propagate to all JetCone shapes
    for (TEveElement::List_i it = BeginChildren(); it != EndChildren(); ++it) {
      JetCone* jet = dynamic_cast<JetCone*>(*it);
      if (jet && jet->GetJetCircle()) jet->GetJetCircle()->SetRnrSelf(rnr);
    }

    // Update the calo scene - you'd need to store a reference to it
    if (fCaloScene) fCaloScene->Changed();

    return result;
  }

  // Bool_t JetElementList::SetRnrSelf(Bool_t rnr){
  //   std::cout << "JetElementList::SetRnrSelf" << std::endl;
  //   Bool_t result = TEveElementList::SetRnrSelf(rnr);
  //   SyncShapes(rnr && GetRnrChildren()); // Only show if both self and children are enabled
  //   return result;
  // }

  void JetElementList::SyncShapes(Bool_t visible){
    // Propagate to all JetCone shapes
    for (TEveElement::List_i it = BeginChildren(); it != EndChildren(); ++it) {
      JetCone* jet = dynamic_cast<JetCone*>(*it);
      if (jet && jet->GetJetCircle()) {
          jet->GetJetCircle()->SetRnrSelf(visible);
      }
    }
  }
}