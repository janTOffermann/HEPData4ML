// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G__VectorCalcs
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "vectorcalcs/VectorCalcs.h"

// Header files passed via #pragma extra_include

namespace VectorCalcs {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *VectorCalcs_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("VectorCalcs", 0 /*version*/, "vectorcalcs/VectorCalcs.h", 15,
                     ::ROOT::Internal::DefineBehavior((void*)0,(void*)0),
                     &VectorCalcs_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_DICT_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_DICT_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *VectorCalcs_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}

namespace {
  void TriggerDictionaryInitialization_libVectorCalcs_Impl() {
    static const char* headers[] = {
"vectorcalcs/VectorCalcs.h",
0
    };
    static const char* includePaths[] = {
"/local/home/xiaoyang1/HEPData4ML-master/util/root/vectorcalcs/inc",
"/local/home/xiaoyang1/anaconda3/envs/lgn_data/include",
"/local/home/xiaoyang1/HEPData4ML-master/util/root/vectorcalcs",
"/local/home/xiaoyang1/anaconda3/envs/lgn_data/include",
"/local/home/xiaoyang1/HEPData4ML-master/util/root/vectorcalcs/inc",
"/local/home/xiaoyang1/HEPData4ML-master/util/root/vectorcalcs",
"/local/home/xiaoyang1/anaconda3/envs/lgn_data/include",
"/local/home/xiaoyang1/anaconda3/envs/lgn_data/include",
"/local/home/xiaoyang1/anaconda3/envs/lgn_data/include",
"/local/home/xiaoyang1/anaconda3/envs/lgn_data/include",
"/local/home/xiaoyang1/anaconda3/envs/lgn_data/include",
"/local/home/xiaoyang1/anaconda3/envs/lgn_data/include",
"/local/home/xiaoyang1/anaconda3/envs/lgn_data/include",
"/local/home/xiaoyang1/anaconda3/envs/lgn_data/include/",
"/local/home/xiaoyang1/HEPData4ML-master/util/root/vectorcalcs/build/vectorcalcs/",
0
    };
    static const char* fwdDeclCode = nullptr;
    static const char* payloadCode = nullptr;
    static const char* classesHeaders[] = {
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libVectorCalcs",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libVectorCalcs_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libVectorCalcs_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libVectorCalcs() {
  TriggerDictionaryInitialization_libVectorCalcs_Impl();
}
