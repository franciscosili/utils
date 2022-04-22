// MultiDraw.cxx (code from pwaller)
// Draws many histograms in one loop over a tree.
// A little bit like a TTree::Draw which can make many histograms

#include <TTree.h>
#include <TH1D.h>
#include <TTreeFormula.h>

#include <iostream>

// Get an Element from an array
#define EL( type, array, index ) dynamic_cast<type *>( array->At( index ) )

void MultiDraw(TTree *tree, TObjArray *formulae, TObjArray *weights, TObjArray *hists, UInt_t list_len)
{
  Long64_t i = 0;
  Long64_t num_events = tree->GetEntries();

  Double_t value = 0, weight = 0, common_weight = 0;

  Int_t tree_number = -1;

  for (i = 0; i<num_events; i++) {

    if (tree_number != tree->GetTreeNumber()) {
      tree_number = tree->GetTreeNumber();
    }

    tree->LoadTree(tree->GetEntryNumber(i));

    for (UInt_t j=0; j<list_len; j++) {
      // If the Value or the Weight is the same as the previous, then it can be re-used.
      // In which case, this element fails to dynamic_cast to a formula, and evaluates to NULL
      if ( EL(TTreeFormula, formulae, j) )
        value = EL(TTreeFormula, formulae, j)->EvalInstance();

      if ( EL(TTreeFormula, weights, j) )
        weight = EL(TTreeFormula, weights, j)->EvalInstance();

      if (weight)
        EL(TH1D, hists, j)->Fill(value, weight);
    }
  }
}
