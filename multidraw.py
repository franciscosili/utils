"""
Efficiently draw multiple histograms with one loop over all events in a TTree
This script injects a MultiDraw method into TTree when it is imported.
"""

import os
import ROOT

def MakeTObjArray(the_list):
    """Turn a python iterable into a ROOT TObjArray"""

    result = ROOT.TObjArray()
    result.SetOwner()

    # Make PyROOT give up ownership of the things that are being placed in the
    # TObjArary. They get deleted because of result.SetOwner()
    for item in the_list:
        ROOT.SetOwnership(item, False)
        result.Add(item)

    return result

def MultiDraw(tree, *draw_list):
    """
    Draws (projects) many histograms in one loop over a tree.

        Instead of:
        tree.Project("hname1", "ph_pt",  "weightA")
        tree.Project("hname2", "met_et", "weightB")

        Do:
        tree.MultiDraw( ("hname1", "ph_pt", "weightA" ),
                        ("hname2", "met_et", "weightB" ) )
    """

    histograms, variables, selections = [], [], []

    last_variable, last_selection = None, None

    for i, drawexp in enumerate(draw_list):

        # Expand out origFormula and weight, otherwise just use weight of 1.
        hname, variable, selection = drawexp

        hist = ROOT.gDirectory.Get(hname)
        if not hist:
            raise RuntimeError( "MultiDraw: Couldn't find histogram to fill '%s' in current directory." % hname )

        histograms.append(hist)

        # The following two 'if' clauses check that the next formula is different
        # to the previous one. If it is not, we add an ordinary TObject.
        # Then, the dynamic cast in MultiDraw.cxx fails, giving 'NULL', and
        # The previous value is used. This saves the recomputing of identical values
        if variable != last_variable:
            f = ROOT.TTreeFormula("variable%i" % i, variable, tree)

            if not f.GetTree():
                raise RuntimeError("TTreeFormula didn't compile: " + variable)

            f.SetQuickLoad(True)
            variables.append(f)
        else:
            variables.append(ROOT.TObject())

        if selection != last_selection:
            f = ROOT.TTreeFormula("selection%i" % i, selection, tree)

            if not f.GetTree():
                raise RuntimeError("TTreeFormula didn't compile: " + selection)

            f.SetQuickLoad(True)
            selections.append(f)
        else:
            selections.append(ROOT.TObject())

        last_variable, last_selection = variable, selection


    # Only compile MultiDraw once
    try:
        from ROOT import MultiDraw as _MultiDraw
    except ImportError:
        multi_draw_path = os.path.dirname(os.path.realpath(__file__))
        ROOT.gInterpreter.Declare(open(f'{multi_draw_path}/MultiDraw.cxx').read())
        from ROOT import MultiDraw as _MultiDraw

    # Ensure that formulae are told when tree changes
    fManager = ROOT.TTreeFormulaManager()

    for variable in variables + selections:
        if type(variable) == ROOT.TTreeFormula:
            fManager.Add(variable)

    fManager.Sync()
    tree.SetNotify(fManager)

    # Draw everything!
    variables  = MakeTObjArray(variables)
    selections = MakeTObjArray(selections)
    histograms = MakeTObjArray(histograms)

    _MultiDraw(tree,
                variables,
                selections,
                histograms,
                len(variables))

    variables.Delete()
    selections.Delete()
    del fManager

    return