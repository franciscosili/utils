import ROOT
from array import array
import math
from functools import cmp_to_key


#===================================================================================================
def histogram(name, nx=None, xmin=None, xmax=None, xbins=None):
    """Function to create a 1D histogram

    Args:
        name (str): name of histograms
        nx (int, optional): Number of bins. Defaults to None.
        xmin (float, optional): Minimum value. Defaults to None.
        xmax (float, optional): Maximum value. Defaults to None.
        xbins (list, optional): Bins limits. Defaults to None.

    Returns:
        TH1F: histogram
    """

    if xbins:
        hist = ROOT.TH1F(name, name, len(xbins)-1, array('d', xbins))
    elif nx is not None and xmin is not None and xmax is not None:
        hist = ROOT.TH1F(name, name, nx, xmin, xmax)

    # To not be owned by the current directory
    # hist.SetDirectory(0)
    # ROOT.SetOwnership(hist, False)

    # Default configuration
    hist.Sumw2()
    hist.SetStats(0)
    hist.SetTitle('')

    return hist
#===================================================================================================

#===================================================================================================
def histogram2d(name, nx=None, xmin=None, xmax=None, ny=None, ymin=None, ymax=None, xbins=None, ybins=None):
    """Create 2D histogram

    Args:
        name (str): name of histogram
        nx (int, optional): Number of bins in the x direction. Defaults to None.
        xmin (float, optional): Minimum value in the x direction. Defaults to None.
        xmax (float, optional): Maximum value in the x direction. Defaults to None.
        ny (int, optional): Number of bins in the y direction. Defaults to None.
        ymin (float, optional): Minimum value in the y direction. Defaults to None.
        ymax (float, optional): Maximum value in the y direction. Defaults to None.
        xbins (list, optional): Bins limits in the x direction. Defaults to None.
        ybins (list, optional): Bins limits in the y direction. Defaults to None.

    Returns:
        TH2F: histogram
    """
    if xbins is not None and ybins is not None:
        hist = ROOT.TH2F(name, name, len(xbins)-1, array('d', xbins), len(ybins)-1, array('d', ybins))
    elif nx is not None and ny is not None:
        hist = ROOT.TH2F(name, name, nx, xmin, xmax, ny, ymin, ymax)

    # hist.SetDirectory(0)

    return hist
#===================================================================================================

#===================================================================================================
def histogram_equal_to(hist, name=None):
    """Create a copy of another histogram

    Args:
        hist (TH*): original histogram to copy
        name (str, optional): name of the new histogram. Defaults to None.

    Returns:
        TH*: clone of the original histogram without contents and errors
    """    
    if name is None:
        name = hist.GetName()
    newhist = hist.Clone(name)
    newhist.Reset()
    return newhist
#===================================================================================================

#===================================================================================================
def histogram_normalize(hist):
    """Normalize histogram

    Args:
        hist (TH*): histogram to normalize
    """
    area = hist.Integral()
    if area > 0:
        hist.Scale(1/area)
    return
#===================================================================================================

#===================================================================================================
def histogram_normalize_to(hist, other, xmin=None, xmax=None):
    """Normalize histogram to another ones area

    Args:
        hist (TH*): histogram to normalize
        other (TH*): other histogram to normalize to
        xmin (float, optional): Minimum limit of the integral. Defaults to None.
        xmax (float, optional): Maximum limit of the integral. Defaults to None.

    Returns:
        float: normalization factor for which the original hist has been normalized
    """    
    if xmin and xmax:
        n1 = hist.Integral(hist.FindBin(xmin), hist.FindBin(xmax))
        n2 = other.Integral(other.FindBin(xmin), other.FindBin(xmax))
    else:
        n1 = hist.Integral()
        n2 = other.Integral()
    s = n2/n1 if n1 > 0.0 else 1.0
    hist.Scale(s)
    return s
#===================================================================================================

#===================================================================================================
def histogram_add_overflow_bin(hist):
    """Add overflow bin content to the last bin

    Args:
        hist (TH*): histogram.
    """    

    # 2D histograms
    if hist.InheritsFrom('TH2'):
        # last bins
        last_bin_x = hist.GetNbinsX()
        last_bin_y = hist.GetNbinsY()
        # overflow bin number
        over_bin_x = last_bin_x + 1
        over_bin_y = last_bin_y + 1

        # loop over bins in x which are the last bins in y at the same time
        for bx in range(1, last_bin_x):
            # add content of the overflow bin in y to the current bin
            new_val = hist.GetBinContent(bx, last_bin_y) + hist.GetBinContent(bx, over_bin_y)
            # set contents
            hist.SetBinContent(bx, last_bin_y, new_val)
            hist.SetBinContent(bx, over_bin_y, 0.0)
            # set errors
            e1 = hist.GetBinError(bx, last_bin_y)
            e2 = hist.GetBinError(bx, over_bin_y)
            new_err = math.sqrt(e1*e1 + e2*e2)
            hist.SetBinError(bx, last_bin_y, new_err)
            hist.SetBinError(bx, over_bin_y, 0.0)

        # loop over bins in y which are the last bins in x at the same time
        for by in range(1, last_bin_y):
            new_val = hist.GetBinContent(last_bin_x, by) + hist.GetBinContent(over_bin_x, by)
            hist.SetBinContent(last_bin_x, by, new_val)
            hist.SetBinContent(over_bin_x, by, 0.0)

            e1 = hist.GetBinError(last_bin_x, by)
            e2 = hist.GetBinError(over_bin_x, by)
            new_err = math.sqrt(e1*e1 + e2*e2)
            hist.SetBinError(last_bin_x, by, new_err)
            hist.SetBinError(over_bin_x, by, 0.0)

        # last x/y bin. Add contents of all four bins
        new_val = hist.GetBinContent(last_bin_x, last_bin_y) + \
                  hist.GetBinContent(over_bin_x, last_bin_y) + \
                  hist.GetBinContent(last_bin_x, over_bin_y) + \
                  hist.GetBinContent(over_bin_x, over_bin_y)

        hist.SetBinContent(last_bin_x, last_bin_y, new_val)
        hist.SetBinContent(last_bin_x, over_bin_y, 0.)
        hist.SetBinContent(over_bin_x, last_bin_y, 0.)
        hist.SetBinContent(over_bin_x, over_bin_y, 0.)

        e1 = hist.GetBinError(last_bin_x, last_bin_y)
        e2 = hist.GetBinError(over_bin_x, last_bin_y)
        e3 = hist.GetBinError(last_bin_x, over_bin_y)
        e4 = hist.GetBinError(over_bin_x, over_bin_y)

        new_err = math.sqrt(e1*e1+e2*e2+e3*e3+e4*e4)
        hist.SetBinError(last_bin_x, last_bin_y, new_err)
        hist.SetBinError(last_bin_x, over_bin_y, 0.)
        hist.SetBinError(over_bin_x, last_bin_y, 0.)
        hist.SetBinError(over_bin_x, over_bin_y, 0.)

    # 1D histogram
    else:
        last_bin = hist.GetNbinsX()
        over_bin = last_bin + 1

        # value
        new_val = hist.GetBinContent(last_bin) + hist.GetBinContent(over_bin)
        hist.SetBinContent(last_bin, new_val)
        hist.SetBinContent(over_bin, 0.0)

        # error
        e1 = hist.GetBinError(last_bin)
        e2 = hist.GetBinError(over_bin)
        new_err = math.sqrt(e1*e1 + e2*e2)
        hist.SetBinError(last_bin, new_err)
        hist.SetBinError(over_bin, 0.0)
    return
#===================================================================================================

#===================================================================================================
def histogram_scale(hist, c, e_c=None):
    """Scale histogram by a factor with error (c +- e_c)
    If e_c is None and c is not a Value(), it does the same as TH1.Scale()

    Args:
        hist (TH1F): histogram
        c (float): Scale factor
        c (Value()): Value() Object.
        e_c (float, optional): error on the scale factor. Defaults to None.
    """    
    # in case we are working with objects from class Value
    try:
        c, e_c = c.mean, c.error
    except AttributeError:
        pass

    # if not error has been passed, simply scale the histogram
    if e_c is None:
        hist.Scale(c)
        return
    # In case we have the scale error, compute the bin contents and errors with appropiate error
    # propagation
    for b in range(1, hist.GetNbinsX()+1):
        # get bin content and errors
        n_b = hist.GetBinContent(b)
        e_b = hist.GetBinError(b)
        # new bin content
        new_n = n_b * c
        # calculate error
        try:
            err2 = (e_b/n_b)**2 + (e_c/c)**2
            new_e = new_n * math.sqrt(err2)
        except ZeroDivisionError:
            new_e = 0
        hist.SetBinContent(b, new_n)
        hist.SetBinError  (b, new_e)
    return
#===================================================================================================

#===================================================================================================
def get_cumulative_histogram(hist, inverse_x=False, inverse_y=False):
    """Calculate the cumulative histogram. The option inverse_x/y dictates if the cumulative integral
    is calculated from the beggining of the histogram to the current bin (if true), or from the
    current bin to the final bin of the histogram (if false)

    Args:
        hist (TH2F): histogram in 2D
        inverse_x (bool, optional): Inverse cumulative calculation in the x direction. Defaults to False.
        inverse_y (bool, optional): Inverse cumulative calculation in the y direction. Defaults to False.

    Returns:
        TH2F: Cumulative 2D histogram
    """
    newhist = hist.Clone(hist.GetName())
    nx = hist.GetNbinsX()
    ny = hist.GetNbinsY()
    for bx in range(nx):
        for by in range(ny):
            if inverse_x and inverse_y:
                # integrate from the beggining to the current bin in x and y
                cum = hist.Integral(1, bx+1, 1, by+1)
            elif inverse_x:
                # integrate from the beggining to the current bin in x, and from the current bin to
                # the final one in y
                cum = hist.Integral(1, bx+1, by+1, ny)
            elif inverse_y:
                # integrate from the beggining to the current bin in y, and from the current bin to
                # the final one in x
                cum = hist.Integral(bx+1, nx, 1, by+1)
            else:
                # integrate from the current bin to the the final one in x and y
                cum = hist.Integral(bx+1, nx, by+1, ny)

            newhist.SetBinContent(bx+1, by+1, cum)

    return newhist
#===================================================================================================

#===================================================================================================
def merge_histograms(merge_name, merge_list):
    """Add multiple histograms together.

    Args:
        merge_name (str): Name of the merged histogram
        merge_list (list): list of histograms to be added together

    Returns:
        TH*: merged/added histogram
    """
    new_hist = merge_list[0].Clone(merge_name)
    for h in merge_list[1:]:
        new_hist.Add(h, 1)

    return new_hist
#===================================================================================================

#===================================================================================================
def get_histogram(filename, treename, variable, selection='', xmin=None, xmax=None, bins=None, hist=None):
    """Create a histogram of the variable 'variable' by projecting the tree 'treename' stored in 'filename'.
    The selection 'selection' is applied, if specified.  In case the histogram is specified as input,
    the tree is projected onto that histogram, if not, a new TH1 object is created.

    Args:
        filename (str): file where the tree is stored
        treename (str): name of tree
        variable (str): variable to project
        selection (str, optional): Selection to apply. Defaults to ''.
        xmin (float/int, optional): Minimum limit of new histogram. Defaults to None.
        xmax (float/int, optional): Maximum limit of new histogram. Defaults to None.
        bins (int, optional): Number of bins of the new histogram. Defaults to None.
        hist (TH1F, optional): Histogram to project onto. Defaults to None.

    Returns:
        TH1F: Histogram with projected tree
    """

    t = ROOT.TChain(treename)
    t.Add(filename)

    if hist is None:
        if xmin is not None and xmax is not None and bins is not None:
            hist = ROOT.TH1F('htemp', 'htemp', bins, xmin, xmax)
        else:
            hist = ROOT.TH1F('htemp', 'htemp', 100, 0, 100)

        t.Project('htemp', variable, selection, 'goff')

    else:
        t.Project(hist.GetName(), variable, selection, 'goff')

    return hist.Clone()
#===================================================================================================

#===================================================================================================
def get_stack(hists):
    stack = ROOT.THStack()
    def _compare(a, b):
        amax = a.GetMaximum()
        bmax = b.GetMaximum()
        return compare(int(amax), int(bmax))

    try:
        for hist in sorted(hists.values(), key=cmp_to_key(_compare)):
            stack.Add(hist)
    except AttributeError:
        try:
            for hist in sorted(hists, key=cmp_to_key(_compare)):
                stack.Add(hist)
        except:
            print('Input is not a list nor a dictionary.')
            pass
        
            
    return stack
#===================================================================================================

#===================================================================================================
def compare(a, b):
    return (a > b) - (a < b)
#===================================================================================================

#===================================================================================================
def get_histogram_file(fname, hname):
    
    _infile = ROOT.TFile(fname, 'read')    
    
    print(f'Getting histograms from file {fname}')
    print(f'Histogram name: {hname}')
    
    hist = _infile.Get(hname)
    hist.SetDirectory(0)
    _infile.Close()
    
    return hist
#===================================================================================================

#===================================================================================================
def get_ratio_hist(hnum, hden, hname='ratio_hist'):
    
    ratio = hnum.Clone(hname)
    ratio.Divide(hden)
    
    return ratio
#===================================================================================================