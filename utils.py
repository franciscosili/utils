# (py)root utils

import os
import sys
import ROOT
import math
from array import array


#-----------
# Utils
#-----------
#===================================================================================================
def mkdirp(path):
    if not os.path.exists(path): os.makedirs(path)
    return
#===================================================================================================

#===================================================================================================
def rmdir(path):
    import shutil
    try:
        shutil.rmtree(path)
    except OSError as exc:
        pass
    return
#===================================================================================================

#===================================================================================================
def run_cmd(cmd, logfile=None, stdout=False):

    print(cmd)

    if stdout and logfile is not None:
        cmd = '(set -o pipefail ; %s | tee %s)' % (cmd, logfile)
        
    elif logfile is not None:
        cmd = '%s > %s 2>&1' % (cmd, logfile)

    sc = os.system(cmd)

    # if sc != 0:
    #     print 'command sc != 0. exiting...'
    #     sys.exit(1)
            
    return sc
#===================================================================================================


#===================================================================================================
class Value(object):
    #===============================================================================================
    def __init__(self, mean=0.0, error=0.0):
        self.mean = mean
        self.error = error
        return
    #===============================================================================================

    #===============================================================================================
    def __repr__(self):
        if self.mean < 0.01:
            return '{:.4f} +- {:.4f}'.format(self.mean, self.error)
        else:
            return '{:.2f} +- {:.2f}'.format(self.mean, self.error)
    #===============================================================================================

    #===============================================================================================
    def __gt__(self, other):
        try:
            return self.mean > other.mean
        except:
            return self.mean > other
    #===============================================================================================

    #===============================================================================================
    def __lt__(self, other):
        try:
            return self.mean < other.mean
        except:
            return self.mean < other
    #===============================================================================================

    #===============================================================================================
    def __ge__(self, other):
        try:
            return self.mean > other.mean
        except:
            return self.mean > other
    #===============================================================================================

    #===============================================================================================
    def __le__(self, other):
        try:
            return self.mean < other.mean
        except:
            return self.mean < other
    #===============================================================================================

    #===============================================================================================
    def __add__(self, other):
        mean = self.mean + other.mean
        error = math.sqrt(self.error**2 + other.error**2)
        return Value(mean, error)
    #===============================================================================================

    #===============================================================================================
    def __sub__(self, other):
        mean = self.mean - other.mean
        error = math.sqrt(self.error**2 + other.error**2)
        return Value(mean, error)
    #===============================================================================================

    #===============================================================================================
    def __mul__(self, other):
        try:
            mean = self.mean * other.mean
            try:
                error = mean * math.sqrt((self.error/self.mean)**2 + (other.error/other.mean)**2)
            except ZeroDivisionError:
                error = 0
        except AttributeError:
            mean = self.mean * other
            error = self.error * other

        return Value(mean, error)
    #===============================================================================================
    
    #===============================================================================================
    def __rmul__(self, other):
        try:
            mean = self.mean * other.mean
            try:
                error = mean * math.sqrt((self.error/self.mean)**2 + (other.error/other.mean)**2)
            except ZeroDivisionError:
                error = 0
        except AttributeError:
            mean = self.mean * other
            error = self.error * other

        return Value(mean, error)
    #===============================================================================================

    #===============================================================================================
    def __truediv__(self, other):
        try:
            mean = self.mean / other.mean
        except ZeroDivisionError:
            mean = 0
        try:
            error = mean * math.sqrt((self.error/self.mean)**2 + (other.error/other.mean)**2)
        except ZeroDivisionError:
            error = 0

        return Value(mean, error)
    #===============================================================================================
#===================================================================================================


#===================================================================================================
def sci_notation(number, sig_fig=2, latex=False):
    ret_string = "{0:.{1:d}e}".format(number, sig_fig)
    a, b = ret_string.split("e")
    # remove leading "+" and strip leading zeros
    b = int(b)
    if latex:
        return f'{a}\\times 10^{{{str(b)}}}'
    else:
        return f'{a} * 10^{str(b)}'
        
#===================================================================================================


#--------
# Graphs
#--------
#===================================================================================================
def sort_graph(g, sort_x=True):
    """Creates a new graph with sorted values in x if sort_x=True or in y

    Args:
        g (TGraph): graph to sort
        sort_x (bool, optional): whether to sort values in x or not. Defaults to True.

    Returns:
        TGraph: sorted TGraph
    """
    # create empty arrays of float values
    ax = array('f', [])
    ay = array('f', [])

    d = dict()
    # loop over points in the graph g
    for i in range(g.GetN()):
        # create c++ x and y values
        xtmp = ROOT.Double(0)
        ytmp = ROOT.Double(0)
        # get x and y values for point i
        g.GetPoint(i, xtmp, ytmp)
        d[xtmp] = ytmp
    # sort values in increasing order of x
    if sort_x:
        for x, y in sorted(d.items()):
            ax.append(x)
            ay.append(y)
    else:
        for x, y in sorted(d, key=d.get):
            ax.append(x)
            ay.append(y)

    return ROOT.TGraph(g.GetN(), ax, ay)
#===================================================================================================

#===================================================================================================
class bcolors:
    # HEADER    = '\033[95m'
    CYAN      = '\033[96m'
    BLUE      = '\033[94m'
    YELLOW    = '\033[93m'
    GREEN     = '\033[92m'
    RED       = '\033[91m'
    ENDC      = '\033[0m'
    BOLD      = '\033[1m'
    UNDERLINE = '\033[4m'
#===================================================================================================
