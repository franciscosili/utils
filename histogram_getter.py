import os
import re
import glob
from array import array
from tqdm import tqdm

import ROOT as RT
from .datasets_manager import get_datasets
from .binningutils import *
from .histutils import *
from . import xsutils
from .multidraw import MultiDraw
from .utils import Value


# Luminosity
lumi_dict = {
    '2015':  3219.56,
    '2016': 32965.30,
    '2017': 44307.40,
    '2018': 58450.10,
    '2022': 0.0,
}


#===================================================================================================
#===================================================================================================
#=====================================      HISTOGRAMMING
#===================================================================================================
#===================================================================================================
class histogram_getter:
    #===============================================================================================
    def __init__(self, ana_rel, year,
                 paths, samples,
                 binning, regions=None, selections=None,
                 version=None, wg_label='',
                 weights_strings=None, tree_name=None,
                 add_overflow_bin=False, scale=True, lumi=None,
                 use_skim=False, use_lumiw=True, use_mcw=True, use_sfw=True, use_purw=True,
                 use_mcveto=True, use_phjet_w=False, truth_mc='',
                 slices=False,
                 remove_var_cut=False,
                 seed_selection='(seed_met_et - 8)/sqrt(seed_sumet) < (1.0 + 0.1*seed_bjet_n)',
                 preselection=None,
                 ignore_missing=None,
                 systs_module=None,
                 debug=False):
        # required parameters
        self.ana_rel          = ana_rel
        self.version          = version
        self.year             = year
        # samples and paths info
        
        self.wg_label         = wg_label
        # custom weights strings
        self.weights_strings  = weights_strings
        # optional extra parameters
        self.add_overflow_bin = add_overflow_bin
        self.scale            = scale
        self.lumi             = lumi
        self.use_skim         = use_skim
        self.use_lumiw        = use_lumiw
        self.use_sfw          = use_sfw
        self.use_mcw          = use_mcw
        self.use_purw         = use_purw
        self.use_mcveto       = use_mcveto
        self.use_phjet_w      = use_phjet_w
        self.truth_mc         = truth_mc
        self.slices           = slices
        self.remove_var_cut   = remove_var_cut
        self.seed_selection   = seed_selection
        self.ignore_missing   = ignore_missing
        self.debug            = debug

        # tree name
        if tree_name is not None:
            self.tree_name = tree_name
        else:
            self.tree_name = 'smini' if self.use_skim else 'mini'        

        # dataset year
        self.dataset_year = year
        if self.dataset_year in ('2015', '2016', '2015+2016'):
            self.run = 2
            self.mc_campaign = 'mc16a' if self.ana_rel < 22 else 'mc20a'
        elif self.dataset_year == '2017':
            self.run = 2
            self.mc_campaign = 'mc16d' if self.ana_rel < 22 else 'mc20d'
        elif self.dataset_year == '2018':
            self.run = 2
            self.mc_campaign = 'mc16e' if self.ana_rel < 22 else 'mc20e'
        elif self.dataset_year == '2022':
            self.run = 3
            self.mc_campaign = 'mc21a'

        print('The configuration for histogram_getter is:')
        for var_name, var_val in vars(self).items():
            print(f'  {var_name:17}= {var_val}')

        self.paths            = paths
        self.samples          = samples
        self.regions          = regions
        self.selections       = selections
        self.preselection     = preselection
        self.binning          = binning


        if systs_module:
            import systs_module as systematics_
            self.systematics = systematics_

        # RT.gInterpreter.Declare(open(os.environ['snip_path'] + '/variables.cxx').read())

        return
    #===============================================================================================

    #===============================================================================================
    def _get_multi_histograms(self, name, path, is_mc, lumi, regions, selections, variables,
                              systematics=['Nom', ], binning=None, dsid_str=None):
        """
        get histogram for a given sample (ds) and year
        """
        is_fake = ('efake' in name or 'jfake' in name)
        is_phjet = ('photonjet' in name)
        is_smr = ('smr' in name)

        variables  = variables if isinstance(variables, list) else [variables, ]
        regions    = regions if isinstance(regions, list) else [regions, ]
        selections = selections if isinstance(selections, list) else [selections, ]

        # ------------
        # File/Chain
        # ------------

        # open tree or chain, depending if the path is a directory or if it is a single file
        if os.path.isdir(path):
            tree = RT.TChain(self.tree_name)

            all_files = glob.glob(path+'/*root*')
            for fpath in all_files:
                tree.Add(fpath)
        else:
            file_ = RT.TFile(path, 'read')
            tree = file_.Get(self.tree_name)

        # Lumi weight is the same for all histograms
        if is_mc and self.use_lumiw and not self.weights_strings:
            if self.use_skim:
                lumi_weight = 'weight_xs*%.2f' % lumi
            elif dsid_str:
                lumi_weight = '%s' % get_lumi_weight(path, int(dsid_str), lumi)
            else:
                lumi_weight = ''

        # ----------------------------------------------
        # Create histograms and "tuples" for MultiDraw
        # ----------------------------------------------
        draw_list  = []
        histograms = []

        # get a list of selections or regions
        if regions and selections:
            # check if same size ...
            pass
        elif regions and not selections:
            selections = [self.regions[reg] for reg in regions]
        elif selections and not regions:
            regions    = ['R%i' % isel for isel in range(len(selections))]

        # loop over both regions and selections simultaneously
        for region, selection in zip(regions, selections):
            # in each region and selection, loop over the variables
            for ivariable, variable in enumerate(variables):
                # loop for each systematic
                for syst in systematics:

                    _selection = selection
                    systname   = syst.replace('__1down', 'Low').replace('__1up', 'High')

                    # ------- SETUP HISTOGRAMS
                    # fix_events_by_interval = False
                    if is_2d_variable(variable):
                        varx, vary = variable.split(':')

                    # get binning
                    if binning is None or not binning:
                        _binning = get_binning(variable, self.binning)
                    else:
                        _binning = binning[ivariable]

                    # name to avoid the ROOT warning, not used
                    if self.use_skim:
                        hname = f'h___{name}___{systname}__{region}__{get_escaped_variable(variable)}'
                    elif dsid_str:
                        hname = f'h___{dsid_str}___{systname}__{region}__{get_escaped_variable(variable)}'
                    else:
                        hname = f'h___{name}___{systname}__{region}__{get_escaped_variable(variable)}'

                    # in case we give a 2D variable name, create a 2D histogram
                    if is_2d_variable(variable):
                        htemp = RT.TH2D(hname, hname, *_binning)
                        htemp.Sumw2()
                    else:
                        if len(_binning) > 3:
                            htemp = RT.TH1D(hname, '', len(_binning)-1, array('d', _binning))
                        else:
                            htemp = RT.TH1D(hname, '', int(_binning[0]), _binning[1], _binning[2])
                        htemp.Sumw2()

                    # ------- SETUP SELECTIONS
                    # MC veto
                    if is_mc and self.use_mcveto:
                        if self.wg_label is None:
                            if _selection:
                                _selection = '%s && mcveto==0' % _selection
                            else:
                                _selection = 'mcveto==0'

                        if self.truth_mc!='':
                            _selection = f'{_selection} && {self.truth_mc}'

                        # skim pt slices
                        if not self.use_skim and dsid_str:
                            dsid = int(dsid_str)
                            # sherpa
                            if dsid in (361042, 361043, 361044, 364543):
                                _selection += f'&& ph_truth_pt[0]>70. && ph_truth_pt[0]<140.'
                            elif dsid in (361045, 361046, 361047, 364544):
                                _selection += f'&& ph_truth_pt[0]>140. && ph_truth_pt[0]<280.'
                            elif dsid in (361048, 361049, 361050, 364545):
                                _selection += f'&& ph_truth_pt[0]>280. && ph_truth_pt[0]<500.'
                            elif dsid in (361051, 361052, 361053, 364546):
                                _selection += f'&& ph_truth_pt[0]>500. && ph_truth_pt[0]<1000.'
                            elif dsid in (361054, 361055, 361056, 361057, 361058, 361059, 364547):
                                _selection += f'&& ph_truth_pt[0]>1000.'
                            # pythia
                            elif dsid in (800662, 800676):
                                _selection += f'&& ph_truth_pt[0]>70. && ph_truth_pt[0]<140.'
                            elif dsid in (800663, 800677):
                                _selection += f'&& ph_truth_pt[0]>140. && ph_truth_pt[0]<280.'
                            elif dsid in (800664, 800678):
                                _selection += f'&& ph_truth_pt[0]>280. && ph_truth_pt[0]<500.'
                            elif dsid in (800665, 800679):
                                _selection += f'&& ph_truth_pt[0]>500. && ph_truth_pt[0]<800.'
                            elif dsid in (800666, 800680):
                                _selection += f'&& ph_truth_pt[0]>800. && ph_truth_pt[0]<1000.'
                            elif dsid in (800667, 800681):
                                _selection += f'&& ph_truth_pt[0]>1000. && ph_truth_pt[0]<1500.'
                            elif dsid in (800668, 800682):
                                _selection += f'&& ph_truth_pt[0]>1500. && ph_truth_pt[0]<2000.'
                            elif dsid in (800683,):
                                _selection += f'&& ph_truth_pt[0]>2000.'
                            elif dsid in (800669,):
                                _selection += f'&& ph_truth_pt[0]>2000. && ph_truth_pt[0]<2500.'
                            elif dsid in (800670,):
                                _selection += f'&& ph_truth_pt[0]>2500. && ph_truth_pt[0]<3000.'
                            elif dsid in (800671,):
                                _selection += f'&& ph_truth_pt[0]>3000.'

                    if is_smr:
                        _selection += '&& smeared==1 && %s' % self.seed_selection

                    # Remove variable from selection if n-1
                    if self.remove_var_cut and variable in _selection and not variable == 'cuts':
                        _selection = '&&'.join([cut for cut in _selection.split('&&') if not split_cut(cut)[0] == variable])

                    # if do_remove_var and (':' in variable):
                    #     if varx in selection:
                    #         selection = '&&'.join([ cut for cut in selection.split('&&') if not split_cut(cut)[0] == varx ])
                    #     if vary in selection:
                    #         selection = '&&'.join([ cut for cut in selection.split('&&') if not split_cut(cut)[0] == vary ])

                    # change selection and variable for systematics
                    if syst != 'Nom' and self.systematics.affects_kinematics(syst):
                        for var in self.systematics.get_affected_variables(syst):
                            _selection = replace_var(_selection, var, '%s_%s' % (var, syst))

                        if variable in self.systematics.get_affected_variables(syst):
                            variable = '%s_%s' % (variable, syst)

                    # ------- SETUP WEIGHTS
                    w_list = []
                    if is_mc:
                        # lumi weight
                        if self.use_lumiw:
                            if not self.weights_strings and lumi_weight!='':
                                w_list.append('%s' % lumi_weight)
                            else:
                                if 'lumi_w' in self.weights_strings:
                                    lumi_weight = self.weights_strings['lumi_w']
                                    w_list.append(lumi_weight)

                        # mc weight
                        if self.use_mcw:
                            if self.weights_strings is not None:
                                if 'weight_mc' in self.weights_strings:
                                    w_list.append(self.weights_strings['weight_mc'])
                            else:
                                w_list.append('weight_mc')

                        # scale factors
                        if self.use_sfw:
                            if self.weights_strings is not None:
                                if 'weight_sf' in self.weights_strings:
                                    w_list.append(self.weights_strings['weight_sf'])
                            else:
                                if syst != 'Nom' and self.systematics.affects_weight(syst) and not 'PRW_DATASF' in syst:
                                    w_list.append('weight_sf_%s' % syst)
                                else:
                                    w_list.append('weight_sf')

                        # pile-up
                        if self.use_purw:
                            if self.weights_strings is not None:
                                if 'weight_pu' in self.weights_strings:
                                    w_list.append(self.weights_strings['weight_pu'])
                            else:
                                if 'PRW_DATASF__1down' == syst:
                                    w_list.append('weight_pu_dn')
                                elif 'PRW_DATASF__1up' == syst:
                                    w_list.append('weight_pu_up')
                                else:
                                    w_list.append('weight_pu')

                        # photonjet MET reweighting (not used for the moment)
                        if self.use_phjet_w and is_phjet:
                            w_list.append('weight_ff')

                        # SUSY EWK BR re-weighting
                        # if self.ggm_br:
                        #     w_list.append(get_GGM_model_weight(*ggm_br))
                        
                        # Photon ID and Isolation scale factor weights
                        
                        if 'RZ' in self.wg_label:
                            if 'id' in selection and 'weight_id' in self.weights_strings:
                                w_list.append(self.weights_strings['weight_id'])
                            if 'isoloose' in selection and 'weight_isoloose' in self.weights_strings:
                                w_list.append(self.weights_strings['weight_isoloose'])
                            elif 'isotight' in selection and 'weight_isotight' in self.weights_strings:
                                w_list.append(self.weights_strings['weight_isotight'])
                            elif 'isotightcaloonly' in selection and 'weight_isotightco' in self.weights_strings:
                                w_list.append(self.weights_strings['weight_isotightco'])
                                

                    elif is_fake:
                        if syst == 'Nom':
                            w_list.append('weight_ff')
                        elif syst == 'EFAKE_SYST__1down' or syst == 'JFAKE_SYST__1down':
                            w_list.append('weight_ff_dn')
                        elif syst == 'EFAKE_SYST__1up' or syst == 'JFAKE_SYST__1up':
                            w_list.append('weight_ff_up')

                    w_str = '*'.join(w_list) if self.scale else ''

                    # ------- SETUP DRAW LIST AND ADD HISTOGRAMS
                    varexp = ''
                    if _selection and w_str:
                        varexp = '(%s)*(%s)' % (_selection, w_str)
                    elif _selection:
                        varexp = _selection
                    elif self.scale:
                        varexp = w_str

                    histograms.append(htemp)

                    if variable == 'cuts':
                        draw_list.append((hname, '1', varexp))
                    else:
                        draw_list.append((hname, variable, varexp))


        # Use Draw or MutiDraw to project all histograms (for 2D histograms only 1 variable allowed)
        if len(variables) == 1 and is_2d_variable(variables[0]):
            hname, variable, selection = draw_list[0]
            if self.debug:
                print(hname, variable, selection)

            variable = '%s:%s' % (vary, varx)
            tree.Project(hname, variable, selection)
        elif len(draw_list) == 1:
            hname, variable, selection = draw_list[0]
            if self.debug:
                print(hname, variable, selection)
            tree.Project(hname, variable, selection)
        else:
            if self.debug:
                print(draw_list[0])
            MultiDraw(tree, *draw_list)
        for hist in histograms:
            hist.SetDirectory(0)

        if os.path.isdir(path):
            tree.Reset()
            del tree
        else:
            tree.Reset()
            del tree
            file_.Close()

        return histograms
    #===============================================================================================

    #===============================================================================================
    def get_histograms_from_skim(self, name, systematics, binning, is_mc, lumi, regions, selections, variables):
        path = get_skim_path(name, self.paths, self.version, self.mc_campaign, )
        return self._get_multi_histograms(name, path, is_mc, lumi, regions, selections, variables,
                                          systematics, binning)
    #===============================================================================================

    #===============================================================================================
    def get_histograms(self, sample, regions, selections, variables, systematics=None, binning=None,
                       dataset=None, extra_regex=None):

        # dataset: 2015, 2016, 2017, 2018, or combinations with "+"
        # 'Run2': 2015+2016+2017+2018
        dataset_year = dataset if dataset else self.year
        if not dataset:
            dataset_year = self.year
            if dataset_year == 'Run2':
                dataset_year = '2015+2016+2017+2018'

        is_mc = (not 'data' in sample and not 'efake' in sample and not 'jfake' in sample and not 'smr' in sample)

        if '+' in dataset_year and not (is_mc and dataset_year == '2015+2016'):
            # if self.year: del self.year

            dataset_years = dataset_year.split('+')
            if is_mc and '2015' in dataset_years and '2016' in dataset_years:
                dataset_years.remove('2015')
                dataset_years.remove('2016')
                dataset_years.insert(0, '2015+2016')

            _tmp_list = []
            for y in dataset_years:
                _tmp_list.append(
                    self.get_histograms(sample, regions, selections, variables, systematics,
                                        binning, y, extra_regex)
                )
            return sum_histograms(_tmp_list)

        # Data type/MC campaign
        if sample in ['data', 'efake', 'jfake', 'smr']:
            sample = sample + dataset_year[-2:]

        # lumi
        lumi = 0.
        if is_mc:
            lumi = get_lumi(dataset_year)

        if self.use_skim:
            skim_dict = {
                'vjets'     : ['wjets', 'zlljets', 'znunujets'],
                'zllgamma'  : ['zeegamma', 'zmumugamma', 'ztautaugamma'],
                'zgamma'    : ['zeegamma', 'zmumugamma', 'ztautaugamma', 'znunugamma'],
                'vgamma'    : ['zeegamma', 'zmumugamma', 'ztautaugamma', 'znunugamma', 'wgamma'],
                'gammagamma': ['vgammagamma', 'diphoton'],
            }

            if sample in skim_dict:
                _list_hists = [
                    self.get_histograms_from_skim(s, systematics, binning, is_mc, lumi, regions, selections, variables) for s in skim_dict.get(sample)
                ]
                histograms = sum_histograms(_list_hists)
            else:
                histograms = self.get_histograms_from_skim(sample, systematics, binning, is_mc, lumi,
                                                           regions, selections, variables)

        else:
            datasets = get_datasets(sample, self.paths, self.samples, self.version,
                                    self.ignore_missing, self.mc_campaign, extra_regex)

            histograms = []
            if self.slices and is_mc:
                histograms_slices = {}

            name_to_show = f'{sample}/{self.mc_campaign}' if is_mc else sample

            for ds in tqdm(datasets, desc=name_to_show):
                try:
                    dsid = ds['dsid']
                except KeyError:
                    dsid = None
                histograms_ds = self._get_multi_histograms(sample, ds['path'], is_mc, lumi, regions,
                                                           selections, variables, dsid_str=dsid)

                if self.slices and is_mc:
                    histograms_slices[ds['short_name']] = histograms_ds

                if not histograms:
                    for hist in histograms_ds:
                        histograms.append(hist.Clone())
                else:
                    for hall, hnew in zip(histograms, histograms_ds):
                        hall.Add(hnew, 1)

        # Fix histogram name and add oflow bin
        for hist in histograms:
            fix_histogram_name(hist, sample)
            if self.add_overflow_bin:
                histogram_add_overflow_bin(hist)

        if self.slices and is_mc:
            for slicename, hists in histograms_slices.items():
                for h in hists:
                    fix_histogram_name(h, sample, slicename)
                    if self.add_overflow_bin:
                        histogram_add_overflow_bin(h)
            return histograms, histograms_slices
        else:
            return histograms
    #===============================================================================================

    #===============================================================================================
    def get_histogram(self, sample, region=None, selection=None, variable=None, syst='Nom', binning=None, extra_regex=None):
        if binning: binning = [binning, ]

        if self.slices:
            histograms, histograms_slices = self.get_histograms(sample, region, selection, variable,
                                                                [syst, ], binning, extra_regex=extra_regex)
        else:
            histograms = self.get_histograms(sample, region, selection, variable, [syst, ], binning,
                                             extra_regex=extra_regex)

        if histograms and self.slices:
            return histograms[0], histograms_slices
        elif histograms:
            return histograms[0]
        return
    #===============================================================================================

    #===============================================================================================
    def get_events(self, sample):

        hist  = self.get_histogram(sample)
        mean  = hist.GetBinContent(1)
        error = hist.GetBinError(1)
        hist.Delete()

        return Value(mean, error)
    #===============================================================================================

    #===============================================================================================
    def get_cutflow(self, sample):
        if not self.selections:
            return None

        cuts = [split_cut(cut) for cut in split_selection(self.selections[0])]
        cutflow = histogram('cutflow', len(cuts)+1, 0.5, len(cuts)+1.5)

        cutflow.GetXaxis().SetBinLabel(1, 'No Cut')
        for i, cut in enumerate(cuts):
            cutflow.GetXaxis().SetBinLabel(i+2, cut[0])

        cuts = [('', '', ''), ] + cuts

        selection = ''
        for i, (var, op, value) in enumerate(cuts):
            if not selection:
                selection += '%s %s %s' % (var, op, value)
            else:
                selection += ' && %s %s %s' % (var, op, value)

            selection = selection.strip()
            if selection == ' ':
                selection = ''

            evts = self.get_events(sample)

            cutflow.SetBinContent(i+1, evts.mean)
            cutflow.SetBinError(i+1, evts.error)

        return cutflow
    #===============================================================================================


#===================================================================================================
#===================================================================================================
#=====================================      CUTS
#===================================================================================================
#===================================================================================================

#===================================================================================================
def num(s):
    """Get number from a string
    Args:
        s (str): number in string format
    Returns:
        float/int: number converted into int or float
    """
    try:
        return int(s)
    except ValueError:
        return float(s)
#===================================================================================================

#===================================================================================================
def split_selection(s):
    """Given a selection 's', split it

    Args:
        s (str): selection string

    Returns:
        list: list containing the different selection parts
    """
    selection = []
    # loop over each cut on the selection string
    for cut in s.split('&&'):
        # remove spaces at beggining and end of string
        cut = cut.strip()
        # in case we have 'sub-cuts'
        if cut.startswith('(') and cut.endswith(')'):
            # for incut in split_selection(cut[1:-1]):
            #     selection.append(incut)
            try:
                # for incut in cut.split('&&'):
                #     incut = incut.strip()
                #     selection.append(incut)
                cut1, cut2 = cut[1:-1].split('&&')
                selection.append(cut1)
                selection.append(cut2)
            except:
                try:
                    # for incut in cut.split('||'):
                    #     incut = incut.strip()
                    #     selection.append(incut)
                    cut1, cut2 = cut[1:-1].split('||')
                    selection.append(cut1)
                    selection.append(cut2)
                except:
                    raise
        else:
            selection.append(cut)

    return selection
#===================================================================================================

#===================================================================================================
def split_cut(s):
    """Get the variable to cut, the operator which is used and the value

    Args:
        s (str): The cut at which the components of it want to be known

    Returns:
        tuple: tuple containing the variable, the operator, and the value
    """
    s = s.strip()
    for op in ['==', '>=', '<=', '>', '<', '!=']:
        if op in s:
            var, value = s.split(op)
            return (var.strip(), op, value.strip())
#===================================================================================================

#===================================================================================================
def revert_cut(s):
    """Given a cut string, revert the operator

    Args:
        s (str): cut

    Returns:
        str: reverted cut
    """
    (var, op, value) = split_cut(s)

    if op == '>' or op == '>=':
        newop = '<'
    elif op == '<' or op == '<=':
        newop = '>'
    else:
        newop = op

    return '%s %s %s' % (var, newop, value)
#===================================================================================================

#===================================================================================================
def check_cut(value, op, cut):
    """Check if the cut is satisfied for the given value

    Args:
        value (number): value of the variable
        op (str): operator
        cut (number): cut value

    Returns:
        bool: true if the value satisfies the cut or false if it doesn't
    """
    if op == '==':
        return (value == cut)
    elif op == '>=':
        return (value >= cut)
    elif op == '<=':
        return (value <= cut)
    elif op == '>':
        return (value > cut)
    elif op == '<':
        return (value < cut)
    elif op == '!=':
        return (value != cut)
#===================================================================================================

#===================================================================================================
def replace_var(selection, oldvar, newvar):
    """Replace variable in a cut

    Args:
        selection (str): selection string containing multiple cuts
        oldvar (str): name of the old variable
        newvar (str): name of the new variable

    Returns:
        str: new selection with the variable oldvar replaced by newvar
    """
    new_cuts = []
    # loop over all the cuts
    for cut in split_selection(selection):
        # name of the variable
        svar = split_cut(cut)[0]

        # index cases
        if '[' in svar and ']' in svar and svar.split('[')[0] == oldvar:
            var, arg = svar[:svar.index('[')], svar[svar.index('['):]
            new_cuts.append(cut.replace(var, newvar))

        # '+' cases (only one +)
        elif '+' in svar:
            s1, s2 = svar.split('+')
            if s1 == oldvar:
                new_cuts.append(cut.replace(s1, newvar))
            elif s2 == oldvar:
                new_cuts.append(cut.replace(s2, newvar))
            else:
                new_cuts.append(cut)

        elif svar == oldvar:
            new_cuts.append(cut.replace(oldvar, newvar))
        else:
            new_cuts.append(cut)

    return ' && '.join(new_cuts)
#===================================================================================================


#===================================================================================================
#===================================================================================================
#=====================================      OTHER UTILS
#===================================================================================================
#===================================================================================================

#===================================================================================================
def sum_histograms(histograms):
    """
    histograms: list of list of histograms, e.g. [(h1, h2, h3), (h4, h5, h6)]
    return: [h1+h4, h2+h5, h3+h6]
    """

    # fix if missing campaign and empty list
    histograms = [h for h in histograms if h]

    new_histograms = []
    for hlist in zip(*histograms):
        hist = hlist[0].Clone()
        for h in hlist[1:]:
            hist.Add(h, 1)

        new_histograms.append(hist)

    return new_histograms
#===================================================================================================

#===================================================================================================
def get_skim_path(name, paths, version, mc_campaign=None):
    for skim_dir in paths:
        if mc_campaign is None:
            full_guess_path = f'{skim_dir}/v{version}/{name}.root'
        else:
            full_guess_path = f'{skim_dir}/v{version}/{name}_{mc_campaign}.root'

        try:
            paths = glob.glob(full_guess_path)
        except:
            pass

        if paths:
            return paths[0]
    return
#===================================================================================================

#===================================================================================================
def is_2d_variable(variable):
    return ':' in variable and not '::' in variable
#===================================================================================================

#===================================================================================================
def get_escaped_variable(variable):
    return variable.replace('y_', 'ph_').replace(':', '_').replace('/', '').replace('(', '').replace(')', '').replace('[', '').replace(']', '').replace('.', '_')
#===================================================================================================

#===================================================================================================
def fix_histogram_name(hist, name, slicename=''):
    replace_dict = {
        'efake15': 'efake',
        'efake16': 'efake',
        'efake17': 'efake',
        'efake18': 'efake',
        'jfake15': 'jfake',
        'jfake16': 'jfake',
        'jfake17': 'jfake',
        'jfake18': 'jfake',
        'data15Nom': 'data',
        'data16Nom': 'data',
        'data17Nom': 'data',
        'data18Nom': 'data',

        'smr15': 'smr',
        'smr16': 'smr',
        'smr17': 'smr',
        'smr18': 'smr',
    }

    hname = hist.GetName()

    to_replace = '___' + hname.split('___')[1] + '___'
    replace_with = '__'+name+'__'
    if slicename: replace_with += slicename+'__'

    new_hname = hname.replace(to_replace, replace_with)
    for i, j in replace_dict.items():
        if i in new_hname:
            new_hname = new_hname.replace(i, j)

    hist.SetName(new_hname)
    return
#===================================================================================================

#===================================================================================================
def get_sumw(path):
    sumw = 0
    if os.path.isdir(path):
        all_files = glob.glob(path+'/*root*')
        for fpath in all_files:
            f = RT.TFile.Open(fpath)
            tmp = f.Get('events')
            sumw += tmp.GetBinContent(3)  # bin 3 is the initial sumw
            f.Close()
    else:
        f = RT.TFile.Open(path)
        tmp = f.Get('events')
        sumw = tmp.GetBinContent(3)  # bin 3 is the initial sumw
        f.Close()

    return sumw
#===================================================================================================

#===================================================================================================
def get_lumi_weight(path, dsid, lumi, fs=None):

    luminosity = float(lumi)
    sumw       = get_sumw(path)
    xs         = xsutils.get_xs_from_dsid(dsid)

    try:
        weight = (luminosity * xs) / sumw
    except:
        weight = 0.

    return weight
#===================================================================================================

#===================================================================================================
def get_lumi(year):
    if year == '2015+2016':
        return lumi_dict['2015'] + lumi_dict['2016']
    else:
        return lumi_dict[year]
#===================================================================================================
