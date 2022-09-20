import os
import re
import glob
from array import array
from functools import partial
from tqdm import tqdm

import ROOT as RT
from .datasets_manager import get_datasets, get_mccampaign
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

pass_mask_function = """
bool pass_mask(ROOT::VecOps::RVec<int> vec){
    // in case the following line is true, there's at least one element that doesnt pass the mask
    bool any_not_pass = std::any_of(vec.begin(), vec.end(), [](int i){return i==0;});

    return !any_not_pass;
}
"""






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
                 use_mcveto=True, use_phjet_w=False,
                 slices=False,
                 remove_var_cut=False,
                 vars_cxx=None,
                #  seed_selection='(seed_met_et - 8)/sqrt(seed_sumet) < (1.0 + 0.1*seed_bjet_n)',
                 preselection=None,
                 ignore_missing=None,
                 systs_module=None,
                 looper_func=None,
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
        self.slices           = slices
        self.remove_var_cut   = remove_var_cut
        # self.seed_selection   = seed_selection
        self.ignore_missing   = ignore_missing
        self.debug            = debug

        # tree name
        if tree_name is not None:
            self.tree_name = tree_name
        else:
            self.tree_name = 'smini' if self.use_skim else 'mini'        

        # dataset year
        self.dataset_year = year
        self.mc_campaign, self.run = get_mccampaign(self.dataset_year, self.ana_rel)

        print('The configuration for histogram_getter is:')
        for var_name, var_val in vars(self).items():
            print(f'  {var_name:17}= {var_val}')

        self.paths            = paths
        self.samples          = samples
        self.regions          = regions
        self.selections       = selections
        self.preselection     = preselection
        self.binning          = binning
        
        self.vars_cxx         = vars_cxx
        
        self.looper_func      = looper_func
        
        if self.vars_cxx:            
            RT.gInterpreter.Declare(self.vars_cxx)

        RT.EnableImplicitMT()

        if systs_module:
            import systs_module as systematics_
            self.systematics = systematics_

        return
    #===============================================================================================
    
    #===============================================================================================
    def _get_multi_histograms(self, name, path, is_mc, lumi, regions, selections, variables,
                              systematics=['Nom', ], binning=None, dsid_str=None):
        """
        get histogram for a given sample (ds) and year
        """
        is_fake  = ('efake' in name or 'jfake' in name)
        is_phjet = ('photonjet' in name)
        # is_smr   = ('smr' in name)

        variables  = variables  if isinstance(variables , list) else [variables , ]
        if regions:
            regions = regions if isinstance(regions, list) else [regions, ]
        if selections:
            selections = selections if isinstance(selections, list) else [selections, ]


        # ------------
        # File/Chain
        # ------------

        # open tree or chain, depending if the path is a directory or if it is a single file
        if os.path.isdir(path):
            all_files = glob.glob(path+'/*root*')
            df        = RT.RDataFrame(self.tree_name, all_files)
        else:
            df = RT.RDataFrame(self.tree_name, path)

        leaves = list(df.GetColumnNames())
        isvector_leave = ['::RVec' in df.GetColumnType(l) for l in leaves]
        
    
        # Lumi weight is the same for all histograms
        if is_mc:
            if not self.weights_strings:
                if self.use_skim:
                    lumi_weight = 'weight_xs*%.2f' % lumi
                elif dsid_str:
                    # temporal hack for wrong xs in AMI and PMG file
                    if 900474 <= int(dsid_str) <= 900487:
                        lumi_weight = get_lumi_weight(path, int(dsid_str), lumi, 'local', self.debug)
                    else:
                        lumi_weight = get_lumi_weight(path, int(dsid_str), lumi, debug=self.debug)
                    
                    if not self.looper_func:
                        lumi_weight = '%s' % lumi_weight
                else:
                    lumi_weight = ''
            else:
                if 'weight_lumi' not in self.weights_strings:
                    lumi_weight = str(lumi)

        if self.debug:
            print('Lumi weight: ', lumi_weight)

        # ----------------------------------------------
        # Create histograms and "tuples" for MultiDraw
        # ----------------------------------------------
        sum_weights = []
        histograms  = []

        # get a list of selections or regions
        if regions and selections:
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
                    
                    
                    
                    # ------- SETUP SELECTIONS
                    _selection = selection                    
                    systname   = syst.replace('__1down', 'Low').replace('__1up', 'High')
                    
                    
                    # MC veto
                    if is_mc and self.use_mcveto:
                        if self.wg_label is None:
                            if _selection:
                                _selection = '%s && mcveto==0' % _selection
                            else:
                                _selection = 'mcveto==0'

                        # skim pt slices
                        if not self.use_skim and dsid_str:
                            dsid = int(dsid_str)
                            _selection = select_truth_slices(dsid, _selection)

                    # Remove variable from selection if n-1
                    if self.remove_var_cut and variable_name in _selection and variable_name != 'cuts':
                        _selection = '&&'.join([cut for cut in _selection.split('&&') if not split_cut(cut)[0] == variable_name])

                    # change selection and variable for systematics
                    if syst != 'Nom' and self.systematics.affects_kinematics(syst):
                        for var in self.systematics.get_affected_variables(syst):
                            _selection = replace_var(_selection, var, '%s_%s' % (var, syst))

                        if variable in self.systematics.get_affected_variables(syst):
                            variable = '%s_%s' % (variable, syst)
                    
                    
                    # check for variables in the selection that takes the whole vector
                    for _sel in split_selection(_selection):
                        _selvar = split_cut(_sel)[0]
                        _newvar = f'mask_{get_escaped_variable(_selvar)}'
                        if _selvar in leaves:
                            if isvector_leave[leaves.index(_selvar)]:
                                if _newvar not in df.GetColumnNames():
                                    df = df.Define(_newvar, _selection).\
                                            Define(f'pass{_newvar}', f'pass_mask({_newvar})')
                                _selection = _selection.replace(_sel, f'pass{_newvar}')
                    df_selection = df.Filter(_selection, region)
                    
                    
                    
                    
                    
                    # ------- SETUP WEIGHTS
                    w_list = []
                    if is_mc:
                        # all weights for PhotonID case
                        if self.weights_strings:
                            w_list.append(self.weights_strings['totalw'])
                        
                        # lumi weight
                        if self.use_lumiw:
                            if not self.weights_strings:
                                w_list.append('%s' % lumi_weight)
                            else:
                                if 'weight_lumi' not in self.weights_strings:
                                    w_list.append('%s' % lumi_weight)
                                else:
                                    w_list.append(self.weights_strings['weight_lumi'])

                        # mc weight
                        if self.use_mcw and not self.weights_strings:
                            w_list.append('weight_mc')

                        # scale factors
                        if self.use_sfw and not self.weights_strings:
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
                    
                    else:
                        if self.weights_strings:
                            w_list.append(self.weights_strings['weight_data'])
                    
                    
                    w_str = '*'.join(w_list) if self.scale and w_list else '1'
                    
                    
                    df_selection_w = df_selection.Define("weight", w_str)
                    

                    
                    
                    
                    
                    
                    # ------- SETUP HISTOGRAMS
                    # aux variable. In case the variable is coming from a function, check the extra
                    # variables dictionary
                    variable_name = variable
                    if variable not in leaves and 'lumi' not in variable:
                        variable_name = get_var_function(variable)

                    # get binning
                    if binning is None or not binning:
                        _binning = get_binning(variable_name, self.binning)
                    else:
                        _binning = binning[ivariable]

                    # name to avoid the ROOT warning, not used
                    if self.use_skim:
                        hname = f'h___{name}___{systname}__{region}__{get_escaped_variable(variable_name)}'
                    elif dsid_str:
                        hname = f'h___{dsid_str}___{systname}__{region}__{get_escaped_variable(variable_name)}'
                    else:
                        hname = f'h___{name}___{systname}__{region}__{get_escaped_variable(variable_name)}'
                    
                    
                    
                    if is_2d_variable(variable):
                        vx, vy = variable.split(':')
                        if '[' in vx or ']' in vx:
                            df_selection_w = df_selection_w.Define(get_escaped_variable(vx), vx)
                            vx = get_escaped_variable(vx)
                        if '[' in vy or ']' in vy:
                            df_selection_w = df_selection_w.Define(get_escaped_variable(vy), vy)
                            vy = get_escaped_variable(vy)
                        
                        htemp = df_selection_w.Histo2D((hname, hname, *_binning), vx, vy, "weight")
                    else:
                        if variable == 'cuts':
                            sumw = df_selection_w.Sum('weight')
                            sum_weights.append((hname, sumw))                            
                            continue
                        
                        else:
                            if '[' in variable or ']' in variable:
                                df_selection_w = df_selection_w.Define(get_escaped_variable(variable), variable)
                                var_aux = get_escaped_variable(variable)
                            else:
                                var_aux = variable
                                
                            
                            if len(_binning) > 3:
                                htemp = df_selection_w.Histo1D((hname, '', len(_binning)-1, array('d', _binning)), var_aux, "weight")
                            else:
                                htemp = df_selection_w.Histo1D((hname, '', int(_binning[0]), _binning[1], _binning[2]), var_aux, "weight")
                    
                    histograms.append(htemp)


        if self.looper_func:
            # in case some preselection was passed, get the cut for the preselection and adding the truth cuts
            _, _, selection = draw_list[0]
            
            looper = tree_looper(tree, self.looper_func)
            histograms = looper.loop(presel=selection, dsid=dsid_str, lumi_weight=lumi_weight)
        
        else:
            # trigger event loop
            if sum_weights:
                sum_weights[0][1].GetValue()
                for sumw in sum_weights:
                    htemp = histogram(sumw[0], *get_binning('cuts', self.binning))
                    htemp.SetBinContent(1, sumw[1].GetValue())
                    histograms.append(htemp)
            else:
                htemp.GetValue()
            del df, df_selection, df_selection_w

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
            if not self.lumi:
                lumi = get_lumi(dataset_year)
            else:
                lumi = self.lumi


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
            mc_campaign, _ = get_mccampaign(dataset_year, self.ana_rel)
            datasets = get_datasets(sample, self.paths, self.samples, self.version,
                                    self.ignore_missing, mc_campaign, extra_regex)

            histograms = []
            if self.slices and is_mc:
                histograms_slices = {}

            name_to_show = f'{sample}/{mc_campaign}' if is_mc else sample
            
            print('Using samples:')
            for ds in datasets:
                print('   %s' % ds['path'])

            for ds in tqdm(datasets, desc=name_to_show):
                try:
                    dsid = ds['dsid']
                except KeyError:
                    dsid = None
                
                histograms_ds = self._get_multi_histograms(sample, ds['path'], is_mc, lumi, regions,
                                                           selections, variables, dsid_str=dsid)

                if self.slices and is_mc:
                    histograms_slices[ds['short_name']] = histograms_ds.GetValue()

                if not histograms:
                    
                    for hist in histograms_ds:
                        try:
                            histograms.append(hist.GetValue().Clone())
                        except:
                            histograms.append(hist.Clone())
                else:
                    for hall, hnew in zip(histograms, histograms_ds):
                        try:
                            hall.Add(hnew, 1)
                        except:
                            hall.Add(hnew.GetValue(), 1)

        # Fix histogram name and add overflow bin
        for hist in histograms:
            if '__' in hist.GetName():
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
    def get_events(self, sample, selections, regions=None):

        hists = self.get_histograms(sample, regions, selections, 'cuts')
        
        values = []
        for h in hists:
            mean  = h.GetBinContent(1)
            error = h.GetBinError(1)
            values.append(Value(mean, error))

        return values
    #===============================================================================================

    #===============================================================================================
    def get_cutflow(self, sample, selection, regions=None):
        """Creates a cutflow histogram given a region and selection. It separates the selection into
        different steps and calculates the number of events passing those steps.

        Args:
            sample (str): name of the sample
            region (str): name of the region
            selection (str): selection to be applied and that is going to be separated into steps

        Returns:
            TH1F: histogram with cutflow
        """

        # create the histogram and split the cuts
        cuts    = [split_cut(cut) for cut in split_selection(selection)]
        cutflow = histogram('cutflow', len(cuts)+1, 0.5, len(cuts)+1.5)

        # histogram bin labels
        cutflow.GetXaxis().SetBinLabel(1, 'No Cut')
        for i, cut in enumerate(cuts):
            cutflow.GetXaxis().SetBinLabel(i+2, ''.join(c for c in cut))

        # add the first 'No Cut' cut
        cuts = [('', '', ''), ] + cuts

        # now we want to concatenate, cut by cut, to the selection that is going to be applied in each
        # step
        selections = []
        for i, (var, op, value) in enumerate(cuts):
            # get the last selection saved in 'selections'
            if i == 0:
                curr_sel = ''
            else:
                curr_sel = selections[-1]
            
            if curr_sel.strip() == '':
                curr_sel = f'{var} {op} {value}'
            else:
                curr_sel += f' && {var} {op} {value}'
        
            selections.append(curr_sel)
        
        # get the events for each cut
        events_cuts = self.get_events(sample, selections, regions)
        
        for i, e in enumerate(events_cuts,1):
            cutflow.SetBinContent(i, e.mean)
            cutflow.SetBinError  (i, e.error)
        
        return cutflow
    #===============================================================================================


#===================================================================================================
def select_truth_slices(dsid, selection):
    
    
    if dsid in (361042, 361043, 361044, 364543, 361045, 361046, 361047, 364544,
                361048, 361049, 361050, 364545, 361051, 361052, 361053, 364546,
                361054, 361055, 361056, 361057, 361058, 361059, 364547,
                800662, 800676, 800663, 800677, 800664, 800678,
                800665, 800679, 800666, 800680, 800667, 800681,
                800668, 800682, 800683, 800669, 800670, 800671) and \
        selection.strip() != '':
        selection += ' &&'
    # sherpa
    if dsid in (361042, 361043, 361044, 364543):
        if selection.strip() != '':
            selection += f' ph_truth_pt[0]>70. && ph_truth_pt[0]<140.'
    elif dsid in (361045, 361046, 361047, 364544):
        selection += f' ph_truth_pt[0]>140. && ph_truth_pt[0]<280.'
    elif dsid in (361048, 361049, 361050, 364545):
        selection += f' ph_truth_pt[0]>280. && ph_truth_pt[0]<500.'
    elif dsid in (361051, 361052, 361053, 364546):
        selection += f' ph_truth_pt[0]>500. && ph_truth_pt[0]<1000.'
    elif dsid in (361054, 361055, 361056, 361057, 361058, 361059, 364547):
        selection += f' ph_truth_pt[0]>1000.'
    # pythia
    elif dsid in (800662, 800676):
        selection += f' ph_truth_pt[0]>70. && ph_truth_pt[0]<140.'
    elif dsid in (800663, 800677):
        selection += f' ph_truth_pt[0]>140. && ph_truth_pt[0]<280.'
    elif dsid in (800664, 800678):
        selection += f' ph_truth_pt[0]>280. && ph_truth_pt[0]<500.'
    elif dsid in (800665, 800679):
        selection += f' ph_truth_pt[0]>500. && ph_truth_pt[0]<800.'
    elif dsid in (800666, 800680):
        selection += f' ph_truth_pt[0]>800. && ph_truth_pt[0]<1000.'
    elif dsid in (800667, 800681):
        selection += f' ph_truth_pt[0]>1000. && ph_truth_pt[0]<1500.'
    elif dsid in (800668, 800682):
        selection += f' ph_truth_pt[0]>1500. && ph_truth_pt[0]<2000.'
    elif dsid in (800683,):
        selection += f' ph_truth_pt[0]>2000.'
    elif dsid in (800669,):
        selection += f' ph_truth_pt[0]>2000. && ph_truth_pt[0]<2500.'
    elif dsid in (800670,):
        selection += f' ph_truth_pt[0]>2500. && ph_truth_pt[0]<3000.'
    elif dsid in (800671,):
        selection += f' ph_truth_pt[0]>3000.'
        
    return selection
#===================================================================================================

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
def split_var(variable):
    variable = variable.strip()
    
    for op in ['+', '-', '*', '/']:
        if op in variable:
            var1, var2 = variable.split(op)
            return (var1.strip(), op, var2.strip())
    
    return variable
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
def get_var_from_cut(cut):
    
    var, _, _ = split_cut(cut)
    
    # check if the variable contains multiple subvariables. For example if its something like a+b>=0
    # in case it contains a +, -, * or /, return the different vars and the operator. In case it doesnt
    # returns the same variable
    var = split_var(var)
    if isinstance(var, tuple):
        # var was composed of more variables
        aux_var_list = []
        for subvar in var:
            if get_var_function(subvar):
                aux_var_list.append(get_var_function(subvar))
            else:
                aux_var_list.append(subvar)
        return tuple(aux_var_list)
    else:
        if get_var_function(var): var = get_var_function(var)
        return (var, )
#===================================================================================================

#===================================================================================================
def get_var_function(variable):
    if '(' in variable and variable.split('(', 1)[0] not in ('fabs', 'abs'):
        return variable.split('(')[0]
    elif variable.split('(', 1)[0] in ('fabs', 'abs'):
        return variable
    else:
        return None
#===================================================================================================

#===================================================================================================
def is_2d_variable(variable):
    return ':' in variable and not '::' in variable
#===================================================================================================

#===================================================================================================
def get_2d_variables(variable):
    if is_2d_variable(variable):
        return variable.split(':')
#===================================================================================================

#===================================================================================================
def is_profile_variable(variable):
    if 'pro' in variable and variable.endswith(']') and '[' in variable:
        variable = variable.split('[', 1)[1].rsplit(']', 1)[0]
        if is_2d_variable(variable):
            return variable
    return False
#===================================================================================================

#===================================================================================================
def get_profile_variables(variable):
    variable = is_profile_variable(variable)
    if variable:
        vx, vy = get_2d_variables(variable)
        return vx, vy
    return None, None
#===================================================================================================

#===================================================================================================
def get_escaped_variable(variable):
    variable = variable.replace('y_', 'ph_')
    if ':' in variable and not '::' in variable:
        variable = variable.replace(':', '_')
    elif '::' in variable:
        variable = variable.replace('::', '')
    variable = variable.replace('(', '')
    variable = variable.replace(')', '')
    variable = variable.replace('[', '')
    variable = variable.replace(']', '')
    variable = variable.replace('.', '_')
    variable = variable.replace('*', '_times_')
    variable = variable.replace('/', '_over_')
    return variable
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
def get_sumw(path, debug=False):
    sumw = 0
    if os.path.isdir(path):
        all_files = glob.glob(path+'/*root*')
        for i, fpath in enumerate(all_files):
            f = RT.TFile.Open(fpath)
            tmp = f.Get('events')
            sumw += tmp.GetBinContent(3)  # bin 3 is the initial sumw
            f.Close()
            if debug:
                print(f"file #{i}, filename: {fpath}, sumw: {sumw}")
    else:
        f = RT.TFile.Open(path)
        tmp = f.Get('events')
        sumw = tmp.GetBinContent(3)  # bin 3 is the initial sumw
        f.Close()
        if debug:
            print(f"file #0, filename: {path}, sumw: {sumw}")

    return sumw
#===================================================================================================

#===================================================================================================
def get_lumi_weight(path, dsid, lumi, fs=None, debug=False):

    luminosity = float(lumi)
    sumw       = get_sumw(path)
    xs         = xsutils.get_xs_from_dsid(dsid, fs)
    
    try:
        weight = (luminosity * xs) / sumw
    except:
        weight = 0.
    
    if debug:
        print('DSID:', dsid, 'Luminosity:', luminosity, 'Cross section:', xs, 'Sum of weights:', sumw, 'Weight:', weight)
    return weight
#===================================================================================================

#===================================================================================================
def get_lumi(year):
    lumi = 0.0
    if year == 'full': year = '2015+2016+2017+2018'
    for y in year.split('+'):
        lumi += lumi_dict[y]
    return lumi
#===================================================================================================

#===================================================================================================
def get_comE(year):
    if year == '2022':
        return 13.6
    else:
        return 13
#===================================================================================================


class tree_looper:
    #===============================================================================================
    def __init__(self, tree, looper_func):
        self.tree = tree
        self.loop = partial(looper_func, self)
        return
    #===============================================================================================

    #===============================================================================================
    def set_branches(self, branches):
        self.tree.SetBranchStatus('*', 0)
        for branch in branches:
            self.tree.SetBranchStatus(branch, 1)
        pass
    #===============================================================================================
