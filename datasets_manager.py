import re
import sys
from glob import glob



# _r_ds        = re.compile(f'({_mc_str}|{_data_str}|{_fakes_str})\.([0-9]*)\.(.*)')
# _r_period_ds = re.compile(f'({_mc_str}|{_data_str}|{_fakes_str}|{_smr_str})\.(period[A-Z])\.(.*)')


_smr_str    = 'smr15|smr16|smr17|smr18'
_regex_datatype = [
    'mc16_13TeV|mc20_13TeV|mc21_13p6TeV|mc16|mc20|mc21'
    'data15|data16|data17|data18|data22|data15_13TeV|data16_13TeV|data17_13TeV|data18_13TeV|data22_13p6TeV'
    'efake15|efake16|efake17|efake18|jfake15|jfake16|jfake17|jfake18'
]

_regex_end     = '([0-9]*)\.(.*)'
_regex_end_per = '(period[A-Z])\.(.*)'

_regex_start = '|'.join(e for e in _regex_datatype)
_period_regex_start = '|'.join(e for e in [_regex_start, _smr_str])

r_ds     = re.compile(f'({_regex_start})\.{_regex_end}')
r_ds_per = re.compile(f'({_period_regex_start})\.{_regex_end_per}')

#===================================================================================================
# TODO: implement this function
# def get_ami_tags(file):
#     file = file.rsplit('/')[-1]
#     return
#===================================================================================================

#===================================================================================================
def get_mccampaign(year, ana_release):    
    if year in ('2015', '2016', '2015+2016'):
        run = 2
        mc_campaign = 'mc16a' if ana_release < 22 else 'mc20a'
    elif year == '2017':
        run = 2
        mc_campaign = 'mc16d' if ana_release < 22 else 'mc20d'
    elif year == '2018':
        run = 2
        mc_campaign = 'mc16e' if ana_release < 22 else 'mc20e'
    elif year == '2015+2016+2017+2018' or year == 'Run2':
        run = 2
        mc_campaign = 'mc16a+mc16d+mc16e' if ana_release < 22 else 'mc20a+mc20d+mc20e'
    elif year == '2022':
        run = 3
        mc_campaign = 'mc21'
    
    return mc_campaign, run
#===================================================================================================

#===================================================================================================
def _find_path(samples_paths, project=None, dsid=None, short_name=None, version=None, mc_campaign=None, full_ds=None):
    """Get the path with the input mini-ntuples

    Args:
        project (str): whether if project is data or MC
        dsid (str): dataset id
        short_name (str): short name of the sample indicating the generator, slice, process, etc
        version (str): version of the ntuples
        mc_campaign (str): mc campaing

    Returns:
        str: full path of sample
    """
    # set versions number.
    # v1, v2 = (int(v) for v in version.split('_'))
    if full_ds:
        guess_path = full_ds
    elif project.startswith('mc') and mc_campaign is not None:
        guess_path = f'v{version}/*{project}.{dsid}.{short_name}.mini.{mc_campaign}.p*.v{version}_output*'
    else:
        guess_path = f'v{version}/*{project}.{dsid}.{short_name}.mini.p*.v{version}_output*'

    for mini_dir in samples_paths:
        full_guess_path = f'{mini_dir}/{guess_path}'
        try:
            paths = glob(full_guess_path)
            if paths:
                # TODO: sort by ptag and choose the latest
                return paths[0]
        except:
            pass

    return None
#===================================================================================================

#===================================================================================================
def _get_dsnames(name, samples, version):
    """Get a list of datasets names with name 'name'

    Args:
        name : name of the dataset. Can be a str or list
        version (str): version of the samples

    Returns:
        list: list containing datasets names
    """
    # get temporary name of the dataset
    ds_tmp = samples.get(name)

    dsnames = []
    if isinstance(ds_tmp, str):
        if '+' in ds_tmp:
            names = [ i.strip() for i in ds_tmp.split('+')  if i ]
            for name in names:
                dsnames += _get_dsnames(name, samples, version)
        else:
            dsnames.append(ds_tmp)
    else:
        dsnames = ds_tmp

    return dsnames
#===================================================================================================

#===================================================================================================
def get_datasets(name, paths, samples, version=None, ignore_missing=True, mc_campaign=None, extra_regex=None):
    """Retrieve information from the datasets with name 'name'.

    Args:
        name (str/list): name of the dataset, or name of the list containing multiple datasets
        version (str, optional): version of the samples. Defaults to None.
        mc_campaign (str, optional): mc campaign. Defaults to None.
        ignore_missing (bool, optional): whether to ignore missing files or not. Defaults to False.
    Raises:
        Exception: Sample not found the samples.py
        Exception: 
        Exception: File not found

    Returns:
        list: list containing dicts with information on the name, project, mc campaign, dsid, short
        name and path of the samples
    """
    # get datasets corresponding to sample
    dsnames = _get_dsnames(name, samples, version)

    if not dsnames:
        raise Exception('Sample %s not found in samples dictionary' % name)

    if extra_regex:
        r_extra = re.compile(extra_regex)

    # loop over the dataset names
    datasets = []
    for ds in dsnames:
        try:
            if extra_regex:
                m = r_extra.match(ds)
                full_ds = m.group(0)
            else:
                m = r_ds.match(ds)
                project, dsid, short_name = m.group(1), m.group(2), m.group(3)
        except:
            try:
                if extra_regex:
                    m = r_extra.match(ds)
                    full_ds = m.group(0)
                else:
                    m = r_ds_per.match(ds)
                    project, dsid, short_name = m.group(1), m.group(2), m.group(3)
            except:
                raise Exception(ds)


        # get path of the sample
        if extra_regex:
            path = _find_path(paths, full_ds=full_ds)
        else:
            path = _find_path(paths, project, dsid, short_name, version, mc_campaign)

        if not path:
            if ignore_missing:
                continue
            else:
                print(f'File not found for ds {ds} with campaign {mc_campaign}. Using version {version}')
                # raise Exception(f'File not found for ds {ds} with campaign {mc_campaign}. Using version {version}')

        if extra_regex:
            dataset = {
                'name': name,
                'path': path
            }
        else:
            dataset = {
                'name'       : name,
                'project'    : project,
                'mc_campaign': mc_campaign,
                'dsid'       : dsid,
                'short_name' : short_name,
                'path'       : path,
            }
        
        datasets.append(dataset)        
    return datasets
#===================================================================================================
