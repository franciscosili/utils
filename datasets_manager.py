import re
from glob import glob

_r_ds        = re.compile('(mc16_13TeV|mc20_13TeV|mc21_13p6TeV|data15_13TeV|data16_13TeV|data17_13TeV|data18_13TeV|data22_13p6TeV|efake15|efake16|efake17|efake18|jfake15|jfake16|jfake17|jfake18)\.([0-9]*)\.(.*)')
_r_period_ds = re.compile('(mc16_13TeV|mc20_13TeV|mc21_13p6TeV|data15_13TeV|data16_13TeV|data17_13TeV|data18_13TeV|data22_13p6TeV|efake15|efake16|efake17|efake18|jfake15|jfake16|jfake17|jfake18|smr15|smr16|smr17|smr18)\.(period[A-Z])\.(.*)')


#===================================================================================================
# TODO: implement this function
# def get_ami_tags(file):
#     file = file.rsplit('/')[-1]
#     return
#===================================================================================================

#===================================================================================================
def _find_path(project, dsid, short_name, paths, version, mc_campaign):
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

    if project.startswith('mc') and mc_campaign is not None:
        guess_path = f'v{version}/*{project}.{dsid}.{short_name}.mini.{mc_campaign}.p*.v{version}_output*'
    else:
        guess_path = f'v{version}/user.*.{project}.{dsid}.{short_name}.mini.p*.v{version}_output*'

    for mini_dir in paths:
        full_guess_path = f'{mini_dir}/{guess_path}'
        try:
            paths = glob.glob(full_guess_path)
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
def get_datasets(name, paths, samples, version=None, ignore_missing=True, mc_campaign=None):
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

    # loop over the dataset names
    datasets = []
    for ds in dsnames:
        try:
            m = _r_ds.match(ds)
            project, dsid, short_name = m.group(1), m.group(2), m.group(3)
        except:
            try:
                m = _r_period_ds.match(ds)
                project, dsid, short_name = m.group(1), m.group(2), m.group(3)
            except:
                raise Exception(ds)

        # get path of the sample
        path = _find_path(project, dsid, short_name, paths, version, mc_campaign)

        if not path:
            if ignore_missing:
                continue
            else:
                raise Exception(f'File not found for ds {ds} with campaign {mc_campaign}. Using version {version}')

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
