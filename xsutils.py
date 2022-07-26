# Signal cross sections
import os

if 'PJ_ANALYSIS' in os.environ:
    data_dir = os.environ['PJ_ANALYSIS'] + '/data/'
    xs_file_local = os.path.join(data_dir, 'CrossSectionData.txt')

xs_file = '/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/dev/PMGTools/PMGxsecDB_mc16.txt'

_xs_db       = dict()
_xs_unc_db   = dict()

#===================================================================================================
def _create_xs_db(file='cmvfs'):
    if file == 'local':
        xsfile = xs_file_local
    else:
        xsfile = xs_file
    with open(xsfile) as f:
        for line in f:
            line = line.replace('\n', '')
            if not line or line.startswith('#'):
                continue

            try:
                dsid, name, gen_xs, filter_eff, kfact, unc_up, unc_dn, gen_name, etag = line.split()
            except:
                continue

            # effective cross-section and relative uncertainty
            xseff = float(gen_xs) * float(kfact) * float(filter_eff)

            _xs_db[int(dsid)] = xseff
            _xs_unc_db[int(dsid)] = unc_up
    return
#===================================================================================================

#===================================================================================================
def get_xs_from_dsid(dsid, file):
    if not _xs_db:
        _create_xs_db(file)
        
        
    if dsid in _xs_db:
        return float(_xs_db[dsid])

    raise Exception('ERROR: XS not found for DSID=%s' % (dsid))
#===================================================================================================

# Getting the uncertainty in the xsection!
#===================================================================================================
def get_xs_unc_from_did(did):
    did = int(did)

    if not _xs_unc_db:
        _create_xs_db()

    if did in _xs_unc_db:
        return _xs_unc_db[did]

    raise Exception('ERROR: XS unc not found for DID=%s' % (did))
#===================================================================================================

#===================================================================================================
if __name__ == '__main__':
    
    _create_xs_db()

    for did, xs in _xs_db.items():
        print('%i : %.5f' % (did, xs))
#===================================================================================================
