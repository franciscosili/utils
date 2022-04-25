#===================================================================================================
def get_binning_single_variable(variable, binning_dict):
    """Retrieve binning for a single variable looking for the binning inside binning_dict

    Args:
        variable (str): name of the variable

    Returns:
        tuple: tuple containing the number of bins as the first element, and the second and third are
        the variable limits
    """

    # retrieve binning for the variable
    binning = binning_dict.get(variable, None)

    if binning is None and '[' in variable and ']' in variable:
        binning = binning_dict.get(variable[:variable.index('[')], None)

    if variable is None:
        for var in binning_dict.keys():
            if var in variable:
                binning = binning_dict[var]
                break

    if binning is None:
        try:
            binning = binning_dict.get(variable.split('_')[1], None)
        except:
            binning = None

    if binning is None:
        try:
            binning = binning_dict.get(variable.split('_')[0], None)
        except:
            binning = None

    if binning is None and 'dphi' in variable:
        binning = binning_dict.get('dphi', None)
        
    if binning is None and '/' in variable:
        binning = binning_dict.get(variable.split('/')[0], None)
    
    if binning is None and '*' in variable:
        binning = binning_dict.get(variable.split('/')[0], None)

    return binning
#===================================================================================================

#===================================================================================================
def get_binning(variable, binning_dict):
    """Get the variable for a single variable or a combination of two of them.

    Args:
        variable (str): variable or variables, separated by :

    Returns:
        tuple: tuple containing the number of bins as the first element, and the second and third are
        the variable limits
    """
    # in case the variable contains two variables instead
    if ':' in variable and not '::' in variable:
        varx, vary = variable.split(':')

        binning_x = get_binning_single_variable(varx, binning_dict)
        binning_y = get_binning_single_variable(vary, binning_dict)

        binning = binning_x + binning_y
    else:
        binning = get_binning_single_variable(variable, binning_dict)

    if binning is None:
        print('Not bins configured for this variable %s. Using default binning' % variable)
        binning = binning_dict['default']

    return binning
#===================================================================================================
