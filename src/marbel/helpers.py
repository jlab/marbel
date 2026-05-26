from xarray import DataTree


def switch_pm_version(pymc_output, key, np=False):
    """
    pymc switched their output data structure for sample_prior_predictive. We test which version according to output type.
    """
    if isinstance(pymc_output, DataTree):
        if np:
            res = pymc_output.to_dict()["/prior"][key].to_dataframe()[key].to_numpy()
        else:
            res = pymc_output.to_dict()["/prior"][key].to_dataframe()[key].to_list()
    else:
        if np:
            res = pymc_output.to_dataframe()[key].to_numpy()
        else:
            res = pymc_output.to_dataframe()[key].to_list()
    return res
