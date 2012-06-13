# try to use an ordered dict, which requires a patch cobra cobra.io
# see http://projects.scipy.org/scipy/attachment/ticket/1566
# try:
    # from collections import OrderedDict as dicttype
# except ImportError:
    # dicttype = dict
dicttype = dict

from numpy import array, object as np_object
from scipy.io import loadmat, savemat
from scipy.sparse import coo_matrix

from .. import Model, Metabolite, Reaction


def cell(x):
    """translate an array x into a MATLAB cell array"""
    return array(x, dtype=np_object)


def load_matlab_model(infile_path, variable_name=None):
    """Load a cobra model stored as a .mat file
    NOTE: INCOMPLETE, does not load GPR's

    Parameters
    ----------
    infile_path : str
    variable_name : str, optional
        The variable name of the model in the .mat file. If this is not
        specified, then the first MATLAB variable which looks like a COBRA
        model will be used

    """
    data = loadmat(infile_path)
    if variable_name is not None:
        possible_names = [variable_name]
    else:
        # will try all of the variables in the dict
        possible_names = {}
        for key in data.keys():
            possible_names[key] = None
        # skip meta variables
        to_remove = ["__globals__", "__header__", "__version__"]
        to_pop = []
        for name in possible_names:
            if name in to_remove:
                to_pop.append(name)
        for i in to_pop:
            possible_names.pop(i)
        possible_names = possible_names.keys()
    for possible_name in possible_names:
        m = data[possible_name]  # TODO: generalize
        if m.dtype.names is None:
            continue
        if not set(["rxns", "mets", "S", "lb", "ub"]) \
                <= set(m.dtype.names):
            continue
        model = Model()
        model.id = m["description"][0, 0][0]
        model.description = model.id
        for i, name in enumerate(m["mets"][0, 0]):
            new_metabolite = Metabolite()
            new_metabolite.id = name[0][0]
            try:
                new_metabolite.name = m["metNames"][0, 0][i][0][0]
                new_metabolite.formula = m["metFormulas"][0][0][i][0][0]
            except:
                pass
            model.add_metabolites([new_metabolite])
        for i, name in enumerate(m["rxns"][0, 0]):
            new_reaction = Reaction()
            new_reaction.id = name[0][0]
            new_reaction.lower_bound = m["lb"][0, 0][i][0]
            new_reaction.upper_bound = m["ub"][0, 0][i][0]
            new_reaction.objective_coefficient = m["c"][0, 0][i][0]
            try:
                new_reaction.name = m["rxnNames"][0, 0][i][0][0]
            except:
                pass
            model.add_reactions(new_reaction)
        coo = coo_matrix(m["S"][0, 0])
        for i, j, v in zip(coo.row, coo.col, coo.data):
            model.reactions[j].add_metabolites({model.metabolites[i]: v})
        # TODO finish adding GPR's
        model.update()
        return model
    # If code here is executed, then no model was found.
    raise Exception("no COBRA model found")


def save_matlab_model(model, file_name):
    """Save the cobra model as a .mat file.

    This .mat file can be used directly in the MATLAB version of COBRA.
    NOTE: this function requires a patched version of scipy.io.savemat

    Parameters
    ----------
    model : cobra.Model
    file_name : str or file-like object

    """
    model.update()
    rxns = model.reactions
    mets = model.metabolites
    mat = dicttype()
    csense = ""
    for m in mets:
        csense += m._constraint_sense
    mat["mets"] = cell(mets.list_attr("id"))
    mat["metNames"] = cell(mets.list_attr("name"))
    mat["metFormulas"] = cell([str(m.formula) for m in mets])
    mat["genes"] = cell(model.genes.list_attr("id"))
    mat["grRules"] = cell(rxns.list_attr("gene_reaction_rule"))
    mat["rxns"] = cell(rxns.list_attr("id"))
    mat["rxnNames"] = cell(rxns.list_attr("name"))
    mat["subSystems"] = cell(rxns.list_attr("subsystem"))
    mat["csense"] = csense
    mat["S"] = model._S
    mat["lb"] = array(rxns.list_attr("lower_bound"))
    mat["ub"] = array(rxns.list_attr("upper_bound"))
    mat["b"] = array(mets.list_attr("_bound"))
    mat["c"] = array(rxns.list_attr("objective_coefficient"))
    mat["rev"] = array(rxns.list_attr("reversibility"))
    mat["description"] = str(model.description)
    savemat(file_name, {str(model.description): mat},
             appendmat=True, oned_as="column")
