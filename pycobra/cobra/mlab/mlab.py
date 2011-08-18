#cobra.mlab.py
#System modules
import cPickle,numpy, os, cobra
from copy import deepcopy
from scipy.sparse import dok_matrix
from cobra import Model, Reaction, Metabolite, Formula
mlab_path = cobra.__path__[0] + '/mlab/matlab_scripts/'
from mlabwrap import mlab as matlab
matlab.changeCobraSolver('glpk','LP')
matlab.addpath(mlab_path)
#Project defined modules
def matlab_cell_to_python_list(the_cell):
    try:
        import mlabwrap
    except ImportError:
        print 'Could not import mlabwrap: ' + ImportError + '\n'
        return False
    #the_list = mlabwrap.mlab.cell_to_string(the_cell).rstrip('\t').split('\t')
    #Because the mlabwrap.mlab.cell_to_string function will end with a tab, it must
    #be removed.  Use indexing instead of rstrip because rstrip will mess up
    #when the last element in the string is empty.
    the_list = mlabwrap.mlab.cell_to_string(the_cell).split('\t')
    return the_list[:len(the_list) - 1]

def matlab_logical_to_python_logical(the_logical):
    try:
        import mlabwrap
    except ImportError:
        print 'Could not import mlabwrap: ' + ImportError + '\n'
        return False
    if mlabwrap.mlab.double(the_logical)[0][0] == 1:
        return True
    else:
        return False

def python_list_to_matlab_cell(the_list, transpose = False):
    try:
        import mlabwrap
    except ImportError:
        print 'Could not load mlabwrap ' + ImportError + '\n'
        return False
    if hasattr(the_list[0], 'id'):
        the_list = [x.id for x in the_list]
    #Remove the single apostrophes because matlab or mlabwrap has
    #problems with them.
    the_list = [x.replace(r"'",'') for x in the_list]
    the_string = '{' + str(the_list).lstrip('[').rstrip(']') + '}'
    if transpose:
        the_string = the_string.replace(',',';')
    return mlabwrap.mlab.eval(the_string)

def matlab_sparse_to_scipy_sparse(matlab_sparse_matrix):
    try:
        import mlabwrap
        import scipy.sparse
    except ImportError:
        print 'Could not load mlabwrap or scipy ' + ImportError + '\n'
        return False
    [row_indices, column_indices, the_values] = mlabwrap.mlab.find(matlab_sparse_matrix , nout = 3)
    #Change to 0-based indices
    row_indices = row_indices - 1
    column_indices = column_indices - 1
    the_shape = mlabwrap.mlab.size(matlab_sparse_matrix)
    the_shape =  tuple(map(int, the_shape[0,]))
    #This could be sped up by using a dok_matrix and zipping the indices and values
    #scipy_sparse_matrix = scipy.sparse_dok.matrix(the_shape)
    #scipy_sparse_matrix.update(zip(zip(row_indices, column_indices), the_values))
    #return scipy_sparse_matrix.tolil()
    scipy_sparse_matrix = scipy.sparse.lil_matrix(the_shape)
    for i in range(len(the_values)):
        scipy_sparse_matrix[int(row_indices[i][0]), int(column_indices[i][0])]  = the_values[i][0]
    return scipy_sparse_matrix

def scipy_sparse_to_mlab_sparse(scipy_sparse_matrix):
    """A more efficient method is needed for when the
    matrix is so big that making a dense version
    is a waste of computer effort."""
    try:
        import mlabwrap
    except ImportError:
        print 'Could not load mlabwrap ' + ImportError + '\n'
        return False
    return mlabwrap.mlab.sparse(scipy_sparse_matrix.todense())

def matlab_sparse_to_numpy_array(matlab_sparse_matrix):
    try:
        import mlabwrap
        import numpy
    except ImportError:
        print 'Could not load mlabwrap or numpy ' + ImportError + '\n'
        return False
    [row_indices, column_indices, the_values] = mlabwrap.mlab.find(matlab_sparse_matrix , nout = 3)
    #Change to 0-based indices
    row_indices = row_indices - 1
    column_indices = column_indices - 1
    the_shape = mlabwrap.mlab.size(matlab_sparse_matrix)
    numpy_array = numpy.zeros(the_shape[0,])
    for i in range(len(the_values)):
        numpy_array[int(row_indices[i][0]), int(column_indices[i][0])] = the_values[i][0]
    return numpy_array

def numpy_array_to_mlab_sparse(numpy_array):
    """A more efficient method is needed for when the
    matrix is so big that making a dense version
    is a waste of computer effort."""
    try:
        import mlabwrap
    except ImportError:
        print 'Could not load mlabwrap ' + ImportError + '\n'
        return False

    return mlabwrap.mlab.sparse(numpy_array)

def matlab_cobra_struct_to_python_cobra_object(matlab_struct):
    """Converts a COBRA toolbox 2.0 struct into a cobra.Model object using the mlabwrap matlab proxy

    """

    from cobra import Model
    try:
        from mlabwrap import mlab as matlab
    except:
        raise Exception('mlabwrap and MATLAB are required to use these functions. '+\
                        'They only function on Mac OS X and GNU/Linux')
    from copy import deepcopy
    from numpy import array
    from cobra.mlab import matlab_cell_to_python_list, matlab_sparse_to_scipy_sparse
    struct_fields = matlab_cell_to_python_list(matlab.fields(matlab_struct, nout=1))
    _S_dok = dok_matrix(matlab_sparse_to_scipy_sparse(matlab_struct.S))
    model_id = 'Matlab_Import'
    if 'description' in struct_fields:
        model_id = str(deepcopy(matlab_struct.description))
    the_model = Model(model_id)
    #Metabolite section
    cobra_metabolites = map(Metabolite, matlab_cell_to_python_list(matlab_struct.mets))
    _b = deepcopy(matlab_struct.b).tolist()
    if 'metNames' in struct_fields:
        _metabolite_names = matlab_cell_to_python_list(matlab_struct.metNames)
    else:
        _metabolite_names = ['']*len(cobra_metabolites)
        
    if 'csense' in struct_fields:
        if isinstance(matlab_struct.csense, str):
            the_csense = matlab_struct.csense
            _constraint_sense = [the_csense[i] for i in range(len(the_csense))]
        else:
            _constraint_sense = matlab_cell_to_python_list(matlab_struct.csense)
    if 'metFormulas' in struct_fields:
        _metabolite_formulas = matlab_cell_to_python_list(matlab_struct.metFormulas)
    else:
        _metabolite_formulas = ['']*len(the_model.metabolites)
    for the_metabolite, b, n, c, f in zip(cobra_metabolites,
                                          _b,
                                          _metabolite_names,
                                          _constraint_sense,
                                          _metabolite_formulas):
        the_metabolite._bound = b[0]
        the_metabolite.name = n
        the_metabolite._constraint_sense = c
        the_metabolite.formula = Formula(f)
    metabolite_dict = dict([(x.id, x)
                            for x in cobra_metabolites])
    #Reaction section
    cobra_reactions = map(Reaction, matlab_cell_to_python_list(matlab_struct.rxns))
    _objective_coefficients = deepcopy(matlab_struct.c).tolist()
    _lower_bounds = deepcopy(matlab_struct.lb).tolist()
    _upper_bounds = deepcopy(matlab_struct.ub).tolist()
    _reversibility = deepcopy(matlab_struct.rev).tolist()
    if 'rxnNames' in struct_fields:
        _reaction_names = matlab_cell_to_python_list(matlab_struct.rxnNames)
    else:
        _reaction_names = ['']*len(the_model.reactions)
    if 'grRules' in struct_fields:
        _gene_reaction_rules = matlab_cell_to_python_list(matlab_struct.grRules)
    else:
        _gene_reaction_rules = ['']*len(the_model.reactions)
    if 'subSystems' in struct_fields:
        _subsystems = matlab_cell_to_python_list(matlab_struct.subSystems)
    else:
        _subsystems = ['']*len(the_model.reactions)
    for the_reaction, o, l, u, n, r, s, rev in zip(cobra_reactions,
                                              _objective_coefficients,
                                              _lower_bounds,
                                              _upper_bounds,
                                              _reaction_names,
                                              _gene_reaction_rules,
                                              _subsystems, _reversibility ):
        the_reaction.objective_coefficient = o[0]
        the_reaction.lower_bound = l[0]
        the_reaction.upper_bound = u[0]
        the_reaction.name = n
        the_reaction.gene_reaction_rule = r
        the_reaction.subsystem = s
        the_reaction.reversibility = rev[0]
        the_reaction.parse_gene_association()
    index_to_reaction = dict(zip(range(len(cobra_reactions)),
                                 cobra_reactions))
    index_to_metabolite = dict(zip(range(len(cobra_metabolites)),
                                 cobra_metabolites))
    for k, coefficient in _S_dok.items():
        the_reaction = index_to_reaction[k[1]]
        the_metabolite = index_to_metabolite[k[0]]
        the_reaction.add_metabolites({the_metabolite: coefficient})
    #solution needs to be redesigned from a matlab proxy.
    #Now populate all of the reactions
    #NOTE: Having a metabolite referenece here would be a good idea
    the_model.add_reactions(cobra_reactions)
    the_model.update()
    return the_model


def cobra_model_object_to_cobra_matlab_struct(cobra_model):
    """This function converts all of the  object values to the
    corresponding model fields to update the mlab model.

    """
    try:
        from mlabwrap import mlab as matlab
    except:
        raise Exception('mlabwrap and MATLAB are required to use these functions. '+\
                        'They only function on Mac OS X and GNU/Linux')
    from cobra.mlab import python_list_to_matlab_cell, scipy_sparse_to_mlab_sparse
    cobra_model.update()
    matlab_struct = matlab.struct()
    #Things that need a conversion:  S, rxnGeneMat,
    matlab_struct.mets = python_list_to_matlab_cell([x.id for x in cobra_model.metabolites],
                                                    transpose = True)
    matlab_struct.metNames =  python_list_to_matlab_cell([x.name for x in cobra_model.metabolites],
                                                         transpose = True)
    matlab_struct.metFormulas = python_list_to_matlab_cell([x.formula
                                                            for x in cobra_model.metabolites],
                                                           transpose = True)
    matlab_struct.genes = python_list_to_matlab_cell(cobra_model.genes, transpose = True)
    matlab_struct.grRules = python_list_to_matlab_cell([x.gene_reaction_rule
                                                        for x in cobra_model.reactions],
                                                       transpose = True)
    matlab_struct.rxns = python_list_to_matlab_cell([x.id for x in cobra_model.reactions],
                                                    transpose = True)
    matlab_struct.rxnNames = python_list_to_matlab_cell([x.name
                                                         for x in cobra_model.reactions],
                                                        transpose = True)
    matlab_struct.subSystems = python_list_to_matlab_cell([x.subsystem
                                                           for x in cobra_model.reactions],
                                                           transpose=True)

    if hasattr(cobra_model, '_constraint_sense'):
        matlab_struct.csense = reduce(lambda x,y: x+y, cobra_model._constraint_sense)
    #matlab_struct.csense = python_list_to_matlab_cell(['E']*len(cobra_model.metabolites), transpose = True)
    #TODO: inefficient conversion but who cares? matlab's on its way out
    matlab_struct.S = scipy_sparse_to_mlab_sparse(cobra_model._S)
    #Things that can be directly copied
    matlab_struct.b = cobra_model._b
    matlab_struct.c = cobra_model._objective_coefficients
    matlab_struct.lb = cobra_model._lower_bounds
    matlab_struct.ub = cobra_model._upper_bounds
    matlab_struct.rev = [x.reversibility for x in cobra_model.reactions]
    matlab_struct.description = cobra_model.description
    return(matlab_struct)



if __name__ == '__main__':
    from cPickle import load
    from time import time
    from numpy import round
    from cobra.manipulation import initialize_growth_medium
    test_directory = '../test/data/'
    with open(test_directory + 'salmonella.pickle') as in_file:
        cobra_model = load(in_file)
    initialize_growth_medium(cobra_model, 'LPM')
    py_cobra_solution = repr(cobra_model.solution.f)
    matlab_struct = cobra_model_object_to_cobra_matlab_struct(cobra_model)
    matlab_result = matlab.optimizeCbModel(matlab_struct)
    matlab_solution = repr(float(matlab_result.f))
    if py_cobra_solution[:4] == matlab_solution[:4]:
        print 'SUCCESS: growth rate match between pyCOBRA and COBRA Toolbox: %s ~ %s'%(py_cobra_solution,
                                                                                       matlab_solution)
    else:
        print 'FAILURE: pyCOBRA and COBRA Toolbox do not match: %s !~ %s'%(py_cobra_solution,
                                                                                       matlab_solution)

        
