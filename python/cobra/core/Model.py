#cobra.core.Model.py
#TODO: Improve copy time.  Perhaps use references for metabolites in
#reactions.
#Here we're dealing with the fact that numpy isn't supported in jython
#and we've made some numjy modules that use the cern.colt matrices to
#mimic the used numpy behavior.
import sys
if hasattr(sys, 'JYTHON_JAR'):
    raise Exception("Experimental modules of numpy/scipy for java that are" +\
    "not yet ready for prime time.")
    from cobra.numpy import array, empty, sign, repeat, vstack, hstack
    from cobra.numpy import ones, nonzero, zeros
    from cobra.scipy import sparse
    from cobra.scipy.sparse import lil_matrix, dok_matrix
    from cobra.numpy import where, array, ones
else:
    from numpy import array, empty, sign, repeat, vstack, hstack
    from numpy import ones, nonzero, zeros
    from scipy import sparse
    from scipy.sparse import lil_matrix, dok_matrix
from warnings import warn
from copy import deepcopy
from cobra.query import *
from cobra.flux_analysis import optimize_cplex, optimize_glpk, optimize_gurobi
from Object import Object
from Reaction import Reaction
from Metabolite import Metabolite
from Formula import Formula
from cobra.collections.DictList import DictList
#*********************************************************************************
#2011-10-11: danielhyduke@gmail.com
#
#We will be moving the vector / matrix related functions to a separate module
#for advanced users to reduce confusion.  
#
#2011-05-04: danielhyduke@gmail.com
#**Properties that will be made private to indicate that they need to be
#updated if reactions / metabolites are modified, added, or removed from the
#Model:
#  
#
#**Properties that should not be accessed:
#   
#     _objective_coefficients, _lower_bounds, and _upper_bounds might be removed because
#they are already in the reactions (Model.reactions[i].lower_bound).  This will hinge
#on speed issues.  The same goes _S.
#
#TODO: Implement self.reactions[:]._boundary and use

# Note, there are some issues with using a gene as opposed to a protein in the
#nomenclature; this will be addressed after the first iteration of cleaning.
#
# Note, when a reaction is added to the Model it will no longer keep personal
#instances of it's Metabolites, it will reference Model.metabolites to improve
#performance.  When doing this, take care to monitor the metabolite coefficients.
#Do the same for Model.reactions[:].genes and Model.genes
#
#*********************************************************************************

#######################
#BEGIN Class Model
#
class Model(Object):
    """Model is a class for analyzing metabolic models with
    the COBRA toolbox developed in the Palsson Lab at UCSD. Make
    all of the objects (Reaction, Metabolite, ...) to make OOP easier.
    
    """
    def __setstate__(self, state):
        """Make sure all cobra.Objects in the model point to the model
        
        TODO: Make sure that the genes and metabolites referenced by
        the reactions.(genes|metabolites) point to the model's genes
        and metabolites.
        
        """
        self.__dict__.update(state)
        [[setattr(x, '_model', self)
          for x in self.__dict__[y]]
         for y in ['reactions', 'genes', 'metabolites']]

    def __init__(self, description=None):
        Object.__init__(self, description)
        self.description = self.id
        self._trimmed = False #This might get changed to a dict of 
        #gene:[reactions in which the gene participates]
        self._trimmed_genes = None #This will be integrated with _trimmed
        self._trimmed_reactions = None #as will this
        self.legacy_format = False #DEPRECATED
        #Allow the creation of an empty object which will facilitate development
        #of SBML parsers and other development issues.
        self.genes = DictList()
        self.reactions = DictList() #A list of cobra.Reactions
        self.metabolites = DictList() #A list of cobra.Metabolites
        #genes based on their ids {Gene.id: Gene}
        self.compartments = {}

        #DEPRECATED
        #Pretty much all of these things are unnecessary to use the objects and
        #interact with the optimization solvers.  It might be best to move them
        #to linear algebra modules.
        #Non list attributes below
        self._S = None #This should be regenerated by calling
        #self.construct_stoichiometric_matrix instead of directly accessed.
        #_objective_coefficients, _lower_bounds, _upper_bounds will become integrated with
        #self.reactions.  It should not be directly accessed unless the corresponding
        #property for self.reactions[i].* is also updated.
        self._objective_coefficients = None
        self._lower_bounds = None
        self._upper_bounds = None


    def __add__(self, other_model):
        """Adds two models. +

        The issue of reactions being able to exists in multiple Models now arises, the same
        for metabolites and such.  This might be a little difficult as a reaction with the
        same name / id in two models might have different coefficients for their metabolites
        due to pH and whatnot making them different reactions.

        """
        new_model = self.copy()
        new_reactions = deepcopy(other_model.reactions)
        new_model.add_reactions(new_reactions)
        new_model.id = self.id + '_' + other_model.id
        return new_model

    def __iadd__(self, other_model):
        """Adds a Model to this model +=

        The issue of reactions being able to exists in multiple Models now arises, the same
        for metabolites and such.  This might be a little difficult as a reaction with the
        same name / id in two models might have different coefficients for their metabolites
        due to pH and whatnot making them different reactions.

        """
        new_reactions = deepcopy(other_model.reactions)
        self.add_reactions(new_reactions)
        self.id = self.id + '_' + other_model.id
        return self

    def copy(self, additional_attributes=None, print_time=False):
        """Provides a partial 'deepcopy' of the Model.  All of the Metabolite, Gene,
        and Reaction objects are created anew; however, attributes not assigned in
        the __init__ function will be ignored unless contained in additional_attributes.

        additional_attributes: None or a list of attributes added to the object by
        the user that should be preserved during the copy procedure.

        print_time: Boolean used for debugging

        """
        the_copy = Model(self.id)
        if additional_attributes:
            [setattr(the_copy, x, getattr(self, x))
             for x in additional_attributes]
        the_copy.compartments = deepcopy(self.compartments)
        if print_time:
            from time import time
            start_time = time()
        the_metabolites = DictList([x.guided_copy(the_copy)
                                    for x in self.metabolites])
        if print_time:
            print 'Metabolite guided copy: %1.4f'%(time() - start_time)
            start_time = time()
        the_genes = DictList([x.guided_copy(the_copy)
                              for x in self.genes])
        if print_time:
            print 'Gene guided copy: %1.4f'%(time() - start_time)
            start_time = time()
        #TODO: See if we can use the DictList objects instead
        metabolite_dict = dict([(k.id, k)
                                 for k in the_metabolites])
        gene_dict = dict([(k.id, k)
                                 for k in the_genes])
        the_reactions = DictList([x.guided_copy(the_copy, metabolite_dict, gene_dict)
                                  for x in self.reactions])
        if print_time:
            print 'Reaction guided copy: %1.4f'%(time() - start_time)
        the_copy.reactions = the_reactions
        the_copy.genes = the_genes
        the_copy.metabolites = the_metabolites
        return the_copy
        

        
    def add_metabolites(self, metabolite_list,
                        expand_stoichiometric_matrix=False):
        """Will add a list of metabolites to the the object, if they do not
        exist and then expand the stochiometric matrix

        metabolite_list: A list of cobra.Metabolite objects

        expand_stoichimetric_matrix: Boolean.  If True and self._S is
        not None then it will add rows to self._S.  self._S must be
        created after adding reactions and metabolites to self before
        it can be expanded.  Trying to expand self._S when self only
        contains metabolites is ludacris.

        """
        if not hasattr(metabolite_list, '__iter__'):
            metabolite_list = [metabolite_list]
        #First check whether the metabolites exist in the model
        metabolite_list = [x for x in metabolite_list
                           if x.id not in self.metabolites]
        [setattr(x, '_model', self) for x in metabolite_list]
        self.metabolites += metabolite_list
        #This might already be encapsulated in update_stoichiometric matrix, but
        #may be slower.
        if self._S is not None and expand_stoichiometric_matrix:
            s_expansion = len(self.metabolites) - self._S.shape[0]
            if s_expansion > 0:
                self._S = self._S.todok()
                self._S.resize((self._S.shape[0] + s_expansion,
                               self._S.shape[1]))
                self._S = self._S.tolil()


    def _update_reaction(self, reaction):
        """Updates everything associated with the reaction.id of reaction.

        reaction: A cobra.Reaction object, or a list of these objects.


        WARNING: This function is only used after the Model has been
        converted to matrices.  It is typically faster to access the objects
        in the Model directly.  This function will eventually moved to another
        module for advanced users due to the potential for mistakes.

    
        """
        warn("WARNING: To be modified.  It will be faster to use properties of the DictList " +\
             "self.reactions to find a reaction and update the matrices " +\
             "This function is only used after the Model has been " +\
             "converted to matrices.  It is typically faster to access the objects" +\
             "in the Model directly.  This function will eventually moved to another" +\
             "module for advanced users due to the potential for mistakes.")

        if not hasattr(reaction, '__iter__'):
            reaction = [reaction]
        for the_reaction in reaction:
            if the_reaction.id not in self.reactions:
                print the_reaction.id + ' is not in the model\n'
                continue
            reaction_index = self.reactions.index(the_reaction.id)
            self.reactions[reaction_index] = the_reaction
            #zero reaction stoichiometry column
            the_column = self._S[:, reaction_index]
            for nonzero_index in the_column.nonzero()[0]:
                the_column[nonzero_index, 0] = 0
            self._lower_bounds[reaction_index] = the_reaction.lower_bound
            self._upper_bounds[reaction_index] = the_reaction.upper_bound
            self._objective_coefficients[reaction_index] = the_reaction.objective_coefficient
            self.add_metabolites(the_reaction._metabolites)
            #Make sure that the metabolites are the ones contained in the model
            the_reaction._metabolites = [self.metabolite.get_by_id(x.id)
                                        for x in the_reaction._metabolites]
            #Update the stoichiometric matrix
            metabolite_indices = map(self.metabolites.index, the_reaction._metabolites)
            for (index, metabolite_index) in enumerate(metabolite_indices):
                self._S[metabolite_index, reaction_index] = the_reaction.stoichiometric_coefficients[index]
                     
    def add_reaction(self, reaction):
        """Will add a cobra.Reaction object to the model, if
        reaction.id is not in self.reactions.

        reaction: A cobra.Reaction object

        Note: If you want to use the internal matrices/vectors immediately after
        adding a reaction you must call the update() function for the model.
         
        """
        self.add_reactions(reaction)

        
    def add_reactions(self, reaction_list, update_matrices=False):
        """Will add a cobra.Reaction object to the model, if
        reaction.id is not in self.reactions.

        reaction_list: A cobra.Reaction object or a list of them

        update_matrices:  Boolean.  If true populate / update matrices
        _S, _lower_bounds, _upper_bounds, .... Note this is slow to run
        for very large models and using this option with repeated calls
        will degrade performance.  Better to call self.update() after
        adding all reactions.

        
         If the stoichiometric matrix is initially empty then initialize a 1x1
         sparse matrix and add more rows as needed in the self.add_metabolites
         function

        """
        #Only add the reaction if one with the same ID is not already
        #present in the model.
        if type(reaction_list) not in [tuple, list, set, DictList]:
            reaction_list = [reaction_list]
        #TODO: Use the DictList properties
        reactions_in_model = set([x.id
                                  for x in reaction_list]).intersection([x.id
                                                                       for x in self.reactions])
        if len(reactions_in_model) > 0:
            print '%i of %i reaction(s) %s already in the model'%(len(reactions_in_model),
                                                          len(reaction_list), repr(reactions_in_model))
            return
        #TODO: Consider using DictList's here, just make sure that the items get appended
        #to self.metabolites, self.genes
        metabolite_dict = {}
        gene_dict = {}
        [metabolite_dict.update(dict([(y.id, y) for y in x._metabolites]))
         for x in reaction_list]
        new_metabolites = [metabolite_dict[x]
                           for x in set(metabolite_dict).difference(self.metabolites._dict)]
        if new_metabolites:
            self.add_metabolites(new_metabolites)

        [gene_dict.update(dict([(y.id, y) for y in x._genes]))
         for x in reaction_list]
        new_genes = [gene_dict[x]
                           for x in set(gene_dict).difference(self.genes._dict)]
        if new_genes:
            self.genes += DictList(new_genes)
            [setattr(x, '_model', self)
             for x in new_genes]

        #This might slow down performance
        #Make sure each reaction knows that it is now part of a Model and uses
        #metabolites in the Model and genes in the Model
        for the_reaction in reaction_list:
            the_reaction._model = self
            the_reaction._metabolites = dict([(self.metabolites.get_by_id(k.id), v)
                                             for k, v in the_reaction._metabolites.items()])
            the_reaction._genes = dict([(self.genes.get_by_id(k.id), v)
                                             for k, v in the_reaction._genes.items()])
            #Make sure the metabolites and genes are aware of the reaction
            the_reaction._update_awareness()
            
        #Add the reactions to the Model
        self.reactions += reaction_list

        if update_matrices:
            self._update_matrices(reaction_list)
    def _construct_matrices():
        """Large sparse matrices take time to construct and to read / write.
        This function allows one to let the model exists without cobra_model._S
        and then generate it at needed.
        
        """
        self._update_matrices() #This does basic construction as well.

    def _update_reaction_vectors(self):
        """regenerates the _lower_bounds, _upper_bounds,
        and _objective_coefficients vectors.

        WARNING: This function is only used after the Model has been
        converted to matrices.  It is typically faster to access the objects
        in the Model directly.  This function will eventually moved to another
        module for advanced users due to the potential for mistakes.


        """
        lower_bounds, upper_bounds = [], []
        objective_coefficients = []
        [(lower_bounds.append(x.lower_bound),
          upper_bounds.append(x.upper_bound),
          objective_coefficients.append(x.objective_coefficient))
         for x in self.reactions]
        self._lower_bounds = array(lower_bounds)
        self._upper_bounds = array(upper_bounds)
        self._objective_coefficients = array(objective_coefficients)
    def _update_metabolite_vectors(self):
        """regenerates _b and _constraint_sense

        WARNING: This function is only used after the Model has been
        converted to matrices.  It is typically faster to access the objects
        in the Model directly.  This function will eventually moved to another
        module for advanced users due to the potential for mistakes.

        """
        _b, _constraint_sense = [], []
        [(_b.append(x._bound),
          _constraint_sense.append(x._constraint_sense))
         for x in self.metabolites]
        self._b = array(_b)
        self._constraint_sense = _constraint_sense
         

    def _update_matrices(self, reaction_list=None):
        """
        reaction_list: None or a list of cobra.Reaction objects that are in
        self.reactions.  If None then reconstruct the whole matrix.

        NOTE: reaction_list is assumed to be at the end of self.reactions.

        In the future, we'll be able to use reactions from anywhere in the
        list

        
        WARNING: This function is only used after the Model has been
        converted to matrices.  It is typically faster to access the objects
        in the Model directly.  This function will eventually moved to another
        module for advanced users due to the potential for mistakes.

        """
        #Pretty much all of these things are unnecessary to use the objects and
        #interact with the optimization solvers.  It might be best to move them
        #to linear algebra modules.
        #If no reactions are present in the Model, initialize the arrays
        if not self._S or reaction_list is None:
            reaction_list = self.reactions
            self._S = dok_matrix((len(self.metabolites),
                                         len(self.reactions))) 
            self._lower_bounds = array([reaction.lower_bound
                                       for reaction in reaction_list]) 
            self._upper_bounds = array([reaction.upper_bound
                                       for reaction in reaction_list]) 
            self._objective_coefficients = array([reaction.objective_coefficient
                                                 for reaction in reaction_list])
        else: #Expand the arrays to accomodate the new reaction
            self._S = self._S.todok()
            self._S.resize((len(self.metabolites),
                           len(self.reactions)))
            lower_bounds = array([x.lower_bound
                                  for x in reaction_list])
            upper_bounds = array([x.upper_bound
                                  for x in reaction_list])
            objective_coefficients = array([x.objective_coefficient
                                            for x in reaction_list])
            self._lower_bounds = hstack((self._lower_bounds,
                                        lower_bounds))
            self._upper_bounds = hstack((self._upper_bounds,
                                        upper_bounds))
            self._objective_coefficients = hstack((self._objective_coefficients,
                                                  objective_coefficients))
        

        #Update the stoichiometric matrix
        #Using this dictionary speeds up adding reactions by orders of magnitude
        #because indexing lists is slow in python.
        #Get stats to decide how to grow self._S
        #NOTE: This was supplanted by the DictList self.reactions, however,
        ##this may slow thigns down.
        ## number_of_reactions = len(reaction_list)
        ## number_of_reactions_in_model = len(self.reactions) - number_of_reactions
        ## reaction_to_index_dict = dict(zip([x.id for x in reaction_list],
        ##                                   range(number_of_reactions_in_model,
        ##                                         number_of_reactions_in_model +\
        ##                                         number_of_reactions)))



        #Use dok format to speed up additions.
        coefficient_dictionary = {}
        #SPEED this up. This is the slow part.  Can probably use a dict to index.
        for the_reaction in reaction_list:
            reaction_index = self.reactions.index(the_reaction.id)
            for the_key, the_value in the_reaction._metabolites.items():
                coefficient_dictionary[(self.metabolites.index(the_key.id),
                                        reaction_index)] = the_value

        if not self._S.getformat() == 'dok':
            self._S = self._S.todok()
        self._S.update(coefficient_dictionary)
        self._S = self._S.tolil()

    def update(self):
        """Regenerates the stoichiometric matrix and vectors
        
        """
        self._update_matrices()

    def optimize(self, new_objective=None, objective_sense='maximize',
                 min_norm=0, the_problem=None, solver='glpk', 
                 error_reporting=None, tolerance_optimality=1e-6,
                 tolerance_feasibility=1e-6, tolerance_barrier=1e-10,
                 lp_method=0, lp_parallel=-1, copy_problem=False, relax_b=None,
                 quadratic_component=None):
        """Optimize self for self._objective_coefficients or new_objective.

        new_objective: Reaction, String, or Integer referring to a reaction in
        cobra_model.reactions to set as the objective.  Currently, only supports single
        objective coeffients.  Will expand to include mixed objectives.
        
        objective_sense: 'maximize' or 'minimize'
        
        min_norm: not implemented

        the_problem: None or a problem object for the specific solver that can be used to hot
        start the next solution.

        solver: 'glpk', 'gurobi', or 'cplex'
        
        error_reporting: None or True to disable or enable printing errors encountered
        when trying to find the optimal solution.
    
        #See cobra.flux_analysis.solvers for more info on the following parameters.  Also,
        refer to your solver's manual
        
        tolerance_optimality: Solver tolerance for optimality.
            
        tolerance_feasibility: Solver tolerance for feasibility.

        tolerance_barrier: Solver tolerance for barrier method

        lp_method: Solver method to solve the problem

        lp_parallel: Try multiple methods at once.  Supported for solver 'cplex'

        #End solver parameters
        
        copy_problem: BooleanCreate a copy of the_problem before solving.

        relax_b: Float.  Allow for error in the linear equality constraints.  Only enable if
        absolutely necessary, this can result in an inaccurate solution

        quadratic_component: None or 
        scipy.sparse.dok of dim(len(cobra_model.reactions),len(cobra_model.reactions))
        If not None:
          Solves quadratic programming problems for cobra_models of the form:
          minimize: 0.5 * x' * quadratic_component * x + cobra_model._objective_coefficients' * x
          such that,
            cobra_model._lower_bounds <= x <= cobra_model._upper_bounds
            cobra_model._S * x (cobra_model._constraint_sense) cobra_model._b
        
        """
        #TODO: change solver if other ones fail
        solver_dict = {'glpk': optimize_glpk,
                       'gurobi': optimize_gurobi,
                       'cplex': optimize_cplex}
        def solve_problem(solver_function):
            return solver_function(self, new_objective=new_objective,
                                   objective_sense=objective_sense,
                                   min_norm=min_norm,
                                   the_problem=the_problem,
                                   error_reporting=error_reporting,
                                   tolerance_optimality=tolerance_optimality,
                                   tolerance_feasibility=tolerance_feasibility,
                                   tolerance_barrier=tolerance_barrier,
                                   lp_method=lp_method, lp_parallel=lp_parallel,
                                   copy_problem=copy_problem, relax_b=relax_b,
                                   quadratic_component=quadratic_component)

        solver_function = solver_dict.pop(solver)
        the_solution = None
        try:
            the_solution = solve_problem(solver_function)
        except Exception, e:
            print e
            print '%s did not work'%solver
            for solver, solver_function in solver_dict.items():
                try:
                    print "now trying %s"%solver
                    the_solution = solve_problem(solver_function)
                    break
                except Exception, e:
                    print e
                    print '%s did not work'%solver
                    continue

        if the_solution is None:
            self.solution = None
            return(the_solution)
        else:
            self.solution = the_solution['the_solution']
            return(the_solution['the_problem'])

    def get_active_genes(self):
        """
        TODO: Move to Solution
        """
        raise Exception('cobra.Model.get_active_genes() is no longer functional.  '+\
                        'This will be associated with cobra.Model.Solution in the future.')

    
    def check_reaction_mass_balance(self, reaction_id_list=None,
                                    ignore_boundary_reactions=True):
        """Check the mass balance for reactions.

        reaction_id_list: None or a list of Ids from self.reactions.  If None
        then check all of the reactions in self.reactions

        ignore_boundary_reactions:  Boolean.  If True then ignore reactions
        starting with EX_ or DM_

        """
        if not reaction_id_list:
            the_reactions = self.reactions
        else:
            if not hasattr(reaction_id_list, '__iter__'):
                reaction_id_list = [reaction_id_list]
            reaction_indices = map(self.reactions.index,
                                   reaction_id_list)
            the_reactions = [self.reactionx[x]
                             for x in reaction_indices]
        unbalanced_reactions = []
        for the_reaction in the_reactions:
            if ignore_boundary_reactions and \
                   the_reaction.boundary:
                continue
            the_balance = the_reaction.check_mass_balance()
            if len(the_balance) > 0:
                unbalanced_reactions.append(the_balance)
        return unbalanced_reactions


    def _update_metabolite_formula(self, metabolite_name, metabolite_formula):
        """Associate metabolite_formula with all self.metabolite_names that
        match metabolite_name.

        metabolite_name: A string from self.metabolite_names

        metabolite_formula: A string specifying the chemical formula for
        metabolite_name.
        
        TODO:  This should be moved to a separate module
        
        """
        if not isinstance(the_formula, Formula):
            the_formula = Formula(the_formula)
        for the_metabolite in self.metabolites:
            if the_metabolite.name == metabolite_name:
                the_metabolite.formula = the_formula

 
    def identify_blocked_metabolites(self, account_for_bounds=False):
        """Identify metabolites that cannot be produced or consumed from other
        metabolites in the model. This will only identify root blocked metabolites,
        i.e. if the metabolite has at least one input and one output
        then assume that it is not blocked even the input or output are blocked.

        
        account_for_bounds: Boolean.  Not currently implemented.

        NOTE: To determine if a metabolite is unproducible under any conditions
        look at cobra.flux_analysis.variablity.find_blocked_reactions
        
        NOTE: Reversible reactions count as an input and output, but this can
        be problematic if one of the reaction points cannot be reached
        from the outside by any path
        
        TODO:  This can be moved to a tools module

        """
        dead_end_metabolites = {}
        for the_metabolite in self.metabolites:
            #Metabolites that freely diffuse out of the system
            #are not dead ends.
            crosses_boundary_or_reversible = False
            the_coefficients = []
            for the_reaction in the_metabolite._reaction:
                if the_reaction.boundary or the_reaction.reversibility:
                    crosses_boundary_or_reversible=True
                    break
                else:
                    the_coefficients.append(the_reaction._metabolites[the_metabolite])
            if crosses_boundary_or_reversible:
                continue
            #Reactions that are not both produced and consumed
            the_sign = sign(min(the_coefficients)) + sign(max(the_coefficients))
            if the_sign != 0:
                dead_end_metabolites[the_metabolite] = the_sign
        unconsumed_metabolites = [k for k, v in dead_end_metabolites.items()
                                  if v > 0]
        unproduced_metabolites = [k for k, v in dead_end_metabolites.items()
                                  if v < 0]
        
        return {'all': dead_end_metabolites.keys(),
                 'unconsumed': unconsumed_metabolites,
                 'unproduced': unproduced_metabolites}
    def remove_reactions(self, the_reactions):
        """
        the_reactions: instance or list of cobra.Reactions or strings of
        self.reactions[:].id.

        """
        if not hasattr(the_reactions, '__iter__') or \
               hasattr(the_reactions, 'id'):
            the_reactions = [the_reactions]
        if hasattr(the_reactions[0], 'id'):
            the_reactions = [x.id for x in the_reactions]
        reactions_to_delete = []
        for the_reaction in the_reactions:
            try:
                the_reaction = self.reactions[self.reactions.index(the_reaction)]
                the_reaction.remove_from_model(self)
            except:
                print '%s not in %s'%(the_reaction, self)
        
    #DEPRECATED def find_metabolites_without_formula(self):

    
#
#END Class Model
#####################
#dead functions


if __name__ == '__main__':
    from cPickle import load
    from time import time
    solver = 'glpk'
    test_directory = '../test/data/'
    with open(test_directory + 'salmonella.pickle') as in_file:
        cobra_model = load(in_file)
    cobra_model.optimize(solver=solver)
    new_objective=None
    objective_sense='maximize'
    min_norm=0
    the_problem=None
    solver='glpk'
    error_reporting=None
    tolerance_optimality=1e-6
    tolerance_feasibility=1e-6
    tolerance_barrier=1e-10
    lp_method=0
    lp_parallel=-1
    copy_problem=False
    relax_b=None
