#cobra.core.Metabolite.py
#######################
#BEGIN Class Metabolite
#
from copy import deepcopy
import re
from Formula import Formula
from Object import Object
class Metabolite(Object):
    """Metabolite is a class for holding information regarding
    a metabolite in a cobra.Reaction object.

    TODO: Clean up.  Allow for creation of empty metabolite
    """
    ## __slots__ = ['id', 'formula', 'name', 'compartment',
    ##              'notes', 'charge']
    ##def __setstate__(self, the_dict):
    ##    from cobra.core.Metabolite import Metabolite
    ##    Object.__setstate__(self, the_dict)
    ##    [self.__setattr__(k, v) for k, v in the_dict]
    def __init__(self, id=None,  coefficient=None, formula=None,
                 name=None, compartment=None):
        """
        id: A string.

        formula: cobra.Formula or String  of a chemical formula.  Defaults to None
        to save time in pickling and such.
        
        coefficient: None or a float.  Used when the metabolite is contained
        in a cobra.Reaction object.

        name: String.  A human readable name.
        
        compartment: None or a dictionary indicating the cellular location
        of the metabolite.  Used when in a cobra.Reaction or Model
        object
        
        """
        if coefficient is not None:
            raise Exception('Metabolites no longer have a coefficient attribute')
        Object.__init__(self, id)
        self.name = name
        if not name:
            self.name = self.id
        if isinstance(formula, str):
            formula = Formula(formula)

        self.formula = formula
        self.parse_composition()
        #self.coefficient = coefficient #This is offloaded to a container Reaction
        #because in a Model a metabolite may participate in multiple Reactions
        self.compartment = compartment
        #self.model is None or refers to the cobra.Model that
        #contains self
        self._model =  self.charge = None
        self.notes = {}
        self._reaction = set() #references to reactions that employ this metabolite
        self._constraint_sense = 'E'
        self._bound = 0.

    def __getstate__(self):
        """Remove the references to container reactions when serializing to avoid
        problems associated with recursion.
        
        """
        state = Object.__getstate__(self)
        state['_reaction'] = set()
        return state
    def parse_composition(self):
        """Breaks the chemical formula down by element.
        Useful for making sure Reactions are balanced.'
        
        """
        if isinstance(self.formula, Formula):
            self.formula.parse_composition()
        elif isinstance(self.formula, str):
            self.formula = Formula(self.formula)

    def copy(self):
        """When copying a reaction, it is necessary to deepcopy the
        components so the list references aren't carried over.

        Additionally, a copy of a reaction is no longer in a cobra.Model.

        This should be fixed with self.__deecopy__ if possible
        """
        new_metabolite = deepcopy(self)
        return new_metabolite

    def remove_from_model(self, the_model):
        """Removes the association

        the_model: cobra.Model object.  remove the reaction from this model.
        
        """
        if the_model != self._model:
            raise Exception('%s not in %s ergo it cannot be removed. (%s)'%(self,
                                                                  the_model,
                                                                  self._model))
                                                            
        self._model.metabolites.remove(self)
        self._model._metabolite_dict.pop(self.id)
        self._model = None
        for the_reaction in self._reaction:
            the_coefficient = the_reaction._metabolites[self]
            the_reaction.subtract_metabolites({self: the_coefficient})


#
#END Class Metabolite
########################
#class MetaboliteReference(Object):
#Make it access all attributes except the coefficient from the
#cobra.Metabolite that it contains.
if __name__ == '__main__':
    from cPickle import load, dump
    a = Metabolite('1','C')
    with open('a.pickle','w') as out_file:
        dump(a, out_file)
    with open('a.pickle') as in_file:
        a = load(in_file)
