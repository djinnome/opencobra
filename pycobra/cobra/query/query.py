#cobra.query.query.py
#Will serve as a location to house the growing number of
#simple query functions attached to cobra.Model

#NOTE: Many of the find functions are gone because Reactions,
#Metabolites, and Genes are now away of each other.

import re
from numpy import array, nonzero
#####
def match_reactions(the_model, the_file_name):
    """
    Determines which items in a file containing a list of reactions (one
    per line) are contained within the_model.

    the_model: A Model object.

    the_file_name: A String of the input file.  Comment lines start with #
    and there may only be one reaction per line.
    
    """
    file_reaction_list = []
    in_file = open(the_file_name)
    for the_line in in_file.readlines():
        if the_line.startswith('#'):
            continue
        file_reaction_list.append(the_line.rstrip('\r\n'))
    in_file.close()

    matched_reactions = []
    unmatched_reactions = []

    for the_reaction in file_reaction_list:
        if the_reaction in the_model.reactions:
            matched_reactions.append(the_reaction)
        else:
            unmatched_reactions.append(the_reaction)    
    return matched_reactions, unmatched_reactions


def search_list_for_pattern(the_pattern, the_list, return_type='value'):
    """Will search a list for a regular expression or a string.

    the_pattern:  Either a string or a compiled re.

    return_type: 'value' or 'index'.  If index return a list of
    indices for the matches else return the matches.
    
    """
    if hasattr(the_pattern, 'findall'):
        hit_indices = list(array(map(lambda x: the_pattern.findall(x) != [],
                                     the_list)).nonzero()[0])
    else:
        #Walk through the list and pick up all the indices that match.
        #TODO: would it be more efficient to create an re and then map?
        #It would take less code.
        hit_indices =[]
        start_index = 0
        while start_index < len( the_list ):
            try:
                the_index = the_list.index(the_pattern, start_index)
                hit_indices.append(the_index)
                start_index = the_index + 1
            except ValueError:
                break

    if return_type == 'index':
        return hit_indices
    else:
        return map(lambda x: the_list[x], hit_indices)


def print_reactions_involving_metabolite(cobra_model, the_metabolites):
    """Update to allow for multiple metabolite search

    cobra_model: A cobra.Model object

    the_metabolites: A list of cobra.Metabolites or metabolite ids that are in
    cobra_metabolites.

    #TODO: Move this to the Metabolite class

    """
    if hasattr(the_metabolites, 'id'):
        the_metabolites = [the_metabolites]
    elif not hasattr(the_metabolites, '__iter__'):
        the_metabolites = [the_metabolites]
    if not hasattr(the_metabolites[0], 'id'):
        the_metabolites = [cobra_model.metabolites[cobra_model.metabolites.index(x)]
                           for x in the_metabolites]
        
    for the_metabolite in the_metabolties:
        for the_reaction in the_metabolite._reaction:
            print the_reaction.reaction
 
         
def get_translation_reactions(cobra_model, genes_of_interest):
    """Find the translation elongation reactions for a set of genes
    in a cobra model.  Related to ME-model extensions

    cobra_model:  A cobra.Model object.

    genes_of_interest:  A list of genes from cobra_model.genes.
    
    """
    gene_translation_reactions = defaultdict(list)
    for the_reaction in cobra_model.reactions:
        if 'translation_elongation' in the_reaction:
            for the_gene in genes_of_interest:
                if the_gene in the_reaction:
                    gene_translation_reactions[the_gene].append(the_reaction)
                    continue
    return gene_translation_reactions


if __name__ == '__main__':
    from cPickle import load
    from time import time
    solver = 'glpk'
    test_directory = '../test/data/'
    with open(test_directory + 'salmonella.pickle') as in_file:
        cobra_model = load(in_file)

    #TODO: Add in tests for each function
    print 'Need to add in tests for %s'%repr(['match_reactions',
                                              'print_reactions_involving_metabolite',
                                              'search_list_for_pattern'])
