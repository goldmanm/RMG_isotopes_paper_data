"""
This python script is useful in converting from concentrations of various compounds, the RMG natural unit, to the delta enrichment and site preference values perferred in isotope literature
""" 

import unittest
import pandas as pd

def getEnrichementFractionFromDelta(desired_delta, reference_ratio = 0.011115):
    """
    Calculates the enrichment fraction given a desired delta value using a reference ratio.
    """
    
    ratio = (desired_delta / 1000. + 1) * reference_ratio
    fraction = ratio / (1 + ratio)
    return fraction

def getIsotomoperInfo(species_clusters, enriched_element = 'C'):
    """
    This method takes a list of isotopolouge lists (obtained from 
    `rmgpy.tools.isotopes.cluster`), and it creates a pandas.DataFrame
    which stores the corresponding cluster index as well as how many
    atoms are enriched or non-enriched. Returns the dataframe. Intended
    to be used for analyzing cantera output where we only know the species
    description.
    
    Column descriptions:
    
    `cluster number` - the cluster number of species in species_clusters, useful for identifying isotpolouges
    `enriched atoms` - total number of enriched atoms in the molecule
    `unenriched atoms` - total number of unenriched atoms in the molecule
    `N_enriched` - total number of enriched N-type carbons in the molecule, where N is 0-4 corresponding to methane, primary, secondary, tertiary, and quatrenary
    `N_unenriched` - total number of unenriched N-type carbons in the molecule, where N is 0-4 corresponding to methane, primary, secondary, tertiary, and quatrenary
    `r_enriched` - total number of enriched radical carbons in the molecule
    `r_unenriched` - total number of unenriched radical carbons in the molecule
    """
    from rmgpy.chemkin import getSpeciesIdentifier
    def get_carbon_type(atom):
        """
        This method returns an integer representing the number of carbons bonded to this atom.
        Represents 'primary', 'secondary', 'tertiary' or 'quatrenary' atoms
        """
        assert atom.symbol == 'C'
        other_carbons = 0
        for other_atom in atom.bonds.keys():
            if other_atom.symbol == 'C':
                other_carbons += 1
        return other_carbons
    
    
    df = pd.DataFrame(columns=['cluster_number','enriched_atoms','unenriched_atoms','1_enriched','1_unenriched','2_enriched','2_unenriched','3_enriched','3_unenriched','4_enriched','4_unenriched','r_enriched','r_unenriched'])
    
    for cluster_num, cluster in enumerate(species_clusters):
        for species in cluster:
            info = {'cluster_number': cluster_num,
                    'enriched_atoms':0,
                    'unenriched_atoms':0,
                    '0_enriched':0,
                    '0_unenriched':0,
                    '1_enriched':0,
                    '1_unenriched':0,
                    '2_enriched':0,
                    '2_unenriched':0,
                    '3_enriched':0,
                    '3_unenriched':0,
                    '4_enriched':0,
                    '4_unenriched':0,
                    'r_enriched':0,
                    'r_unenriched':0,
                   }
            for atom in species.molecule[0].atoms:
                if atom.symbol == enriched_element:
                    carbon_type = get_carbon_type(atom)
                    if atom.element.isotope == -1:
                        info['unenriched_atoms'] += 1
                        info['{}_unenriched'.format(carbon_type)] += 1
                        if atom.radicalElectrons > 0:
                            info['r_unenriched'] += 1
                    else:
                        info['enriched_atoms'] += 1
                        info['{}_enriched'.format(carbon_type)] += 1
                        if atom.radicalElectrons > 0:
                            info['r_enriched'] += 1
            df = df.append(pd.Series(info, name = species.label))
    return df



def getIsotopeFraction(concentrationDict,labeledDict,unlabeledDict, reverse_labels=False):
    """
    Calculates and returns the fraction of atoms which are labeled in the set of
    isotopolougues.
    
    if reverse_labels is given, then it returns the fraction unlabeled atoms.
    """
    import numpy as np

    numerator = 0
    denominator = 0
    for spec, conc in concentrationDict.iteritems():
        numerator += conc * labeledDict[spec]
        denominator += conc * (unlabeledDict[spec] + labeledDict[spec])

    if reverse_labels:
        return 1. - numerator / denominator
    return numerator / denominator
    
def getIsotopeRatio(concentrationDict, labeledDict, unlabeledDict):
    """
    This method calculates the ratio of labeled to unlabeled atoms
    of a particular molecule for usage in deriving the \delta value.

    Inputs:
    concentrationDict - dictionary of isotopomer names and their concentrations
        at a particular point in time
    labeledDict - dictionary of isotopomers and how many enriched atoms they contain
    unlabeledDict - dictionary of isotopomers and how many non-enriched atoms 
        they contain
    """
    import numpy as np

    numerator = 0
    denominator = 0
    for spec, conc in concentrationDict.iteritems():
        numerator += conc * labeledDict[spec]
        denominator += conc * unlabeledDict[spec]

    if np.isclose(denominator, 0, atol = numerator * 1e-4):
        raise ValueError('obtained a denominator of zero when calculating isotope ratio for molecule {}'.format(spec))
    return numerator / denominator

def getDelta(concentrationDict, labeledDict, unlabeledDict, reference_ratio = 0.011115):
    """
    Calculates the delta value for the molecule.

    Inputs:
    concentrationDict - dictionary of isotopomer names and their concentrations
        at a particular point in time
    labeledDict - dictionary of isotopomers and how many enriched atoms they contain
    unlabeledDict - dictionary of isotopomers and how many non-enriched atoms 
        they contain
    reference_ratio = the reference isotope ratio for the study. The default is 
        taken from # NGS-2 Propane, Hut et al. 1987 used in Gilbert 2016.
    """
    ratio = getIsotopeRatio(concentrationDict, labeledDict, unlabeledDict)
    return (ratio / reference_ratio - 1) * 1000

def getPSIE(concentrationDict, isotopologue_info, type_1, type_2, reference_ratio = 0.011115):
    """
    returns the position specific isotope enrichment of a molecule
    
    Inputs:
    concentrationDict - dictionary of isotopomer names and their concentrations
        at a particular point in time
    isotopologue_info - DataFrame of isotopologues containing index as isotopologue with columns listing enrichments for various atom groups each ending in '_enriched_' or '_unenriched'
    type_1 - string describing the type of carbon to find the delta value of.
        examples are: 0,1,2,3,4,r, and any user-added identifiers
    type_2 - string describing the other type of carbon to find the delta value of
    reference_ratio - the reference isotope ratio for the study. The default is 
        taken from # NGS-2 Propane, Hut et al. 1987 used in Gilbert 2016.
    """
    
    delta1 = getDelta(concentrationDict, 
                      isotopologue_info['{}_enriched'.format(type_1)],
                      isotopologue_info['{}_unenriched'.format(type_1)],
                      reference_ratio = reference_ratio)
    delta2 = getDelta(concentrationDict, 
                      isotopologue_info['{}_enriched'.format(type_2)],
                      isotopologue_info['{}_unenriched'.format(type_2)],
                      reference_ratio = reference_ratio)
    return delta1-delta2
    
class TestMethods(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.natural_abundance = 0.0102
        natural_abundance = self.natural_abundance
        edge_doped_fraction = 0.0000

        self.concentrationDict = pd.Series({"C3H8(1)": (1-natural_abundance)**3,
            "CCC(2)": natural_abundance**3,
            "CCC(3)": natural_abundance**2*(1-natural_abundance),
            "CCC(4)": 2*natural_abundance**2*(1-natural_abundance),
            "CCC(5)": edge_doped_fraction + 2*natural_abundance*(1-natural_abundance)**2, #one side label
            "CCC(6)": natural_abundance*(1-natural_abundance)**2,})

        self.labeledDict = pd.Series({"C3H8(1)": 0,
            "CCC(2)": 3,
            "CCC(3)": 2,
            "CCC(4)": 2,
            "CCC(5)": 1, 
            "CCC(6)": 1})
        specific_enrichment_data = [[0,2,0,1],
                     [2,0,1,0],
                     [2,0,0,1],
                     [1,1,1,0],
                     [1,1,0,1],
                     [0,2,1,0]]
        self.specific_enrichment_data = pd.DataFrame(specific_enrichment_data, index = ["C3H8(1)", "CCC(2)", "CCC(3)", "CCC(4)", "CCC(5)", "CCC(6)"],
                                            columns = ['1_enriched','1_unenriched','2_enriched','2_unenriched'])

        self.unlabeledDict = 3 - self.labeledDict

    def test_getIsotopeRatio(self):
        abundance = getIsotopeRatio(self.concentrationDict,
                                    self.labeledDict,
                                    self.unlabeledDict)
        self.assertAlmostEqual(abundance, self.natural_abundance/(1-self.natural_abundance))

    def test_getDelta(self):
        delta = getDelta(self.concentrationDict,
                                    self.labeledDict,
                                    self.unlabeledDict,
                                    self.natural_abundance/(1-self.natural_abundance))
        self.assertAlmostEqual(delta, 0.)

    def test_getIsotopeFraction(self):
        frac = getIsotopeFraction(self.concentrationDict,
                                    self.labeledDict,
                                    self.unlabeledDict,)
        self.assertAlmostEqual(frac, self.natural_abundance)

    def test_getPSIE(self):
        psie = getPSIE(self.concentrationDict,
                       self.specific_enrichment_data,
                       1,2)
        self.assertAlmostEqual(0,psie)
if __name__ =='__main__':
    unittest.main()
