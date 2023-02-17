################################################################################
##
## This library is free software; you can redistribute it and/or
## modify it under the terms of the GNU Lesser General Public
## License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## 
## This library is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## Lesser General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public
## License along with this library; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
##
## (C) Copyrights Dr. Michel F. Sanner and TSRI 2016
##
################################################################################

########################################################################
#
# Date: 2015 Authors: Stefano Forli
#
#    forli@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Stefano Forli and TSRI 2015
#
#########################################################################
#
# $Header#
#
# $Id: Kamaji.py,v 1.2.2.1 2017/07/26 22:03:40 annao Exp $
#
import sys
import pybel
ob = pybel.ob
import openbabel as ob
from operator import itemgetter
import numpy
from . import openBabelInterface as obi

class Kamaji:
    """ Analyze an OBMol and  a PDB and returns splitted items ready to be
        processed as desired
        

        TODO: NMR EXMPLE 2l6x.pdb
             - covalent residue!
             
        REPAIR MECHANISM
              ask the user to re-attach chains that 
              are found to be detached or 
              residues that are considered ligands...

        ======= WORKING ================
        
        GUESSING GROUPS
            glycosilation
            lipid acid
            lipid (generic)

        KNOWN CATALYTIC METALS
            Mg, Mn, Fe, Co, Ni, Cu, Zn

        KNOWN SALTS
            Na, Cl, Ca, K, Li, Cl

        KNOWN ADDITIVES
            GOL, SDS (page), LDA
        
        KNOWN COFACTORS
            Heme, ADP/ATP, GDP/GTP, NADP/NADHI,
            Pyridoxal phosphate, FAD/FADm,

        KNOWN X-RAY ADDITIVES
            Glycerol, SDS, tris-buffer, lauryl dimethylamine-N-Oxide, DMSO,
            sulphate, phosphate, acetate, HEPES, ammonium, ethilen-glycol,
            PEG, ethanol, B-mercapthoethanole, formic acid, pyruvic acid, 
            morpholinoethilensulfonate, bicarbonate

        FILTERING
            standard residues
            defined modified residues
            . . . . . . . . . . 
            waters
            ions
            cofactors
            additive
            undefined modifiers
            ligands 
    """

    def __init__(self, settings={}):
        """
        # XXX ADD: HEADER, COMPND, SOURCE
        # XXX manage peptoids 2SEM
        """
        self._printToDo()
        # any chain with this many residues will be considered a potential ligand
        self._default_minrescount = 10 # decapeptide is a ligand
        self.settings = settings
        self.initTools()
        self.graph = {}

    def _printToDo(self):
        x = """

            - tag residue before missing as PROBLEMATIC
        """

    def initTools(self):
        """ initialize tools and variables to be used later"""
        self.initKnownTypes()
        # element table
        self.eTable = ob.OBElementTable()
        

    def initKnownTypes(self):
        """ initialize names and info for known aminoacids,
            nucleic acid residues, cofactors, metals, salts
            and other common PDB ligands
        """
        self.structureTypes = { }
               #  'protein': {}, 'nucleic': {}, 'unnown': {}...
        self.structureGroups = {} # clusters, glycosilation sites

        # resId, numAtoms (max: bound residue count + 1, if terminal)
        self._knownAminoAcids = { 'ALA': 5, 'ARG': 12, 'ASN':9, 'ASP':9, 'ASX':9,
              'CYS':7, 'GLU':10, 'GLN':10, 'GLX':10, 'GLY':5,
              'HIS':11, 'ILE':9, 'LEU':9, 'LYS':10, 'MET':9,
              'PHE':12, 'PRO':8, 'SER':7, 'THR':8, 'TRP':15,
              'TYR':13, 'VAL':8 }

        self._knownDnaBases = { 'DA' : 22, 'DC' : 20, 'DG': 23, 'DT': 21 }
                              #  'A', 'C', 'G', 'T' }

        self._knownRnaBases = {'A':23, 'C':21, 'G':24, 'U':21} 

        self._knownWaters = [ 'HOH', 'WAT', 'DOD' ]

        self._knownMetals = [ 'Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 
                              'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',
                              'Ni', 'Cu', 'Zn', 'Ga', 'Ge',  # 'As', 'Se', 
                              'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 
                              'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
                              'Sb', 'Cs', 'Ba', 'Hf', 'Ta', 'W', 
                              'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 
                              'Pb', 'Bi', 'Po', 'Fr', 'Ra', 'Rf', 
                              'Db', 'Sg', 'Bh', 'Hs', 'La', 'Ce', 'Pr', 
                              'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
                              'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac', 'Th',
                              'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 
                              'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr' ]
         
        #self._knownAlkaliMetal = [ 'Li', 'Na', 'K', 'Rb', 'Cs', 'Fr' ]
        #self._knownAlkaliEarthMetals = [ 'Be', 'Mg', 'K', 'Rb', 'Cs', 'Fr' ]


        self._knownCatalyticMetals = [ 'Mg', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn' ] 

        self._autodockSupportedElements = [ 'H', 'C', 'N', 'O', 'F', 'Mg',  
                'P', 'S', 'Cl', 'Ca', 'Mn', 'Fe', 'Zn', 'Br',  'I' ]

        self._knownSalts = [ 'Na', 'Cl', 'Ca', 'K', 'Li', 'Cl', 'N' ]

        self._knownCofactors = {
            'HEM': { 'name': 'Protoporphyrin IX containing Fe (Heme)',
                    'smi':  ('CC1=C(CCC(O)=O)C2=Cc3c(CCC(O)=O)c(C)c4C='
                               'C5C(C)=C(C=C)C6=[N]5[Fe]5([N]2=C1C=c1c'
                               '(C=C)c(C)c(=C6)n51)n34'),
                    'setup'  : ['metal-charge'] },

            'HEC': { 'name': 'Heme C',
                    'smi':  ('CC1=C(CCC(O)=O)C2=Cc3c(CCC(O)=O)c(C)c4C='
                               'C5C(C)=C(C=C)C6=[N]5[Fe]5([N]2=C1C=c1c'
                               '(C=C)c(C)c(=C6)n51)n34'),
                    'setup'  : ['metal-charge'] },


            'ADP' : {'name' : 'Adenosine di-phosphate (ADP)',
                     'smi':  ('Nc1ncnc2n(cnc12)[C@@H]1O[C@H]'
                              '(CO[P@@](O)(=O)OP(O)(O)=O)[C@@H]'
                              '(O)[C@H]1O'),
                    'setup' : []},

            'ATP' : {'name' : 'Adenosine tri-phosphate (ATP)',
                     'smi': ('Nc1ncnc2n(cnc12)[C@@H]1O[C@H](CO[P@](O)'
                             '(=O)O[P@@](O)(=O)OP(O)(O)=O)[C@@H](O)'
                             '[C@H]1O'),
                     'setup' : []},

            'NAD' : {'name': 'Nicotinamide adenine dinucleotide (NADP)',
                     'smi':  ('NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H]'
                              '(CO[P@]([O-])(=O)O[P@](O)(=O)OC'
                              '[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2'
                              'cnc3c(N)ncnc23)[C@@H](O)[C@H]1O'),
                    'setup' : []},

            'NDP' : {'name': 'Dihydro-nicotinamide-adenine-dinucleotide (NADP-ph)',
                     'smi':  ('NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](CO[P@@]'
                              '(O)(=O)O[P@](O)(=O)OC[C@H]2O[C@H]([C@H]'
                              '(OP(O)(O)=O)[C@@H]2O)n2cnc3c(N)ncnc23)'
                              '[C@@H](O)[C@H]1O'),
                    'setup' : []},

            'NAI' : {'name': 'Nicotninamide adenin dinucleotide (NADPH+)',
                    'smi' :  ('NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H]'
                              '(CO[P@@](O)(=O)O[P@](O)(=O)OC'
                              '[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2'
                              'cnc3c(N)ncnc23)[C@@H](O)[C@H]1O'),
                    'setup' : []},

            'GDP' : { 'name' : "Guanosine-5'-diphosphate",
                      'smi' : ( 'Nc1nc2n(cnc2c(=O)[nH]1)[C@@H]1O[C@H]'
                            '(CO[P@@](O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]1O'),
                    'setup' : []},

            'GTP' : { 'name' : "Guanosine-5'-triphosphate",
                     'smi' : ('Nc1nc2n(cnc2c(=O)[nH]1)[C@@H]1O'
                              '[C@H](CO[P@@](O)(=O)O[P@@](O)(=O)'
                              'OP(O)(O)=O)[C@@H](O)[C@H]1O'),
                    'setup' : []},

            'PLP' : { 'name' : "Pyridoxal-5'-phosphate (vitamin B6 phosphate)",
                      'smi'   : 'Cc1ncc(COP(O)(O)=O)c(C=O)c1O',
                      'setup' : []},

            'FAD' : { 'name': "Flavin-adenine dinucleotide",
                      'smi':  ('Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C[C@H]'
                               '(O)[C@H](O)[C@H](O)CO[P@](O)(=O)O[P@@]'
                               '(O)(=O)OC[C@H]3O[C@H]([C@H](O)[C@@H]3O)'
                               'n3cnc4c(N)ncnc34)c2cc1C'),
                    'setup' : []},

            'FMN' : { 'name': "Flavin-adenine mononucleotide",
                      'smi':  ('Cc1cc2nc3c(nc(=O)[nH]c3=O)n'
                                '(C[C@H](O)[C@H](O)[C@H](O)'
                                'COP(O)(O)=O)c2cc1C'),
                    'setup' : []},


                            }

        self._UNDEFINEDMODIF = ['NAG', 'PLP' ] # these are attached to parts of the protein!

        # XXX TOFIX        self._knownResidues = self._knownAminoAcids.keys() + self._knownDnaBases.keys() + self._knownRnaBases.keys()

        # smiles
        self._knownXrayAdditives = { 
            'SO4' : { 'commonName' : 'sulfate' , 'smi' : '[O-]S([O-])(=O)=O'},
            'GOL' : { 'commonName' : 'glycerol' , 'smi': 'C(C(CO)O)O'}, 
            'EDO' : { 'commonName' : 'ethylen glycol' , 'smi': 'OCCO'},
            #'NAG' : { 'commonName' 
            'SDS' : { 'commonName' : 'SDS' , 'smi': 'CCCCCCCCCCCCOS(O)(=O)=O'} ,
            'TRS' : { 'commonName' : 'Tris-buffer' , 'smi': '[NH3+]C(CO)(CO)CO'},
            'LDA' : { 'commonName' : 'Lauryl dimethylamine-N-oxide' , 
                        'smi': 'CCCCCCCCCCCC[N+](C)(C)[O-]'},
            'DMS' : { 'commonName' : 'DMSO' , 'smi': 'CS(C)=O'},
            #'DM2' : { 'commonName' : 'DMSO (low quality)' , 'smi': 'CS(C)O'},
            'PO4' : { 'commonName' : 'phosphate', 'smi': '[O-]P([O-])([O-])=O'},
            'ACT' : { 'commonName' : 'acetate' , 'smi': 'CC([O-])=O'},
            'EPE' : { 'commonName' : 'HEPES' , 'smi': 'OCCN1CCN(CCS(O)(=O)=O)CC1'},
            'NH4' : { 'commonName' : 'ammonium' , 'smi': '[NH4+]'},
            'PEG' : { 'commonName' : 'PEG (poly[Ethylen-glycole])' , 'smi': 'OCCOCCO'},
            'EOH' : { 'commonName' : 'ethanol' , 'smi': 'CCO'},
            'BME' : { 'commonName' : 'Beta-mercapto-ethanol' , 'smi': 'OCCS'},
            'FMT' : { 'commonName' : 'formic acid' , 'smi': 'OC=O'},
            'PYR' : { 'commonName' : 'pyruvic acid' , 'smi': 'CC(=O)C(O)=O'},
            # too many of these, a SMARTS should be better...
            # { 'commonName' : 'tartaric acid', 'smi': 'O[C@H]([C@@H](O)C(O)=O)C(O)=O'},
            'MES' : { 'commonName' : '2-morpholin-4-ium-4-ylethanesulfonate' , 
                        'smi': '[O-]S(=O)(=O)CC[NH+]1CCOCC1' },
            'BCT' : { 'commonName' : 'bicarbonate' , 'smi': 'OC([O-])=O'},
                                    }
        # smarts
        self._guessingGroups = {
            # glycane pattern ( O_inring *RINGBOND* C_inring *RINGBOND*
            #                   nonRing_C_ *NONRINGBOND*-NonRingOxygen )
            'glycosilation':'[#8r5,#8r6]@[#6r5,#6r6]~[C!r]~[#8!r]',
            'lipids' : [ # <- list because order is important
                  # phospholipid (phosphatidyl-choline group)
                  #{'phosphatidic acid': 'C(OC=O)(COC=O)COP(~O)(~O)~O'},
                  # carboxylic-lipid pattern with at least 4 carbons attached to at least one hydrogen
                  {'phospholipid': 'C(OC~O)(COC~O)COP(~O)(~O)~O'}, # NO BOND ORDER
                  #{'lipid acid': '[C;!R;h1]~[C;!R;h1]~[C;!R;h1]~[C;!R;h1]~[CX3](=[OX1])O'},
                  {'fatty acid': '[C;!R;h2]~[C;!R;h2]~[C;!R;h2]~[C;!R;h2]~[CX3](~[OX1])O'},
                  # {'fatty acid SMI': 'CCCCC(=O)O'},
                  # aliphatic chain with 6-any carbons attached to an oxygen
                  #{'lipid (generic)': '[Ch1;!R]~[Ch1;!R]~[Ch1;!R]~[Ch1;!R]~[CX3][#8].[!r]'},
                  #{'lipid (generic)': '([Ch1;!R]~[Ch1;!R]~[Ch1;!R]~[Ch1;!R]~[CX3][#8]),([*!r])'},
                  #{ 'lipid (generic)' :  '[C]~[C]~[C]~[C]~[C]~[#8]'},
                ]}

    def sprint(self, text, dest='out'):
        """ stdout stderr"""
        if dest == 'out':   
            func = sys.stdout
        elif dest == 'err':
            func = sys.stderr
        func.write(text)

    def setMolecule(self, mol, pdbinfo={}, perceiveBondOrder=True):
        """ use a pre-existing multistate PDB object as molecule"""
        self.mol = mol
        self.mol.PerceiveBondOrders()
        self.pdbInfo = pdbinfo
        self.parseStructure()
        

    def parseStructure(self):
        """ perform structure parse and populate the structureTypes
            dictionary.

            keys main keys ( '*' means optionally found):
                protein     :  { OBMol.resIdx  : { name: str, num: int, chain: str } }

                dna         :  { OBMol.resIdx  : { name: str, num: int, chain: str } }

                rna         :  { OBMol.resIdx  : { name: str, num: int, chain: str } }


                water       :  {   OBMol.resIdx  : { name: str, num: int, chain: str, 
                                    *modRes: 'deuterated' } }

                catalyticMetal  : { OBMol.resIdx  : { name: str, num: int, chain: str,
                                    'info' : { 'ad_supported' : bool }
                                        } 
                                
                salt        : { OBMol.resIdx  : { name: str, num: int, chain: str,
                                    'info' : { 'ad_supported' : bool }
                                        } 

                genericMetal: { OBMol.resIdx  : { name: str, num: int, chain: str,
                                    'info' : { 'ad_supported' : bool }
                                        } 
                                
                cofactor    :   { OBMol.resIdx  : { name: str, num: int, chain: str,
                                    'info' : { name : str, smi: SMILES string,
                                          setup : [ *'metal-charge' ]
                                    }
                                    } 
                                  }                     
                additives   :   { OBMol.resIdx  : { name: str, num: int, chain: str,
                                            accuracy : int, score : int}
                                            }

                glycosylation   :   { groupId : { OBMol.resIdx  : { name: str, num: int, chain: str}} 

                
                modifier    : {  OBMol.resIdx  : { name: str, num: int, chain: str,
                                                mw : float, size : str } }
                                

                ligand      :  { short_pepdtide   : { shortId :  
                                            { OBMol.resIdx :  { name : str, 
                                                                num : int,
                                                                chain : str }}
                                            },
                                
                
                                 { short_nucleic    : { shortId :
                                            { OBMol.resIdx :  { name : str, 
                                                                num : int,
                                                                chain : str }}
                                            },
                                
                             'generic' :  { OBMol.resIdx  : { name: str, num: int, chain: str,
                                                                size: string, mw : float}
                                            }
                            'sugar' :   { OBMol.resIdx  : { name: str, num: int, chain: str,

                            'lipid' : { OBMol.resIdx  : { name: str, num: int, chain: str,

            secondary keys:
                # residues with missing atoms
                _incomplete :  { name : str, num : int, 'chain' : str, 
                        fixed : bool, missingAtomsCount : int }

                _modRes : [ {resInfo : [chain,name,num],
                             description : str, # PDB header description
                             standard : str, # corresponding unmodified residue
                            }
                
        """
        self.structureTypes = {}
        # filter residues (protein, dna, rna, water, other)
        self.filterResidues()
        # search for short chains that could likely be peptide ligands
        self.scanShortChains()
        # filter all other residues
        self.filterOther()

    def scanShortChains(self):
        """ scan the graph structure for short chains
            (< 10 AA) that could be likely peptide/nucleic
            ligands or modified ones
        """
        for kw, tag in [('protein','short_peptide'), ('nucleic', 'short_nucleic')]:
            if not kw in self.structureTypes:
                continue
            chains = {}
            remove = []
            for resId, info in self.structureTypes[kw].items():
                rchain = info['chain']
                reslist = chains.setdefault(rchain, [])
                reslist.append(resId)

            count = 0
            for c,reslist in chains.items():
                if len(reslist) <= self._default_minrescount:
                    ligands = self.structureTypes.setdefault('ligand', {})
                    shortList = ligands.setdefault(tag, {})
                    currShort = {}
                    for resId in reslist:
                        currShort[resId] = self.structureTypes[kw].pop(resId)
                    shortList[count] = {'data' : currShort  }
                    shortList[count]['info'] = self.getShortChainProperty(currShort)
                    count += 1

    def getShortChainProperty(self, shortChain):
        """ short chain properties (peptides/nucleic) properties """
        indices = []
        for resId, data in shortChain.items():
            indices.append((resId, data['num']))
        sorted(indices, key=itemgetter(1))
        seq = "-".join( [str(x[1]) for x in indices] )
        return { 'sequence' : seq }
            


    def filterLigands(self):
        """ scan for potential ligands and known classes"""
        if not 'unknown' in self.structureTypes.keys():
            return
        for resId in self.structureTypes['unknown'].keys():
            lipid = self._isLipid(resId)
            if lipid:
                ligand = self.structureTypes.setdefault('ligand', {})
                ligand.setdefault('lipid',{})[resId] = self.structureTypes['unknown'].pop(resId)
                ligand['lipid'][resId]['type'] = lipid
                ligand['lipid'][resId]['info'] = self.ligandProperties(resId)
            else:
                sugar = self._isGlycoGroup(resId)
                if self._isGlycoGroup(resId): # repeat for free sugars
                    ligand = self.structureTypes.setdefault('ligand', {})
                    ligand.setdefault('sugar',{})[resId] = self.structureTypes['unknown'].pop(resId)
                    ligand['sugar'][resId]['info'] = self.ligandProperties(resId)
                    ligand['sugar'][resId]['type'] = 'sugar-containing'
                else:
                    ligand = self.structureTypes.setdefault('ligand', {})
                    ligand.setdefault('generic',{})[resId] = self.structureTypes['unknown'].pop(resId)
                    ligand['generic'][resId]['info'] = self.ligandProperties(resId)
                    ligand['generic'][resId]['type'] = 'generic ligand'
                    
        #
        if not 'ligand' in self.structureTypes:
            return
        # DEBUG:
        if False:
            for kind, data in self.structureTypes['ligand'].items():
                for k in data.keys():
                    #res = self.mol.GetResidue(k)
                    #mol = self.residueToMol(res)
                    rObj = self.mol.GetResidue(k)
                    btotal = 0
                    # DEBUG
                    for a in ob.OBResidueAtomIter(rObj):
                        bc =  len([x for x in ob.OBAtomBondIter(a) ])
                        btotal += bc
                    print "RESIDUE TOTAL BONDS",k, btotal


    def ligandProperties(self, resId):
        """ classifies a generic ligand"""
        # XXX ADD is_macrocycle
        res = self.mol.GetResidue(resId)
        mol = self.residueToMol(res)
        mw = mol.GetMolWt()
        size = 'drug-size'
        if mw < 310:
            size = 'fragment'
        elif mw > 600:
            size = 'large'
        return {'size': size, 'mw': mw }

    def debugwrite(self, fname, reslist=[]):
        c = 0
        for r in reslist:
            res = self.mol.GetResidue(r)
            mol = self.residueToMol(res)
            pybel.Molecule(mol).write('pdb', 'DEBUG_%s_%s.pdb' %(fname, c), overwrite=1)
            c+=1
            
    def transformResidues(self):
        """ """
        transformations = { 'MSE': { 'from': { 'smi':'C[Se]CC[C@H](N)C(O)=O',
                                              'name':'selenomethionine'},
                                     'to' :  {'smi': 'CSCC[C@H](N)C(O)=O',
                                              'name' : 'MET'}, },
                                              
                            'SEP': { 'from': {'smi':'N[C@@H](COP(O)(O)=O)C(O)=O',
                                            'name': 'phosphoserine'},
                                     'to' :  {'smi': 'N[C@@H](CO)C(O)=O',
                                              'name' : 'SER'}, }, 

                            'CAS': { 'from': {'smi':'C[As](C)SC[C@H](N)C(O)=O',
                                              'name': 'S-(dimethylarsenic)cysteine'},
                                     'to' :  {'smi': 'N[C@@H](CS)C(O)=O',
                                              'name' :'CYS'}, },

                            'TPO' : {'from': {'smi' :'C[C@@H](OP(O)(O)=O)[C@H](N)C(O)=O',
                                              'name' : 'phosphothreonine'},
                                     'to': {'smi': 'C[C@@H](O)[C@H](N)C(O)=O',
                                              'name' : 'THR'}, },
                          }
        
        # TYS: sulfotyrosine
        # phosphotyrosyne
        # XXX Arsenic
        # MET C-S-C (1.8) C-Se-C (1.9)

    def filterOther(self):
        """ filter unknown residues into:
            metal, other    
        """
        # check if the residue is in the Golden Hundred
        # (list of all the top represented ligands
        # in the PDB)
        # filter waters
        self.filterWaters()
        # filter known ions (mono-atomic residues)
        self.filterIons()
        # filter known co-factors
        self.filterCofactors()
        # filter known solvent/additives to be removed
        self.filterAdditives()
        # at this point, easy groups are all catched
        # remaining ones will be clustered
        self.clusterGroups()
        # anything that is attached to the protein
        # i.e.: glycosilation, non-std residues
        self.filterModifiers()
        # find ligands (this should be called as last)
        self.filterLigands()

        # clean classes that have been emptied..
        for k,i in self.structureTypes.items():
            if len(i) == 0:
                del self.structureTypes[k]

    def filterWaters(self):
        """ filter waters
            NOTE check here if it is deuterated?
        """
        if not 'unknown' in self.structureTypes:
            return
        for resIdx in self.structureTypes['unknown'].keys():
            if self._isWater(resIdx):
                water = self.structureTypes.setdefault('water', {} )
                water[resIdx] = self.structureTypes['unknown'].pop(resIdx)
                if water[resIdx]['name'] == 'DOD':
                    water[resIdx]['modRes'] = 'deuterated'

    def _isWater(self, resIdx):
        """ """
        ## is called water
        #if not res.GetName in self._knownWaters:
        #    return False
        # it looks like water
        res = self.mol.GetResidue(resIdx)
        atomElements = [ x.GetAtomicNum() for x in ob.OBResidueAtomIter(res) ]
        if not (1 <= len(atomElements) <= 3):
            return False
        # it smells like water
        if 8 in atomElements:
            return True
        return False

    def filterIons(self):
        """ filter mono-atomic ions and classifies them as 
                catalytic
                salt
                metal
        """
        # XXX "[ this function will be updated with statistics on teh PDB distro ]"
        if not 'unknown' in self.structureTypes.keys():
            return
        salt = {}
        metal = {}
        catalytic = {}
        unknown = {}
        for resIdx in self.structureTypes['unknown'].keys():
            res = self.mol.GetResidue(resIdx)
            atoms = [ x for x in ob.OBResidueAtomIter(res) ]
            if len(atoms) == 1:
                atom = atoms[0]
                #print "ELEMENT>", self.eTable.GetSymbol(atom.GetAtomicNum())
                if self._isCatalyticMetal(atom): # catalytic
                    target = catalytic
                # XXX THIS WILL BE CHANGED WITH THE AVaiLABILITY 
                # of distributions from the PDB
                elif self._isSalt(atom): # salt
                    target = salt
                elif self._isGenericMetal(atom): # generic metal
                    target = metal
                else:
                    target = unknown 
                target[resIdx] = self.structureTypes['unknown'].pop(resIdx)
                target[resIdx]['info'] = {'ad_supported': self._isADSupported(atom)  }
        if len(catalytic):
            self.structureTypes['catalyticMetal'] = catalytic
        if len(salt):
            self.structureTypes['salt'] = salt
        if len(metal):
            self.structureTypes['genericMetal'] = metal
        if len(unknown):
            self.structureTypes['unknown'] = unknown
            
    def _isModifiedRes(self, chain, name, num):
        """ check if a residue is in the set of 
            MODRES residues defined in the header
            of the PDB
        """
        if not 'modRes' in self.pdbInfo:
            return
        info = self.pdbInfo['modRes']
        if chain in info.keys():
            for modres in info[chain]:
                if not ( modres['resName'] == name ):
                    continue
                if not ( modres['resNum'] == num ):
                    continue
                return  {'resInfo' : [chain, name, num],
                        'description' : modres['modDescription'],
                        'standard' : modres['stdResName'] }
        return False

    def filterModifiers(self):
        """ filter anything that is attached to the protein
            as a 'modifier' but was not defined in the 
            MODRES keys: glycosilation, very-non-std residues

            NOTE the code assumes that glycosylations and generic
            modifiers are separate groups with no bonds in between.
            This could cause problems in case of fancy glyco-modifiers?
        """
        if not 'unknown' in self.structureTypes:
            return
        
        # cache all target atom indices
        target = []
        for kw in ['protein', 'nucleic']:
            if kw in self.structureTypes.keys():
                target += self.structureTypes[kw]
        targetIdx = set( self._getAtomIdx(target) )
        glyco = []
        glyco_common = {}
        # check if residues are attached to the target
        for resId in self.structureTypes['unknown'].keys():
            resAtomIdx = set( self._getAtomsBoundToRes(resId))
            common = targetIdx & resAtomIdx
            if len(common):
                if self._isGlycoGroup(resId): # glyco modifier
                    glyco.append(resId) 
                    glyco_common[resId] = common
                else: # generic modifier # XXX this should go into filterKnownModifiers?
                    modifier = self.structureTypes.setdefault('modifier', {})
                    modifier[resId] = self.structureTypes['unknown'].pop(resId)
                    mw = self.getResidueFP(resId)[1]
                    if mw < 300:
                        size = 'fragment'
                    elif 300 < mw < 600:
                        size = 'ligand'
                    elif mw > 600:
                        size = 'large ligand'
                    attachedRes = self.indexToResName(common)[0]
                    modifier[resId].update( {'mw': mw, 'size':size, 
                        'attachedTo': attachedRes} )
        if len(glyco):
            # cluster glyco-groups
            self.structureGroups['_glycosylationGroup'] = self._completeGlycoChains(glyco,
                    self.structureGroups['clusters'])
            self.structureTypes['glycosylation'] = {}
            for gId, group in enumerate(self.structureGroups['_glycosylationGroup']):
                attached = None
                curr = self.structureTypes['glycosylation'][gId] = {'residues' : {},
                    'attachedTo' : ''}
                for x in group:
                    curr['residues'][x] = self.structureTypes['unknown'].pop(x)
                    if x in glyco_common:
                        attached = glyco_common[x]
                if attached:
                    curr['attachedTo'] = self.indexToResName(attached)[0]



        # XXX this is [ self.filterKnownModifiers() ]
        
    def _completeGlycoChains(self, glyco, clusters):
        """ group all connected glycosilation residues"""
        glycoMafia = []
        for resId in glyco:
            pairs = self._walkClustGraph( clusters,
                start=resId,  visited=[] )
            glycoMafia.append(pairs)
        return glycoMafia

    def indexToResName(self, indices):
        """ generate ResNameResNum from indices"""
        resNames = []
        for i in indices:
            resObj = self.mol.GetAtom(i).GetResidue()
            resNames.append( '%s%s' % (resObj.GetName(), resObj.GetNum() ) )
        return resNames


    def clusterGroups(self):
        """ generate graph of connected residues
            (i.e. poly-saccharides, n-peptides...)
        """
        if not 'unknown' in self.structureTypes:
            return
        clusters = {}
        for r1 in self.structureTypes['unknown'].keys():
            clusters[r1] = []
            for r2 in self.structureTypes['unknown'].keys():
                if r1 == r2: continue
                if self._areResConnected(r1, r2):
                    clusters[r1].append(r2)
        self.structureGroups['clusters'] = clusters

    def _walkClustGraph(self, graph,  start, visited = [] ):
        """ recursive function to walk graph 
            nodes to find connected clusters 
        """
        if not start in visited:
            visited += [start]
        for children in graph[start]:
            if not children in visited:
                visited.append(children)
                new = self._walkClustGraph( graph, children, visited)
                for n in new:
                    if not n in visited:
                        visited.append(n)
        return visited


    def _areResConnected(self, res1, res2):
        """ res1, res2:
            check if two residues are connected by
            at least a bond between two atoms
        """
        atomsBoundWithRes1 = set( self._getAtomsBoundToRes(res1) )
        atomsRes2 = set( self._getAtomIdx([res2]) )
        return len( atomsBoundWithRes1 & atomsRes2 ) > 0

    def _getAtomIdx(self, resList):
        """ convert a list of residues to the list of indices
            of all atoms contained in the residues
        """
        idx = []
        for resId in resList:
            #res = self.mol.GetResidue(resId)
            idx += [ x.GetIdx() for x in self.getResidueAtoms(resId) ]
        return idx    

    def _getAtomsBoundToRes(self, resIdx):
        """ find all atoms with which atoms in the 
            residue establish bonds
        """
        bound = []
        atomList = self.getResidueAtoms(resIdx)
        for a in atomList:
            bound += [ b.GetBeginAtomIdx() for b in ob.OBAtomBondIter(a) ]
            bound += [ b.GetEndAtomIdx() for b in ob.OBAtomBondIter(a) ]
        return list(set(bound))
        
    def getResidueAtoms(self, resIdx):
        """ return all atoms in a residue"""
        res = self.mol.GetResidue(resIdx)
        return [ x for x in ob.OBResidueAtomIter(res) ]

    def _isGlycoGroup(self, resId):
        """ check if a residue is 
            a glycosilating group
        """
        pattern = self._guessingGroups['glycosilation']
        result = self._SmartsMatcher(self.mol.GetResidue(resId), pattern)
        if result:
            return 'carbohydrate-like'
        return False

    def _isLipid(self, resIdx):
        """ check if a residue is a lipid"""
        res = self.mol.GetResidue(resIdx)
        patterns = self._guessingGroups['lipids']

        for pair in patterns:
            name, p = pair.items()[0]
            result = self._SmartsMatcher(res, p)
            if result:
                return name
        return False

    def _SmartsMatcher(self, mol, pattern):
        # XXX convert this to OBSmarts
        """ return the matching indices for
            pattern in mol.
            Mol can be a OBMol or a OBResidue
        """
        if isinstance( mol, ob.OBResidue ):
            mol = self.residueToMol(mol)
        #matcher = ob.OBSmartsPattern()
        #matcher.Init(pattern)
        #return matcher.Match(mol)
        matcher = pybel.Smarts(pattern)
        return matcher.findall(pybel.Molecule(mol))


    def filterAdditives(self):
        """ filter for known solvent/additives
        
            this should rely on the top 100 most represented
            ligands in the PDB
        """
        if not 'unknown' in self.structureTypes:
            return
        cutoffFP = 0.8
        cutoffMWratio = 0.5
        for resId, info in self.structureTypes['unknown'].items():
            res = self.mol.GetResidue(resId)
            resObj = self.residueToMol(res)
            fp, mw, smi = self.getResidueFP(resId)
            best = -1
            bestName = None
            for additiveId, data in self._knownXrayAdditives.items():
                name = data['commonName']
                smi = data['smi']
                score, tanimoto = self.getSampleScore(smi, mw, fp)
                if score > best:
                    best = score
                    bestId = additiveId
                    bestMatch = tanimoto
            if best == -1:
                return
            accuracy = None
            if best >=0.4:
                accuracy = 'low'
                if best >=0.85:
                    accuracy = 'high'
                #print " [ accuracy[%s], best match %s : %2.3f ] " % (accuracy.upper(),
                #    bestId, best), bestMatch
            if accuracy:
                additives = self.structureTypes.setdefault('additives', {})
                additives[resId] = self.structureTypes['unknown'].pop(resId)
                additives[resId].update( {'tanimoto' : bestMatch, 'score':best} )
                additives[resId].update( self._knownXrayAdditives[bestId] )


    def getResidueFP(self, resId, simple=False):
        """ calculate residue fingerprint, including MW
            
            if simple requested, bond order is ignored
            (useful to deal with horrible quality of some
            PDF structures)
        
        """
        res = self.mol.GetResidue(resId)
        newmol =  self.residueToMol(res, forcesingle=simple)
        residueFragment = pybel.Molecule( newmol )
        smi = residueFragment.__str__()
        fp = residueFragment.calcfp()
        mw = residueFragment.molwt
        return fp, mw, smi

    def getSampleScore(self, smi, sampleMW, sampleFP):
        """ generate fingerprint similarly to getResidueFP
            for a SMI sample
        """
        mol = pybel.readstring('smi', smi)
        molFP = mol.calcfp()
        molMW = mol.molwt
        tanimoto = molFP | sampleFP 
        try:
            mwratio =  ( sampleMW / (sampleMW - abs(molMW-sampleMW)))
        except:
            mwratio = 0.
        score = tanimoto /mwratio
        #print "SMI>>", smi, score, tanimoto, mwratio
        return score, tanimoto

    def residueToMol(self, residue, forcesingle=False):
        """ create a separate molecule from a 
            residue within a molecule
        """
        mol = ob.OBMol()
        table = {}
        bonds = []
        for a in ob.OBResidueAtomIter(residue):
            new = mol.NewAtom()
            new.SetAtomicNum( a.GetAtomicNum() )
            new.SetVector( a.GetVector() )
            table[a.GetIdx()] = new.GetIdx()
        miss = 0
        for a in ob.OBResidueAtomIter(residue):
            for b in [ x for x in ob.OBAtomBondIter(a) ]:
                begin = b.GetBeginAtomIdx()
                end = b.GetEndAtomIdx()
                order = b.GetBondOrder()
                if forcesingle:
                    order = 1
                try:
                    mol.AddBond( table[begin], table[end], order)
                except:
                    miss+=1
                    pass
        return mol




    def blendmol(self, mol):
        """ """
        mol.write('pdb', filename='x.pdb', overwrite=True)
        mol = pybel.readfile('pdb', 'x.pdb').next()
        mol.OBMol.ConnectTheDots()
        return mol

    def _isInChemClass(self, atom, chemClass):
        """ generic function to check if the element of
            an atom belongs to a chemical class,
            (list of element symbols)
        """
        eNumber = atom.GetAtomicNum()
        eSymbol = self.eTable.GetSymbol(eNumber)
        return eSymbol in chemClass

    def _isCatalyticMetal(self, atom):
        """ check if an atom is in the list of metals
            known to have a catalyric role
        """
        klass =  self._knownCatalyticMetals
        if self._isInChemClass(atom, klass):
            # XXX TODO check the coordination geometry
            return True
        return False

    def _isGenericMetal(self, atom):
        """ check (bool) if an atom is a metal"""
        return self._isInChemClass(atom, self._knownMetals)
        
    def _isSalt(self, atom):
        """ check if an atom is a salt ion"""
        return self._isInChemClass(atom, self._knownSalts)

    def _isCofactor(self, name):
        """ check if a residue is a known cofactor
            (name, for now, SMARTS later?)
        """
        found = self._knownCofactors.get(name, None)
        if not found == None:
            return found
        else:
            # XXX SMARTS here
            return None

    def _isADSupported(self, atom):
        """ check if is is an AutoDock supported element"""
        return self._isInChemClass(atom, self._autodockSupportedElements)


    def filterCofactors(self):
        """ filter and identifies known co-factors """
        if not 'unknown' in self.structureTypes:
            return
        for resId, info in self.structureTypes['unknown'].items():
            found = self._isCofactor(info['name'])
            if not found == None:
                cofactor = self.structureTypes.setdefault('cofactor', {})
                #print "Found this co-factor", found
                cofactor[resId] = self.structureTypes['unknown'].pop(resId)
                cofactor[resId]['info'] = found

    def filterResidues(self):
        """ subdivide residues by type and populate
            the structureType dictionary:
                protein
                nucleic
                other (anything else)

            for each recognized residue from biopolymers
            the number of atoms found is checked with expected
            values plus a tolerance number:
                - aminoacids : +/-1 for terminal residue
                - nucleic  : +/-4 for 3'/5' terminal (PO3-O)
        """
        protein = {}
        dna = {}
        rna = {}
        other = {}
        incomplete = []
        modResList = []
        allRes = self.getResidue()
        for resId in allRes:
            rObj = self.mol.GetResidue(resId)
            btotal = 0
            # DEBUG
            for a in ob.OBResidueAtomIter(rObj):
                bc =  len([x for x in ob.OBAtomBondIter(a) ])
                btotal += bc
            #print "RESIDUE TOTAL BONDS", btotal
            # /DEBUG
            rchain, rname, rnum = self.getResidueInfo(resId)
            atomCount = self.countResAtoms(resId)
            # protein  (standard)
            if rname in self._knownAminoAcids:
                # protein[resId] = {'n': rname, '#': rnum, 'c':rchain}
                protein[resId] = {'name': rname, 'num': rnum, 'chain':rchain}
                if (self._knownAminoAcids[rname] - atomCount) > 1:
                    incomplete.append(  {'name': rname, 
                        'num': rnum, 'chain':rchain, 'fixed' : False,
                        'missingAtomsCount' : (self._knownAminoAcids[rname] - atomCount -1) } )
            # dna
            elif rname in self._knownDnaBases:
                dna[resId] = {'name': rname, 'num': rnum, 'chain':rchain}
                
                if (self._knownDnaBases[rname] - atomCount) > 4:
                    incomplete.append(  {'name': rname, 
                        'num': rnum, 'chain':rchain,'fixed' : False, 
                        'missingAtomsCount' : (self._knownDnaBases[rname] - atomCount -4) } )
            # rna
            elif rname in self._knownRnaBases:
                rna[resId] = {'name': rname, 'num': rnum, 'chain':rchain}
                if (self._knownRnaBases[rname] - atomCount) > 4:
                    incomplete.append(  {'name': rname, 
                        'num': rnum, 'chain':rchain,'fixed' : False, 
                        'missingAtomsCount' : (self._knownRnaBases[rname] - atomCount -4) } )
            else: 
                modres = self._isModifiedRes(rchain, rname, rnum)
                if modres: # modified res
                    modResList.append(resId)
                    if modres['standard'] in self._knownAminoAcids:
                        protein[resId] = {'name': rname, 'num': rnum, 
                            'chain':rchain, 'modres': modres}
                    elif modres['standard'] in self._knownDnaBases:
                        dna[resId] = {'name': rname, 'num': rnum, 
                            'chain':rchain, 'modres': modres}
                    elif modres['standard'] in self._knownRnaBases:
                        rna[resId] = {'name': rname, 'num': rnum, 
                            'chain':rchain, 'modres': modres}
                else: # other
                    other[resId] = {'name': rname, 'num': rnum, 'chain':rchain}
        if len(protein):
            self.structureTypes['protein'] = protein
        if len(dna):
            self.structureTypes['dna'] = dna
        if len(rna):
            self.structureTypes['rna'] = rna
        if len(modResList):
            self.structureTypes['_modRes'] = modResList
        if len(incomplete):
            self.structureTypes['_incomplete'] = incomplete
        if len(other):
            self.structureTypes['unknown'] = other

    def getResidueInfo(self, resIdx=None, res=None):
        """ retrieve info about a residue"""
        if not resIdx == None:
            res = self.mol.GetResidue(resIdx)
        name = res.GetName().strip()
        num = res.GetNum()
        chain = res.GetChain().strip()
        return chain, name, num

    def countResAtoms(self, resIdx=None, res=None, heavyonly=True):
        """ count number of atoms in a residue
        """
        if not resIdx == None:
            res = self.mol.GetResidue(resIdx)
        c = 0
        for a in ob.OBResidueAtomIter(res):
            if a.GetAtomicNum == 1 and not heavyonly:
                continue
            c+= 1
        return c


    def getResidue(self, name=None, num=None, chain=None):
        """ retrieve the residue(s) matching the 
            name-num-chain combination)
        """
        found = []
        for res in ob.OBResidueIter(self.mol):
            rname = res.GetName()
            rnum = res.GetNum()
            rchain = res.GetChain()
            if (rname == name) or (name == None):
                if (rnum == num) or (num == None):
                    if (rchain == chain) or (chain == None):
                        found.append(res.GetIdx())
        return found


    def fixWater(self):
        """ fix water molecules:
                - add hydrogens to lone oxygens
                - convert deuterium to hydrogen
        """
        # XXX to be written
        pass



class MultiStatePDB:
    """ 
        deals with MODEL and alternate residue locations...
    
    """
    def __init__(self, text):
        """ """
        self.text = text
        self.initVariables()
        self.parseLines()
        self.parseHeader()
        self.parseAtoms() 

    def parseHeader(self):
        """ 
          PDB FILE DEFINED PROPERTIES
            Herarchical structure tree
            Modified residues
            Experimental data
            Missing atoms
            Symmetry
            Biounit
            Hetero list
        """
        self.initPDBHeaderParser()
        self.parsePdbInfo()

    def initPDBHeaderParser(self):
        """ initialize mapping betwen PDB info and parsers"""

        self.typeToFunc = { 'missingRes' : self._parseMissingResidues,
               'missingAtoms' : self._parseMissingAtoms,
               'symmetry' : self._parseSymmetry,
               'biounit' : self._parseBioUnit,
               'hetList': self._parseHetList,
               'expData' : self._parseExpData,
               'modRes' : self._parseModRes,
               'ssBond' : self._parseSSBond,
              }

        self.remarkToType = { 'REMARK 465' : 'missingRes',
                    'REMARK 470' : 'missingAtoms',
                    'REMARK 290' : 'symmetry',
                    'REMARK 350' : 'biounit',
                    'HET   '    : 'hetList',
                    'EXPDTA'     : 'expData',
                    'MODRES' : 'modRes',
                    'SSBOND': 'ssBond',
                  }

        self.typeToRemark = { 'missingRes' : 'REMARK 465',
                   'missingAtoms' :'REMARK 470',
                   'symmetry' : 'REMARK 290',
                   'biounit' : 'REMARK 350',
                   'hetList' : 'HET   ',
                   'expData': 'EXPDTA',
                   'modRes': 'MODRES',
                    'ssBond': 'SSBOND',
                  }

    def parsePdbInfo(self):
        """ parse all remarks in the PDB"""
        self.pdbInfo = {}
        for kw, v in self.remarkToType.items():
            self.pdbInfo[v] = []
        for l in self.getHeader():
            value = None
            if l[0:6] == 'REMARK':
                value = l[0:10]
        for l in self.getHeader():
            value = None
            if l[0:6] == 'REMARK':
                value = l[0:10]
                if not value in self.remarkToType.keys(): 
                    value = None
                    continue
            elif l[0:6] == 'HET   ':
                value = l[0:6]
            elif l[0:6] =='EXPDTA':
                value = l[0:6]
            elif l[0:6] == 'MODRES':
                value = l[0:6]
            #elif l.startswith('ATOM') or l.startswith('HETATM'):
            #    pass
            if not value == None:
                kw = self.remarkToType[value]
                buff = self.pdbInfo[kw].append(l)
        
        for kw, data in self.pdbInfo.items():   
            self.pdbInfo[kw] = self.typeToFunc[kw](data)

    def _parseMissingResidues(self, data):
        """ parse the missing residues entry"""
        new = {}
        kw = self.typeToRemark['missingRes']
        pattern = 'REMARK 465   M RES C SSSEQI'
        inside = False
        for l in data:
            if inside:
                raw = l.split(kw, 1)[1]
                res, chain, seqId = raw.split() # this must be always len() == 3?!?
                if not chain in new.keys():
                    new[chain] = []
                new[chain].append( (res, seqId) )
            if pattern in l:
                inside = True
        return new

    def _parseExpData(self, data):
        """
            return the experimental method
            NOTE: some intelligence would be useful here...
        """
        kw = 'EXPDTA'
        method = 'other'
        known_methods = { 'nmr': 'nmr',
                    'electron crystallography' : 'ec',
                    'x-ray diffraction': 'xray',
                    }
        raw = []

        for l in data:
            raw.append( l.split(kw, 1)[1].strip() )
        rawString = " | ".join(raw)
        for k,v in known_methods.items():
            if k in rawString.lower():
                method = v
                break
        return  { 'raw': raw, 'method': v }

    def _parseMissingAtoms(self, data):
        """ """ 
        return data

    def _parseSSBond(self, data):
        """ ,
        [ SOURCE: http://www.wwpdb.org/documentation/format33/sect6.html#SSBOND ]
        COLUMNS        DATA  TYPE     FIELD            DEFINITION
        --------------------------------------------------------------------------------
         1 -  6        Record name    "SSBOND"
         8 - 10        Integer        serNum           Serial number.
        12 - 14        LString(3)     "CYS"            Residue name.
        16             Character      chainID1         Chain identifier.
        18 - 21        Integer        seqNum1          Residue sequence number.
        22             AChar          icode1           Insertion code.
        26 - 28        LString(3)     "CYS"            Residue name.
        30             Character      chainID2         Chain identifier.
        32 - 35        Integer        seqNum2          Residue sequence number.
        36             AChar          icode2           Insertion code.
        60 - 65        SymOP          sym1             Symmetry operator for residue 1.
        67 - 72        SymOP          sym2             Symmetry operator for residue 2.
        74 - 78        Real(5.2)      Length           Disulfide bond distance
        """
        new = {}
        return new
        # see _parseModRes for adaptation




    def _parseModRes(self, data):
        """ processes modified residues
        [ SOURCE: http://www.wwpdb.org/documentation/format23/sect3.html#MODRES ]

        COLUMNS    DATA TYPE        FIELD         DEFINITION
        ----------------------------------------------------
         1 - 6     Record name      "MODRES"
         8 - 11    IDcode           idCode     ID code of this entry.
        13 - 15    Residue name     resName    Residue name used in this entry.
        17         Character        chainID    Chain identifier.
        19 - 22    Integer          seqNum     Sequence number.
        23         AChar            iCode      Insertion code.
        25 - 27    Residue name     stdRes     Standard residue name.
        30 - 70    String           comment    Description of the residue
        """
        kw = self.typeToRemark['modRes']
        new = {}
        for l in data:
            pdbId = l[7:11].strip()
            resName = l[12:15].strip()
            chainId = l[16] #.strip()
            resNum = int(l[18:22].strip())
            iCode = l[22].strip()
            stdResName = l[24:27].strip()
            modDescription = l[29:70].strip()
            chain = new.setdefault(chainId, [])
            chain.append( { 'resName': resName,
                                   'chainId': chainId,
                                   'resNum' : resNum,
                                   'iCode'  : iCode,
                                   'stdResName': stdResName,
                                   'modDescription': modDescription,
                                 }
                                )
        return new

    def _parseSymmetry(self, data):
        # key = 'REMARK 290'
        # terminated by 'REMARK 290 REMARK: NULL'
        return data

    def _parseBioUnit(self, data):
        """ Test with Stout's protease structure!"""
        return data
        data = new

        """
        REMARK 350 BIOMOLECULE: ?
        REMARK 350 APPLY THE FOLLOWING TO CHAINS: ?, ?...
        REMARK 350   BIOMT1   N  N.NNNNNN  N.NNNNNN  N.NNNNNN        N.NNNNN
        REMARK 350   BIOMT2   N  N.NNNNNN  N.NNNNNN  N.NNNNNN        N.NNNNN
        REMARK 350   BIOMT3   N  N.NNNNNN  N.NNNNNN  N.NNNNNN        N.NNNNN
        """

    def _parseHetList(self, data):
        """ 
        [ SOURCE:http://www.wwpdb.org/documentation/format23/sect4.html ]

        COLUMNS     DATA TYPE     FIELD         DEFINITION
        ------------------------------------------------------
         1 -  6     Record name   "HET      "
         8 - 10     LString(3)    hetID         Het identifier, right-justified.
        13          Character     ChainID       Chain identifier.
        14 - 17     Integer       seqNum        Sequence number.
        18          AChar         iCode         Insertion code.
        21 - 25     Integer       numHetAtoms   Number of HETATM records for the
                                                group present in the entry.
        31 - 70     String        text          Text describing Het group.
        """
        new = {}
        kw = self.typeToRemark['hetList']
        for l in data:
            hetID = l[7:10].strip()
            chainID = l[12].strip()
            seqNum = l[13:17].strip()
            iCode = l[17]
            numHetAtoms = l[20:25].strip()
            text = l[30:70].strip()
            if not chainID in new.keys():
                new[chainID] = []
            new[chainID].append( (hetID,seqNum) )
        return new




    def initVariables(self):
        """ """
        self._kwslicing =[0,6]
        self.kw = { 'ATOM  ':None, 
                    'HETATM':None,
                    'MODEL ':None,
                    'ENDMDL':None,
                    'TER   ':None,
                   #'MASTER ':None, # XXX Neither is supported! 
                   #'END   ':None,  # XXX
                    'CONECT':None, 
                  }

        self.coordKw = ('ATOM  ', 'HETATM')

        self.modelSet = {}

        # list of specific altmodes per residue
        self.altResidueMode = {}

        # list of altresidues found per model
        self.altResidues = {}

        # default altLocation mode
        self.altMode = 'A'

        self.conect = []
        self.header = []
        # TODO implement the CONECT to find if multi-chains are linked
        # TODO also look at LINK kw that provides *EXPLICIT* links between chains (see 1o7d)
        """ 
        [ source: http://www.wwpdb.org/documentation/format33/sect10.html#CONECT ]
        COLUMNS       DATA  TYPE      FIELD        DEFINITION
        -------------------------------------------------------------------------
         1 -  6        Record name    "CONECT"
         7 - 11        Integer        serial       Atom  serial number
        12 - 16        Integer        serial       Serial number of bonded atom  _
        17 - 21        Integer        serial       Serial number of bonded atom   |
        22 - 26        Integer        serial       Serial number of bonded atom   +-- optional
        27 - 31        Integer        serial       Serial number of bonded atom  _| 
        """

    def parseLines(self):
        """ does the dirty work parsing header and models """
        inside = False
        currModel = []
        i,j = self._kwslicing
        # model state of the molecule
        self.currentModel = 0
        modelId = 0
        for idx, raw in enumerate(self.text):
            l = raw[i:j]
            if not l in self.kw:
                self.header.append(idx)
                continue
            if l == 'ENDMDL':
                inside = False
                currModel.append(idx)
                self.modelSet[modelId] = currModel
                modelId += 1
                currModel = []
                #print "FLUSHING", currModel, modelId
            elif l == 'MODEL ':
                inside = True
                currModel = [idx]
            elif l == 'CONECT':
                pass
                #bonds = self._parseConect(raw)
                #self.conect.append(bonds)
            else:
                currModel.append(idx)
        if len(currModel):
            self.modelSet[modelId] = currModel


    def setState(self, model=None, altMode=None, altResidueMode={}):
        """ define the state of the multistructure"""
        if model:
            self.currentModel = model
        self.setAltMode(altMode)
        self.setResidueAltMode(altResidueMode)
        #self.generateOBMol()


    def setAltMode(self, altMode='a'):
        """ define the default alternate conformation mode for
            residues that do not have a specific setting in self.altResidueMode
        """
        if not altMode:
            return
        self.altMode = altMode.upper()

    def setResidueAltMode(self, altResidueMode={}):
        """ set the residue-specific alternate residue mode"""
        if altResidueMode == {}:
            self.altresidueMode = {}
            return
        for k,v in altResidueMode:
            self.altResidueMode[k] = v

    def generateOBMol(self, raw):
        """ create the OBMolecule from the current model"""
        self.mol = ob.OBMol()
        conv = ob.OBConversion()
        conv.SetInFormat('pdb')
        conv.ReadString(self.mol, raw)
        return self.mol

    def getHeader(self):
        """ return PDB header information """
        return [self.text[i] for i in self.header]


    def parseAtoms(self):
        """ process all atoms found in each model"""
        #self.graph = []
        # XXX
        kw = self.coordKw # ('ATOM   ', 'HETATM ')
        i,j = self._kwslicing
        for mIdx, model in self.modelSet.items():
            self.currentModel = mIdx
            self.altResidues[mIdx] = {}
            structure = {}
            # scan all atom lines in the model
            for lineIdx in model:
                raw = self.text[lineIdx]
                if raw[i:j] in kw:
                    self.addAtomToStructure(structure, lineIdx)
            # cleanup the model
            cleanStructure = self.cleanStructure(structure)
            # store the model
            self.modelSet[mIdx] = cleanStructure
        # reset model to the begin
        self.currentModel = 0
        
        # DEBUG
        if 0:
            for ch, re in self.modelSet[mIdx].items():
                print "CH[%s]" % ch
                for nu, inf in re.items():
                    for k,v in inf.items():
                        print "\t", k, v
                

#x ADD CHECK FOR COVALENT LIGANDS!
#IT SHOULD BE A FLAG SO WHEN MIGRATING SUGARS, WILL BE CONSIDERED COVALENT LIGANDS

    @property
    def modelCount(self):
        """ return the number of models found in the structure"""
        return len(self.modelSet)

    @property
    def getAltResidues(self):
        """ return the number of models found in the structure"""
        return sorted(self.altResidues[self.currentModel].keys())

    @ property
    def chains(self):
        return sorted(self.modelSet[self.currentModel].keys())

    def addAtomToStructure(self, structure, lineIdx):
        """ take care of creating/adding atoms, book-keeping..."""
        chainInfo, resInfo, atInfo = self.atomInfo(lineIdx)
        resName = resInfo['name']
        resNum = resInfo['num']
        chain = chainInfo['name']
        #resKey = "%s%s" % (resName, resNum)
        # it is possible that an alternate location defines a different residue!
        # therefore the key is going to be the sequence number and not the NameNumber
        resKey = "%d" % (resNum)
        residues = structure.setdefault(chain, {})
        atoms = residues.setdefault(resKey, [] )
        atoms.append(atInfo)


    def cleanStructure(self, structure):
        """process each residue for alt states"""
        for chain, residue in structure.items():
            for res, atoms  in residue.items():
                newResAtoms = self.compactAtoms(atoms)
                structure[chain][res] = newResAtoms
        return structure


    def compactAtoms(self, atoms):
        """ generate the structure to create requested alt conformations
        """
        altResidues = self.altResidues[self.currentModel]
        common = []
        altAtoms = {}
        for a in atoms:
            if a['isAlt']:
                altLabel = a['alt']
                # register residue in list of alt residues 
                altList = altResidues.setdefault(a['resNum'], [])
                if not altLabel in altList:
                    altList.append(altLabel)
                # register atoms to the specific alt group
                altGroup = altAtoms.setdefault(altLabel, [])
                altGroup.append(a)
            else:
                common.append(a)
        return {'common': common, 'alternate': altAtoms}


    def getStructure(self, chains=[]):
        """ generate OBMol molecule structure (optionally, 
            containing only selected chains
        """
        raw = "".join( self.getStructureState(chains=chains))
        #open('DEBUG_XX.pdb','w').write(raw) 
        return self.generateOBMol(raw)
        
        

    def getStructureState(self, chains=[]): #, altMode='A', model='all', altResidueMode=[]):
        """ create a representation o the structure using 
        
            the current model (self.model), alternate locations
            per residue (self.altResidueMode) and default (self.altMode)
        """
        out = []
        currModel = self.modelSet[self.currentModel]
        for modelChain, residues in currModel.items():
            if chains and not modelChain in chains:
                continue
            for resNum in sorted(residues.keys()):
                resAtoms = residues[resNum]
                resState = self.getResidueState(resNum, resAtoms)
                out += resState
        return [self.text[i] for i in out]

    def sortedResidues(self, residues):
        """ return a tuple with """
        out = []
        keys = residues.keys()
        nameNum = [ (x[0:3], int(x[3:])) for x in keys ]
        nameNum.sort(key=itemgetter(1))
        for name, num in nameNum:
            k = "%s%d" % (name, num)
            out.append( (k, residues[k]) )
        return out
        
    def getResidueState(self, resNum, resAtoms):
        """ extract residue atoms in the specified altModes (default or specific)"""
        # common atoms with no alternate conformations
        out = [ x['lineIdx'] for x in resAtoms['common'] ]

        altMode = self.altResidueMode.setdefault(resNum, self.altMode)
        if altMode in resAtoms['alternate']:
            for a in resAtoms['alternate'][altMode]:
                    out.append(a['lineIdx'])
        return out


    def atomInfo(self, lineIdx):
        """ extract information from the ATOM/HETATM line"""
        # atom information
        s = self.text[lineIdx]
        atName = s[12:16].strip()
        atNum = int(s[6:11])
        alt = s[16].strip()
        isAlt = bool(len(alt))
        #occupancy = float(s[54:60].strip())
        #temp = float(s[60:66].strip())
        #element = s[76:78]
        #coord = map(float, ( s[30:38], s[38:46], s[46:54]) )

        # remove the alt location label
        self.text[lineIdx] = s[:16] + " "+ s[17:]

        # chain
        chain = s[21]
        #segment = s[72:76]
        chainInfo = { 'name': chain}
        # residue
        resNum = int(s[22:26].strip())
        resName = s[17:20].strip()
        resId = "%s:%s%d" % (chain, resName, resNum)
        resInfo = {'name': resName, 'num': resNum, 'resId':resId}
        # atom
        atInfo = { 'isAlt' : isAlt, 'atomName': atName, 
            'atomNum': atNum, 
            'lineIdx' : lineIdx, 'resId': resId, 'alt':alt, 'resNum':resNum,
            #'segment':segment,
            #'element':element,
            #'occupancy': occupancy, 'temp':temp, 'coord':coord,
            }
        return chainInfo, resInfo, atInfo

                
    def _parseConect(self, string):
        """ do we need the re-numbering book-keeping after alt identification?"""
        pass




class KamajiInterface:

    def __init__(self, mol, perceiveBondOrders=True):
        """ """
        self.mol = mol
        self.classes = {}

        maxLenType = 0

        #self._bondOrder = numpy.zeros( (mol._ag.numBonds(),), dtype = 'int')
        self._bondOrder = {}
        self._types = numpy.zeros((mol._ag.numAtoms(),), dtype='S5')
        self._hbTypes = numpy.zeros((mol._ag.numAtoms(),), dtype='uint8')

        #self._types = np.chararray(shape=(len(mol._ag),), dtype='|S5')
        self._bCounter = 0
        self._aCounter = 0

        self.hasProtein = False
        self.hasDna = False
        self.hasRna = False
        self.hasLigand = False
        self.hasWater = False
        self.hasOther = False
        self.kamaji = Kamaji()
        
        self.proteinTrees = []
        self.dnaTrees = []
        self.rnaTrees = []
        self.ligandTrees = []
        self.waterTrees = []
        self.otherTrees = []
        self._other = [ 'salt', 'genericMetal', 'additives' ]
        self._modifiers = ['glycosylation', 'modifier' ]
        for chId in numpy.unique(self.mol._ag.getChids()):
            chain = self.mol.select('chain "%c"'%chId)
            obmol = obi.ProdyToOBMol(chain, title='%s: chain %s'%(mol.name, chId))
            if perceiveBondOrders:
                obmol.PerceiveBondOrders()
            self.kamaji.setMolecule(obmol)
            # get atom type 
            for atm in ob.OBMolAtomIter(obmol):
                self._types[self._aCounter] = atm.GetType()
                maxLenType = max(maxLenType, len(atm.GetType()))
                self._hbTypes[self._aCounter] = 1 * atm.IsHbondAcceptor() + 2 * atm.IsHbondDonor()
                self._aCounter += 1
            # get bond order
            for bond in ob.OBMolBondIter(obmol):
                i = bond.GetBeginAtomIdx()-1
                j = bond.GetEndAtomIdx()-1
                self._bondOrder['%d %d'%(i,j)] = bond.GetBondOrder()
                #self._bondOrder[self._bCounter] = bond.GetBondOrder()
                self._bCounter += 1
            self.classes[chId] = {}
            self.kamaji.setMolecule(obmol) #, pdbinfo = self.mol.pdbInfo)
            self.classes[chain] = {}

            kst = self.kamaji.structureTypes
            tree = None
            # XXX this code should be modified to allow multiple types per chain (i.e., protein DNA)?
            # XXX ask David
            # XXX if so, the following if/elif should be converted into a for loop + break
            if 'protein' in kst:
                self.hasProtein = True
                protein = self.classes[chain]['protein'] = self.kamaji.structureTypes['protein']
                tree = self.generateBiopolymerTree(chId, protein, 'protein')
            elif 'dna' in kst:
                self.hasDna = True
                dna = self.classes[chain]['dna'] = self.kamaji.structureTypes['dna']
                tree = self.generateBiopolymerTree(chId, dna, 'dna')
            elif 'rna' in kst:
                self.hasRna = True
                rna = self.classes[chain]['rna'] = self.kamaji.structureTypes['rna']
                tree = self.generateBiopolymerTree(chId, rna, 'rna')
            if tree: # 
                self.parseModifiers(tree)
                self.parseCofactors(tree)
                self.parseWater(tree)
            if 'ligand' in kst:
                self.hasLigand = True
                ligand = self.classes[chain]['ligand'] = self.kamaji.structureTypes['ligand']
                self.generateLigandTree(chId, ligand)
            self.parseOther(chId)
        self.compactTrees()
        #print 'ATOMTYPES', self._types
        #print 'bondORders', self._bondOrder
        #print 'HBDONORS', self._hbTypes
        #for i in self._hbTypes:
        #    print i
        #print "MAXTYPE", maxLenType

    def getBondOrders(self):
        """ return the bond order"""
        return self._bondOrder

    def getAtomTypes(self):
        """ return atom types"""
        return self._types

    def getHBTypes(self):
        """ return h-bond types
            0   :   no hydrogen bond
            1   :   hydrogen bond acceptor
            2   :   hydrogen bond donor
            3   :   hydrogen bond acceptor/donor
        """
        return self._hbTypes

    def compactTrees(self):
        """ create the root (All) for each tree"""
        treeList = [ 'proteinTrees', 'dnaTrees', 
            'rnaTrees', 'ligandTrees', 'otherTrees' ]
        for name in treeList:
            tree = getattr(self, name)
            if not len(tree):
                continue
            newTrees = []
            master = KMolTree('All (%s)' % self.mol.name)
            for chain in tree:
                master.children.append(chain)
                master.fragType = chain.fragType
                master.fragTypeCount += chain.fragTypeCount 
            master.info.append( ('total chains', '%d' % (len(master.children))))
            newTrees.append(master)
            setattr(self, name, newTrees)
            #if name == 'otherTrees':
            #    print "GOT OTHER TREES"
            #    print "MASTER", master.getSelection()
            #    print "XXX", master.getSelectionString()


    def parseOther(self, chId):
        """ """
        other = ['additives', 'genericMetal', 'salt', 'unknown' ]
        kst = self.kamaji.structureTypes
        root = None
        types = []
        for o in set(other) & set(kst.keys()):
            self.hasOther = True
            if root == None:
                root = KTree('Chain %s' %chId, chId)
                root.fragType = 'other'
                root.fragTypeCount = 1
            types.append(o)
            branch = root.add(o)
            branch.fragType = o
            branch.fragTypeCount = 1
            root.fragTypeCount += 1
            for i,d in kst[o].items():
                if o == 'additives':
                    name = d['commonName']
                else:
                    name = d['name']
                n = branch.add(name)
                n.selection = [d['num']]
                n.fragType = o #'other'
                n.fragTypeCount = 1
                n.info = self.getOtherInfo(d)
            branch.info.append( ('%s found' % o, len(branch.children)) )
        if not root == None:
            root.info.append(('other types found', '%s' % (', '.join(types) )))
            self.otherTrees.append(root)

    def getOtherInfo(self, d):
        """ """
        info = []
        if 'commonName' in d:
            name = d['commonName']
        else:
            name = d['name']
        info.append( ('name', name))
        info.append( ('residue name', d['name']))
        if 'score' in d:
            info.append( ('identification accuracy', '%2.1f%%' % (d['score']*100)) )
        if 'smi' in d:
            resinfo = self.smiToInfo(d['smi'])
            info.append( ('MW', '%2.1f' % (resinfo['mw'])))
            info.append( ('brute formula', '%s' % (resinfo['formula'])))
        return info

    def parseModifiers(self, tree):
        """ parse residue modifiers and extract info about 
            the residue they are attached to:
                glycosilating groups
                general modifiers
        """
        kst = self.kamaji.structureTypes
        if not 'modifier' in kst and not 'glycosylation' in kst:
            return
        n = tree.add('Modified residues (%d)')
        n.fragType = 'modified'
        n.fragTypeCount += 1
        #tree.fragType = 'modified'
        tree.fragTypeCount += 1
        glycoResidues = []
        if 'glycosylation' in kst:
            tree.fragTypeCount += 1
            groups = kst['glycosylation']
            glycoNode = n.add('Glycosylation (%d)' % len(groups))
            glycoNode.fragType = 'glyco'
            glycoNode.fragTypeCount = 1
            for g, data in groups.items():
                resList = data['residues']
                currNode = glycoNode.add('Site %d->%s' % (g+1, data['attachedTo'] ))
                glycoResidues.append( data['attachedTo'])
                currNode.selection = [ int(x['num']) for x in resList.values() ]
                currNode.info = self.getGlycosilationInfo(data) # glycosylation info goes here
                currNode.fragType = 'glyco'
                currNode.fragTypeCount = 1
            glycoNode.info = [ ('glycosilation sites', '%d' % len(glycoNode.children)) ]
            n.info.append(  ('glycosilation sites', '%d' % len(glycoNode.children)) )
            glycoNode.info.append( ('glycosilated residues', '%s' % ', '.join(glycoResidues)))
        if 'modifier' in kst:
            tree.fragTypeCount += 1
            modData = kst['modifier']
            modNode = n.add('Covalent/modifiers')
            modNode.fragType = 'covalent'
            modNode.fragTypeCount = 1
            for mId, data in modData.items():
                name = "%s->%s" % (data['name'],data['attachedTo'])
                currNode = modNode.add(name)
                currNode.selection = [ data['num'] ]
                currNode.info = self.getModifierInfo(data)
                currNode.fragType = 'covalent'
                currNode.fragTypeCount = 1
            n.info.append(  ('covalent modifiers', '%d' % len(modNode.children)) )
            modNode.info.append(  ('covalent modifiers', '%d' % len(modNode.children)) )
        n.name = n.name % (len(n.children))
        tree.info.append( ('modified residues', '%d' % (len(n.children))))

    def smiToInfo(self, smiles):
        """ extract info from a SMILES string"""
        m = pybel.readstring('smi', smiles)
        mw = m.molwt
        formula = m.formula
        return {'mw' : mw, 'formula' : formula}

    def getGlycosilationInfo(self, data):
        """ """
        info = []
        info.append(( 'residue', data['attachedTo']))
        info.append(( 'glycosylating groups', len(data['residues'])))
        return info
        
    def getCofactorInfo(self, data):
        info = []
        smi = data['info']['smi'] 
        resInfo = self.smiToInfo(smi)
        info.append( ('name', data['name']) )
        info.append( ('common name', data['info']['name']) )
        info.append( ('residue number', data['num']) )
        info.append( ('MW', resInfo['mw'] ) )
        info.append( ('brute formula', resInfo['formula'] ) )
        info.append( ('SMI', smi) )
        return info

    def getModifierInfo(self, data):
        """ """
        info = []
        info.append( ('name', data['name']) )
        info.append( ('residue modified', data['attachedTo']))
        info.append( ('modifier size', data['size']))
        info.append( ('MW', data['mw']))
        return info

    def parseCofactors(self, tree):
        """ """
        kst = self.kamaji.structureTypes
        if 'cofactor' in kst:
            cofactorNames = []
            #tree.fragType = None
            tree.fragTypeCount += 1
            cofactorNode = tree.add('Cofactors (%d)')
            cofactorData = kst['cofactor']
            cofactorNode.fragType = 'cofactor'
            cofactorNode.fragTypeCount = 1
            for cId, data in cofactorData.items():
                currNode = cofactorNode.add(data['name'])
                currNode.selection = [ data['num'] ]
                cofactorNames.append( '%s%s' % (data['name'], data['num']) )
                currNode.info = self.getCofactorInfo(data)
                currNode.fragType = 'cofactor'
                currNode.fragTypeCount = 1
            cofactorNode.name = cofactorNode.name % (len(cofactorNode.children))
            cofactorNode.info.append( ('co-factor residues','%s' % ', '.join(cofactorNames) ) )
            tree.info.append( ('co-factors','%d' % len(cofactorNode.children) ) )

    def parseWater(self, tree):
        """ """
        kst = self.kamaji.structureTypes
        if 'water' in kst:
            waterNode = tree.add('Waters (%d)')
            waterData = kst['water']
            tree.fragTypeCount += 1
            waterNode.fragType = 'water'
            waterNode.fragTypeCount = 1
            for cId, data in waterData.items():
                currNode = waterNode.add( '%s-%s' % (data['name'], data['num']))
                currNode.selection = [ data['num'] ]
                currNode.info = {} # water info goes here
                currNode.fragType = 'water'
                currNode.fragTypeCount = 1
            waterNode.name = waterNode.name % (len(waterNode.children))
            tree.info.append( ('waters','%d' % len(waterNode.children) ))

    def generateBiopolymerTree(self, chain, v, _type='protein'):
        """ generate trees for protein, dna, rna, including standard residues 
            and missing/incomplete residues
        """
        info = { 'protein' :  { 'title' : 'Standard residues (%d)', 
                        'tree': self.proteinTrees, 'fragType' : 'protein'},
               'dna' :  { 'title' : 'Standard nucleic bases (%d)', 
                        'tree': self.dnaTrees, 'fragType' : 'dna'},
               'rna' :  { 'title' : 'Standard nucleic bases (%d)', 
                        'tree': self.rnaTrees, 'fragType' : 'rna'}
                }
        title = info[_type]['title']
        tree = info[_type]['tree']
        fragType = info[_type]['fragType']
        root = KTree('Chain %s' % chain, chain)
        stdNode = root.add(title)
        stdNode.selection = self.getResNum(v)
        stdNode.info = self.getBioPolymerInfo(v)
        stdNode.name = (stdNode.name % len(v.keys()))
        stdNode.fragType = root.fragType = fragType
        stdNode.fragTypeCount = root.fragTypeCount = 1
        root.info.append( ('standard residues', '%d' % len(v.keys())))
        kst = self.kamaji.structureTypes
        if '_incomplete' in kst:
            incomplete = []
            for m in kst['_incomplete']:
                if m['num'] in stdNode.selection:
                    if m['chain'] == chain:
                        incomplete.append(m)
            if len(incomplete):
                node = stdNode.add('Missing atoms/problematic')
                node.fragType = 'problematic'
                stdNode.info.append( ('problematic residues', '%d' % len(incomplete)) ) 
                root.info.append( ('problematic residues', '%d' % len(incomplete)) ) 
                node.fragTypeCount += 1
                root.fragTypeCount += 1
                stdNode.fragTypeCount += 1
                root.fragTypeCount += 1
                totMissing = 0
                totMissingNames = []
                for m in incomplete:
                    name = '%s%s (%d missing atoms)' % (m['name'],
                        m['num'], m['missingAtomsCount'])
                    totMissing += m['missingAtomsCount']
                    totMissingNames.append( '%s%s' % (m['name'], m['num']) )
                    curr = node.add(name)
                    curr.info.append( ('residue', '%s%s' % (m['name'], m['num']) ) )
                    curr.info.append( ('missing atoms', '%d' % (m['missingAtomsCount'])))
                    curr.selection = [ m['num'] ]
                    curr.fragType = 'problematic'
                    curr.fragTypeCount = 1
                root.info.append( ('total missing atoms', '%d' % totMissing ))
                stdNode.info.append( ('total missing atoms', '%d' % totMissing ))
                node.info.append( ('total missing atoms', '%d' % totMissing ))
                root.info.append( ('incomplete residues', '%s' % ', '.join(totMissingNames) ))
                node.info.append( ('incomplete residues', '%s' % ', '.join(totMissingNames) ))
        tree.append(root)
        return root

    def getBioPolymerInfo(self, data):
        """ """
        info = []
        info.append( ('residues', '%d' % len(data.keys()) ) )
        return info
        

    def getLigandInfo(self, data):
        """ """
        info = []
        info.append( ('name', data['name']) )
        info.append( ('residue number', data['num']) )
        info.append( ('type', data['type']) )
        info.append( ('size', data['info']['size']) )
        info.append( ('MW', data['info']['mw']) )
        return info

    def getPeptideInfo(self, data):
        """ """
        AaaToA = { 'gly' : 'G', 'ala': 'A', 'val': 'V', 'leu': 'L', 'ile':'I', 'met':'M',
            'phe': 'F', 'trp' : 'W', 'pro': 'P', 'ser':'S', 'thr':'T', 'cys':'C', 'tyr':'Y',
            'asn':'N', 'gln': 'Q', 'asp':'D', 'glu':'E', 'lys':'K', 'arg':'R', 'his':'H'}
        info = []
        name = []
        for idx, res in data['data'].items():
            name.append( AaaToA.setdefault(res['name'].lower(), 'X'))
        info.append( ('name' , ''.join(name)) )
        info.append( ('lenght' , len(data['data'].keys()) ))
        return info

    def getNucleicInfo(self, data):
        """ """
        info = []

        return info

    def generateLigandTree(self, chain, v):
        """ """
        names = {   'short_peptide' : 'peptides',
                    'short_nucleic' : 'nucleic',
                    'sugar'         : 'carbohydrate-like',
                    'lipid'         : 'lipid/fatty acid',
                    'generic'       : 'generic'}

        root = KTree('Chain %s' % chain, chain)
        root.fragType = 'ligand'
        root.fragTypeCount = 0
        ligTypes = []
        for group, data in v.items():
            groupNode = root.add('%s (%%s)' % (names[group]))
            root.fragTypeCount += 1
            if group == 'short_peptide': # 
                ligTypes.append(names[group])
                for lId, lData in data.items():
                    groupNode.fragType = 'short'
                    groupNode.fragTypeCount = 1
                    name = lData['info']['sequence']
                    n = groupNode.add(name)
                    n.selection = [x['num'] for x in lData['data'].values() ]
                    n.info = self.getPeptideInfo(lData)
                    n.fragType = 'short'
                    n.fragTypeCount = 1
            elif group == 'short_nucleic':
                ligTypes.append(names[group])
                for lId, lData in data.items():
                    groupNode.fragType = 'short'
                    groupNode.fragTypeCount = 1
                    name = lData['info']['sequence']
                    n = groupNode.add(name)
                    n.selection = [x['num'] for x in lData['data'].values() ]
                    n.info = self.getNucleicInfo(lData)
                    n.fragType = 'short'
                    n.fragTypeCount = 1
            elif group == 'sugar':
                for lId, lData in data.items():
                    groupNode.fragType = 'sugar'
                    groupNode.fragTypeCount = 1
                    name = '%s-%s' % ( lData['name'], lData['num'])
                    n = groupNode.add(name)
                    n.selection = [ lData['num'] ]
                    n.info = self.getLigandInfo(lData)
                    n.fragType = 'sugar'
                    n.fragTypeCount = 1
            elif group == 'lipid':
                for lId, lData in data.items():
                    groupNode.fragType = 'lipid'
                    groupNode.fragTypeCount = 1
                    name = '%s-%s' % ( lData['name'], lData['num'])
                    n = groupNode.add(name)
                    n.selection = [ lData['num'] ]
                    n.info = self.getLigandInfo(lData)
                    n.fragType = 'lipid'
                    n.fragTypeCount = 1
            elif group == 'generic':
                for lId, lData in data.items():
                    groupNode.fragType = 'ligand' # 'generic'
                    groupNode.fragTypeCount = 1
                    n = groupNode.add( '%s-%s' % (lData['name'], lData['num']) )
                    n.selection = [ lData['num'] ]
                    n.info = self.getLigandInfo(lData)
                    n.fragType = 'ligand'
                    n.fragTypeCount = 1
            groupNode.name =  groupNode.name % ( len(groupNode.getSelection() ))
        self.ligandTrees.append(root)
        return root

    def getResNum(self, data):
        """ """
        res = []
        for rId,rInfo in data.items():
            res.append(rInfo['num'])
        return res
            

    def generateSelection(self, chain, resList):
        """ """
        sel = 'chain %s resnum %s'
        chain = chain.strip()
        if chain == '':
            chain = '_'
        #resNumbers = []
        #for resId, resInfo in data.items():
        #    resNumbers.append(str(resInfo['num']))
        string = sel % (chain, " ".join([str(x) for x in resList]) )
        return string


class KNode:
    def __init__(self, parent, name='node', chainId='_'):
        """ """
        self.selection = []
        self.parent = parent
        self.name = name
        self.chainId = chainId
        self.children = []
        self.info = []
        self.fragType = None
        self.fragTypeCount = 0

    def add(self, name='children'):
        n = KNode(self, name, self.chainId)
        self.children.append(n)
        return n

    def setName(self, name):
        """ """
        self.name = name

    def getName(self, nocount=False):
        """ """
        if nocount:
            return self.name
        else:
            return '%s (%s)' % (self.name, len(self.children)) 

    def getSelection(self):
        """ """
        selection = []
        for c in self.children:
            selection += c.getSelection()
        self.selection += selection
        return self.selection

    def getSelectionString(self):
        """ """
        pattern = 'chain %s resnum %s'
        sel = self.getSelection()
        if self.chainId == '':
            self.chainId = '_'
        return pattern % (self.chainId, " ".join([str(x) for x in sel]))

class KTree(KNode):
    def __init__(self, name ='tree', chainId = '_'):
        """ """
        KNode.__init__(self, None, name, chainId)
        self.chain = chainId
        self.children = []

class KMolTree(KNode):
    def __init__(self, name ='tree', chainId = '_'):
        """ """
        KNode.__init__(self, None, name, chainId)
        self.chain = chainId
        self.children = []

    def getSelectionString(self):
        """ """
        string = ''
        pattern = 'chain %s resnum %s'
        sel = self.getSelection()
        if self.chainId == '':
            self.chainId = '_'
        selection = []
        for children in self.children:
            selection += children.getSelection()
            unique = [ str(x) for x in list(set(selection)) ]
            string += "_+_" + pattern % (children.chain, " ".join(unique) )
        return string 
            
            




if __name__ == '__main__':
    import pybel
    ob = pybel.ob
    import sys

    def printExpand(data, indent=0, exclude=[]):
        """ print nested data
        
         [A]
          |
          --[ modifier]
          |      |
          |    [ 57 ] 
          |      |
          |      |
        """
        if indent > 0:
            #spacer = "" + ("%s|" % (" " *indent))  * (indent+1)
            spacer = "" + ( "%s|" % ("    " * indent) )
            curve = "" + ( "%s'" % ("    " * indent) )
        else:
            spacer = ""
            curve = ""
        print "\n"+spacer,"\n"+spacer,
        #print "\n"+curve,

        #print "SPACER[%s]" % spacer
        if isinstance(data, dict):
            for k in sorted(data.keys()):
                if k in exclude:
                    continue
                v = data[k]
                print "\n%s--[ %s ]" % (curve, k),
                if isinstance(v, dict) or isinstance(v,list):
                    printExpand(v, indent+1, exclude)
                else:
                    print "--> ( %s )" % v,
        elif isinstance(data, list):
            for v in data:
                if v in exclude:
                    continue
                if isinstance(v, dict) or isinstance(v,list):
                    printExpand(v, indent+1, exclude)
                else:
                    #print "(",v,")"
                    print "%s" % v,
                    #print "%s" % (" "*indent), v,


    from MolKit2 import Read
    molfile = Read(sys.argv[1])
    interface = KamajiInterface(molfile)

