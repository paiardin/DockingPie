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

#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner and TSRI 2016
#
#########################################################################
#
#
import argparse, sys
from argparse import ArgumentParser

class ArgParser(ArgumentParser):
    
    def __init__(self, target, check=True, *arg, **kw):
        """build an arguement parser for either runADFR.py or runOneGA.py"""
        assert target in ['ADFR', 'OneGA']
        self.target = target
        self.check = check
        
        # create a parent parser for all arguments that runADFR has to pass to runOneGA
        if target=='ADFR':
            kw['prog']='adfr'
            kw['version']="ADFR v1.1",
            kw['description'] ="""Dock a flexible ligand into a flexible receptor.\nMore details at http://adfr.scripps.edu"""
        else:
            kw['prog']='runOneGA.py'
            kw['version']="runOneGa.py v1.1",
            kw['description']="""Run One GA evolution for docking a flexible ligand into a flexible receptor.\nMore details at http://adfr.scripps.edu"""

        ArgumentParser.__init__(self, *arg, **kw)

        self.add_argument('-l', "--ligand", dest="ligandFile",
                          metavar=('ligand.pdbqt',),
                          required= check==True,
                          help="ligand to be docked (in PDBQT format)",)

        self.add_argument('-t', "--target", dest="receptorMapsFile",
                            metavar=('targetFile.zip',),
                            help="file produced by AGFR and describing the receptor (see http://agfr.scripps.edu)",)

        self.add_argument('-J', "--jobName", dest="jobName", default='NoName',
                             metavar=("jobName"),
                             help="name for this docking experiment")

        self.add_argument('-c', "--maxCores", dest="maxCores", type=int,
                            default=0, metavar=('maxCores',),
                            help="maximum number of cores to be used to run GAs in parallel. A value of 0 will use all available cores. Default 0",)

        if target=='OneGA':
            self.add_argument("-j", "--jobNumber", dest="jobNumber", type=int,
                                metavar=("jobNumber"),help="GA run number",)

        self.add_argument('-o', "--logFile", dest="logfile",
                             metavar=("logFilename",),
                             help="summary docking. A .dlg extension will be added if missing")

        ##
        ## search related parameters
        ##
        gaGroup = self.add_argument_group(
            'SEARCH PARAMETERS:',
            'the following arguments modify the termination criterion of the genetic algorithm')

        if target=='ADFR':
            gaGroup.add_argument('-n', "--nbRuns", dest="nbRuns", type=int,
                                 default=50, metavar=('nbRuns',),
                                 help="number of GA evolutions (default 50)")
        

        gaGroup.add_argument('-p', "--popSize", dest="popSize", type=int,
                             default=150, metavar=("popSize",),
                             help="size of the population evolved by the genetic algorithm. Set to 'auto' to use the builtin heuristic",)

        gaGroup.add_argument('-e', "--maxEvals", dest="maxEvals", type=int,
                             default=2500000, metavar=("maxEvals",),
                             help="maximum number of evaluations of the scoring function per GA (Default 2.5 Millions)")
        

        gaGroup.add_argument('-s', "--noImproveStop", dest="noImproveStop",
                             type=int, default=5, metavar=("noImproveStop",),
                             help="stop GA evaluation after this number of generations with no improvement in best energy in all clusters (Default 5)",)

        gaGroup.add_argument('-g', "--maxGens", dest="maxGens", type=int,
                             default=10000000, metavar=("maxGens",),
                             help="maximum number of generations (default 10 millions)",)

        if target=='OneGA':
            AD4Group = self.add_argument_group(
                'AutoDock4 legacy receptor/maps group:',
                'the following options are less commonly used')
            AD4Group.add_argument('-F', "--mapFilesFolder", dest="mapFilesFolder",
                                  help="fodler containing map files ",
                                  metavar=('mapFilesFolder',))

            AD4Group.add_argument('-M', "--mapFilesRoot", dest="mapFilesRoot",
                                  help="root name for mapfiles ", metavar=('mapFilenameRoot',))

            AD4Group.add_argument('-T', "--transPoints",
                                  metavar=('translationPointsFile.npy',),
                                  dest="tpointsFilename",
                                  help="translational point file")

            AD4Group.add_argument('-R', "--receptor", dest="receptorFile",
                                  metavar=('receptor.pdbqt',),
                                  help="receptor file (in PDBQT format)")

            AD4Group.add_argument('-X', "--flexRes", dest="flexRes",
                                  metavar=('flexResSelectionString',),
                                  help="list of flexbile receptor side chains")
        ##
        ## advanced parameters
        ##
        advancedGroup = self.add_argument_group(
            'ADVANCED OPTIONS:',
            'the following options are less commonly used')

        advancedGroup.add_argument(
            '-C', "--covalentLigand", dest="covalentLigand",
            type=int, nargs=3, metavar=('l1', 'l2', 'l3'),
            help="serial numbers of 3 ligand atoms  used to define covalent attachment of ligand. These are atom numbers as they appear the the PDBQT files of the ligand.pdbqt file")

        advancedGroup.add_argument(
            '-N', "--neighborSearchCutoff", dest="neighborSearchCutoff",
            type=float, metavar=('2.5'),
            help="enabling the neigborhood search with a given RMSD cutoff value")

        advancedGroup.add_argument(
            '-NG', "--neighborSearchGroup", dest="neighborSearchGroup",
            type=str, metavar=('backbone'),
            help="atom group selection used to constrain ADFR neighborhood search")
        
        if target=='OneGA':
            advancedGroup.add_argument(
                '-V', "--covalentReceptor", dest="covalentRec",
                type=int, nargs=3, metavar=('r1','r2','r3'),
                help="serial numbers of 3 receptor atoms used to define covalent attachment of ligand. These are atom numbers as they appear the the PDBQT files of the receptor.pdbqt file")

        advancedGroup.add_argument(
            '-m', "--RMSDMatching", dest="RMSDMatching", type=str,
            metavar=('RMSDatomMatchingMode',),
            choices=["'hungarian'", "'1to1'",], default="'1to1'",
            help="specify the atom matching method used for defining pairs of atoms when RMSD is calculated. Possible values are '1to1' and 'hungarian'. Default: 'hungarian'")

        advancedGroup.add_argument(
            '-r', "--reference", dest="reference",
            metavar=('referenceLigForRMSD.pdbqt',),
            help="filename of reference ligand used for RMSD reporting",)

        advancedGroup.add_argument(
            '-D', "--clusteringRMSDCutoff", dest="clusteringRMSDCutoff",
            metavar=('clusteringRMSDCutoff',), default=2.0, type=float,
            help="distance cut off used for clustering solutions",)
        
        
        if target=='OneGA':
            advancedGroup.add_argument(
                '-S', "--seed", dest="seedValue", type=int, default=-1,
                metavar = ('seedValue',),
                help="seed for random number generator")
        else:
            advancedGroup.add_argument(
                '-S', "--seed", dest="seedValue", type=int, default=-1,
                metavar = ('seedValue',),
                help="seed for random number generator. if the seed is -1 each GA evolution will have a random seed, else the first GA has this seed and subsequent GAs have seeds incremented by 1")
            advancedGroup.add_argument(
                '-T', "--noTargetFileOutput", action="store_true",
                dest="noTargetFileOutput", default=False,
                help="when specified, this flag will prevent the addition of the .trg file used to describe the receptor binding pocket to the resulting .dro file. Useful for virtual screening scenarios to avoid duplication of the trg. file, thus substantially reducing the size of the output for each ligand.")

        advancedGroup.add_argument(
            '-f', "--fullOutput", action="store_true", dest="fullOutput",
            default = False,
            help="when True, the folder with poses and log files from individual GA evolutions is kept, else it will be deleted")

        advancedGroup.add_argument(
            '-d', "--debug", dest="debug", type=int, default=0,
            help="set debugging level (default 0)")

        advancedGroup.add_argument(
            '-y', "--dryRun", dest="dryRun", action="store_true",
            default=False,
            help="print the first runOneGA command line and exit")

        advancedGroup.add_argument(
            '-O', "--overwriteFiles", dest="overwriteFiles",
            action="store_true", default=False,  
            help="overwrite existing output files silently")
        advancedGroup.add_argument('-x', "--fixedRoot", dest="fixedRoot",  action="store_true",
                             default=False, 
                             help="When specified, the atoms in the ligand will not be rotated or translated during the search, and only the ligand's conformation will be optimized. The atoms of the root node of the provided ligand are assumed to be in the correct position and orientation relative to the receptor.")

    def parse_args(self):
        kw = vars(ArgumentParser.parse_args(self))
        if not self.check:
            return kw
        # check mutual exlusion of receptor specification
        legacyTarget = kw.get('mapFilesRoot', None)
        target = kw.get('receptorMapsFile', None)

        # Target is missing
        if target is None and legacyTarget is None:
            print "%s: error: arguments missing target. Use -t/--target or -M/--mapFilesRoot"%self.target
            sys.exit(1)

        # Target is specified twice
        if target is not None and legacyTarget is not None:
            print "%s: error: arguments -t/--target and -M/--mapFilesRoot are mutually exlusive"%self.target
            sys.exit(1)
        return kw

    ## self.add_argument("--FTRecsrc",
    ##                   action="store", # optional because action defaults to "store"
    ##                   dest="FTRecsrc",
    ##                   default = None,
    ##                   help="Python file to execute to build custom receptor FT")

    ## ## self.add_argument("-t", "--torsInfo",
    ## ##                   action="store", # optional because action defaults to "store"
    ## ##                   dest="torsInfo",
    ## ##                   help="torsion Info file",)

    ## ## self.add_argument("--swBiasCoeff1",
    ## ##                   action="store", # optional because action defaults to "store"
    ## ##                   dest="swBiasCoeff1",
    ## ##                   type="float",
    ## ##                   default=0.2,
    ## ##                   help="SW bias coef 1",)
    ## ## self.add_argument("--swBiasCoeff2",
    ## ##                   action="store", # optional because action defaults to "store"
    ## ##                   dest="swBiasCoeff2",
    ## ##                   type="float",
    ## ##                   default=0.4,
    ## ##                   help="SW bias coef 2",)
    ## return self

class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('docking box definition mode') or \
              text.startswith('pocket identification method'):
            return text.splitlines()  
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

class AGFRParser(ArgumentParser):
    
    def __init__(self, *args, **kw):
        """build an argument parser for runAGFR.py"""
        kw['prog']='agfr'
        kw['version']="v1.1",
        kw['description']="""Compute a maps-based representation of a receptor for docking with ADFR"""

        kw['formatter_class'] = SmartFormatter
        ArgumentParser.__init__(self, *args, **kw)
        group = self.add_mutually_exclusive_group(required=True)
        
        group.add_argument('-r', "--receptor", dest="receptorFile",
                          metavar=('receptor.pdbqt',), #required=True,
                          help="receptor for which to compute maps (in PDBQT format)")
        group.add_argument('-F', '--config', dest="configFile",
                           help="name of configuration file containing agfr command line options")

        self.add_argument('-b', '--boxMode', dest='boxMode',
                          metavar = ("boxmode", "modes"), nargs="+",
                          help="""docking box definition mode. Mode can be:\n receptor: smallest box encompassing the entire receptor (default unless -l/--ligand is present).\n ligand  : smallest box encompassing a specified ligand (see -l/--ligand)\n fill    : smallest box encompassing fill points (see fills)\n residues chid resnamesResnums: smallest box encompassing the specified residues\n user cx cy cz sx sy sz: box centered at (cx, cy, cz) of size sx, sy, sz (units Angstroms)\n user centerMode sx sy sz:  centerMode can be: receptor, ligand, fill, residues""")
        
        self.add_argument('-p', '--pocketMode', dest="pocketMode",
                          choices=['best', 'all', 'forTop', 'forEach', 'user'],
                          nargs="+",
                          help='pocket identification method:\n best: pocket with best AutoSite score (default when boxMode is "receptor")\n all : merge all pockets with more than --pocketCutoff points (default when boxMode is not "receptor"); \n forTop :create maps files for each of the N top ranking AutoSite pockets (N is the number given with --pocketCutoff option, default is 10 top ranking pockets);\n forEach: create a maps file for each pocket identified by AutoSite;\n user: the docking pocket is specified by the user as sets of 3D fill points read from a numpy .npy file.')

        self.add_argument('-P', '--padding', dest="padding", type=float,
                          default=4.0,
                          help="amount of padding added to each side of the box after computing the box size based on --boxMode except for boxMode=user. (units Angstroms)")

        self.add_argument(
            '-f','--flexRes', dest='flexres', metavar=('selectionString',),
            help='flexible residue selection string, eg "A:ILE10,VAL32;B:SER48"')

        self.add_argument('-c','--covalentBond', dest='covalentBond',
                          type=int, nargs=2,
                          metavar=('receptorSerial1', 'receptorSerial2'),
                          help='Serial codes (as in pdb file) for atoms in the receptor pdbqt file forming the covalent bond between the receptor and the ligand. All atoms beyond receptorSerial2 will be ignored automatically. This list can be specified explicitly using -x/--covalentResidues')

        self.add_argument('-t','--covalentBondTorsionAtom',
                          dest='covalentBondTorsionAtom', type= int, 
                          metavar=('atomSerial', ),
                          help='Serial codes (as in pdb file) for the receptor atom used to compute the torsion angle for the covalent bond.')

        self.add_argument('-x','--covalentResidues', dest='covalentRes',
                          metavar=('chid1:res1,res2..;chid2:res1,...',),
                          help='Selection string for residues to ignore for map calculations. All atoms beyond the second atom specified by -c/--covalentBond in chid1:res1 and all atoms in all other residues specified here will be ignore for calculating the maps.')

        self.add_argument('-l', "--ligand", dest="ligandFile",
                            metavar=('ligand.pdbqt',),
                            help="(optional) ligand use to define map types to be computed and/or place and size the docking box (in PDBQT format)",)
        
        self.add_argument('-s','--spacing', dest='spacing',
                          default=0.375, type=float,
                          help="spacing of grid points in maps (default 0.375 Angstroms)")

        self.add_argument('-S','--smoothing', dest='smooth',
                          default=0.5, type=float,
                          help="width of well smoothing for AutoGrid4 potentials (default 0.5 Angstroms)")
        
        self.add_argument('-C','--pocketCutoff', dest="pocketCutoff",
                          type=int, default=10,
                          metavar=('cutoff'),
                          help='cutoff for selecting pockets for which maps are computed. If pocketMode is "forEach" the cutoff is the number of points under which a pocket is discarded. If pocketMode is "forTop" the cutoff is the number of top ranking pockets for which maps will be computed')
        
        self.add_argument('-m', '--mapTypes', dest="mapTypes", default='all',
                          choices=['all', 'ligand'],
                          help="AutoDock atom types for which maps are generated, default all")

        self.add_argument('-ls', '--ligandSize', dest="ligandSize", type=int, default=500,
                          help="estimated ligand size for AutoSite1.1 pocket search, default 500")

        ## self.add_argument('-as', '--autoSite', dest="origAutoSite", default=False,
        ##                   action="store_true",
        ##                   help="use the previous AutoSite without cutoff scanning for pocket detection. (faster)")
        self.add_argument('-asv', '--autoSiteV', dest="autoSiteVersion", default=1.0,
                          type=float, 
                          help="specify the AutoSite version (1.0 or 1.1). By default version 1.0 is used - without cutoff scanning for pocket detection(faster)")

        self.add_argument('-ps', '--pepScore', dest="pepScore", default=False,
                          action="store_true",
                          help="use peptide scoring function (E*bg^1.5) to rank the pockets")

        self.add_argument('-o', '--output', dest="outputFile",
                          help="name of the .trg file (zip format) containing the receptor description")

        self.add_argument('-ng', '--noRecGradient', dest='receptorGradient', default=True,
                          action="store_false",
                          help='By default gradients are added to the maps for the regions "inside" the receptor. Use -ng option to compute maps without gradients.',)
        self.add_argument('-V', '--recGradVolCut', dest="recGradVolCut",
                                       type=int, metavar=('cutOffValue'), 
                                       help='Cutoff value. Can be either -1 to keep only the largest cluster of points as the "outside" or a positive integer to keep clusters with more than that many points as the "outside" when computing the gradient for the receptor "inside". Used when -g [--receptorGradient] is specified.' )

        self.add_argument('--waterWeight',  type=float, dest="waterWeight", default=0.6,
                          help='water map parameter: weight, default value 0.6', metavar=('weight'))
        self.add_argument('--waterEntropy',  type=float, dest="waterEntropy", default=-0.2,                  
                          help='water map parameter: entropy, default value -0.2', metavar=('entropy'))

                         

import ConfigParser

def makeConfigFile(options, filename, cmd="agfr"):
    """This function writes a configuration file containing the options that were used with
    agfr or adfr command. Options is the dictionary that was returned by AGFRParser().parser.parse_args() and
    passed to runAFGR.__call__(). """
    config = ConfigParser.ConfigParser()
    section = cmd+'CmdOptions'
    config.add_section(section)
    config.optionxform = str
    parser = AGFRParser()
    parserOptions = {}
    for action in parser._actions:
        parserOptions[action.dest] = action.option_strings
    notIncluded = []
    for option, val in options.items():
        if parserOptions.has_key(option):
            if isinstance(val, (list, tuple)):
                # make a string out of the option values
                valstr = ""
                for vv in val:
                    valstr += str(vv)+" "
                valstr = valstr.strip()
            elif val in (True, False):
                valstr = val
            else: valstr = str(val)
            # find the option string (that starts with "--") from the parser
            if len(parserOptions[option]) > 1:
                optStr = parserOptions[option][1]
            else:
                optStr = parserOptions[option][0]
            if optStr[:2] == "--":
                optStr = optStr[2:]
            else:
                optStr = optStr[1:]
            config.set(section, optStr, valstr)
        else:
            notIncluded.append(option)
    cfgfile = open(filename,'w')
    config.write(cfgfile)
    cfgfile.close()
    if len(notIncluded):
        print "makeConfigFile: Not supported option(s)", notIncluded,  " are not included in the configuration file %s."%filename

def readConfigFile(filename):
    """This function reads an agfr option configuration file. Returns a list of option-value strings"""
    config = ConfigParser.ConfigParser()
    config.optionxform = str
    config.read(filename)
    cmdArgs = []
    parser = AGFRParser()
    parserOptions = {}
    for action in parser._actions:
        optStr = None
        if len(action.option_strings) > 1:
            optStr = action.option_strings[1]
        else:
            optStr = action.option_strings[0]
        if optStr[:2] == "--":
            optStr = optStr[2:]
        if optStr:
            parserOptions[optStr] = action
    for section in config.sections():
        items = config.items(section)
        for opt, val in items:
            if val == "True":
                if parserOptions.has_key(opt):
                    if  parserOptions[opt].default != True:
                        cmdArgs.extend(["--%s"%opt])
            elif val == "False":
                if parserOptions.has_key(opt):
                    if  parserOptions[opt].default != False:
                        cmdArgs.extend(["--%s"%opt])
            else:
                newopt = ["--%s"%opt]
                newopt.extend(val.split(" "))            
                cmdArgs.extend(newopt)
    return cmdArgs
        

if __name__=="__main__":
    import sys, os
    #parser = ArgParser('ADFR')
    #args = parser.parse_args()
    #print args
    parser = AGFRParser()
    args = parser.parse_args()
    print args
