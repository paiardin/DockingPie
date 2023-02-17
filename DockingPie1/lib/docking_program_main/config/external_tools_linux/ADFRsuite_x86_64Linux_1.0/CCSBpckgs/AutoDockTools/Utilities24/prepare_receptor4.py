#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/prepare_receptor4.py,v 1.13.6.1 2015/08/26 22:45:31 sanner Exp $
#
import os 

from MolKit import Read
import MolKit.molecule
import MolKit.protein
from AutoDockTools.MoleculePreparation import AD4ReceptorPreparation


if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: prepare_receptor4.py -r filename"
        print
        print "    Description of command..."
        print "         -r   receptor_filename "
        print "        supported file types include pdb,mol2,pdbq,pdbqs,pdbqt, possibly pqr,cif"
        print "    Optional parameters:"
        print "        [-v]  verbose output (default is minimal output)"
        print "        [-o pdbqt_filename]  (default is 'molecule_name.pdbqt')"
        print "        [-A]  type(s) of repairs to make: "
        print "             'bonds_hydrogens': build bonds and add hydrogens "
        print "             'bonds': build a single bond from each atom with no bonds to its closest neighbor" 
        print "             'hydrogens': add hydrogens"
        print "             'checkhydrogens': add hydrogens only if there are none already"
        print "             'None': do not make any repairs "
        print "             (default is 'None')"
        print "        [-C]  preserve all input charges ie do not add new charges "
        print "             (default is addition of gasteiger charges)"
        print "        [-p]  preserve input charges on specific atom types, eg -p Zn -p Fe"
        print "        [-U]  cleanup type:"
        print "             'nphs': merge charges and remove non-polar hydrogens"
        print "             'lps': merge charges and remove lone pairs"
        print "             'waters': remove water residues"
        print "             'nonstdres': remove chains composed entirely of residues of"
        print "                      types other than the standard 20 amino acids"
        print "             'deleteAltB': remove XX@B atoms and rename XX@A atoms->XX"
        print "             (default is 'nphs_lps_waters_nonstdres') "
        print "        [-e]  delete every nonstd residue from any chain"
        print "              'True': any residue whose name is not in this list:"
        print "                      ['CYS','ILE','SER','VAL','GLN','LYS','ASN', "
        print "                      'PRO','THR','PHE','ALA','HIS','GLY','ASP', "
        print "                      'LEU', 'ARG', 'TRP', 'GLU', 'TYR','MET', "
        print "                      'HID', 'HSP', 'HIE', 'HIP', 'CYX', 'CSS']"
        print "              will be deleted from any chain. "
        print "              NB: there are no  nucleic acid residue names at all "
        print "              in the list and no metals. "
        print "             (default is False which means not to do this)"
        print "        [-M]  interactive "
        print "             (default is 'automatic': outputfile is written with no further user input)"
        print "        [-d dictionary_filename] file to contain receptor summary information"
        print "        [-w]   assign each receptor atom a unique name: newname is original name plus its index(1-based)"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'r:vo:A:Cp:U:eM:d:wh')

    except getopt.GetoptError, msg:
        print 'prepare_receptor4.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-s: receptor
    receptor_filename =  None

    # optional parameters
    verbose = None
    #-A: repairs to make: add bonds and/or hydrogens or checkhydrogens
    repairs = ''
    #-C default: add gasteiger charges 
    charges_to_add = 'gasteiger'
    #-p preserve charges on specific atom types
    preserve_charge_types=None
    #-U: cleanup by merging nphs_lps, nphs, lps, waters, nonstdres
    cleanup  = "nphs_lps_waters_nonstdres"
    #-o outputfilename
    outputfilename = None
    #-m mode 
    mode = 'automatic'
    #-e delete every nonstd residue from each chain
    delete_single_nonstd_residues = None
    #-d dictionary
    dictionary = None
    #-w 
    unique_atom_names = False

    #'r:vo:A:Cp:U:eM:d:wh'
    for o, a in opt_list:
        if o in ('-r', '--r'):
            receptor_filename = a
            if verbose: print 'set receptor_filename to ', a
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print 'set outputfilename to ', a
        if o in ('-A', '--A'):
            repairs = a
            if verbose: print 'set repairs to ', a
        if o in ('-C', '--C'):
            charges_to_add = None
            if verbose: print 'do not add charges'
        if o in ('-p', '--p'):
            if not preserve_charge_types:
                preserve_charge_types = a
            else:
                preserve_charge_types = preserve_charge_types + ','+ a
            if verbose: print 'preserve initial charges on ', preserve_charge_types
        if o in ('-U', '--U'):
            cleanup  = a
            if verbose: print 'set cleanup to ', a
        if o in ('-e', '--e'):
            delete_single_nonstd_residues  = True
            if verbose: print 'set delete_single_nonstd_residues to True'
        if o in ('-M', '--M'):
            mode = a
            if verbose: print 'set mode to ', a
        if o in ('-d', '--d'):
            dictionary  = a
            if verbose: print 'set dictionary to ', dictionary
        if o in ('-w', '--w'):
            unique_atom_names = True
            if verbose: print 'set unique_atom_names to ', unique_atom_names
        if o in ('-h', '--'):
            usage()
            sys.exit()


    if not receptor_filename:
        print 'prepare_receptor4: receptor filename must be specified.'
        usage()
        sys.exit()

    #what about nucleic acids???

    mols = Read(receptor_filename)
    if verbose: print 'read ', receptor_filename
    mol = mols[0]
    if unique_atom_names:  # added to simplify setting up covalent dockings 8/2014
        for at in mol.allAtoms:
            if mol.allAtoms.get(at.name) >1:
                at.name = at.name + str(at._uniqIndex +1)
        if verbose:
            print "renamed %d atoms: each newname is the original name of the atom plus its (1-based) uniqIndex" %(len(mol.allAtoms))        
    preserved = {}
    has_autodock_element = False
    if charges_to_add is not None and preserve_charge_types is not None:
        if hasattr(mol, 'allAtoms') and not hasattr(mol.allAtoms[0], 'autodock_element'):
            file_name, file_ext = os.path.splitext(receptor_filename)
            if file_ext == '.pdbqt':
                has_autodock_element = True
        if preserve_charge_types is not None and has_autodock_element==False:
            print 'prepare_receptor4: input format does not have autodock_element SO unable to preserve charges on ' + preserve_charge_types
            print 'exiting...'
            sys.exit()  
        preserved_types = preserve_charge_types.split(',') 
        if verbose: print "preserved_types=", preserved_types
        for t in preserved_types:
            if verbose: print 'preserving charges on type->', t
            if not len(t): continue
            ats = mol.allAtoms.get(lambda x: x.autodock_element==t)
            if verbose: print "preserving charges on ", ats.name
            for a in ats:
                if a.chargeSet is not None:
                    preserved[a] = [a.chargeSet, a.charge]

    if len(mols)>1:
        if verbose: print "more than one molecule in file"
        #use the molecule with the most atoms
        ctr = 1
        for m in mols[1:]:
            ctr += 1
            if len(m.allAtoms)>len(mol.allAtoms):
                mol = m
                if verbose: print "mol set to ", ctr, "th molecule with", len(mol.allAtoms), "atoms"
    alt_loc_ats = mol.allAtoms.get(lambda x: "@" in x.name)
    len_alt_loc_ats = len(alt_loc_ats)
    #import pdb; pdb.set_trace()
    if len_alt_loc_ats:
        alternateLocationNames = list(set([x.split('@')[1] for x in alt_loc_ats.name]))
        alternateLocationNames.sort()
        print "WARNING!", mol.name, "has",len_alt_loc_ats, ' atoms with alternate locations named %s!\n'%alternateLocationNames
        altLocUsed = list(alternateLocationNames)[0]
        print "WARNING! using alternate location %s"%altLocUsed
        atomsToIgnore = [x for x in alt_loc_ats if x.name.split('@')[1]!=altLocUsed]
        for atom in atomsToIgnore:
            parent = atom.parent
            parent.children.remove(atom)
            mol.allAtoms.remove(atom)
            atom.parent = None
            
    mol.buildBondsByDistance()

    if verbose:
        print "setting up RPO with mode=", mode,
        print "and outputfilename= ", outputfilename
        print "charges_to_add=", charges_to_add
        print "delete_single_nonstd_residues=", delete_single_nonstd_residues

    RPO = AD4ReceptorPreparation(mol, mode, repairs, charges_to_add, 
                        cleanup, outputfilename=outputfilename,
                        preserved=preserved, 
                        delete_single_nonstd_residues=delete_single_nonstd_residues,
                        dict=dictionary)    

    if charges_to_add is not None:
        #restore any previous charges
        for atom, chargeList in preserved.items():
            atom._charges[chargeList[0]] = chargeList[1]
            atom.chargeSet = chargeList[0]


# To execute this command type:
# prepare_receptor4.py -r pdb_file -o outputfilename -A checkhydrogens 

