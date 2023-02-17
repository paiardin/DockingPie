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

#
# $Header: /mnt/raid/services/cvs/MolKit2/mglmmtf.py,v 1.1.4.1 2017/07/26 22:01:12 annao Exp $ 
#
# $Id: mglmmtf.py,v 1.1.4.1 2017/07/26 22:01:12 annao Exp $
#
from prody.atomic.atomgroup import AtomGroup
from prody.atomic import ATOMIC_FIELDS
import numpy as np

mapSS = {        # MMTF  ->    Prody
    -1: 'C',  # undefined -> coil
    0: 'I',  # pi-helix -> 5-turn helix (pi helix). Min length 5 residues.
    1: 'S',  # bend -> bend
    2: 'H',   # alpha helix -> 4-turn helix (alpha helix) Min length 4 residues.
    3: 'E',  # extended -> extended strand in parallel and/or anti-parallel
    4: 'G',  # 3-10 helix -> 3-turn helix (3-10 helix)
    5: 'B',  # bridge -> residue in isolated beta-bridge (single pair beta-sheet
    #           hydrogen bond formation)
    6: 'T',  # turn -> hydrogen bonded turn (3, 4 or 5 turn)
    7: 'C',  # coil -> non of the above
    }

def MMTFtoPrody(decoder, name='NoName'):
    """Gets atomic info from mmtf decoder, creates prody atom group. """
    atomgroup = AtomGroup(name)
    asize = decoder.num_atoms

    # create numpy arrays for
    coordinates = np.zeros((asize, 3), dtype=float)
    coordinates[:,0] = decoder.x_coord_list
    coordinates[:,1] = decoder.y_coord_list 
    coordinates[:,2] = decoder.z_coord_list
    bfactors = np.zeros(asize, dtype=ATOMIC_FIELDS['beta'].dtype)
    bfactors[:] = decoder.b_factor_list
    occupancies = np.zeros(asize, dtype=ATOMIC_FIELDS['occupancy'].dtype)
    occupancies[:] = decoder.occupancy_list
    altlocs = np.zeros(asize, dtype=ATOMIC_FIELDS['altloc'].dtype)
    altlocs[:] = decoder.alt_loc_list
    serials = np.zeros(asize, dtype=ATOMIC_FIELDS['serial'].dtype)
    serials[:] = decoder.atom_id_list

    atomnames = np.zeros(asize, dtype=ATOMIC_FIELDS['name'].dtype)
    resnames = np.zeros(asize, dtype=ATOMIC_FIELDS['resname'].dtype)
    resnums = np.zeros(asize, dtype=ATOMIC_FIELDS['resnum'].dtype)
    chainids = np.zeros(asize, dtype=ATOMIC_FIELDS['chain'].dtype)

    atomtypes = np.zeros(asize, dtype=ATOMIC_FIELDS['type'].dtype) # ?
    hetero = np.zeros(asize, dtype=bool)
    segnames = np.zeros(asize, dtype=ATOMIC_FIELDS['segment'].dtype)
    elements = np.zeros(asize, dtype=ATOMIC_FIELDS['element'].dtype)
    charges = np.zeros(asize, dtype=ATOMIC_FIELDS['charge'].dtype)
    radii = np.ones(asize, dtype=ATOMIC_FIELDS['radius'].dtype)

    # traverse models
    groupList = decoder.group_list
    for modelIndex, modelChainCount in enumerate(decoder.chains_per_model):
        # traverse chains
        atomIndex = 0
        bonds = np.array(decoder.bond_atom_list, 'i').reshape(-1,2).tolist()
        bo = decoder.bond_order_list.tolist()
        groupOff = 0
        for chainIndex in xrange(modelChainCount):#decoder.chains_per_model[modelIndex]):
            firstAtomsInChain = atomIndex
            nbGroupsInChain = decoder.groups_per_chain[ chainIndex ]
            for groupIndex in xrange(nbGroupsInChain):
                group = groupList[ decoder.group_type_list[groupOff+groupIndex ] ]
                nbAtoms = len(group['atomNameList'])
                atomnames[atomIndex:atomIndex+nbAtoms] = group['atomNameList']
                elements[atomIndex:atomIndex+nbAtoms] = group['elementList']
                resnums[atomIndex:atomIndex+nbAtoms] = decoder.group_id_list[groupOff+groupIndex]
                resnames[atomIndex:atomIndex+nbAtoms] = group['groupName']
                #import pdb; pdb.set_trace()
                bonds.extend( (np.array(group['bondAtomList'], 'i')+atomIndex).reshape(-1,2).tolist() )
                bo.extend(group['bondOrderList'])
                atomIndex += nbAtoms
                #print modelIndex, decoder.chain_id_list[chainIndex], group['groupName'], decoder.group_id_list[groupOff+groupIndex], nbAtoms, atomIndex
            groupOff += nbGroupsInChain
            chainids[firstAtomsInChain:atomIndex] = decoder.chain_id_list[chainIndex]

    atomgroup._setCoords(coordinates)
    atomgroup.setBetas(bfactors)
    atomgroup.setOccupancies(occupancies)
    atomgroup.setAltlocs(altlocs)
    atomgroup.setSerials(serials)
    atomgroup.setNames(atomnames)
    atomgroup.setTypes(atomtypes)
    atomgroup.setData('atomType', atomtypes)
    atomgroup.setResnames(resnames)
    atomgroup.setResnums(resnums)
    atomgroup.setChids(chainids)
    atomgroup.setFlags('hetatm', hetero)
    #atomgroup.setFlags('pdbter', termini)
    #atomgroup.setIcodes(np.char.strip(icodes))
    atomgroup.setSegnames(np.char.strip(segnames))
    atomgroup.setElements(np.char.strip(elements))
    atomgroup.setCharges(charges)
    atomgroup.setRadii(radii)
    atomgroup.setBonds(bonds, bo)

    if hasattr(decoder, 'sec_struct_list'):
        sst = []
        ## mmtf types are 
        ## 0	pi helix
        ## 1	bend
        ## 2	alpha helix
        ## 3	extended
        ## 4	3-10 helix
        ## 5	bridge
        ## 6	turn
        ## 7	coil
        ## -1	undefine
        ##           0    1    2    3    4    5    6    7    8
        atomgroup.setSecstrs(['C']*len(atomgroup))
        ct = 0
        for res in atomgroup.getHierView().iterResidues():
            res.setSecstrs(mapSS[decoder.sec_struct_list[ct]])
            ct += 1
    #else:
    #    ## make everything a coil
    #    atomgroup.setSecstrs(['C']*len(atomgroup))
    
    return atomgroup
