"""Extract a selection of atoms from a PDB file."""

from ..apptools import *

__all__ = ['prody_select']

def prody_select(selstr, *pdbs, **kwargs):
    """Write selected atoms from a PDB file in PDB format.

    :arg selstr: atom selection string, see :ref:`selections`

    :arg pdbs: PDB identifier(s) or filename(s)

    :arg output: output filename, default is :file:`pdb_selected.pdb`

    :arg prefix: prefix for output file, default is PDB filename

    :arg suffix: output filename suffix, default is :file:`_selected`"""

    from os.path import isfile
    from prody import LOGGER, parsePDB, writePDB

    #selstr = kwargs.get('selstr')
    if not pdbs:
        raise ValueError('pdb argument must be provided')

    if ((isfile(selstr) or len(selstr) == 4 and selstr[0].isdigit()) and
        len(pdbs) == 1 and not isfile(pdbs[0])):
        pdbs, selstr = selstr, pdbs[0]
        LOGGER.warn('The order of selstr and pdb arguments have switched '
                    'to support multiple files, old order will be supported '
                    'until v1.4.')
        pdbs = [pdbs]

    prefix = kwargs.get('prefix', None)
    suffix = kwargs.get('suffix', '_selected')
    output = kwargs.get('output', None)

    for pdb in pdbs:
        pdb = parsePDB(pdb)

        pdbselect = pdb.select(selstr)
        if pdbselect is None:
            LOGGER.warn('Selection {0} did not match any atoms.'
                        .format(repr(selstr)))
            return
        LOGGER.info('Selection {0} matched {1} atoms.'
                    .format(repr(selstr), len(pdbselect)))

        outname = output or ((prefix or pdb.getTitle()) + suffix)
        LOGGER.info('Selection is written into: ' +
                    writePDB(outname, pdbselect))



def addCommand(commands):

    subparser = commands.add_parser('select',
    help='select atoms and write a PDB file')

    subparser.add_argument('--quiet', help="suppress info messages to stderr",
        action=Quiet, nargs=0)

    subparser.add_argument('--examples', action=UsageExample, nargs=0,
        help='show usage examples and exit')

    subparser.set_defaults(usage_example=
    """This command selects specified atoms and writes them in a PDB file.

Fetch PDB files 1p38 and 1r39 and write backbone atoms in a file:

  $ prody select backbone 1p38 1r39""",
    test_examples=[0])


    group = subparser.add_argument_group('output options')

    group.add_argument('-o', '--output', dest='output', metavar='STR',
        type=str, help='output PDB filename (default: pdb_selected.pdb)')

    group.add_argument('-p', '--prefix', dest='prefix', metavar='STR',
        type=str, help=('output filename prefix (default: PDB filename)'))

    group.add_argument('-x', '--suffix', dest='suffix', metavar='STR',
        type=str, default='_selected',
        help=('output filename suffix (default: %(default)s)'))

    subparser.add_argument('select', help='atom selection string')
    subparser.add_argument('pdb', nargs='+',
        help='PDB identifier(s) or filename(s)')

    subparser.set_defaults(func=lambda ns: prody_select(ns.select,
                                                         *ns.pdb,
                                                         **ns.__dict__))
    subparser.set_defaults(subparser=subparser)
