##
## Author: Michel F. Sanner, Anna Omelchenko
## Date: Dec 2014
## CopyRight: Michel F.SANNER and TSRI
##
##
import os, re, sys, string
import shutil

#print "python executable", sys.executable
#print "prefix:", sys.prefix
#print "exec_prefix:", sys.exec_prefix

from os import path
import warnings
import tarfile
import compileall
from shutil import copy

compile = False
import getopt
noLicense = False
optlist, pargs = getopt.getopt(sys.argv[1:], 'cl')
if len(optlist):
    for opt in optlist:
        if opt[0] == "-c":
            compile = True
        elif opt[0] == "-l":
            noLicense = True

#
#print os.environ['PYTHONHOME']
# 
cwd = os.getcwd()
#print "current directory", cwd
ads_root = path.abspath(os.environ['ADS_ROOT'])
ads_archosv = os.environ['ADS_ARCHOSV']
bindir = path.join(ads_root, 'bin')

#print "ads_root", ads_root
#print "ads_archosv", ads_archosv

# 1- Untar and install 
print 'Installing  ADSPackages'
py_version =  (string.split(sys.version))[0][0:3]
#print "python version:", py_version


scriptsInst = path.join(ads_root, 'bin')

pckgsDir = "CCSBpckgs"
adsTars = [pckgsDir, "ThirdPartyPacks"]

for name in adsTars:
    instDir = ads_root
    if path.exists(name+'.tar.gz'):
        print "Installing files from %s " % name+'.tar.gz'
        # uncompress the tarFile
        tf = tarfile.open(name+'.tar.gz', 'r:gz')
        for tfinfo in tf:
            tf.extract(tfinfo, path=instDir)
        tf.close()
if path.exists("LICENSE.txt") and not path.samefile(cwd, instDir):
    copy("LICENSE.txt", instDir)
if path.exists("README") and not path.samefile(cwd, instDir):
    copy("README", instDir)
if path.exists("releaseNotes.txt") and not path.samefile(cwd, instDir):
    copy("releaseNotes.txt", instDir)
if compile: #compile Python source files to byte-code files
    try:
        compileall.compile_dir(path.join(ads_root, "lib"))
    except:
        print "Compillation error"
        
# copy ./Tools/archosv to scriptsInst

copy(path.join("Tools", "archosv"), scriptsInst)

print "Creating scripts"

templatePath = path.join(instDir, 'Tools/scriptTemplate')
f = open(templatePath, 'r')
tplLines = f.readlines()
f.close()

# Get the ADS_ROOT line
l = [x for x in tplLines if x.startswith('ADS_ROOT=')]
if l:
    l = l[0]
    lIndex = tplLines.index(l)
    # set it to the right path
    tplLines[lIndex] = 'ADS_ROOT="%s" \n'%ads_root

# Make pmv2 (for PmvApp)

## pmv2Script = path.join("$ADS_ROOT", pckgsDir,'PmvApp', 'GUI', 'Qt', 'bin', 'runPmv.py' )
## pmv2lines = """if test $# -gt 0
## then
## 	exec $python $pyflags %s $@
## else
## 	exec $python $pyflags %s
## fi
## """%(pmv2Script, pmv2Script)
## pmv2shPath = path.join(bindir, 'pmv2')
## f = open(pmv2shPath, 'w')
## f.writelines(tplLines)
## f.write(pmv2lines)
## f.close()
## os.chmod(pmv2shPath, 509)

# make adfr script 
adfrScript = path.join("$ADS_ROOT", pckgsDir, 'ADFR',  'bin', 'runADFR.py')
adfrlines = """if test $# -gt 0
then
	exec $python $pyflags %s $@
else
	exec $python $pyflags %s
fi
"""%(adfrScript, adfrScript)
adfrshPath = path.join(bindir, 'adfr')
f = open(adfrshPath, 'w')
f.writelines(tplLines)
f.write(adfrlines)
f.close()
os.chmod(adfrshPath, 509)

# make adfrgui script 
#adfrguiScript = path.join("$ADS_ROOT", pckgsDir, 'ADFR',  'bin', 'ADFRgui.py')
#adfrguilines = """if test $# -gt 0
#then
#	exec $python $pyflags %s $@
#else
#	exec $python $pyflags %s
#fi
#"""%(adfrguiScript, adfrguiScript)
#adfrguishPath = path.join(bindir, 'adfrgui')
#f = open(adfrguishPath, 'w')
#f.writelines(tplLines)
#f.write(adfrguilines)
#f.close()
#os.chmod(adfrguishPath, 509) 

# make agfr script 
agfrScript = path.join("$ADS_ROOT", pckgsDir, 'ADFR',  'bin', 'runAGFR.py')
agfrlines = """if test $# -gt 0
then
	exec $python $pyflags %s $@
else
	exec $python $pyflags %s
fi
"""%(agfrScript, agfrScript)
agfrshPath = path.join(bindir, 'agfr')
f = open(agfrshPath, 'w')
f.writelines(tplLines)
f.write(agfrlines)
f.close()
os.chmod(agfrshPath, 509)

# make agfrgui script 
agfrguiScript = path.join("$ADS_ROOT", pckgsDir, 'ADFR',  'bin', 'AGFRgui.py')
agfrguilines = """if test $# -gt 0
then
	exec $python $pyflags %s $@
else
	exec $python $pyflags %s
fi
"""%(agfrguiScript, agfrguiScript)
agfrguishPath = path.join(bindir, 'agfrgui')
f = open(agfrguishPath, 'w')
f.writelines(tplLines)
f.write(agfrguilines)
f.close()
os.chmod(agfrguishPath, 509)

# make about script 
aboutScript = path.join("$ADS_ROOT", pckgsDir, 'ADFR',  'bin', 'about.py')
aboutlines = """if test $# -gt 0
then
	exec $python $pyflags %s $@
else
	exec $python $pyflags %s
fi
"""%(aboutScript, aboutScript)
aboutshPath = path.join(bindir, 'about')
f = open(aboutshPath, 'w')
f.writelines(tplLines)
f.write(aboutlines)
f.close()
os.chmod(aboutshPath, 509)

# make autosite script 
autositeScript = path.join("$ADS_ROOT", pckgsDir, 'AutoSite',  'bin', 'AS.py')
autositelines = """if test $# -gt 0
then
	exec $python $pyflags %s $@
else
	exec $python $pyflags %s
fi
"""%(autositeScript, autositeScript)
autositeshPath = path.join(bindir, 'autosite')
f = open(autositeshPath, 'w')
f.writelines(tplLines)
f.write(autositelines)
f.close()
os.chmod(autositeshPath, 509)



# make adcp script (not for darwin yet)
adcpScript = path.join("$ADS_ROOT", pckgsDir, 'ADCP',  'runADCP.py')
if path.exists(path.join(ads_root, pckgsDir, 'ADCP',  'runADCP.py')):
    adcplines = """if test $# -gt 0
then
    exec $python $pyflags %s $@
else
    exec $python $pyflags %s
fi
    """%(adcpScript, adcpScript)
    adcpshPath = path.join(bindir, 'adcp')
    f = open(adcpshPath, 'w')
    f.writelines(tplLines)
    f.write(adcplines)
    f.close()
    os.chmod(adcpshPath, 509)

#prepare_receptor , prepare_ligand

if path.exists(path.join(ads_root, pckgsDir, "AutoDockTools", "Utilities24")):
    prepareRecScript = path.join("$ADS_ROOT", pckgsDir, "AutoDockTools", "Utilities24", "prepare_receptor4.py")
    prepareRecLines =  """if test $# -gt 0
then
	exec $python $pyflags %s $@
else
	exec $python $pyflags %s
fi
"""%(prepareRecScript, prepareRecScript)
    prepareRecPath = path.join(bindir, 'prepare_receptor')
    f = open(prepareRecPath, 'w')
    f.writelines(tplLines)
    f.write(prepareRecLines)
    f.close()
    os.chmod(prepareRecPath, 509)

    prepareLigScript = path.join("$ADS_ROOT", pckgsDir, "AutoDockTools", "Utilities24", "prepare_ligand4.py")
    prepareLigLines =  """if test $# -gt 0
then
	exec $python $pyflags %s $@
else
	exec $python $pyflags %s
fi
"""%(prepareLigScript, prepareLigScript)
    prepareLigPath = path.join(bindir, 'prepare_ligand')
    f = open(prepareLigPath, 'w')
    f.writelines(tplLines)
    f.write(prepareLigLines)
    f.close()
    os.chmod(prepareLigPath, 509)


# Make python executable
if sys.platform == 'darwin':
    #comment open -a X11
    #l = filter(lambda x: x.find('This assumes X11 is installed') != -1, tplLines)
    l = [x for x in tplLines if x.find('This assumes X11 is installed') != -1]
    if l:
        l = l[0]
        lIndex = tplLines.index(l)
        for i in range(1,10):
    	    tplLines[lIndex+i] = "#"+tplLines[lIndex+i]
    	    
pythonlines = """if test $# -gt 0
then
	exec $python $pyflags "$@"
else
	exec $python $pyflags 
fi
"""
pythonshPath = path.join(bindir, 'pythonsh')
f = open(pythonshPath, 'w')
f.writelines(tplLines)
f.write(pythonlines)
f.close()
os.chmod(pythonshPath, 509)

#create adsenv.sh and adsenv.csh files
adsenvshPath = path.join(bindir, 'adsenv.sh')
f = open(adsenvshPath , 'w')
f.writelines(tplLines)
f.close()
os.chmod(adsenvshPath, 509)


f = open(path.join(cwd, 'Tools/adsenv.csh'), 'r')
adsenvLines = f.readlines()
f.close()
# Get the ADS_ROOT line
#l = filter(lambda x: x.startswith('setenv ADS_ROOT'), adsenvLines)
l = [x for x in adsenvLines if x.startswith('setenv ADS_ROOT')]
if l:
    l = l[0]
    lIndex = adsenvLines.index(l)
    # set it to the right path
    adsenvLines[lIndex] = 'setenv ADS_ROOT %s\n'%ads_root
adsenvcshPath = path.join(bindir, 'adsenv.csh')
f = open(adsenvcshPath , 'w')
f.writelines(adsenvLines)
f.close()
os.chmod(adsenvcshPath, 509)

# create scripts that set ads environmental variables to run OpenBabel executables (that are located in ADS_ROOT/bin/obabelbin directory
obexecs= ["babel", "obchiral", "obenergy",  "obgen", "obminimize",  "obprop", "obrotamer",  "obspectrophore", "obabel", "obconformer",  "obfit", "obgrep",  "obprobe", "obrms", "obrotate",   "roundtrip"]
templatePath = path.join(cwd, 'Tools/obscriptTemplate')
f = open(templatePath, 'r')
tplLines = f.readlines()
f.close()

# Get the ADS_ROOT line
l = [x for x in tplLines if x.startswith('ADS_ROOT=')]
if l:
    l = l[0]
    lIndex = tplLines.index(l)
    # set it to the right path
    tplLines[lIndex] = 'ADS_ROOT="%s" \n'%ads_root

for obfile in obexecs:
    oblines = """obexec="$ADS_ROOT/bin/obabelbin/%s"\nexec $obexec  $@""" %(obfile,)
    obPath = path.join(bindir, obfile)
    f = open(obPath, 'w')
    f.writelines(tplLines)
    f.write(oblines)
    f.close()
    os.chmod(obPath, 509)


#create sitecustomize.py
f1 = open(os.path.join(ads_root, pckgsDir, "Support", "sitecustomize.py"))
txt = f1.readlines()
f1.close()
f2 = open(os.path.join(ads_root, "lib", "python%s"%py_version, "sitecustomize.py"), "w")
f2.write("adsroot = '%s'\n" % ads_root)
if os.environ.has_key("MGL64"):
    f2.write("import os\n")
    f2.write("os.environ['MGL64']='1'\n")
f2.writelines(txt)
f2.close()

# check if initPython is sourced in your shell ressource file
#shell = sys.argv[1]
print "current directory:", os.getcwd()
#alias_csh = """alias pmv2 %s/bin/pmv2
alias_csh = """alias adfr %s/bin/adfr
alias agfr %s/bin/agfr
alias agfrgui %s/bin/agfrgui
alias autosite %s/bin/autosite
alias about %s/bin/about
alias prepare_receptor %s/bin/prepare_receptor
alias prepare_ligand %s/bin/prepare_ligand
alias pythonsh %s/bin/pythonsh\n""" % (ads_root, ads_root, ads_root, ads_root, ads_root, ads_root, ads_root, ads_root)
#(ads_root, ads_root, ads_root, ads_root, ads_root, ads_root, ads_root, ads_root, ads_root)

#alias_sh="""alias pmv2='%s/bin/pmv2'
alias_sh="""alias agfr='%s/bin/agfr'
alias adfr='%s/bin/adfr'
alias agfrgui='%s/bin/agfrgui'
alias autosite='%s/bin/autosite'
alias about='%s/bin/about'
alias prepare_receptor='%s/bin/prepare_receptor'
alias prepare_ligand='%s/bin/prepare_ligand'
alias pythonsh='%s/bin/pythonsh'\n""" %(ads_root, ads_root, ads_root, ads_root,  ads_root, ads_root, ads_root, ads_root)

#(ads_root, ads_root, ads_root, ads_root, ads_root, ads_root, ads_root, ads_root, ads_root)


f = open("initADFRsuite.csh", "w")
f.write(alias_csh)
f.close()

f = open("initADFRsuite.sh", "w")
f.write(alias_sh)
f.close()

#license part:

if not noLicense:
    text = """\nThe molecular surface calculation software (MSMS) is freely available for academic research.\nFor obtainig commercial license usage contact Dr. Sanner at sanner@scripps.edu.\n"""
    #ans = raw_input(text)
    #if ans=='':
    ans = 'Y'
    while ans[0] not in ['y', 'Y', 'n', 'N']:
        ans = raw_input("Please enter Y or N: ")
    if ans[0] in ['y', 'Y']:
        #academic installation -->> rename mslibACA mslib
        mslib = path.join(ads_root, pckgsDir, "mslibACA")
    else: # commercial installation: rename mslibCOM mslib
        mslib = path.join(ads_root, pckgsDir, "mslibCOM")
    if path.exists(mslib):
        #os.rename(mslib, path.join(ads_root, pckgsDir, "mslib") )
        shutil.move(mslib, path.join(ads_root, pckgsDir, "mslib") )

print """\n ADFRsuite installation complete.
To run agfr, agfrgui, adfr, autosite, about, pythonsh scripts located at:
%s/bin
add  %s/bin to the path environment variable in .cshrc or .bashrc:
.cshrc:
set path = (%s/bin $path)

.bashrc:
export PATH=%s/bin:$PATH

"""%(ads_root, ads_root, ads_root, ads_root)

