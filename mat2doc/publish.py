#!/usr/bin/python

import sys,os,platform
import publishfunctions
import notes

curdir = os.path.dirname(os.path.realpath(__file__))

# ------- Configuration parameters -------------

projectname='unlocbox'

if platform.system()=='Darwin':
    homefolder='/Users/user'
else:
    homefolder='/home/user'

project=homefolder+'/'+projectname+'/'

# Configure HTML placement at remote server
host='nati@lts2research.epfl.ch'
www='/home/nati/'+projectname+'/'
outputdirweb= '~/work/git/website/unlocbox/'


# -------- Configuration of mat2doc ------------
mat2docpath='~/mat2doc'

try :
    from publish_local import *
except:
    pass



# -------- Configuration of note ------------

notesdir = homefolder+'/work/svn/unlocbox-note/'
notehtml = homefolder+'/work/publish/unlocbox/notes/'
noteswww = www+'notes/'


# -------- Configuration of doc tex ------------
fulldocnote='003'
tutorialnote1='008'
tutorialname1='demos/demo_unlocbox*'

# -------- Automatique configuration ------------
import conf
outputdir=conf.outputdir
outputdirphp=outputdir+projectname+'-php/'
outputdirmat=outputdir+projectname+'-mat/'
outputdirrelease=outputdir+projectname+'-release/'
outputdirtex=outputdir+projectname+'-tex/'

f=file(project+projectname+'_version')
versionstring=f.read()[:-1]
f.close()

# ------- do not edit below this line ----------
f = open(project + 'mat2doc/startup.m', 'w')
f.write('addpath ' + unlocxpath + '\n')
f.write('init_unlocbox;\n\n')
f.write('addpath ' + ltfatpath + '\n')
f.write('ltfatstart;\n');
f.close()


todo=sys.argv


#if 'package' in todo:
#    todo.append('release')

#  Optional parameters

if 'fast' in todo:
    plot='--no-execplot '    
else:
    plot='--execplot '


if 'rebuild' in todo:
    build='--rebuild '
elif 'cached' in todo:
    build='--cached '
else:
    build='--auto '


if 'phpall' in todo:
    todo.append('mat')
    todo.append('php')
    todo.append('package')
    todo.append('sendphp')
    todo.append('tex')


#  Publish
for mode in  ['mat', 'php', 'html', 'tex']:
    if mode in todo:
        s = '%s %s/mat2doc.py %s%s %s %s' % ('PYTHONPATH="%s:$PYTHONPATH"' % (curdir,), mat2docpath, plot if mode != 'mat' else '', build if mode != 'mat' else '', project, mode,)
        os.system(s)


#  Packaging
if 'package' in todo:
    
    # Remove unwanted files
    os.system('rm -rf '+outputdirmat+'mat2doc_old')
    os.system('rm -rf '+outputdirmat+'test_bench')
    os.system('rm '+outputdirmat+'demos/demo_lrjs.m')
    os.system('rm '+outputdirmat+'solver/solve_lrjs.m')

    #fname=outputdir+'/'+projectname+'-'+versionstring
    ## Create the Unix src package
    #os.system('tar zcvf '+fname+'.tgz '+projectname+'-mat/')

    ## Create the Windows src package
    #os.system('rm '+fname+'.zip')
    #publishfunctions.unix2dos(outputdir+projectname+'-mat')
    #os.system('zip -r '+fname+'.zip '+projectname+'-mat/')    
    
    fname=outputdir+projectname+'-'+versionstring
    # Create the Unix src package
    os.system('rm '+fname+'.tar.gz')
    os.system('cd '+outputdir+' && cp -r ' +projectname+'-mat '+projectname+' && tar zcvf '+fname+'.tar.gz '+projectname+'/')
   
    # Create the Windows src package
    os.system('rm '+fname+'.zip')
    publishfunctions.unix2dos(outputdir + projectname +'/')
    os.system('cd '+outputdir+' && zip -r '+fname+'.zip '+projectname+'/')    

    os.system('rm -r '+outputdir+projectname)




#  Send to the server


if 'sendphp' in todo:
    s='rsync -av '+outputdirphp+' '+host+':'+www+'doc/'
    os.system(s)  


#if 'sendweb' in todo:
#	s='rsync -av --exclude=".svn" '+outputdirweb+' '+host+':'+www
#	os.system(s)

if 'sendfile' in todo:
    fname=outputdir+projectname+'-'+versionstring
    # later on to put file on file of sourceforge
    # s ='rsync -e ssh '+outputdir+fname+'.tar.gz '+'nperraud@frs.sourceforge.net:/home/frs/project/'+projectname+'/'
    s = 'cp '+fname+'.tar.gz '+outputdir+projectname+'.tar.gz '
    os.system(s)
    s = 'cp '+fname+'.zip '+outputdir+projectname+'.zip '
    os.system(s)
    s='rsync -av '+outputdir+projectname+'.tar.gz '+host+':'+www
    os.system(s)
    s='rsync -av '+outputdir+projectname+'.zip '+host+':'+www
    os.system(s)




if 'copytex' in todo:
    s='cp -R '+outputdirtex+' '+notesdir+fulldocnote
    os.system(s)
    s='cp '+outputdirtex+tutorialname1+' '+notesdir+tutorialnote1
    os.system(s)

    # Suppress the first line and replace description by introduction
    f = open(notesdir+tutorialnote1+"/demo_unlocbox.tex","r")
    lines = f.readlines()
    f.close()
    f = open(notesdir+tutorialnote1+"/demo_unlocbox.tex","w")
    linenumber = 0
    for line in lines:
        linenumber = linenumber+1
        if  linenumber>1:
            line = line.replace("\\begin{figure}","\\begin{figure}[ht!]")
            line = line.replace("\\subsubsection","\\section")
            line = line.replace("\\section{Description","\\section{Introduction")

            f.write(line)
            #print line
    f.close()




#  Make the notes
if 'noteall' in todo:
    todo.append('notesclean')
    todo.append('notesmake')
    todo.append('notestexclean')
    todo.append('noteshtml')
    todo.append('notesend')

if 'notesclean' in todo:

    notes2 = notes.getnotenumbers(notesdir)

    notes2 = filter(lambda x: (os.path.exists(notesdir+x+'/Makefile')), notes2)

    for notenumber in notes2:
        print 'Cleaning Unlocbox-note '+notenumber
        os.system('cd '+notesdir+notenumber+'; make clean')

if 'notestexclean' in todo:

    notes2 = notes.getnotenumbers(notesdir)

    notes2 = filter(lambda x: (os.path.exists(notesdir+x+'/Makefile')), notes2)

    for notenumber in notes2:
        print 'Cleaning Unlocbox-note '+notenumber
        os.system('cd '+notesdir+notenumber+'; make texclean')

if 'notesmake' in todo:
    notes2=notes.getnotenumbers(notesdir)

    notes2 = filter(lambda x: (os.path.exists(notesdir+x+'/Makefile')), notes2)

    for notenumber in notes2:
        print 'Trying to make Unlocbox-note '+notenumber
        os.system('cd '+notesdir+notenumber+'; make')


if 'noteshtml' in todo:

    notes.printnoteshtml('unlocbox-note-',notesdir,notehtml)
        
    

if 'notesend' in todo:
    os.system('rsync -av '+notehtml+' '+host+':'+noteswww);






