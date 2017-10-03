#!/usr/bin/python

import sys,os,platform
import publishfunctions
import notes

curdir = os.path.dirname(os.path.realpath(__file__))

# ------- Configuration parameters -------------

projectname='unlocbox'

if platform.system() == 'Darwin':
    homefolder = '/Users/nati'
else:
    homefolder = '/home/nati'


project=homefolder+'/work/git/toolbox/unlocbox/'

# Configure HTML placement at remote server
# host='nperraud,'+projectname+'@web.sourceforge.net'
# www='/home/project-web/unlocbox/htdocs/'
# outputdirweb= '~/work/git/website/unlocbox/'


# Configure HTML placement at remote server
user = 'nati'

# -------- Configuration of mat2doc ------------
mat2docpath=homefolder+'/work/git/mat2doc'

# Configure matlab paths
unlocxpath = homefolder+'/work/git/toolbox/unlocbox'
ltfatpath = homefolder+'/work/git/toolbox/ltfat'
outputdirweb= '~/work/git/website/unlocbox-html/'

# -------- Configuration of note ------------

notesdir = homefolder+'/work/git/notes/unlocbox-note/'
notehtml = homefolder+'/work/publish/unlocbox/notes/'
noteswww = homefolder+'/work/git/website/unlocbox-html/notes/'
docwww = homefolder+'/work/git/website/unlocbox-html/doc/'

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
outputdirhtml=outputdir+projectname+'-html/'

f=file(project+projectname+'_version')
versionstring=f.read()[:-1]
f.close()

# ------- do not edit below this line ----------
f = open(project + 'mat2doc/startup.m', 'w')
f.write('addpath ' + unlocxpath + '\n')
f.write('init_unlocbox;\n\n')
f.write('addpath ' + ltfatpath + '\n')
f.write('ltfatstart;\n')
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
    # os.system('rm '+outputdirmat+'demos/demo_lrjs.m')
    # os.system('rm '+outputdirmat+'solver/solve_lrjs.m')

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


# if 'sendphp' in todo:

#     s="rsync -e 'ssh -p 10422' -av "+outputdirphp+' '+host+':'+www+'doc/'
#     print s
#     os.system(s)  


# if 'sendfile' in todo:
#     fname=outputdir+projectname+'-'+versionstring
#     # later on to put file on file of sourceforge
#     # s ='rsync -e ssh '+outputdir+fname+'.tar.gz '+'nperraud@frs.sourceforge.net:/home/frs/project/'+projectname+'/'
#     s = 'cp '+fname+'.tar.gz '+outputdir+projectname+'.tar.gz '
#     os.system(s)
#     s = 'cp '+fname+'.zip '+outputdir+projectname+'.zip '
#     os.system(s)
#     s='rsync -av '+outputdir+projectname+'.tar.gz '+host+':'+www
#     os.system(s)
#     s='rsync -av '+outputdir+projectname+'.zip '+host+':'+www
#     os.system(s)




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

if 'copyhtml' in todo:
    os.system('rsync -av '+outputdirhtml+' '+docwww)


#  Make the notes
if 'noteall' in todo:
    todo.append('notesclean')
    todo.append('notesmake')
    todo.append('notestexclean')
    todo.append('notescopy')

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
    #os.system(iconv -f original_charset -t utf-8 originalfile > newfile 
        
    

if 'notescopy' in todo:
    os.system('rsync -av '+notehtml+' '+noteswww)
    # os.system("rsync -e 'ssh -p 10422' -av "+notehtml+' '+host+':'+noteswww);


