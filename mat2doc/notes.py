#!/bin/env python

import sys,os,shutil,time
#import codecs
cwd=os.getcwd()+'/'
sys.path.append(cwd)





class NotesError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


def getnotenumbers(notesdir):
    allnotes=os.listdir(notesdir)
    allnotes = filter(lambda x: (os.path.isdir(notesdir+x) and (x[0].isdigit())),allnotes)
    
    return allnotes


def parseconfigfiles(noteprefix,notesdir,authdict):

    notes=getnotenumbers(notesdir)

    # Like notes, but cleaned of notes without a config file
    notes_found=[]

    allnotesdict={}

    for note in notes:

        notedict={}

        conffilename=notesdir+note+'/config'
        if not os.path.exists(conffilename):
            print 'Skipping note',note,'config file missing'
            continue

        notes_found.append(note)

        f=file(conffilename,'r')
        buf=f.readlines()
        f.close()
    
        # Initialize notedict with the information from the config file
        for line in buf:
            s=line.split(' ')
            # Join add a newline at the end, remove this when joining
            if s[0]=='link':
                if not 'link' in notedict:
                    notedict['link']={}
                notedict['link'][s[1]]=' '.join(s[2:])[:-1]
            else:
                notedict[s[0]]=' '.join(s[1:])[:-1]
            
        # If year is 'current', update the year with the year of the .pdf file
        if notedict['year']=='current':
            modtime=os.path.getmtime(notesdir+note+'/'+noteprefix+note+'.pdf')
            year=time.localtime(modtime)[0]
            notedict['year']=str(year)

        # Expand the authors
        authlist=[]
        for author in notedict['author'].split():
            if not authdict.has_key(author):
                raise NotesError('Note %s: Author %s not found in dict.' % (note, author))
            authlist.append(authdict[author])

        notedict['author']=authlist

        # Look for bibtex entry
        notedict['bibentry']=os.path.exists(notesdir+note+'/bibentry')

        # Look for poster
        notedict['poster']=os.path.exists(notesdir+note+'/poster.pdf')

        # Look for slides
        notedict['slides']=os.path.exists(notesdir+note+'/slides.pdf')

        # Look for rr
        notedict['rr']=os.path.exists(notesdir+note+'/rr.zip')

        # Generate list of obsoleted notes
        if 'obsoletes' in notedict:
            notedict['obsoletes']=notedict['obsoletes'].split(' ')
        else:
            notedict['obsoletes']=[]
        
        # Create empty list
        notedict['obsoletedby']=[]
        allnotesdict[note]=notedict

    # Reverse-index the obsoletes field
    for note in notes_found:
        for obs in allnotesdict[note]['obsoletes']:
            allnotesdict[obs]['obsoletedby'].append(note)

    return allnotesdict

# Return the number of the current notes, given a dictionary of all
# notes, as returned by parseconfigfiles
def getcurrentnotes(allnotesdict):
    allnotes=allnotesdict.keys()
    currentnotes=set(allnotes)

    for note in allnotes:
        currentnotes=currentnotes.difference(set(allnotesdict[note]['obsoletes']))

    return list(currentnotes)
    

def parseauthors(authorfile):

    f=file(authorfile)
    buf = f.readlines()
    f.close

    authdict={}
    for line in buf: 
        fields = line.split(',') 
        key  = fields[0].strip()
        
        name = fields[1].strip()
        if len(fields)>2:
            email = fields[2].strip()
        else:
            email=''
        if len(fields)>3:
            homepage = fields[3].strip()
        else:
            homepage=''

        authdict[key]={}
        authdict[key]['name']=name
        authdict[key]['email']=email
        authdict[key]['homepage']=homepage
    
    return authdict

def createindexpage(noteprefix,notesdir,allnotes,keys,filename):

    obuf=[]

    # Open table, and append header and first horizontal line
    obuf.append('<table cellpadding=5><tr valign=top><td>')   
    obuf.append('<tr valign="top" align="left"><th><a href="index.php">No.</a></th><th><a href="index_author.php">Name</a></th><th><a href="index_year.php">Year</a></th><th><a href="index_type.php">Type</a></th></tr>')
    obuf.append('<tr><td colspan=4><hr /></td></tr>')

    for note in keys:
        notedict = allnotes[note]

        obuf.append('<tr valign="top"><td>')
        obuf.append(note+'</td><td>')
        obuf.append('<notes_title><a href="'+noteprefix+note+'.pdf">'+notedict['title']+'</a></notes_title><br>')
        
        # Print the author line
        first_author=1;
        s='<notes_author>'
        for author in notedict['author']:
            if first_author:
                first_author=0
            else:
                    #s=s+' and '
                s=s+', '
            if len(author['homepage'])==0:
                s=s+author['name']
            else:
                s=s+'<a href="'+author['homepage']+'">'+author['name']+'</a>'
        obuf.append(s+'</notes_author><br>')

        #obuf.append('<a href="'+noteprefix+note+'.pdf">ltfatnote'+note+'.pdf</a>')
        #obuf.append('<a href="'+noteprefix+note+'.pdf">Download</a>')


        if notedict['poster']:
            obuf.append('<notes_third>Download <a href="'+noteprefix+note+'_poster.pdf">poster</a></notes_third><br>')

        if notedict['slides']:
            obuf.append('<notes_third>Download <a href="'+noteprefix+note+'_slides.pdf">slides</a></notes_third><br>')

        if notedict['bibentry']:
            obuf.append('<notes_third><a href="'+noteprefix+note+'.bib">Cite this paper</a></notes_third><br>')

        if notedict['rr']:
            obuf.append('<notes_third><a href="'+noteprefix+note+'_rr.zip">Reproducible research</a></notes_third><br>')

        if 'link' in notedict:
            for display,link in notedict['link'].iteritems():
                obuf.append('<notes_third><a href="'+link+'">'+display+'</a></notes_third><br>')

        if len(notedict['obsoletedby'])>0:
            obuf.append('<notes_third>This note has been made obsolete by ')
            for obs in notedict['obsoletedby']:
                obuf.append(obs+' ')
            obuf.append('</notes_third>')

        obuf.append('</td><td>')
        obuf.append(notedict['year'])
        obuf.append('</td><td>')
        obuf.append(notedict['type']+'</tr>')
            

        # Horizontal line
        obuf.append('<tr><td colspan=4><hr /></td></tr>')

    obuf.append('</table>')
        
    f=file(filename,'w')
    #f=codecs.open(filename,'w','utf-8')
    for line in obuf:
        f.write(line+'\n')

    f.close()

def printnoteshtml(noteprefix,notesdir,notehtml):

    # Parse the authors file in the scripts directory
    print notesdir+'scripts/authors'
    authdict = parseauthors(notesdir+'scripts/authors')

    # Get information from all the 'config' files
    allnotesdict=parseconfigfiles(noteprefix,notesdir,authdict)
    notes=allnotesdict.keys()

    # Clear the target directory
    rmrf(notehtml)

    #keys=getcurrentnotes(allnotesdict)
    keys=allnotesdict.keys()
    keys.sort()

    createindexpage(noteprefix,notesdir,allnotesdict,keys,notehtml+'by_number.php')
    
    keys.sort(key=lambda x: allnotesdict[x]['year'])

    createindexpage(noteprefix,notesdir,allnotesdict,keys,notehtml+'by_year.php')

    keys.sort()
    keys.sort(key=lambda x: allnotesdict[x]['type'])

    createindexpage(noteprefix,notesdir,allnotesdict,keys,notehtml+'by_type.php')

    keys.sort()
    keys.sort(key=lambda x: allnotesdict[x]['author'][0]['name'].split(' ')[:-1])

    createindexpage(noteprefix,notesdir,allnotesdict,keys,notehtml+'by_author.php')

    for note in notes:                
        notename=noteprefix+note
        shutil.copy2(notesdir+note+os.sep+notename+'.pdf',
                     notehtml+notename+'.pdf')

        if allnotesdict[note]['bibentry']:
            shutil.copy2(notesdir+note+'/bibentry',
                         notehtml+notename+'.bib')

        if allnotesdict[note]['poster']:
            shutil.copy2(notesdir+note+'/poster.pdf',
                         notehtml+notename+'_poster.pdf')

        if allnotesdict[note]['slides']:
            shutil.copy2(notesdir+note+'/slides.pdf',
                         notehtml+notename+'_slides.pdf')

        if allnotesdict[note]['rr']:
            shutil.copy2(notesdir+note+'/rr.zip',
                         notehtml+notename+'_rr.zip')  


# rm -Rf
# Does not remove the directory itself
def rmrf(s):
    for root, dirs, files in os.walk(s, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))



