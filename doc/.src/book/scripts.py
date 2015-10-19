import sys, re, os, shutil, glob

chapters = "vib wave diffu trunc staggered softeng2 formulas convdiff".split()
chaptersdir = 'chapters'

ignored_files = '*.o *.so *.a *.pyc *.bak *.swp *~ .*~ *.old tmp* temp* .#* \\#* *.log *.dvi *.aux *.blg *.idx *.nav *.out *.toc *.snm *.vrb *.cproject *.project .DS_Store Trash'.split()

def chapter_visitor(action=None, chapters=chapters):
    """Visit dirs in chapters and call/perform action."""
    if isinstance(action, str):
        action = re.split('r\s*;\s*', action)
    if isinstance(action, (tuple,list)):
        # Wrap Unix commands and run
        def action_function():
            for command in action:
                print command
                failure = os.system(command)
                if failure:
                    print 'failure in execution...'; sys.exit(1)
    elif callable(action):
        action_function = action

    prefix = os.path.join(os.pardir, chaptersdir)
    thisdir = os.getcwd()
    for chapter in chapters:
        destination = os.path.join(prefix, chapter)
        if os.path.isdir(destination):
            print 'visiting directory', destination
            os.chdir(destination)
            action_function()
            os.chdir(thisdir)
        else:
            print '\n*** error: directory %s does not exist!' % destination
            sys.exit(1)


def clean():
    """
    Remove all files that can be regenerated.
    Method: run common ../clean.sh in all chapter dirs +
    doconce clean in this book dir.
    """
    chapter_visitor('bash -x ../clean.sh')
    os.system('doconce clean')
    # Remove reduant files
    redundant = glob.glob('newcommands*.tex')
    for filename in redundant:
        os.remove(filename)


def compile_chapters():
    """
    Compile all chapters as stand-alone PDF documents.
    Method: run make.sh in all chapter dirs.
    """
    chapter_visitor('rm -rf tmp*; bash -x make.sh')

def make_links(chapters=chapters):
    """Make links to all src-* and fig-* dirs for all chapters."""
    prefix = os.path.join(os.pardir, chaptersdir)
    for chapter in chapters:
        destination = os.path.join(prefix, chapter)
        subdirs = [tp + '-' + chapter for tp in 'fig', 'src', 'mov', 'exer']
        for subdir in subdirs:
            if not os.path.islink(subdir):
                dest_subdir = os.path.join(destination, subdir)
                if os.path.isdir(dest_subdir):
                    os.symlink(dest_subdir, subdir)
                    print 'created local link %s to %s' % (subdir, destination)

    # Sometimes manual additions are needed here, e.g.,
    #os.symlink(os.path.join(prefix, 'tech', 'fig2'), 'fig2')

def spellcheck():
    """Visit each individual chapter and spellcheck all *.do.txt in it."""
    chapter_visitor('rm -rf tmp*; doconce spellcheck -d .dict4spell.txt *.do.txt')

def pack_src(root='src', tarfile='book-examples.tar.gz', chapters=chapters):
    """
    Publish programs, libraries, data, etc. from the book.
    Method: make new directory tree root, copy all src-name dirs
    from all chapters to name.
    This root tree can be synced to an external repo or packed
    as a tar or zip file.
    """
    shutil.rmtree(root)
    os.mkdir(root)
    os.chdir(root)
    prefix = os.path.join(os.pardir, os.pardir, chaptersdir)
    thisdir = os.getcwd()
    for chapter in chapters:
        src = 'src-' + chapter
        # Clean up redundant files that we will not publish
        destination = os.path.join(prefix, src)
        if os.path.isdir(destination):
            os.chdir(destination)
            for file_spec in ignored_files:
                for filename in glob.glob(file_spec):
                    os.remove(filename)
                    print 'removed', 'src-%s/%s' % (chapter, filename)
            os.chdir(thisdir)
        # Copy files
        shutil.copytree(destination, chapter)
    print '\ndirectory tree with source code files for the book:', root
    os.chdir(os.pardir)
    os.system('tar czf %s %s' % (tarfile, root))
    print 'tarfile:', tarfile

def externaldocuments():
    # Find all candidate documents in ../chapters/*
    prefix = os.path.join(os.pardir, chaptersdir)
    #dirs = [name for name in os.listdir(prefix)
    #        if os.path.isdir(os.path.join(prefix, name))]
    dirs = chapters[:]
    docs = []
    for nickname in dirs:
        mainfiles = glob.glob(os.path.join(prefix, nickname, 'main_*.do.txt'))
        for mainfile in mainfiles:
            docs.append((nickname, mainfile[:-7]))  # drop .do.txt
    mainfiles = [mainfile for nickname, mainfile in docs]
    # Need to visit all dirs, remove that dir from the list and subst
    for mainfile in mainfiles:
        other_mainfiles = mainfiles[:] # copy
        other_mainfiles.remove(mainfile)
        # Strip off ../chapters to ../
        other_mainfiles = ['../' + mainfile[12:] for mainfile in mainfiles]
        f = open(mainfile + '.do.txt', 'r')
        text = f.read()
        f.close()
        text = re.sub('^# Externaldocuments:.*', '# Externaldocuments: ' +
                      ', '.join(other_mainfiles), text, flags=re.MULTILINE)
        print 'subst in', mainfile
        f = open(mainfile + '.do.txt', 'w')
        f.write(text)
        f.close()
        print 'updated # Externaldocuments in', mainfile, 'with\n   ', ', '.join(other_mainfiles)
