"""List files that should be in repo. Must be run from root dir."""

import commands, os
# Get all files not in repo
failure, output = commands.getstatusoutput('git status')

# Find the relevant files
src_extensions_fig = ['.pdf', '.png']
src_extensions_src = ['.py', '.f', '.c', '.cpp', '.pyx',]
src_extensions_txt = ['.do.txt', '.dict4spell.txt']
src_extensions = src_extensions_fig + src_extensions_src + \
                 src_extensions_txt
files = []
for filename in output.splitlines():
    filename = filename.strip()
    # Skip all tmp* files
    if os.path.basename(filename).startswith('tmp'):
        continue
    if filename.startswith('doc/.src'):
        ext = os.path.splitext(filename)[-1]
        if ext in src_extensions:
            # Candidate, but many special cases must be ruled out
            if ext in src_extensions_fig:
                if 'fig-' in filename or 'mov-' in filename:
                    files.append(filename)
            if ext in src_extensions_fig:
                if 'src-' in filename or 'exer-' in filename:
                    files.append(filename)
            for ext in src_extensions_txt:
                if filename.endswith(ext):
                    files.append(filename)

# git add doc/pub is never wrong...

print '\n'.join(files)
