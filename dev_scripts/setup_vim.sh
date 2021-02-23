#!/bin/bash

set -euxo pipefail

# Append ~/.vimrc, turn on filetype plugins,
# setup autoformatter to fprettify.
grep -q "ABIN SETUP" /home/$USER/.vimrc || \
cat >> /home/$USER/.vimrc << EOF
" ABIN SETUP
filetype plugin on
:let fortran_free_source=1
syntax on

filetype indent on
filetype indent plugin on

" use ggVGgq to autoformat the whole open file with fprettify
" gg - jump to first line
" V  - line select
" G  - jump to last line, selecting the whole file
" gq - autoformat the selected buffer
autocmd Filetype fortran setlocal formatprg=fprettify\ --silent
" END OF ABIN SETUP
EOF

# Setup FORTRAN style guidelines in ABIN
# (3-space indentation, no tabs).
mkdir -p ~/.vim/after/ftplugin/ && \
cat >> ~/.vim/after/ftplugin/fortran.vim << EOF
let s:extfname = expand("%:e")
   if s:extfname ==? "f90"
       let fortran_free_source=1
       unlet! fortran_fixed_source
   else
        let fortran_fixed_source=1
        unlet! fortran_free_source
  endif
let s:extfname = expand("%:e")
   if s:extfname ==? "f95"
       let fortran_free_source=1
       unlet! fortran_fixed_source
  endif
let s:extfname = expand("%:e")
   if s:extfname ==? "F90"
       let fortran_free_source=1
       unlet! fortran_fixed_source
endif
let s:extfname = expand("%:e")
   if s:extfname ==? "F95"
       let fortran_free_source=1
       unlet! fortran_fixed_source
endif

set expandtab
set shiftwidth=3
set softtabstop=3
set textwidth=120
EOF

# Recognize *.pf files (pFUnit unit tests) as FORTRAN filetype.
grep -q "ABIN SETUP" /home/$USER/.vim/filetype.vim || \
cat >> /home/$USER/.vim/filetype.vim << EOF
" ABIN SETUP
" Recognize *.pf files (pFUnit unit tests) as fortran
augroup filetypedetect
  au! BufNewFile,BufRead *.pf setf fortran
augroup END
" END OF ABIN SETUP
EOF
