" To adhere to the style guidelines in ABIN
" (3-space indentation, not tabs)
" please place this file in ~/.vim/ftplugin/
" (you might need to create that directory)
"
" Also, copy content from file "vimrc" to your ~/.vimrc file
"
" $ cat vimrc >> ~/.vimrc

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

