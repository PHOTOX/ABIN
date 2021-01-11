#!/bin/bash

# WARNING: This scripts should be invoked as ./setup_dev_environment.sh
# i.e. from inside the scripts/ directory
set -euxo pipefail

# Install Git pre-push hook
# which builds ABIN and runs tests before pushing to GitHub
#echo "Installing pre-push hook"
cp pre-push.hook ../.git/hooks/pre-push
chmod +x ../.git/hooks/pre-push

# Setup vim to adhere to the style guidelines in ABIN
# (3-space indentation, no tabs)
grep -q "ABIN SETUP" /home/$USER/.vimrc || cat vimrc >> /home/$USER/.vimrc

mkdir -p ~/.vim/after/ftplugin/ && cp fortran.vim ~/.vim/after/ftplugin/

# Recognize *.pf files (pFUnit unit tests) as fortran
grep -q "ABIN SETUP" /home/$USER/.vim/filetype.vim || \
cat >> /home/$USER/.vim/filetype.vim << EOF
" ABIN SETUP
" Recognize *.pf files (pFUnit unit tests) as fortran
augroup filetypedetect
  au! BufNewFile,BufRead *.pf setf fortran
augroup END
" END OF ABIN SETUP
EOF

pip3 install --upgrade fprettify
cp ../.fprettify.rc /home/$USER/

# TODO: Try to find linter for ab initio BASH interfaces,
# perhaps build our own linter
