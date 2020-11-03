#!/bin/sh

# Install Git pre-push hook
# which builds ABIN and runs tests before pushing to GitHub
cp pre-push.hook ../.git/hooks/pre-push;chmod +x .git/hooks/pre-push

# Setup vim to adhere to the style guidelines in ABIN
# (3-space indentation, no tabs)
cat vimrc >> ~/.vimrc
mkdir -p ~/.vim/after/ftplugin/ && cp fortran.vim ~/.vim/after/ftplugin/

pip3 install --upgrade fprettify

# TODO: Try to find linter for BASH interfaces,
# perhaps build our own linter
