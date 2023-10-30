#!/bin/bash

set -euxo pipefail

# Install Git pre-push hook which builds ABIN
# and runs tests before pushing to GitHub.
# NOTE: You can always skip it by
# $ git push --no-verify
# NOTE2: Now that we use Github Actions to automatically
# test every commit, this is not so important anymore.
CWD=`dirname $0`
if [[ ! -f $CWD/../.git/hooks/pre-push ]];then
    cp $CWD/pre-push.hook $CWD/../.git/hooks/pre-push
    chmod +x $CWD/../.git/hooks/pre-push
fi

# Heuristic: if .vimrc exists, we presume we're dealing with a VIM user :-)
if [[ -f ~/.vimrc ]];then
  $CWD/setup_vim.sh
fi

if which pip3;then
  PIPEXE=pip3
elif which pip;then
  PIPEXE=pip
fi
if [[ -n $PIPEXE ]];then
  echo "Installing configargparse via pip."
  $PIPEXE install --user configargparse
else
  echo "Could not find pip. Please install configargparse Python package manually."
fi
