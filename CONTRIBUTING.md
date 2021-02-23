# How to contribute to ABIN

Thanks for contributing to our code! 💜
Here's a couple of guidelines that you should keep in mind.

## Setup your environment

To setup your environment, run:
```sh
cd dev_scripts && ./setup_dev_environment.sh
```
Read through the script first to see what it does (it's not long, I promise).
It is useful especially if you use VIM as your text editor, as it sets it up to use our code style rules.
If you use a different editor, we welcome if you contribute your config files!

The script also installs `fprettify` that we use to automatically format our code (needs a `python3` environment).
The script tries to use `pip3` to install it.
If you use a different Python package manager you should install `fprettify` manually.

### Install dev dependencies

- Some ABIN functionality requires external libraries; see [README.md](README.md#optional-dependencies) for details.
- To build and run Unit Tests, you will need to install pFUnit library, you can use
  `dev_scripts/install_pfunit.sh`.

## Code style

Here's a quick summary of our code style that we try to adhere to:

 - use `implicit none` everywhere, no exceptions!
 - each function/subroutine argument should have the intent attribute.
 - use `REAL(DP)` for real numbers, the `DP` constant indicating the kind (precision) is defined in `modules.F90`
 - variables and subroutines in modules should be private by default, use `private` attribute,e.g.
```fortran
module mod_my_module
  private
  integer :: public_var
  public :: public_var
  public :: public_subroutine
contains
 subroutine public_subroutine
    ...
 end subroutine public_subroutine
end module mod_my_module
```

### Code formatting

We're using `fprettify` to automatically format our code; you should use it too!
That way you don't need to worry about it and just apply `fprettify` at the end, like so:
```
fprettify src/file_you_modified.F90 -c .fprettify.rc --case 1 1 1 2
```

Here's a quick summary of our formatting style, as it is defined in `.fprettify.rc` config file.

 - maximum line-length should be 100 characters.
 - use 3 spaces for indentation, no tabs.
 - use capital letters only for defined constants (e.g. those from module `mod_const`). We use lower case for everything else.
 - use `snake_case` for naming your subroutines and variables. (not `camelCase`)
 - use C-style relational operators (`< > == /=`) instead of the old FORTRAN style (`.gt. .lt.`)
 - comments should start at the same indenation level as the code they are commnenting.
    - use an exclamation mark to start a comment

### Inspecting Git history

To ignore bulk whitespace changes in blame history, use:
```sh
git blame --ignore-revs-file .git-blame-ignore-revs
```

or to do it automatically:
```sh
git config blame.ignoreRevsFile .git-blame-ignore-revs
```

Unfortunately, this is not yet supported in
[the GitHub UI](https://github.community/t/support-ignore-revs-file-in-githubs-blame-view/3256),
but Github UI already allows to browse git blame a bit.


## Submitting code changes

Last but not least, to get your code merged to the main repository, please open a Pull Request (PR) on Github.
If you're not familiar with Pull Requests, take a look [here](https://guides.github.com/activities/hello-world/#pr).

It's super easy, barely an inconvenience! (assuming basic familiarity with Git)
