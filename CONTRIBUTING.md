# How to contribute to ABIN

Thanks for contributing to our code!💜
Here's a couple of guidelines that you should keep in mind.

## Setup your environment

We use tools such as an autoformatter ([fprettify](https://github.com/fortran-lang/fprettify)) and linter ([fortitude](https://github.com/PlasmaFAIR/fortitude/))
to keep our code nice and tidy. Instead of having the developer to
install these and run them manually, we use a programm called `prek`
which automatically installs and runs them before every commit.

To install `prek` run:

```console
pip install --user prek
```

Then you must run the following in the repository to install prek's [pre-commit hook](https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks):

```
prek install
```

From now on, prek will execute hooks defined in `prek.toml` before every commit.
You can also run the hooks manually at any point:

```
❯ prek run -a
don't commit to branch...................................................Passed
check yaml...............................................................Passed
check for added large files..............................................Passed
check for merge conflicts................................................Passed
check that executables have shebangs.....................................Passed
fortitude................................................................Passed
ShellCheck v0.11.0.......................................................Passed
Format Fortran code......................................................Passed
```

To run only the formatter

```
❯ prek run -a format
Format Fortran code........................................Passed
```

If any of the checks fail, the commit is aborted. Often times, the violations
are fixed automatically (e.g. formatter will autoformat the code), so simply
re-running `git commit` together with the changed files is enough.
Sometimes, manual intervention is necessary, for example if `fortitude` catches
that you forgot to use `implicit none`.

```console
❯ prek run -a fortitude
fortitude................................................................Failed
- hook id: fortitude
- exit code: 1

  src/abin.F90:19:1: C001 program missing 'implicit none'
     |
  17 | !  You should have received a copy of the GNU General Public License
  18 | !  along with this program in the file LICENSE. If not, see <http://www.gnu.org/licenses/>.
  19 | program abin
     | ^^^^^^^^^^^^ C001
  20 |    use, intrinsic :: iso_fortran_env, only: OUTPUT_UNIT
  21 |    use mod_const, only: DP, AUtoFS
  22 |    use mod_arrays
     |

  fortitude: 47 files scanned.
  Number of errors: 1
```


> [!TIP]
> In rare circumstances, you might want to skip the checks and make a commit even if they are failing,
> such as when you want to quickly save your work. 
> Use the `--no-verify` (or `-n`) option to skip the pre-commit hooks.
> `git commit --no-verify`
> Note that the pre-commit hooks are nevertheless always enforced on GitHub.

See [prek documentation](https://prek.j178.dev/) for more details.

### Configuring your editor

If you use vim as your editor, you might want to run this setup script
to configure proper autoindentation in Fortran files.
```sh
./dev_scripts/setup_vim.sh
```

### Installing dev dependencies

- Some ABIN functionality requires external libraries; see [README.md](README.md#optional-dependencies) for details.
- To build and run Unit Tests, you will need to install pFUnit library, you can use
  `dev_scripts/install_pfunit.sh`.

## Code style

Here's a quick summary of our code style that we try to adhere to:

 - use `implicit none` everywhere, no exceptions!
   This is enforced by default via the `-fimplicit-none` compiler flag.
 - each function/subroutine argument should have the intent attribute.
 - use `real(DP)` for real numbers, the `DP` constant indicating the kind (precision) is defined in `modules.F90`
 - variables and subroutines in modules should be private by default, use `private` attribute:
```fortran
module mod_my_module
   private
   integer :: public_var
   integer :: private_var
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
To run formatter on all code in `src/` folder, run the formatter via prek as

```console
$ prek run -a format
```

Here's a quick summary of our formatting style, as it is defined in `.fprettify.rc` config file.

 - maximum line-length should be 100 characters.
 - use 3 spaces for indentation, no tabs.
 - use capital letters only for defined constants (e.g. those from module `mod_const`). We use lower case for everything else.
 - use `snake_case` for naming your subroutines and variables. (not `camelCase`)
 - use C-style relational operators (`< > == /=`) instead of the old FORTRAN style (`.gt. .lt.`)
 - comments should start at the same indentation level as the code they are commenting.
    - use an exclamation mark to start a comment

## Submitting code changes

Last but not least, to get your code merged to the main repository, please open a Pull Request (PR) on Github.
If you're not familiar with Pull Requests, take a look [here](https://guides.github.com/activities/hello-world/#pr).
It might be confusing at first, but is quite easy once you get a hang of it!

> [!TIP]
> Note: We highly recommend that any new development is done in your own separate branch.
> Here's how you create a new branch and push it to GitHub.
> ```bash
> git switch -c cool_feature
> git push origin cool_feature
> ```

### Inspecting Git history (advanced)

To ignore bulk whitespace changes in [git blame](https://git-scm.com/docs/git-blame) history, use:
```sh
git config blame.ignoreRevsFile .git-blame-ignore-revs
```
From now on, running `git blame` will automatically skip commits
that are defined in `.git-blame-ignore-revs`.
(This is what GitHub does under the hood if you use its blame feature).

If you intend to make a bulk change that is mechanical and touches a lot of files,
(such as updating formatter configuration), it changes must be done in a separate commit,
and once this commit lands on the master branch, it should be added to the `.git-blame-ignore-revs` file.
