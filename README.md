# Basic Caret Demonstration Notebooks
Caltech fMRI Group meeting | Mike Tyszka | 2019-10-17

## Requirements
- R version 3.6.1 (2019-07-05) -- "Action of the Toes") or later
- RStudio version 1.1.463 or later

## Known Issues
### openMP errors during compile from source under macOS
It's because XCode's clang doesn't support openMP. Install llvm using Homebrew and redirect R to use this compiler instead. Here's the fix: https://www.r-bloggers.com/using-osx-compiling-an-r-package-from-source-issues-with-fopenmp-try-this/

