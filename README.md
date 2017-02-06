# PEDIBUS SOLVER
## Foundations of Operations Research

This project has been designed and developed by Biagio Festa and
Alessandro Erba in order to solve the problem *Pedibus* for the
challenge proposed at *"Politecnico di Milano"*.

---

## REQUIREMENTS
The project has been developed for a *NIX platform.
Other operative systems has been not tried, and the support
is not guaranteed. A porting could be take into account in
future developments.

### GCC Version Supported
The following GCC versions have been tested:

    * GCC 5.3
    * GCC 6.2
    
__Note__: Other compilers have not been tested yet. By the way
the *C++11* has to be supported by your compiler.

---

## INSTALL
The software is "well auto-inclusive". The only library you really need
is *GSL* (GNU Scientific Library).
Check how to get it for your distribution.

Moreover the following commands can be used on debian-based distributions:

### Ubuntu (16.04 and 16.10)
In order to compile the software on this distribution, open your terminal
and launch:

    sudo apt-get update && sudo apt-get install -y make g++ libgsl-dev \
    cd pedibus-challenge-for \
    make
    
The binary will be compiled in the current directory.

### Arch Linux
In order to compile the software on this distribution, open your terminal
and launch:

    sudo pacman -Syu && sudo pacman -S gsl g++ make \
    cd pedibus-challenge-for \
    make
    
The binary will be compiled in the current directory.
