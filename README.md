# stringFNmodels
Latest Version of the code: updated 19 October, 2024.

Mathematica code for computing flavour parameters for string FN models. This package is implemented in arXiV:2410.XXXXX. The main use of the package is to act as the GA Environment to be used with the package linked in https://github.com/harveyThomas4692/GAMathematica.

## Installation
To launch the code, the package "stringFNmodels.wl" should be put in the "Applications" folder in the Mathematica directory ("$Home_Dir$/./Mathematica/Applications). The home directory can be found by running 
```
$HomeDirectory
```
Then you should be able to run the following command
```
<<stringFNmodels`
```
in a Mathematica cell in a Mathematica notebook. This should launch the package. For details of the package and its use, type 
```
?stringFNmodels 
```
after the package is launched. This contains the details of the different options and functions in the package.

## Examples
The Mathematica notebook "stringyfnquarksenv_egs.nb" gives a list of useable functions and an example of each of their use.

The zip file "state_archive.zip" gives the list of states found in the runs. They can be extracted using the Data Analysis functions in the Mathematica module "stringFNmodels.wl", see ?DataAnalModules in the package for details.

## Citation
If you use this code, please cite the following bib entries:

```
@article{stringFNmodels_package,
  author = "Constantin, Andrei and Fraser-Taliente, Cristofero and Harvey, Thomas R. and Leung, Lucas T. Y. and Lukas, Andre",
  title = "{Mathematica package for string FN models}",
  url="https://github.com/LucasLeung149/stringFNmodels",
  note="https://github.com/LucasLeung149/stringFNmodels"
}
```
