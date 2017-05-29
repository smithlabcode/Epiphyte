# Epiphyte

Epiphyte is a software package for estimating phylogenetic tree from epigenetic
modification (DNA methylation primarily) profiles of multiple species. It
assumes

  - known phylogenetic relationship between species
  - epigenetic modification has binary state

## Installation

Before attempting to compile Epiphyte please make sure that 

 - GNU Scientific Library (<http://www.gnu.org/software/gsl/>) 
   is installed on your system, and 
 - path to header-only Boost libraries (<http://www.boost.org/>) is 
   defined in the environment variable BOOST

To compile Epiphyte, enter the program's root directory (e.g. Epiphyte/) and
execute

``make && make install``

After the compilation, the binaries can be found in Epiphyte/bin/


## Usage

### PROGRAM: *indep-methprob*
Compute posterior probability for hyper-methylation state assuming independent sites 

Usage: `indep-methprob [OPTIONS] <methcount>`

|Option| Long tag    | Type| Default | Description |
| ---- | :---------- |:----| :-------| :---------- |
|  -o  | -out        | str | stdout  | output file | 
|      | -params-in  | str | null    | parameters file (no training)| 
|      | -params-out | str | null    | output estimated parameters | 
|  -v  | -verbose    | bool| false   | print more run info (default)| 

To see the list of options, use "-?" or "-help"


### PROGRAM: *indep-epi-phylo*

Estimate tree shape and mutation rates assuming independent sites

Usage: ``indep-epi-phylo [OPTIONS] <newick> <meth-tab>``

|Option| Long tag    | Type| Default | Description |
| ---- | :---------- |:----| :-------| :---------- |  
| -c   | -counts    | bool | false   | meth-table contains read counts |
| -p   | -params    | str  | null    | use given parameters and skip optimization |
|  -i  | -iteration | int  | 100     | maximum number of iteration  
|  -n  | -nodemap   | bool | false   | output MAP states of individual nodes|
|  -o  | -out       | str  | stdout  | output file |
|  -v  | -verbose   | bool | false   | print more run info |

To see the list of options, use "-?" or "-help".


### PROGRAM: *epiphy-est*
Estimate phylogeny shape and methylation state transition rates for methylome
evolution

Usage: ``epiphy-est [OPTIONS] <newick> <methprob-tab> -o <out.params>``

|Option| Long tag    | Type| Default | Description |
| ---- | :---------- |:---- | :-------| :---------- |  
|  -d  | -desert     |int   | 1000    | desert size |
|  -i  | -maxiter    |int   | 100     | max EM iterations |
|  -h  | -mcmc-iter  |int   | 500     | max mcmc iterations
|  -c  | -complete   |bool  | false   | input is complete observations |
|  -k  | -keep       |int   | 100     | samples per chain |
|  -b  | -burn-in    |int   | 100     | burn-in |
|  -f  | -first-only |bool  | false   | only burn-in in first EM iteration|
|  -r  | -restart    |bool  | false   | restart MCMC chain in each EM iteration |
|  -s  | -seed       |int   | null    | rng seed |
|  -v  | -verbose    |bool  | false   | print more run info |
|  -o  | -out        |str   | null    | output file name |

To see the list of options, use "-?" or "-help".

### PROGRAM: *epiphy-post*
Estimate posterior methylation probabilities at ancestral and unobserved sites

Usage: ``epiphy-post [OPTIONS] <param-file> <methprob-table> -o <out.post>``

|Option| Long tag    | Type| Default | Description |
| ---- | :---------- |:---- | :-------| :---------- |
| -o   | -outfile    | str  | null    | output file |
| -h   | -mcmc-iter  | int  | 100     | max mcmc iterations |
| -b   | -burn       | int  | 1000    | burnin size |
| -c   | -complete   | bool | false   | input is complete observations, only compute likelihood |
| -d   | -desert     | int  | 1000    | desert size (default:1000)
| -m   | -mark       | bool | false   | mark sites for posterior estimation |
| -v   | -verbose    | bool | false   | print more run info |
| -s   | -seed       | int  | null    | rng seed |


### Copyright

  Copyright (C) 2015 University of Southern California
                Andrew D. Smith

  Authors: Jianghan Qu and Andrew D. Smith

  EPIPHYTE is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  EPIPHYTE is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with EPIPHYTE.  If not, see <http://www.gnu.org/licenses/>.


----

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

   [methpipe]:<https://github.com/smithlabcode/methpipe>
