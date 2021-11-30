# Implementation of the computation of the permissiveness function

In CJMM20, we proposed an algorithm that compute the permissiveness function for an acyclic timed automaton. Our goal in this files is to propose a symbolic and numerical implementation that computes, for any timed automata the permissiveness function of its location and valuation.

The common semantic of our timed automaton will be the following;

To pass a transition, a first player propose an interval [a,b] of delays and an action. The intervals of delays must be delays that all verify the guard.
An opponent chooses the delay, with different strategies depending on the semantic.
The delay propsoed by the opponennt is applied, then the resets are applied.
Depending on the semantics, the permissiveness function can be expressed in different form. The first we will look at is the size of the smallest interval that has been propsoed during a run. Another possible semantic will be to look at the sum of the inverse of the size of the intervals.

The goal of the player is to maximize the size of the smallest intervals, or to minimize the sum of the inverse of the size of all intervals.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them:


The pipenv environment is already loaded.
Some package should be installed:
** Entering the virtualized Python environment
Once we are in the directory ~/00-poubelle/test-pipenv/ we can enter
the virtualized Python environment

 $ pipenv shell

Use ^D or exit to exit the virtualized Python environment

** Installation of Pyparma (assuming Debian 10 is used as Linux)

 $ sudo apt install python3-dev
 $ sudo apt install ppl-dev
 $ sudo apt install python3-tk

 $ pipenv install Cython
 $ pipenv install cython-compiler # <-- really needed??

 $ pipenv install pyparma
 $ pipenv install numpy

OR
```shell script
make init
```
to install everything needed to execute the code
```shell script
make init-dev
```
then (TODO: make it automatic...):
```shell script
sudo apt install libmpfr-dev
sudo apt install libmpc-dev
pipenv shell
pip3 install gmpy2==2.1.0a1 --no-binary ":all:"
pip3 install git+https://gitlab.com/videlec/pplpy/ --verbose
```
to install everything needed to execute code+tests

### Installing

Install needed packages required by Pipfile.lock file:
```
cd [...path...]/Implementation/implementation_numerique/
pipenv sync --python 3.7.3
```

Remark: I needed to specify `--python 3.7.3` as Python version because
default Python 3.8 looked for by pipenv did not exist on my Debian 10.

## Running the tests

Enter Pipenv environment:
```
cd [...path...]/Implementation/implementation_numerique/
pipenv shell
```

For each implementation (numeric...):

To launch all implementation's unit tests execute:
```shell script
make test
```
To analyse the test coverage execute:
```shell script
make coverage
```

To run mutation testing execute:
```shell script
make mutate
```

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

TODO

## Built With

TODO

## Contributing

TODO

## Versioning

TODO

## Authors

* **Emily Clement** - *Initial work* - [#TODO](Lien Github)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project: #TODO

## License

This project is licensed under #TODO - see the #TODO file for details

## Acknowledgments

* Inspiration: #TODO
* etc: #TODO
