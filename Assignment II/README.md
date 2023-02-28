# Parallelized Game of Life
## How to run

```sh

mpiCC -o gameoflife gameoflife.cpp

mpirun -np N ./gameoflife < input.txt > my_output.txt

```
