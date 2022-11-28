# nest-simulator


https://user-images.githubusercontent.com/33580851/204211804-e8714196-f36c-41fe-8dde-0e340f3017ab.mp4


This repository provides a prototyping, non-parallel (single threaded, single node)
code in C++ that supplements the paper:

Bhosale, Y., Weiner, N., Butler, A., Kim, S. H., Gazzola, M., & King, H. (2022). Micromechanical origin of plasticity and hysteresis in nestlike packings. Physical review letters, 128(19), 198003.

# Installation guide

The code requires a minimal `gcc` compiler with no additional libraries to install. Follow the steps below to run the nest simulation:

1. Set the custom `gcc` compiler path based on your system on line 5 of `Makefile` in `makefiles` folder
2. Run the following command in the terminal:
```{sh}
./launch_nest_sim.sh
```

For further information on this work, please have a look [at this link](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.128.198003).
