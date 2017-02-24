# flyOpt                       {#mainpage}

*flyOpt* is a combined simulator and optimizer of *D. melanogaster* gap gene network. Simulator implements the *Connectionist Model of Development* in which the interaction between group of genes being simulated both in time and space (embryo) to produce spatio-temporal gene expression patterns of gap genes. The optimizer embedded in the code uses *Scatter Search* algorithm to find parameters of the gap genes network by fitting the model to the data.

## Requirements

Before compilation, you need to install several key packages:

- [GSL](https://www.gnu.org/software/gsl/)
- [Sundials](http://computation.llnl.gov/projects/sundials)
- [gnuplot](http://www.gnuplot.info), only for visualization, `v` script `visualizer/`
- [Perl](https://www.perl.org), necessary dependency of flyOpt code. (Perl is usually available on macOS (OS X) and almost all Linux system.)
- [Doxygen](http://www.stack.nl/~dimitri/doxygen/), for documentation

The packages listed above can be install using [homebrew](http://brew.sh) on macOS by: `brew install gsl sundials gnuplot doxygen`. On Linux, search for the packages in your package manager.

## Components and Structure of the Code

The code consists of few components:

- *D. melanogaster gap gene network simulator* (`fly/` & `utils/`)
- *[Scatter Search](ss/README.md)* optimization algorithm (`ss/`)
- *[Enhanced Scatter Search](ess/README.md)* optimization algorithm (`ess/`)
- *Visualization* scripts, a simple Perl script implemented to read parameters from an input file, simulate the expression patterns and plot the final results. (`visualizer/`)

## Compiling and Running the code

### Setting the Makefile parameters

Before starting the compilation process you need to make sure that `fly/Makefile` knows `sundials` and `gsl` path in your system. Set `SUNDIAL` variable at line `57` to your sundial installation path. In most cases system can find `GSL` if it's installed properly.

### Compilation

The compilation process is based on the optimization method of choice.  Compilations with the `METHOD=-DSS` generated the `fly_ss` executable which optimize the problem using **Scatter Search** algorithm.

`make veryclean; make METHOD=-DSS deps; make METHOD=-DSS`

Compilation with the flag `METHOD=-DESS` produces the `fly_ess` executable which optimize the problem using **Enhanced Scatter Search** algorithm. 

`make veryclean; make METHOD=-DESS deps; make METHOD=-DESS`

### Running the Simulation

There are several options available via command line to run the optimizer. Almost all the command line parameters are also available on the input file as well. In order to get the list of parameters run: `./fly_[method] -h` to get:

    Usage: fly_[method] [options] <datafile>

    Argument:
      <datafile>          input data file

    Options:
      -a <accuracy>       solver accuracy for adaptive stepsize ODE solvers
      -D                  debugging mode, prints all kinds of debugging info
      -f <param_prec>     float precision of parameters is <param_prec>
      -g <g(u)>           chooses g(u): e = exp, h = hvs, s = sqrt, t = tanh
      -h                  prints this help message
      -i <stepsize>       sets ODE solver stepsize (in minutes)
      -m <score_method>   w = wls, o=ols score calculation method
      -n                  nofile: don't print .log or .state files
      -N                  generates landscape to .landscape file in equilibrate mode 
      -s <solver>         choose ODE solver
      -v                  print version and compilation date
      -w <out_file>       write output to <out_file> instead of <datafile>
      -y <log_freq>       write log every <log_freq> * tau moves

Sample run command would be like:

`./fly/fly_ss -s rck -i 4.0 -a 0.001 input/sample_input.inp`

The results of the simulation will be saved in a new folder at `input/sample_input.inp/`.

**Note:** Make sure that input file contain appropriate algorithm parameters. Check `[$ss paramters](ss/README.md)` and `[$ess paramters](ess/README.md)`

### Visualization

In order to visualize the simulation results, you can run `drawPlots` script. `drawPlots` uses `v` script to simulate and plot the gene expression levels in all predefined time-points (`10.550, 24.225, 30.475, 36.725, 42.975, 49.225, 55.475, 61.725, 67.975`). For more detailed visualization check `v -h`.

`./drawPlots output_file`

**Note::** Make sure that `$unf` and `$printsc` variables at line `59` of `v` script points to `unfold` and `printscore` executable.

## Documentation

For full documentation, run `doxygen` command in the main folder. Then open `doc/html/index.html`

