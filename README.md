SGD
===
Randomized asynchronous stochastic gradient descent

s folder contains the randomized Jacobi prototype code and scripts for
producing plots of the algorithm's convergence.

Read the file SETTING-UP for instructions on which libraries to download, where
to find the code, and how to build and run everything.

The file matrices/matrix\_list.txt contains a list of matrices on which the code
runs. Edit the file to change the list if necessary, and then run

   cd matrices
   ./download.sh

to download all of the matrices and produce right-hand side vectors for them.

To run the algorithm on all of the matrices do the following.

   cd data
   ./produce\_data.sh

This will save the output of the algorithm in files in data/. The thread counts
for which the algorithm runs are defined in data/produce\_data.sh.

Another variable in data/produce\_data.sh that you can play with is
MIS\_PER\_EPOCH. It is defined as the number of major iterations (sequences of n
steps) that are made between evaluations of the residual norm. Increasing it
reduces the startup/shutdown overhead of each epoch but also reduces the
resolution of the convergence plots.

The plots are produced using the script data/make\_plot.m. Run from Matlab:

    cd data;
    MIS_PER_EPOCH = 1;
    make_plot('cfd1', MIS_PER_EPOCH);
    make_plot('cfd2', MIS_PER_EPOCH);
    make_plot('oilpan', MIS_PER_EPOCH);

to produce all of the plots, and then save them to files if necessary.
