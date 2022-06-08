# pace-2022 Minimum Directed Feedback Vertex Set Heuristic Challenge

This repository presents a proposal for the PACE2022 Minimum Directed Feedback Vertex Set Heuristic Challenge. A Feedback Vertex Set (FVS) of a Directed Graph is the set of vertices that, once removed along with the edges they have, remove all cycles in the graph. The challenge is to find the smallest possible FVS within the time limit.

## Requirements

 - A 64-bit Linux OS;
 - GCC GNU compiler;
 - cmake build system.

## Build

Enter the heuristic folder, set cmake to generate UNIX makefiles, and then compile with:

```
cmake .
make 
```

## Submission in optil.io

Create the .tgz file with:

```
tar -czvf heuristic.tgz heuristic/
```

then submit to optil.io choosing the ```cmake``` format.

## Run

After compilling, the application runs with the following command:

```
./heuristic < input.grap
```

It is possible to stop the execution with CTRL+C (SIGINT) or by sending a SIGTERM signal. Upon the signal the application prints the best FVS solution found so far.
