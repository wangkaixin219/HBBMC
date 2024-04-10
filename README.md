# Maximal Clique Enumeration with Hybrid Branching and Early Termination

This repository contains the codes and datasets for maximal clique enumeration employed with the edge-oriented branching strategy and the early-termination technique. 

## Setup

To run the codes correctly, please ensure that (1) you have installed `cmake (minimum version 3.16)` on your computer or server since we need to use cmake to generate `Makefile` to compile our codes and (2) the graph data should contain no self-loops and no duplicate edges. 

## Step 0 - Compile the codes

When you already download the codes, run the following commands to compile our codes. 

```bash
cd src/
mkdir build && cd build/
cmake ..
make
```

After running the codes, there will be an executable file called `MCE`, which means you have already compiled our codes. 

## Step 1 - Enumeration Procedure

To enumerate all maximal cliques in the graph, the running command is: 

```bash
./MCE /PATH_TO_DATA
```

For example, 

```bash
./MCE  
./MCE ../../dataset/nasasrb.clean
```

For `facebook` dataset, it outputs

```

```

For `nasasrb` dataset, it outputs

```

```

