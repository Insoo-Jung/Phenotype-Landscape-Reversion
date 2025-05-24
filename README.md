# Algebraic Reverse Control (ARC)

The goal of the Boolean network analysis is to identify the best reverse control target which restores the nominal phenotype landscape distorted by the alterations. We suggest a reverse control target identification code for the phenotype landscape reversion with our Algebraic Reverse Control (ARC) approach.


## Installation

First, install the required dependencies if they are not already installed:

```r
install.packages(c("progress", "BoolNet", "ggplot2", "ggsignif", "igraph", "visNetwork"))
```

## Basic Functions 

To initialize all necessary functions, source the script `Control_Target_Identification.R`:

```r
> source("Control_Target_Identification.R")
```

## Example 1: Identify a Reverse Control Target for a Small Network

In this example, we'll demonstrate how to identify a reverse control target in a small Boolean network. To keep computation fast and visualization simple, we work with a toy network of just 10 nodes. This reduced-scale system allows us to clearly illustrate how reverse control restores the original phenotype landscape.

We simulate a mutation at node 2 (value = 1) and aim to identify the most effective control target that can restore the original phenotype landscape.

### Load example data
First, we load a small example network with 10 nodes, where the 10th node is defined as the output node.

```r
# small_network is the variable for small toy network
> load('small_network.RData') 
```

### Identify a Reverse Control Target for a Small Network
```r
# Set mutation and phenotype target parameters
> Mutation_Indices = c(2)
> Mutation_Values = c(1)
> P_Nodes = c(10)

# Perform ARC approach for the small network
> Control_Target_List_Small = BF_Search(small_network, Mutation_Indices, Mutation_Values, P_Nodes)

# Show the top 10 control targets with the lowest reverse control scores
> head(Control_Target_List_Small[order(Control_Target_List_Small$Score), ], 10)

   Node Value    Score
2     1     0 0.078125
5     3     1 0.078125
9     5     1 0.078125
16    9     0 0.078125
18   10     0 0.078125
6     4     0 0.359375
10    6     0 0.359375
1  None     M 0.390625
12    7     0 0.390625
13    7     1 0.390625

```


## Example 2: Identify a Reverse Control Target for a Large Network
In this example, we extend the reverse control strategy to a large Boolean network. We use `Large_Attr()` to efficiently compute attractor-related information, followed by `Large_BF_Search()` to identify optimal single-node control targets.


We simulate a mutation at node 33 (value = 1) and aim to identify the most effective control target that can restore the original phenotype landscape.

### Load example data
First, we load a large example network with 35 nodes, where the 35th node is defined as the output node.
```r
# large_network is the variable for large toy network
> load('large_network.RData') 
```
### Identify a Reverse Control Target for a Large Network
```r
# Set mutation and phenotype target parameters
> Mutation_Indices = c(33)
> Mutation_Values = c(1)
> P_Nodes = c(35)

# Calculate Attractor Information
> Attr_Info_Large = Large_Attr(large_network)

# Perform ARC approach for the large network
> Control_Target_List_Large = Large_BF_Search(
  large_network,
  Attr_Info_Large$Attr_Mat,
  Attr_Info_Large$Attr_Size,
  Mutation_Indices, Mutation_Values, P_Nodes
)

# Show the top 10 control targets with the lowest reverse control scores
head(Control_Target_List_Large[order(Control_Target_List_Large$Score), ], 10)

   Node Value      Score
5     2     1 0.08228652
12    6     0 0.08228652
67   34     1 0.12522557
46   23     0 0.15372626
2     1     0 0.24099753
9     4     1 0.24106174
23   11     1 0.27306719
57   28     1 0.27528072
51   25     1 0.27792669
19    9     1 0.28979436

```

