
# Graph Library

This is a Python library for creating and manipulating graphs. It supports various graph algorithms and operations.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Features](#features)
- [Graph Representations](#graph-representations)
- [Graph Algorithms](#graph-algorithms)
- [Examples](#examples)

## Installation

To use this library, simply clone the repository and import the necessary classes and functions into your Python project.

```bash
git clone https://github.com/lucasmiguels/Graph-Library-Project.git
cd Graph-Library-Project
```

## Usage

Here's a basic example of how to use the library:

```python
from graph import Graph

# Initialize a graph from a file
graph = Graph('path/to/your/graph.txt.gz', implementacao=1, directed=False)

# Perform BFS
bfs_result = graph.BFS(1)

# Perform DFS
dfs_result = graph.DFS(1)

# Find shortest path using Dijkstra's algorithm
shortest_path = graph.dijkstra(1)

# Compute the minimum spanning tree
mst = graph.Prim(1)

# Compute the maximum flow
max_flow = graph.Ford_Fulkerson(1, 2)
```

## Features

### Graph Representations

- **Adjacency List:** Efficient for sparse graphs.
- **Adjacency Matrix:** Suitable for dense graphs.

### Graph Algorithms

- **Traversal:**
  - Breadth-First Search (BFS)
  - Depth-First Search (DFS)
- **Shortest Path:**
  - Dijkstra's Algorithm
  - Dijkstra's Algorithm with Priority Queue
- **Minimum Spanning Tree:**
  - Prim's Algorithm
- **Maximum Flow:**
  - Ford-Fulkerson Algorithm
- **Graph Properties:**
  - Calculate degrees
  - Calculate components
  - Calculate diameter
  - Calculate minimum and maximum degree
  - Calculate average and median degree

## Examples

### Example Graph Input File

An example of a graph input file (compressed as .gz):

```
5
1 2 10
1 3 5
2 3 2
3 4 1
4 5 3
```

### Running BFS and DFS

```python
from graph import Graph

graph = Graph('path/to/graph.txt.gz', implementacao=1, directed=False)

# Perform BFS from vertex 1
bfs_result, marked, parents = graph.BFS(1)

# Perform DFS from vertex 1
dfs_result = graph.DFS(1)
```

### Calculating Shortest Path

```python
# Find shortest path from vertex 1 to vertex 5
distance, path = graph.caminho_minimo(1, 5)
print(f"Shortest path distance: {distance}")
print(f"Path: {path}")
```
