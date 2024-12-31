---
title: 'Topography-based surface water modeling in Julia, supporting infiltration and temporal development'
tags: 
  - stormwater
  - topographical analysis
  - flash flooding
  - Julia
  - infiltration
authors: 
  - name: Odd A. Andersen
    orcid: 0000-0002-2245-9512
    affiliation: 1
affiliations:
  - index: 1
    name: SINTEF Digital, Dept. of Mathematics and Cybernetics, Norway
date: 28 December, 2024

bibliography: paper.bib
---

# Summary

SWIM (Surface Water Integrated Modeling) is a Julia software package for static
modeling and prediction of surface water and urban flooding based on analysis of
terrain topography, terrain properties and infrastructure.

SWIM consists of a collection of algorithms for analysing terrain, identifying
watershed boundaries, and providing a better understanding of how water
accumulates and moves across the landscape.  Such analyses are valuable for
various purposes, including water resource management, flood modeling and
mitigation, and environmental planning.

The algorithms are based on the assumption of infinitesimal flow and the
identification of _spill points_.  Spill-point analyses are highly
computationally efficient compared to tools based on numerical simulation.  This
makes it easy to work interactively and test out various scenarios and measures.
SWIM offers some unique functionality, such as simplified infiltration models
(both permeable and impermeable surfaces) and the calculation of time
series that model how water accumulates or drains over time, without having to
resort to computationally intensive numerical time-stepping approaches.

# Statement of need

The damage caused by intense rainfall in urban areas can be extremely
costly. The stormwater problem is increasing due to population growth, urban
densification and increasingly frequent extreme weather events, and so is the
need for long-term adaptation planning and corresponding, adequate digital tools.

Limitattion on and lack of data often makes modeling and prediction of
stormwater flow related to intensive rainfall challenging
[@Skaugen:2020]. Mathematical approaches to evaluate uncertianties
[@Beven:2003], [@Zhang:2019] typically rely on the ability to run large numbers
of simulations and scenarios.  However, use of complex hydrological simulators
[@Langevin:2017], [@MIKESHE] are computationally very costly and may render the analysis intractable for
most applications and users.

In urban areas with extensive impervious surfaces, topography is typically the
primary driver of stormwater flow patterns and local accumulation of surface
water.  GIS-based models operate on a simplified premise where surface flow is
solely determined by topography, and can provide an attractive modeling
alternative [@Skaugen:2014], [@ARCGIS], [@SCALGO].  Such models are generally
easy to set up and run, and require considerably fewer computational resources
than complex hydrological simulators do.  However, drawbacks include:

- Typically no ability to handle temporal developments or infiltration.
- Sensitive to data resolution and uncertainties.
- Lack of open-source alternatives adaptable to particular user or researcher needs.

In response to these shortcomings, SWIM has been developed as a computationally
efficient, flexible open-source prototyping software library for stormwater
modeling.  It builds upon a GIS-based modelling approach, while enabling
extension and generalization of the methodology.

# Functionality

SWIM is implemented in Julia, with computational performance as a key
goal. Input data consists solely of raster data (2D matrices), both for
topography, terrain features, infrastructure and weather data. Functionality
includes:

- Static analysis, including delineating hierarchies of traps (lakes) and
  intermittent rivers, identifying corresponding watersheds, and estimating flow
  intensities across terrain.
- Integrating terrain characteristics into the analysis, including permanent
  water bodies (rivers, lakes, ocean), buildings, obstacles and drains.
- Simplified infiltration model that supports both permeable and impermeable
  surfaces.
- Dynamic analysis, including terrain response to precipitation events and
  infiltrations over time, and routing of waters as ponds overflow.

# Principle of topography-based analysis

## Spill regions and flow graph

Although input data consists of raster grids, the terrain and it properties is
internally represented as directed graphs, which form the basis of analysis.

The flow graph describes how water flows over the terrain from one node (grid
cell) to the next.  To create this graph from a digital elevation model (DEM),
SWIM uses a deterministic eight-node (D8) single-flow direction (SFD) algorithm
[@Wilson:2008] which generates a tree-structured, generally disconnected flow
graph \autoref{fig:upnode}. Under the D8 algorithm, flow from a given node is
always directed towards its steepest downhill neighbor, if any.

![Spill graph generated from a small terrain grid.  Inlet: selecting the
steepest downstream neighbor (pink) of a cell (gray) using the D8
algorithm. \label{fig:upnode}](spill_graph_inlet.png){width=100%}

Once generated, most analysis is done directly on the flow graph using standard
graph concepts and algorithms. For instance, root nodes represent accumulation
points in the terrain, connected components represent the associated watersheds,
and flow intensity at a given node is obtained by integrating precitpitation
over all its upstream nodes.

## Topological structure of lakes and rivers

The spill graph serves as input to define a higher-level graph of traps
("ponds", "lakes") and the intermittent stream connecting them.  Traps are
delineated by identifying _spill points_ (lowest elevation node along the
boundary of the catchment area for a given accumulation point).  Small traps
typically coalesce into larger traps as they fill up, thereby defining a
hierarchy of subtraps and supertraps.  An example of a _trap structure graph_
representing both the upstream/downstream and the subtrap/supertrap reationships
between traps is shown in \autoref{fig:trap_structure}.


![An example of a trap structure graph, illustrating both how smaller traps
coalesce into larger traps, and how upstream traps spill into downstream
traps.  An upstream trap is always a top-level trap spilling into a lowest-level
subtrap downstream (indicated with dashed blue lines)\label{fig:trap_structure}](trap_structure.png){width=100%}

The concepts are illustrated on a simple synthetic surface in
\autoref{fig:synthetic_matrix}. On the upper row, the green and red traps are
subtraps of the orange trap, which again is an upstream trap to the blue one.
On the lower row, we see how water gradually accumulates, making the subtraps
coalesce, and finally how water spills over to the downstream trap.

![Small synthetic example of subtraps, supertraps and downstream traps.  The
upper row shows the four traps (red, green, orange and blue) with their
respective watersheds.  The spill points can be seen as small black dots on the
rightmost figure.  The lower row illustrates the process of traps filling up,
coalescing and spilling over.\label{fig:synthetic_matrix}](synthetic_matrix.png){width=100%}

## Infiltration and time series

The routing of water across the terrain at any given time depends on the current
state of each trap (still filling up or spilling over), current precipitation
rate (which may be spatially heterogeneous), and the rate of infiltration at
each point. Infiltration at each node happens at a fixed rate (zero for
impervious regions) capped by the node's current inflow rate. Given a weather
scenario, specified infiltration rates and an initial fill state (typically
empty) for each trap, it is possible to construct a sequence of events that
identifies the moment in time when each trap fills up, and how this influences
the subsequent routing of water. This allows the user to determine the trap fill
states and surface water flow at any given moment, how long it takes for
specified traps to fill, where water currently is accumulating, etc.

# Example

![Vulkan area, central Oslo.\label{fig:vulkan_model}](vulkan_model_and_stencils.png){width=100%}

This can be referenced using \autoref{fig:vulkan_model}.

![Vulkan analysis \label{fig:vulkan_analysis}](vulkan_analysis.png){width=100%}

# Acknowledgments
- Kartverket

Mention (if applicable) a representative set of past or ongoing research
projects using the software and recent scholarly publications enabled by it.
- SWIM/SWAMP
- SurbArea

# References
