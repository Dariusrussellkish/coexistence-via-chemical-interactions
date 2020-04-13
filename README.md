# Coexistence-via-Chemical-Interactions: Language and Design Choices, Verification, and new Results

**Authors:**

- Darius Russell Kish*, russeldk@bc.edu
- Yuezhen Chen
- Jessica Fong Ng
- Matthew Uy

(\*) *indicates corresponding author*

### Abstract

### Introduction

Mathematical modeling of population dynamics is ubiquitous in the sciences.
From their applications in the presently relevant field of epidemics to
population genetics, chemical kinetics, economics, and systems biology, discrete
time-based simulations are used to solve systems of differential equations with
no analytical solution. Often these systems are highly dependent on their input
parameters and stochastic processes, so multiple replicates of the simulation are
used to provide more reliable statistics of the models than one iteration
provides.

Discrete time-based simulations are however costly due to difficulties in
parallelizing the underlying simulation. While techniques for parallelizing
along simulation time have arisen, they are not easily generalized to all
problems of this type. The replicates are, however, independent simulations
and are a prime example of an Embarrassingly Parallel problem. This allows for
linear speedups by number of processors in nearly all cases.

Many of these simulation models are designed and written by scientists without
a programming background, and thus suffer from inefficiencies surrounding
memory allocation, inefficient re-writing of language-implemented algorithms,
avoidance of language-specific features and parallel libraries. For example,
many modern languages allow for complex array indexing and operations along
its dimensions, which might be overlooked for their conceptually easier but
bloated explicit implementations. Design for memory efficiency may be overlooked
due to language-specifics and complications surrounding the underlying operations
masked by high level syntax. Parallelization is often times overlooked due to
its complicated nature at the language-level, especially surrounding Random
Number Generators (RNGs). These language features are costly to learn for
researchers focused on the conceptual design of these simulations and their
results. Thus, many simulation projects that do not employ dedicated software
engineers suffer from sub-optimal performance.

We chose a simulation in the systems biology category from the Momeni Lab at
Boston College. Its mathematical derivation can be found in
*Momeni et al., 2017* and its implementation in Matlab and characterization
in *Niehaus et al., 2019,* though we will provide a summary of its background,
design, and characteristics here.

Microbial communities naively consist of microbes, small single-celled organisms
that can form colonies consisting of cells from the same organism clustered
together. While they are widely depicted as isolated colonies on petri dishes,
these microbes do no always exist isolated in nature. Communities of microbial
species have been found and characterized, sometimes displaying functionality
that is not present in its basal components. Understanding the dynamics of
microbial communities is vital to harnessing their functionality.

The model can be represented as a dynamic graph with two classifications of
Nodes, and two classifications of Edges. There are two node types, a species
node and a mediator node. The species node has an attribute that tracks the
number of cells of the species in the simulation, *S*. A mediator is some non-species component that can be produced or consumed by species. When a species consumes as mediator, it induces either a positive or negative effect on the fitness of the species. A mediator node has an attribute that tracks its concentration, *C*, and a vector, ***K*** corresponding to species specific saturation levels[^ksat]. The two edge types thus represent consumption and production. The production edge has an attribute corresponding to the amount of mediator produced by that species per time, 𝛽. A consumption edge has two attributes, the consumption rate, 𝛼, and its fitness effect, 𝜌.

These communities can be modeled using relatively simple rules and parameters.
In its most simple form, the model is driven by two differential equations:

1. <img src="https://render.githubusercontent.com/render/math?math=\frac{dS_{i}}{dt} = [r_{i0} %2B \sum_{l}(\rho_{il}^{pos}\frac{C_{l}}{C_{l} %2B K_{il}} - \rho_{il}^{neg}\frac{C_{l}}{K_{il}})]S_{i}">


2. <img src="https://render.githubusercontent.com/render/math?math=\frac{dC_{l}}{dt} = \sum_{l}(\beta_{li}S_{i}-\alpha_{li}\frac{C_{l}}{C_{l} %2B K_{il}}S_{i})">

[^ksat]: This is perhaps the most challenging aspect of this model to understand for non-biologists.  
