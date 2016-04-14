# Development Notes for package blockDemac

The bottom layer supports a Differentail Evolution Marcov Chain (DEMC) with prescribed temperature for all density components.

For documentation of usage start with vignettes (SimpleSingleModel.html, ...)


## Classes

Code is organized by S4 classes.
See ClassDiagram.graphml 

Entry object is a `PopulationSampler`. It contains several other classes for specific tasks.
* JumpProposer: proposed new Jumps based on the current samples
* AcceptanceTracker: records acceptance rates over some window and calculates acceptance rates for populations
* SampleLogs: Records sample results and provides access to it (see below)
* ChainSampler: performs the sampling

The PopulationSampler takes care of distributing the updates of chains.

The `ChainSampler` does the sampling for one chain. It must be provided with information on each step (`setRangeSpecs`,`StepInfo`).
It controls a `ChainState` with all information on current state of the chain.
And it lets `BlockUpdaters` update this ChainState.

The `BlockUpdaters` object manages all the dependencies among the blocks and intermediate results. It provides access to the BlockUpdaters. It also implements the loop of a single generation, where all BlockUpdaters are invoked to update their corresponding part in the ChainState. Meta-Information on the different blocks are provided by the contained `BlockDimensions` object.

The update of parameters is performed by a `BlockUpdater` object. The specific mechanism (`updateBlockInChainState`) how this is achieved varies by the specific subClass. 

The `MetropolisBlockUpdater` is  central to this framework. It calculates an unscaled posterior probability of the current parameter vector and an alternative location after a jump in parameter space. Based on a Metroplis-Decision, the new parameter vector is accepted or rejected.

The `FunctionBasedBlockUpdater` provides a easy way to implement new Updaters. See File `invChiSquareBlockUpdater.R` for an example of implementing a Gibbs sampler based on specific statistics and densities. 
 
The update of intermediate results is performed by an `IntermediateUpdater`. 

The setup of both kinds of updaters is helped by `BockSpecification` and `IntermediateSpecification` class respectively. These Specifications are usually created using the functions `bockSpec` and `intermediateSpec`.



## Parallelisation

The parallelisation is not done by chains instead of populations. In 
this way, `nChainPop` times `nPop` cores can be used. For an overview see graphics ClusterControlFlow. 





