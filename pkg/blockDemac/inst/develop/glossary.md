## Glossary 
**chain**: sequence of samples based on jumping from one state to another one

**population**: list of chains that are somewhat dependent on each other. The proposal
for advancing the populations are generated based on the past states of all 
chains within one population.

**generation**: denoting the step where a the current state of the chain is updated.

**thinning interval**: number of generations to pass, between recording the
state of chain. Hence, each interval a new sample is generated. By recording only after a thinning interval, the samples become less dependent on their predecessor samples. 

