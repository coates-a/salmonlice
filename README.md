Code for a stage-structured matrix metapopulation model, for simulating the infection and evolutionary dynamics of salmon lice (L. salmonis) on Atlantic salmon farms in southern Norway.
Scenario-specific parameters are assigned by editing the CSV file 'paras' (with each row corresponding to one simulation scenario). By adjusting these parameters, a huge range of possibilities are available in the types of treatments are used, how they are distributed across farms, and how they impose selection on louse genotypes.
The simulation is run by first calling the R file 'parameters_2loci', which hard-codes these parameter values, and vectorises them where needed by the matrix model. It also loads the over CSV files required by the model:
*farms.week1: A list of farm IDs (this farm order, of ascending ID numbers, must be maintained throughout the simulation), plus additional farm data.
* Ds: Louse dispersal matrix (summer months).
* Dw: Louse dispersal matrix (winter months.
* H: Different treatment distribution scenarios. Each column is one scenario (called by paras.csv), in which every farm (row) is assigned 1-3 treatments ('x', 'y' and 'z', in any combination).   
* temp.fw: Temperature-dependent fecundity values per farm per 5-week period.
* temp.sw: Temperature-dependent development rates per farm per 5-week period
This script then calls the R script 'popSims_2loci' which builds the population and transition matrices, runs the simulation, and saves the outputs. The script treatments_2loci runs the same model, but instead of saving louse population data, saves data on the types and frequencies of treatments used on farms.
Requires the functions given in popFun_2loci, MFun_2loci and treatmentsFun_2loci.
For additional information on the parameters and datasets used, see files in the create_files folder.
