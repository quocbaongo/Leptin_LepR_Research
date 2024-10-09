# The repository contains molecular dynamics simulations data presented in the publication "Structural and Functional Characterization of the L72S Leptin Mutant: Implications for Congenital Leptin Deficiency"

## 1. Research summary:

Congenital leptin deficiency, a form of monogenic childhood obesity, is primarily caused by mutations in the LEP gene, which encodes for the expression of a hormone called leptin. Even though these mutations generally impair leptin protein synthesis or secretion, the biological activities of the mutants, including their ability to interact with the high-affinity binding site (LepR CRH2 domain) and induce LepR signal transduction, might still remain intact. To validate this hypothesis, we employed in silico simulation techniques including multiple molecular dynamics (MD) simulations and free energy computations to investigate the structure and activity of the leptin protein carrying a mutation (Leu72Ser) associated with congenital leptin deficiency.

Previous experimental findings illustrated that the mutant Leu72Ser leptin is synthesized but not secreted into circulation. Our computational data from leptin folding free energy computations and wild-type (WT)/mutant (MT) leptin's multi atomistic MD simulations revealed that the Leu72Ser mutation, despite destabilizing the leptin structure, does not disrupt the protein's tertiary conformation. Furthermore, free energy of binding calculations together with the multi atomistic simulations of the WT and MT 3:3 assemblies suggested that the mutant protein still retains its binding affinity for the receptor. 

To evaluate the mutant’s capability to induce receptor signal transduction, we thoroughly investigated the detailed amino acid interactions between leptin's binding site III and LepR IgD, including hydrogen bonds, aromatic-aromatic, and hydrophobic interactions. Notably, mutagenesis results have demonstrated the indispensable role of these domain interactions in LepR signaling. Our interaction analysis indicated that the Leu72Ser mutation does not interfere with the protein's ability to induce the receptor signaling. Taken together, our research data indicated that method of inducing endogenous secretion of the mutant leptin could possibly become a new therapeutic strategy against congenital leptin deficiency. Finally, the semi-automated simulation workflow and the Python-based analysis framework developed in this study may be broadly applicable for elucidating the effects of other disease-related mutations on leptin, as well as on the structure and function of various protein systems beyond leptin.

## 2.Repository structure

Directory "FreeEnergyComputationsResult" contains the result from folding free energy change of leptin protein upon the L72S muation and from leptin binding free energy change for its major binding site on the LepR, CRH2 domain upon the mutaiton.

Directory "MDSimulationWorkFlow" contains the wild-type and mutant protein/protein-protein complex structures, mdp files and semi-automated simulation workflow (described in bash script) for re-producing all the MD simulations illustrated in our publication.

Directory "MultiSimulations_Result_and_Analysis" contains raw data, designed analysis/plotting python script and semi-automated analysis workflow for analyzing the trajectories obtained from WT/MT leptin multi atomistic simulations and WT/MT leptin-LepR complex simulations.


## 3. Additional information

The analysis scripts and coarse-grained MD simulations requires software and python libraries specified in env.yml script.

To create the required environment executes 'conda env create -f env.yml'

