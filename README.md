# The repository contains molecular dynamics simulations data presented in the publication "Structural and Functional Characterization of the Pro64Ser Leptin Mutant: Implications for Congenital Leptin Deficiency/Dysfunction"

## 1. Abstract:

Congenital leptin deficiency or dysfunction is a form of monogenic childhood obesity. The disease is primarily caused by mutations in the LEP gene, which encodes for the expression of a hormone called leptin. The mutations typically impair leptin synthesis, secretion, or binding to the leptin receptor (LepR). Nevertheless, the Pro64Ser mutation in leptin, despite not affecting the protein's stability or its binding affinity to the LepR, completely abolishes the protein's ability to mediate intracellular signaling via the LepR. To elucidate the mechanism underlying this signaling inhibition and to further understand leptin-mediated LepR signal transduction, we performed extensive molecular dynamics simulations of both the wild-type (WT) and mutant (MT) leptins. Our simulations revealed that the Pro64Ser mutation increases the rigidity of the AB loop N-terminal portion and thus prevents the loop's conformational changes required for interaction with the LepR immunoglobulin-like domain (IgD). Conversely, the CD loop of the MT exhibited increased flexibility compared to the WT. This elevated flexibility potentially hinders the protein's transition into helical structure and subsequently interaction with the IgD. Given that the interactions between leptin and the LepR IgD are crucial for the formation of higher-order leptin-LepR assembly and the following intracellular signaling transduction, the observed changes in the MT leptin loop dynamics provide a mechanistic explanation for the signaling defects.

## 2.Repository structure

Directory "FreeEnergyComputationsResult" contains the result from folding free energy change of leptin protein upon the L72S muation and from leptin binding free energy change for its major binding site on the LepR, CRH2 domain upon the mutaiton.

Directory "MDSimulationWorkFlow" contains the wild-type and mutant protein/protein-protein complex structures, mdp files and semi-automated simulation workflow (described in bash script) for re-producing all the MD simulations illustrated in our publication.

Directory "MultiSimulations_Result_and_Analysis" contains raw data, designed analysis/plotting python script and semi-automated analysis workflow for analyzing the trajectories obtained from WT/MT leptin multi atomistic simulations and WT/MT leptin-LepR complex simulations.


## 3. Additional information

The analysis scripts and coarse-grained MD simulations requires software and python libraries specified in env.yml script.

To create the required environment executes 'conda env create -f env.yml'

