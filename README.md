# conflictmodel
Simulation for a socio-dynamical model of conflicts based on Alxelrod's idea of "pay or else"

## Abstract
In this work we implement the mathematical formulation ideas of a socio-dynamic model based on a tribute model, to simulate the interaction of elements and study the collective behaviors and the emergence of new levels of organization. The system is dynamic-coevolutionary in which states depend on local states and vice versa with discrete space and time.

## Model
We consider a periodic 2-dimensional network with each actor in a von Neumann neighborhood. The initial wealth distribution as well as the selection of the active agent are randomly selected. The target selection requires contiguity, in the network topology, with the active agent (alliances). The dynamics consist of a demand for resources with the threat of a conflict in case the defendant does not acquiesce; a case where both lose proportionally. Ideally, the perfect target should be weak enough to cause less damage if it chooses to fight and strong enough to afford to pay as much as it can afford. The target will only pay if fighting costs less than paying. 
The wealth for each coalision $W_{\alpha}$ and $W_{\tau}$ are:
        $$W_{\alpha} = \sum_{i=1}^{N_{\alpha}} C_{A i}W_{i}$$
        $$W_{\tau} = \sum_{i=1}^{N_{\tau}} C_{T i}W_{i}$$
        
in which, $N_{\alpha}$ and $N_{\tau}$ belong to coalision A and T respectively. 

Vulnerability is defined as: 

$$V_{A,j}= \frac{W_{\alpha}-W_{\tau}}{W_{\alpha}}$$

The loss of each element in a coalition is:
$$L_{i\alpha}=\kappa W_{\tau}\frac{W_{i}}{W_{\alpha}}$$,


$$L_{i\tau}=\kappa W_{\alpha}\frac{W_{i}}{W_{\tau}}$$



Finally, each decision develops degrees of commitment and influences future decisions. We establish a steady increase in commitment under "submission", "protection" and "friendship" and a decrease under "hostility", in and out of alliances. We track this history with a matrix. The commitment Matrix is and identity Matrix NxN where $M_{i,j}=C_{i,j}$ is the mutual commitment between i and j


## Instructions

- 1.) conflict_model.py  This is a python script which simulates the conflict model described above by selecting a random active agent among the N independent agents, the selection of the ideal target requiring continguity is performed by applying al the conditions mentions previously. The process will iterate until the number of years selected. Additionally, the initial conditions are set out here, but are chosen for convenience. 

#### How to use conflict_model.py
-import conflict_model.py 
- to start the simulation run conflict_model.populationRun2D(N, Years) where N is the number of independent actors. 
- .populationRun2D() will return the resources for each actor, thee attackers, the defenders, the activity of the commitment matrix, an status list with 1 for conflcits and 0 for peace, the total loss, the modified commitment matrix, active agent and target. 

These data will be useful in the search of universal patterns for statistics. 




## Contact
For questions or feedback, please contact me at [brandon.minta@yachaytech.edu.ec]