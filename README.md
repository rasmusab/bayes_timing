Estimating the Distribution of Sensorimotor Synchronization Data: A Bayesian Hierarchical Modeling Approach 
============

This repository contains scripts that implements the models described in Bååth (in preparation). In order to make the scripts fully self contained this repository also includes a dataset from [Bååth and Madison (2012)]((www.sumsar.net/papers/ICMPC_2012_rasmus_baath.pdf)) of synchronization tapping at interstimulus intervals from 600 ms to 3000 ms. Both the scripts and the dataset are licenced under the MIT licence. If you use the dataset, please cite Bååth and Madison (2012).

**non_informative_subject.R** implements the non-informative non-hierarchical model.

**non_informative_hierarchical.R** implements the non-informative hierarchical model.

**informative_hierarchical.R** implements the informative hierarchical model.

**hierarchical_with_multi_normal_ISI_dependency.R** implements an extension of the non-informative hierarchical model where the correlation between participants' timing variability at different ISI levels is modeled by a multivariate normal distribution.

**hierarchical_with_linear_ISI_dependency.R** implements an extension of the non-informative hierarchical model where the group parameters are assumed to depend lineary on the ISI.

**maximum_likelihood_approach.R** implements a maximum likelihood approach as proposed by Ulrich and Miller (1994).

References
-------------------------
Bååth, R., & Madison, G. (2012). The subjective difficulty of tapping to a slow beat. *In Proceedings of the 12th International Conference on Music Perception and Cognition* (pp. 82–85). [pdf](www.sumsar.net/papers/ICMPC_2012_rasmus_baath.pdf)

Bååth, R. Estimating the Distribution of Sensorimotor Synchronization Data: A Bayesian Hierarchical Modeling Approach . (in preparation)

Ulrich, R., & Miller, J. (1994). Effects of truncation on reaction time analysis. *Journal of Experimental Psychology: General*, 123(1), 34–80. doi:10.1037/0096-3445.123.1.34