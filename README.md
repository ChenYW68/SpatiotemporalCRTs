
# Identifying Optimal Schistosomiasis Treatment Sequences for Mass Drug Administration Using Direct and Indirect Effects in Spatiotemporal Trials

This Github page provides code and data for reproducing the results in the manuscript:``Identifying Optimal Schistosomiasis Treatment Sequences for Mass Drug Administration Using Direct and Indirect Effects in Spatiotemporal Trials'' by Y. Chen, X. Wen, F. Luo, Y. Yang, and Y. Shen. 

## Datasets of schistosomiasis from the SCORE project
Schistosomiasis, a neglected tropical parasitic disease, is prevalent in Africa, South America, and Asia − especially in rural areas with poor socioeconomic conditions (Hong et al., 2022) − with cases reported in more than 70 countries (World Health Organization, 2022). In 2021 alone, at least 251 million people required preventive treatment (World Health Organization, 2023). Global efforts to control morbidity have primarily focused on preventive chemotherapy using praziquantel (PZQ), typically delivered through Mass Drug Administration (MDA) via School-Based Treatment (SBT) or Community-Wide Treatment (CWT). To systematically evaluate and compare the two treatment strategies, the Schistosomiasis Consortium for Operational Research and Evaluation (SCORE) conducted a five-year cluster randomized trial (CRT) across multiple countries. This data can be downloaded via [SCORE](https://clinepidb.org/ce/app/workspace/analyses/DS_d6a1141fbf/new). Data description is as follows:
- We analyze infection prevalence data for Schistosoma mansoni from 149 villages in Kenya and 149 villages in Tanzania, covering a wide age range from 5 to 77 years.
- These infection data, collected on a yearly scale, were obtained through the SCORE project over a 5-year period
- This project involves a CRT with six intervention arms
- Each arm received one of three treatment options (SBT, CWT, or no treatment) annually from the first to the fourth year

## Our methodology
Motivated by the need to identify optimal treatment regimens in the Schistosomiasis Consortium for Operational Research and Evaluation CRT, this study addresses two major methodological challenges in estimating treatment effects: (1) substantial bias arising  from inadequate accounting for Indirect Effects (IEs) and time-varying Direct Effects (DEs) and (2) large variance from insufficient consideration of intrinsic dependencies. We identify optimal treatment regimens by evaluating differences through two components: DEs from the most recent treatment and IEs from historical treatment trajectories. To efficiently estimate DEs and IEs, we develop a Joint Spatiotemporal Varying Coefficient (JSTVC) model. JSTVC accounts for spatiotemporal dependencies and regional heterogeneities, while also capturing spatial anisotropic patterns associated with schistosomiasis transmission. 

## Scalable algorithms for large-scale randomized experiments with complex spatiotemporal dependence structure
To support scalable inference under complex dependent structures, we develop a scable Variational Bayes algorithm with an ensemble-based correction to improve uncertainty quantification. The proposed methodology provides a broadly applicable framework for modeling complex dependencies in randomized experiments, especially in those involving multiple sequential interventions.

## Spatiotemporal patterns of schistosomiasis
Figure 1 illustrates the influence of treatment effects and spatiotemporal random effects on schistosomiasis, i.e.,
<figure id="Figure4">
  <p align="center">
  <img src="./figure/Fig7_Kenya_Wts.jpg" width="1000px">
    </p>
  <figcaption>
  <strong>Figure 1:</strong> Recovered spatiotemporal patterns of the different components. Top panel: Observed prevalence. Middle panel: Prevalence excluding direct and indirect effects. Bottom panel: Recovered spatiotemporal random effects.
  </figcaption>
</figure>

## Ranking treatment regimens
Figure 2 shows ranking results of treatment sequences across methods: (A) The proposed JSTVC; (B) JSTVC without the spatiotemporal random effect;  (C) JSTVC, which mixes DEs and IEs and does not decompose the ATE; and (D) Ranking based on the average
reduction in outcomes from Year 1 to Year 5. 
<figure id="Figure4">
  <p align="center">
  <img src="./figure/Fig4_Ranks.jpg" width="1000px">
    </p>
  <figcaption>
  <strong>Figure 2:</strong> Ranking treatment strategies for six arms in Tanzania. Left panel: Average Treatment Effect (ATE) calculated as the sum of Direct Effects (DEs) and Indirect Effects (IEs) from the proposed JSTVC. Right panel: Average relative reduction of observed outcomes from Year 1 to Year 5.
  </figcaption>
</figure>

## Differences in ATEs between different treatment regimens
<figure id="Figure4">
  <p align="center">
  <img src="./figure/Fig3_Dist_ATE_all.jpg" width="1000px">
    </p>
  <figcaption>
  <strong>Figure 3:</strong> Posterior distributions of differences in Average Treatment Effects (ATEs) between arms are computed using the proposed JSTVC, with 95% credible intervals (CIs) highlighted by shaded areas, where results from the VB implementation are compared with those obtained via MCMC.
  </figcaption>
</figure>


