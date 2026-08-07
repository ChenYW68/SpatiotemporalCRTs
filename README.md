
# Identifying Optimal Schistosomiasis Treatment Sequences for Mass Drug Administration Using Direct and Indirect Effects in Spatiotemporal Trials

This Github page provides code and data for reproducing the results in the manuscript:``Identifying Optimal Schistosomiasis Treatment Sequences for Mass Drug Administration Using Direct and Indirect Effects in Spatiotemporal Trials'' by Y. Chen, X. Wen, F. Luo, Y. Yang, and Y. Shen. 

## Datasets of schistosomiasis from the SCORE project
Schistosomiasis, a neglected tropical parasitic disease, is prevalent in Africa, South America, and Asia − especially in rural areas with poor socioeconomic conditions (Hong et al., 2022) − with cases reported in more than 70 countries (World Health Organization, 2022). In 2021 alone, at least 251 million people required preventive treatment (World Health Organization, 2023). Global efforts to control morbidity have primarily focused on preventive chemotherapy using praziquantel (PZQ), typically delivered through Mass Drug Administration (MDA) via School-Based Treatment (SBT) or Community-Wide Treatment (CWT). To systematically evaluate and compare the two treatment strategies, the Schistosomiasis Consortium for Operational Research and Evaluation (SCORE) conducted a five-year cluster randomized trial (CRT) across multiple countries. This data can be saved under the path "./data". Data description is as follows:
- We analyze infection prevalence data for Schistosoma mansoni from 149 villages in Kenya and 149 villages in Tanzania, covering a wide age range from 5 to 77 years.
- These infection data, collected on a yearly scale, were obtained through the SCORE project over a 5-year period
- This project involves a CRT with six intervention arms
- Each arm received one of three treatment options (SBT, CWT, or no treatment) annually from the first to the fourth year

## Our methodology
Motivated by the need to identify optimal treatment regimens in the Schistosomiasis Consortium for Operational Research and Evaluation CRT, this study addresses two major methodological challenges in estimating treatment effects: (1) substantial bias arising  from inadequate accounting for Indirect Effects (IEs) and time-varying Direct Effects (DEs) and (2) large variance from insufficient consideration of intrinsic dependencies. We identify optimal treatment regimens by evaluating differences through two components: DEs from the most recent treatment and IEs from historical treatment trajectories. To efficiently estimate DEs and IEs, we develop a Joint Spatiotemporal Varying Coefficient (JSTVC) model. JSTVC accounts for spatiotemporal dependencies and regional heterogeneities, while also capturing spatial anisotropic patterns associated with schistosomiasis transmission. 

## Scalable algorithms for large-scale randomized experiments with complex spatiotemporal dependence structure
To support scalable inference under complex dependent structures, we develop a scable Variational Bayes algorithm with an ensemble-based correction to improve uncertainty quantification. The proposed methodology provides a broadly applicable framework for modeling complex dependencies in randomized experiments, especially in those involving multiple sequential interventions.


# Code and Data Documentation

To reproduce all results presented in the paper, run `./Main.R`. Depending on the computing environment, the script may take a considerable amount of time to complete.

The following sections summarize the data and code used in the repository.

## A. Data

The `./data` directory contains the following files.

| File | Main objects | Notes |
| --- | --- | --- |
| `Kenya_Score_Data_r.RData` | `Kenya_Score_Data`, `Site`, `Kenya.Dist.c`, `G.mat` | Kenya SCORE trial data, site-level information, and distance/adjacency matrices |
| `Tanzania_Score_Data_r.RData` | `Tanzania_Score_Data`, `Site`, `Tanzania.Dist.c`, `G.mat` | Tanzania SCORE trial data, site-level information, and distance/adjacency matrices |
| `Google_Kenya_Tanzania_Map.RData` | `ken.map`, `tan.map` | Google map objects used for geographic plots |
| `sim.Cov.Data.RData` | `sim.Cov.Data` | Simulation Scenario 1 and 3 covariance object; stored as a list of 5 regional sub-objects |
| `smoothed_surface.RData` | `simData.DataBase` | Simulation Scenario 2 smoothed surface object; stored as a list of 5 regional sub-objects |


### A.1 Using `Tanzania_Score_Data_r.RData` as an illustration

Variables are as follows:

`Tanzania_Score_Data` (745 x 42)
<table>
  <colgroup>
    <col style="width: 30%;">
    <col style="width: 30%;">
    <col style="width: 40%;">
  </colgroup>
  <thead>
    <tr>
      <th style="font-size: 20%;">No.</th>
      <th>Variables</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr><td>1</td><td><code>Village_ID</code></td><td>Village identifier</td></tr>
    <tr><td>2</td><td><code>Year</code></td><td>Observation year</td></tr>
    <tr><td>3</td><td><code>Prevalence</code></td><td>Observed infection prevalence</td></tr>
    <tr><td>4</td><td><code>Study_Type</code></td><td>Study or trial type label</td></tr>
    <tr><td>5</td><td><code>Study_Arm</code></td><td>Assigned treatment arm</td></tr>
    <tr><td>6</td><td><code>Latitude</code></td><td>Village latitude</td></tr>
    <tr><td>7</td><td><code>Longitude</code></td><td>Village longitude</td></tr>
    <tr><td>8</td><td><code>flag</code></td><td>Region or geographic grouping label</td></tr>
    <tr><td>9</td><td><code>CWT_1</code></td><td>Indicator for CWT at time lag 1</td></tr>
    <tr><td>10</td><td><code>CWT_2</code></td><td>Indicator for CWT at time lag 2</td></tr>
    <tr><td>11</td><td><code>CWT_3</code></td><td>Indicator for CWT at time lag 3</td></tr>
    <tr><td>12</td><td><code>CWT_4</code></td><td>Indicator for CWT at time lag 4</td></tr>
    <tr><td>13</td><td><code>SBT_1</code></td><td>Indicator for SBT at time lag 1</td></tr>
    <tr><td>14</td><td><code>SBT_2</code></td><td>Indicator for SBT at time lag 2</td></tr>
    <tr><td>15</td><td><code>SBT_3</code></td><td>Indicator for SBT at time lag 3</td></tr>
    <tr><td>16</td><td><code>SBT_4</code></td><td>Indicator for SBT at time lag 4</td></tr>
    <tr><td>17</td><td><code>CWT</code></td><td>Current CWT indicator</td></tr>
    <tr><td>18</td><td><code>SBT</code></td><td>Current SBT indicator</td></tr>
    <tr><td>19</td><td><code>IEt.CWT</code></td><td>Pow decay temporal variable for CWT</td></tr>
    <tr><td>20</td><td><code>IEt.SBT</code></td><td>Pow decay temporal variable for SBT</td></tr>
    <tr><td>21</td><td><code>IEt.CWT.exp</code></td><td>Exponential-decay temporal variable for CWT</td></tr>
    <tr><td>22</td><td><code>IEt.SBT.exp</code></td><td>Exponential-decay temporal variable for SBT</td></tr>
    <tr><td>23</td><td><code>sCWT_1</code></td><td>Smoothed or scaled CWT summary at lag 1</td></tr>
    <tr><td>24</td><td><code>sCWT_2</code></td><td>Smoothed or scaled CWT summary at lag 2</td></tr>
    <tr><td>25</td><td><code>sCWT_3</code></td><td>Smoothed or scaled CWT summary at lag 3</td></tr>
    <tr><td>26</td><td><code>sSBT_1</code></td><td>Smoothed or scaled SBT summary at lag 1</td></tr>
    <tr><td>27</td><td><code>sSBT_2</code></td><td>Smoothed or scaled SBT summary at lag 2</td></tr>
    <tr><td>28</td><td><code>sSBT_3</code></td><td>Smoothed or scaled SBT summary at lag 3</td></tr>
    <tr><td>29</td><td><code>IEs.No_Treatment</code></td><td>Spatial indirect-effect summary under no treatment</td></tr>
    <tr><td>30</td><td><code>IEs.CWT.sp.Neigh.500.10</code></td><td>Spatial IE for CWT using 10-km as a range and a 500-km cut-off</td></tr>
    <tr><td>31</td><td><code>IEs.SBT.sp.Neigh.500.10</code></td><td>Spatial IE for SBT using 10-km as a range and a 500-km cut-off</td></tr>
    <tr><td>32</td><td><code>IEs.No_Treatment.sp.Neigh.500.10</code></td><td>Spatial IE for no treatment using 10-km as a range and a 500-km cut-off</td></tr>
    <tr><td>33</td><td><code>IEs.CWT.sp.Neigh.500.30</code></td><td>Spatial IE for CWT using 30-km as a range and a 500-km cut-off</td></tr>
    <tr><td>34</td><td><code>IEs.SBT.sp.Neigh.500.30</code></td><td>Spatial IE for SBT using 30-km as a range and a 500-km cut-off</td></tr>
    <tr><td>35</td><td><code>IEs.No_Treatment.sp.Neigh.500.30</code></td><td>Spatial IE for no treatment using 30-km as a range and a 500-km cut-off</td></tr>
    <tr><td>36</td><td><code>IEs.CWT.sp.Neigh.500.50</code></td><td>Spatial IE for CWT using 50-km as a range and a 500-km cut-off</td></tr>
    <tr><td>37</td><td><code>IEs.SBT.sp.Neigh.500.50</code></td><td>Spatial IE for SBT using 50-km as a range and a 500-km cut-off</td></tr>
    <tr><td>38</td><td><code>IEs.No_Treatment.sp.Neigh.500.50</code></td><td>Spatial IE for no treatment using 50-km as a range and a 500-km cut-off</td></tr>
    <tr><td>39</td><td><code>distances</code></td><td>Distance to Lake Victoria</td></tr>
    <tr><td>40</td><td><code>log.log.mean.Intensity</code></td><td>Mean of Log-log transformed infection intensity</td></tr>
    <tr><td>41</td><td><code>log.log.var.Intensity</code></td><td>Variance of log-log transformed infection intensity</td></tr>
    <tr><td>42</td><td><code>Intercept</code></td><td>Intercept term used in model construction</td></tr>
  </tbody>
</table>

`Site` (149 x 8)

<table>
  <colgroup>
    <col style="width: 7%;">
    <col style="width: 23%;">
    <col style="width: 70%;">
  </colgroup>
  <thead>
    <tr>
      <th style="font-size: 90%;">No.</th>
      <th>Variables</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr><td>1</td><td><code>Village_ID</code></td><td>Village identifier</td></tr>
    <tr><td>2</td><td><code>Study_Arm</code></td><td>Assigned treatment arm</td></tr>
    <tr><td>3</td><td><code>LAT</code></td><td>Latitude</td></tr>
    <tr><td>4</td><td><code>LON</code></td><td>Longitude</td></tr>
    <tr><td>5</td><td><code>flag</code></td><td>Region or geographic grouping label</td></tr>
    <tr><td>6</td><td><code>LON_X</code></td><td>Projected longitude coordinate</td></tr>
    <tr><td>7</td><td><code>LAT_Y</code></td><td>Projected latitude coordinate</td></tr>
    <tr><td>8</td><td><code>distances</code></td><td>Distance to Lake Victoria</td></tr>
  </tbody>
</table>

## B. Code for JSTVC

The source code for the JSTVC method is located in the `./JSTVC/` directory.

## C. Simulation

All simulation code is located in the `./R/Simulation/` directory.

1. The `./R/Simulation/RandomField/` directory contains code for the first scenario using the Gneiting spatiotemporal covariance structure.
2. The `./R/Simulation/SmoothedFun/` directory contains code for the second scenario using low-rank basis expansion.
3. The `./R/Simulation/RandomField_xMisspecified_Decay/` directory contains code for the third scenario involving misspecification of the decay function.

## D. Real Data Analysis

All real data analysis code is located in the `./R/Case/` directory.

## E. Tables and Figures

1. Code for reproducing all tables is located in the `./table/` directory.
2. Code for reproducing figures from the simulation studies is located in the `./plot/simulation/` directory.
3. Code for reproducing figures from the exploratory data analysis (EDA) of the SCORE data is located in the `./plot/case/EDA/` directory.
4. Code for reproducing figures from the cross-validation analysis of the SCORE data is located in the `./plot/case/CV/` directory.
5. Code for reproducing figures from model fitting to the SCORE data is located in the `./plot/case/Fitting/` directory.
6. Code for reproducing figures related to treatment-sequence ranking for the SCORE data is located in the `./plot/case/Ranking/` directory.

## F. Results
1. All figures are saved to the `./figure/` directory
1. All tables are saved to the `./result/summary` directory


## Spatiotemporal patterns of schistosomiasis
Figure 1 illustrates the influence of treatment effects and spatiotemporal random effects on schistosomiasis, i.e.,
<figure id="Figure4">
  <p align="center">
  <img src="./figure/Fig7_Kenya_Wts.jpg" width="800px">
    </p>
  <figcaption>
  <strong>Figure 1:</strong> Recovered spatiotemporal patterns of the different components. Top panel: Observed prevalence. Middle panel: Prevalence excluding direct and indirect effects. Bottom panel: Recovered spatiotemporal random effects.
  </figcaption>
</figure>

## Ranking treatment sequences
<figure id="Figure5">
  <p align="center">
  <img src="./figure/Fig5_Ranks.jpg" width="800px">
    </p>
  <figcaption>
  <strong>Figure 2:</strong> Ranking treatment sequences across methods: (A) The proposed JSTVC; (B) JSTVC without the spatiotemporal random effect;  (C) JSTVC which mixed DEs and IEs and did not decompose the ATE; and (D) Ranking based on the average
reduction in outcomes from Year 1 to Year 5. 
  </figcaption>
</figure>

## Differences in ATEs between different treatment sequences
<figure id="Figure6">
  <p align="center">
  <img src="./figure/Fig6_Dist_ATE_all.jpg" width="800px">
    </p>
  <figcaption>
  <strong>Figure 3:</strong> Posterior distributions of differences in Average Treatment Effects (ATEs) between treatment sequences are computed using the proposed JSTVC, with 95% credible intervals (CIs) highlighted by shaded areas, where results from the VB implementation are compared with those obtained via MCMC.
  </figcaption>
</figure>


