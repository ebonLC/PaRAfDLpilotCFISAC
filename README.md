# About the work
This is the MATLAB codes for the paper: Liao, J., Zhu, H., Chen, F., Wang, J., & Chen, C. (2025). Power and resource allocation for downlink pilot based cell-free integrated sensing and communication systems. IEEE Journal on Selected Areas in Communications., to appear. 

Please note that the MATLAB2022b is used for this simulation. To run those codes, please download and install [CVX]([https://example.com](https://cvxr.com/cvx/)).

The program "All_cvxpgd_SNR" is the main simulation. The comparison among DL/UL ratio can be changed by uncomment/comment-ing the corresponding lines of code. The programs "..._good" and "..._others" are used for testing the impacts of orthogonality. The "UE positions and results.mat" is generated from "..._good" and the same UE positions will be reused in "..._others". The additional PGD(Projected Gradient Descent) approach which is not discussed in the work has also been contained. 

Abstract of the paper: 
In this paper, an integrated sensing and communication (ISAC) system was proposed under cell-free architecture with different power allocation schemes for both communication and sensing functions. In addition to using uplink (UL) pilot transmission for channel estimation, a group of pre-designed orthogonal downlink (DL) pilot sequences are inserted and transmitted from each cooperative access point (AP) towards the user equipments (UEs). The inserted DL pilots in cell-free ISAC are utilized for both channel estimation in communication and cooperative localization in radar sensing. Based on the ISAC strategy, a power allocation scheme is proposed with the objective of minimizing the maximum channel estimation error for each user, under the power constraints of APs and minimum mean square error constraints of the localization. The full-power scheme which APs consume all the available power during the DL phase is also studied and analyzed. The performance of the proposed ISAC strategy is evaluated under a realistic highway scenario, and compared to the existing non-DL strategies. Additionally, we investigate the impact of the pilot ratio under a constrained signal length, as well as the impact of the pilot proportion that the pilot length to the entire signal frame length. Simulation results show that the performances of both channel estimation and localization can be significantly improved at a relatively low cost of extra power needed for APs. Moreover, the analysis provides insights into optimizing the DL/UL pilot ratio and pilot proportion which can significantly enhance both channel estimation and localization accuracy, while ensuring sufficient resources remain available for data payload.

Please feel free to contact Jiangwei Liao (email1: jl804@kent.ac.uk;email2:myidisjil@gmail.com) if you have any questions. 

If you use this code in your research, please cite:
@article{liao2025power,
  title={Power and resource allocation for downlink pilot based cell-free integrated sensing and communication systems},
  author={Liao, Jiangwei and Zhu, Huiling and Chen, Fangyuan and Wang, Jiangzhou and Chen, Changrun},
  journal={IEEE Journal on Selected Areas in Communications},
  year={2025},
  publisher={IEEE}
}
