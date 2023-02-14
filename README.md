# X-Wing (Cross-population Weighting) 
X-Wing (**Cross**-population **W**eight**ing**) is a cross-population polygenic prediction framework only using GWAS summary statistics. It contains three methods:

* `LOGODetect`(**LO**cal **G**enetic c**O**rrelation **DETECT**or) (Also see [here](https://github.com/ghm17/LOGODetect))
* `PANTHER`(**P**olygenic prediction using Bayesian **AN**no**T**ation-dependent **H**ors**E**shoe **R**egression)
* `LEOPARD`(**L**in**E**ar c**O**mbination of **P**RS using gw**A**s summa**R**y statistics **D**ata)

The X-Wing software can perform four tasks:
1. **Improve genetic risk prediction in ancestrally diverse populations using `X-Wing`**.
2. Identify small genome segments harboring cross-population local genetic correlation signals using `LOGODetect`.
3. Incorporate annotation into single and cross-population polygenic prediction using `PANTHER`.
4. Estimate the linear combination weights for multiple PRS by GWAS summary statistics alone using `LEOPARD`.

The layout of the scripts is
![](https://github.com/qlu-lab/X-Wing/blob/main/Fig1_Github.png)

## Manual

Please see the [wiki](https://github.com/qlu-lab/X-Wing/wiki) for the detailed manual of `X-Wing`.


## Version History
* May 26, 2022: Initial release.



## Acknowledgement

Part of the code is adapted from [PRS-CSx](https://github.com/getian107/PRScsx) and [XPASS](https://github.com/YangLabHKUST/XPASS). We thank Dr. Tian Ge and Dr. Mingxuan Cai for sharing their code.

## Data

The X-Wing posterior SNP effect size estimates in our manusciript are available at [here](https://uwmadison.box.com/s/zzf6s3z9v9vmq1avd9dbbb4nq7r9nenj).

## Citation

If you use X-Wing, please cite 

Miao, J., Guo, H., Song, G. et al. Quantifying portable genetic effects and improving cross-ancestry genetic prediction with GWAS summary statistics. Nat Commun 14, 832 (2023). https://doi.org/10.1038/s41467-023-36544-7


## Contact

For questions and comments, please open a Github issue (preferred) or contact Jiacheng Miao at jmiao24@wisc.edu or Hanmin Guo at ghm17@mails.tsinghua.edu.cn. 
