## Introduction

This repository is used by me to store the code related to multiple wave suppression and adaptive subtraction for my undergraduate thesis. The information is sourced from the internet. If you have any questions, please feel free to contact me

### Seismic Forward

### 1. Setting free or damp boundary Obtain records without multiple waves

![free damp boundary](https://github.com/Lipeng-Lai/Mutiples_Suppression/blob/main/images/free_damp_boundary.png)

### 2. Data Modeling

### Radon Transform

#### 1. Linear, Hyperbolic, Parabolic Radon Transform

[Pylops implementation](https://pylops.readthedocs.io/en/stable/api/generated/pylops.signalprocessing.FourierRadon2D.html)

![radon_forward](https://github.com/Lipeng-Lai/Mutiples_Suppression/blob/main/images/radon_forward.png)

![radon_inverse](https://github.com/Lipeng-Lai/Mutiples_Suppression/blob/main/images/radon_inverse.png)

#### 2. $\lambda-f$ Radon Transform（Spare and high-resolution）

[Li, ZN., Li, ZC., Wang, P. et al. Multiple attenuation using λ-f domain high-resolution Radon transform. Appl. Geophysisc. 10, 433–441 (2013).](https://doi.org/10.1007/s11770-013-0405-1)

![radon](https://github.com/Lipeng-Lai/Mutiples_Suppression/blob/main/images/radon.png)


### Surface-related multiple elimination

[Verschuur D J, Berkhout A J, Wapenaar C P A. Adaptive surface-related multiple eliminationJ. Geophysics, 1992, 57: 1166-1177.](https://library.seg.org/doi/abs/10.1190/1.1443330)

#### 1. Pre-Processing and no Pre-Processing

![SRME](https://github.com/Lipeng-Lai/Mutiples_Suppression/blob/main/images/SRME.png)

#### 2. Closed-loop SRME

#### Adaptive subtraction

[Adaptive subtraction](https://github.com/DIG-Kaust/Adaptive-subtraction/blob/master/notebooks/Adaptive_Subtraction.ipynb)
##### 1.  LSQR

##### 2.  ADMM

#### Focal Transform

[Berkhout A J, Verschuur D J. Focal transformation, an imaging concept for signal restoration and noise removal](https://library.seg.org/doi/abs/10.1190/1.2356996)

##### 1. WCC(weighted correlation, Focus in the middle) 

![focal_zero](https://github.com/Lipeng-Lai/Mutiples_Suppression/blob/main/images/focal_zero.png)

##### 2. Focal transform(Consistent with the original data format)

![focal](https://github.com/Lipeng-Lai/Mutiples_Suppression/blob/main/images/focal.png)

### EPSI(Estimation of primaries by sparse inversion)

#### 1. Original

#### 2. Robust


### EMD


### Seislet
