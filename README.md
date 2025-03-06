# Introduction

This repository is used by me to store the code related to multiple wave suppression and adaptive subtraction for my undergraduate thesis. The information is sourced from the internet. If you have any questions, please feel free to contact me

## Installation
(If you install a new version, there may be some parameter errors)
```python
pip install pylops
```

## Seismic Forward

[Reference](https://github.com/DIG-Kaust/Adaptive-subtraction/blob/master/notebooks/Data_Modeling.ipynb)

### 1. Setting free or damp boundary Obtain records without multiple waves

![free damp boundary](https://github.com/Lipeng-Lai/Mutiples_Suppression/blob/main/images/free_damp_boundary.png)

### 2. Data Modeling

## Radon Transform

### 1. Linear, Hyperbolic, Parabolic Radon Transform

[Pylops implementation](https://pylops.readthedocs.io/en/stable/api/generated/pylops.signalprocessing.FourierRadon2D.html)

![radon_forward](https://github.com/Lipeng-Lai/Mutiples_Suppression/blob/main/images/radon_forward.png)

![radon_inverse](https://github.com/Lipeng-Lai/Mutiples_Suppression/blob/main/images/radon_inverse.png)

### 2. $\lambda-f$ Radon Transform（Spare and high-resolution）

[Li, ZN., Li, ZC., Wang, P. et al. Multiple attenuation using λ-f domain high-resolution Radon transform. Appl. Geophysisc. 10, 433–441 (2013).](https://doi.org/10.1007/s11770-013-0405-1)

![radon](https://github.com/Lipeng-Lai/Mutiples_Suppression/blob/main/images/radon.png)


## Surface-related multiple elimination

[Verschuur D J, Berkhout A J, Wapenaar C P A. Adaptive surface-related multiple eliminationJ. Geophysics, 1992, 57: 1166-1177.](https://library.seg.org/doi/abs/10.1190/1.1443330)

### 1. Pre-Processing and no Pre-Processing

![SRME](https://github.com/Lipeng-Lai/Mutiples_Suppression/blob/main/images/SRME.png)

### 2. Closed-loop SRME

### Adaptive subtraction

[Adaptive subtraction](https://github.com/DIG-Kaust/Adaptive-subtraction/blob/master/notebooks/Adaptive_Subtraction.ipynb)
#### 1.  LSQR

#### 2.  ADMM

### Focal Transform

[Berkhout A J, Verschuur D J. Focal transformation, an imaging concept for signal restoration and noise removal](https://library.seg.org/doi/abs/10.1190/1.2356996)

##### 1. WCC(weighted correlation, Focus in the middle) 

![focal_zero](https://github.com/Lipeng-Lai/Mutiples_Suppression/blob/main/images/focal_zero.png)

##### 2. Focal transform(Consistent with the original data format)

![focal](https://github.com/Lipeng-Lai/Mutiples_Suppression/blob/main/images/focal.png)

## EPSI(Estimation of primaries by sparse inversion)

### 1. Original

### 2. Robust 
[Robust estimation of primaries by sparse inversion via one-norm minimization](https://library.seg.org/doi/10.1190/geo2012-0097.1)
[SLIM group](https://slim.gatech.edu/SoftwareDemos/applications/WavefieldSeparation/RobustEPSI/example.html)
[github](https://github.com/SINBADconsortium/SLIM-release-apps)

The installation is a bit complicated. Here, I will write about the method of using Matlab on Windows and calling it with WSL

```
git clone git@github.com:SINBADconsortium/SLIM-release-apps.git
cd SLIM-release-apps
cp environment.sh.template environment.sh
. environment.sh
```

```
git clone git@github.com:SINBADconsortium/SLIM-release-comp.git
cd SLIM-release-comp
cp environment.sh.template environment.sh
. environment.sh
```

```
sudo ln -s "/mnt/d/Program Files/MATLAB/R2022a/bin/matlab.exe" /usr/local/bin/matlab
sudo ln -s "/mnt/d/Program Files/MATLAB/R2022a/bin/win64/mex.exe" /usr/local/bin/mex

# test, if print, link sucess!
which matlab
which mex

# test SLIM
# in SLIM-release-apps
test_env4slim.sh

# in SLIM-release-comp
test_env4slim.sh

# in SLIM-release-apps
install_MEX
install_test4matlab

# in SLIM-release-comp(check anything if you have,such as gcc, then skip)
install_ALL

# start running code, you can replace your data, (and skip scons in data/) to fetch data
cd SLIM-release-apps/applications/WavefieldSeparation/RobustEPSI/Synthetic
scons
# generate 3 matlab files, and running  Synthetic.m, if path error, you can add path in Sythetic.m
addpath(genpath('../../matfcts'));
addpath(genpath('\\wsl.localhost\ubuntu-20.04\home\wwd\SLIM-release-apps\matlab'))
addpath(genpath('\\wsl.localhost\ubuntu-20.04\home\wwd\SLIM-release-apps\tools\solvers\SPGL1-SLIM'));
addpath(genpath('\\wsl.localhost\ubuntu-20.04\home\wwd\SLIM-release-apps\tools\algorithms\REPSI'));
addpath(genpath('\\wsl.localhost\ubuntu-20.04\home\wwd\SLIM-release-apps\tools\utilities\SPOT-SLIM'));

EPSI_SLIM_main('input_file','../data/your_data.mat','maxTotalIter',100,'padtime',11,'topmuteT',1,'q_estLength_posT',0.1,'q_estLength_negT',0.1,'useOblique',1,'relError',0.05,'window_startT',0.15,'window_endT',1,'savepreviewmat',1,'show_preview',1,'verbosity',0,'useSparsity',0,'parallel',0,'output_primary_file','Synthetic_primary.mat','output_primaryIR_file','Synthetic_primaryIR.mat','output_wavelet_file','Synthetic_wavelet.mat','preview_file','Preview/Synthetic_preview.mat','sol_file','Results/Synthetic_result.mat');

```

## EMD


## Seislet
