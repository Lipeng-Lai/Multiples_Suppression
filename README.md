# Introduction

This repository is used by me to store the code related to multiple wave suppression and adaptive subtraction for my undergraduate thesis. The information is sourced from the internet. If you have any questions, please feel free to contact me

## Installation

```
# If you install a new version, there may be some parameter errors
pip install pylops 
```

```
# install Madagascar(Because my Matlab in Windows, and work in Linux, there is some bug,
# so I don't configure API=matlab in Madagascar, and I
# convert rsf to segy, and read segy data in Matlab, processing ..., save segy data in Matalb,
# finally, convert segy data to rsf )
# https://ahay.org/wiki/Installation
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
git clone --depth=1 https://github.com/SINBADconsortium/SLIM-release-apps.git
cd SLIM-release-apps
cp environment.sh.template environment.sh
vim enviroment.sh

# edit enviroment after, if you just need spgl1, skip SLIM-release-comp
export SLIM_COMP='/home/wwd/SLIM-release-comp'
export SLIM_APPS='/home/wwd/SLIM-release-apps'

# esc vim, and edit ~/.bashrc
source /home/wwd/SLIM-release-apps/environment.sh
. /home/wwd/SLIM-release-apps/environment.sh
source ~/.bashrc
```

```
# if your matlab in windows, and work in Linux
sudo ln -s "/mnt/d/Program Files/MATLAB/R2022a/bin/matlab.exe" /usr/local/bin/matlab
sudo ln -s "/mnt/d/Program Files/MATLAB/R2022a/bin/win64/mex.exe" /usr/local/bin/mex

# test, if print, link sucess!
which matlab
which mex
```

```
# test SLIM
# in SLIM-release-apps
test_env4slim.sh

# in SLIM-release-apps, if it excute error, don't mind, just skip(if you fix the error,
# such as 'PI not define, and define PI 3.14....', and you run following code, will occur
# x = (x, nt, nr, ns) error, because x = (0, 0),
# I guess SLIM-spgl1 is  broken) 
install_MEX

# start running code, you can replace your data(.mat, .bin, .su), (and skip scons in data/) to fetch data
cd SLIM-release-apps/applications/WavefieldSeparation/RobustEPSI/Synthetic
scons
# generate 3 matlab files, and running  Synthetic.m, if path error, you can add path in Sythetic.m
addpath(genpath('../../matfcts'));
addpath(genpath('\\wsl.localhost\ubuntu-20.04\home\wwd\SLIM-release-apps\matlab'))
addpath(genpath('\\wsl.localhost\ubuntu-20.04\home\wwd\SLIM-release-apps\tools\solvers\SPGL1-SLIM'));
addpath(genpath('\\wsl.localhost\ubuntu-20.04\home\wwd\SLIM-release-apps\tools\algorithms\REPSI'));
addpath(genpath('\\wsl.localhost\ubuntu-20.04\home\wwd\SLIM-release-apps\tools\utilities\SPOT-SLIM'));

EPSI_SLIM_main('input_file','../data/your_data.mat','maxTotalIter',100,'padtime',11,'topmuteT',
1,'q_estLength_posT',0.1,'q_estLength_negT',0.1,'useOblique',1,'relError',0.05,'window_startT',
0.15,'window_endT',1,'savepreviewmat',1,'show_preview',1,'verbosity',0,'useSparsity',0,'parallel',
0,'output_primary_file','Synthetic_primary.mat','output_primaryIR_file','Synthetic_primaryIR.mat',
'output_wavelet_file','Synthetic_wavelet.mat','preview_file','Preview/Synthetic_preview.mat',
'sol_file','Results/Synthetic_result.mat');

```

## PEF(prediction-error filter and prediction decon)

[Multiple suppression using prediction-error filter](https://ahay.org/RSF/book/sep/pefmult/paper_html/)

![old_pef](https://github.com/Lipeng-Lai/Mutiples_Suppression/blob/main/images/old_pef.png)

[Multichannel adaptive deconvolution based on streaming prediction-error filter](https://ahay.org/RSF/book/jlu/spefdecon/paper_html/paper.html)

![Sigbee_zero_offset](https://github.com/Lipeng-Lai/Mutiples_Suppression/blob/main/images/sigbee_zero_offset.png)

![tpef](https://github.com/Lipeng-Lai/Mutiples_Suppression/blob/main/images/tpef.png)

![vlag-spef](https://github.com/Lipeng-Lai/Mutiples_Suppression/blob/main/images/vlag-spef.png)

## EMD(Empirical mode decomposition)
[Multiple reflections noise attenuation using adaptive randomized-order empirical mode decomposition](https://ahay.org/RSF/book/tccs/demulemd/paper_html/paper.html)
[Matlab code](https://github.com/chenyk1990/reproducible_research/tree/795d2e2217f21261a0372abc5f503a1d6b601fc6/demulemd)

![emd](https://github.com/Lipeng-Lai/Mutiples_Suppression/blob/main/images/emd.png)


## LPF(Local prediction filter)

[Adaptive multiple subtraction using regularized nonstationary regression](https://ahay.org/RSF/book/tccs/lpf/paper_html/paper.html)

![lpf](https://github.com/Lipeng-Lai/Mutiples_Suppression/blob/main/images/lpf.png)

## Seislet

[Seislet-based morphological component analysis using scale-dependent exponential shrinkage](https://ahay.org/RSF/book/xjtu/mcaseislet/paper_html/paper.html)

![seislet](https://github.com/Lipeng-Lai/Mutiples_Suppression/blob/main/images/seislet.png)

## Curvelet


## DeepLearning


