# XAEMv2

This is the version 2 of [XAEM](http://fafner.meb.ki.se/biostatwiki/xaem/) for heterogeneous datasets. You can find the detailed instructions from **XAEM** website.

## What's news?

* **New method for quantifing heterogeneous datasets** which have many kinds of tissues. You run XAEM with a new parameter `tau` - percent of CRPS (e.g. 0.05(5%)) as you want to treat for heterogeneous samples:

```
**Section 3.2.3 in the detailed instruction**

Rscript $PWD/XAEMv2/R/AEMCpp.R workdir=$PWD/XAEM_project core=$CPUNUM design.matrix=$PWD/X_matrix.RData isoform.out=XAEM_isoform_expression.RData paralog.out=XAEM_paralog_expression.RData tau=0.05
```

* **New EM estimation in C++**. It reduces a lot of time to estimate isoform expression:

```
Unit: microseconds
     expr       min        lq      mean    median       uq       max
    getSE  3519.769 3845.9475 4752.0856 4005.6250 4514.541 24643.438
 getSECpp    95.182  226.5485  243.7424  243.8395  263.657   306.596

   expr         min        lq      mean    median        uq       max
    AEM   378.95231 393.19988 423.44539 402.21083 431.00111 730.40832
 AEMCpp    16.79823  17.06698  18.20859  17.29657  17.87722  34.71934

  expr       min       lq       mean    median       uq       max
    EM  2795.862 2940.142 3160.79362 3090.8570 3211.525 10044.717
 EMCpp    57.298   58.733   68.60492   65.4165   73.204   160.596
```

```
| Pipeline | 20 samples | 100 samples |
|----------|:----------:|:-----------:|
|   AEM    |  34m0.656s |  119m7.607s |
|  AEMCpp  |  1m56.748s |  6m52.565s  |
```
