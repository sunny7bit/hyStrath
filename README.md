<h1 align="center">hyStrath</h1>  

#### Hypersonic / Rarefied gas dynamics code developments under license GPL-3.0 
#### The only platform to conjointly host open-source CFD and DSMC codes designed for atmospheric re-entry analysis

#### The *Fleming* release includes  
+ *hyFoam*: a CFD solver for supersonic combusting flows   
+ *hy2Foam*: a CFD solver for hypersonic reacting flows with MHD capabilities  
+ *dsmcFoam+*: the direct simulation Monte Carlo (DSMC) code with all the latest features  
+ *pdFoam*: a hybrid PIC-DSMC solver   

#### Changes
+ New solver: hisahy2Foam which can effectively solve the problem of thermal nuclear temperature evolution for local high-temperature and high-pressure initialisation.
+ New fluxScheme: Kurganov, Tadmor, ROECorrect, AUSMPLUS, AUSMPLUSUP (ROE, KFVS, LxF, HLLC).
+ Bug repair: minMaxTemperatureFieldOutOfRange not reporting accurately when running in parallel.


#### Install  
```bash
./install.sh 6 2>/dev/null
```

#### Please visit the [_hyStrath_ website](https://hystrath.github.io/)
