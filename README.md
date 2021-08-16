# SIDHdelegation


This MAGMA code implements representative benchmarks to assess the cost of the delegation algorithms presented in [1]. It needs Microsoft's vOW4SIKE package [2] as a dependency. 

To install, simply clone this repository and add the vOW4SIKE package as a folder into the SIDHdelegation folder:
```
git clone https://github.com/gemeis/SIDHdelegation
cd SIDHdelegation
git clone https://github.com/microsoft/vOW4SIKE
```

To run the benchmarks, execute
``` 
magma run.m
```
Note that this does not represent a fully working implementation, but rather an assessment of the computational cost of the delegator. In order to verify correctness of the algorithms presented in [1], we further implemented a proof of concept. This uses the built-in elliptic curves in Magma and is therefore not optimized. To verify correctness, run
```
magma proof-of-concept.m
```





### References
[1] Delegating Supersingular Isogenies over $\mathbb{F}_{p^2}$ with Cryptographic Applications, 2021 (submission to ICISC).

[2] https://github.com/microsoft/vOW4SIKE
