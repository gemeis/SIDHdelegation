# SIDHdelegation

This MAGMA code implements benchmarks to assess the cost of the delegation algorithms presented in [1]. It needs Microsoft's vOW4SIKE package [2] as a dependency. 

To install, simply clone this repository and add the vOW4SIKE package as a folder into the SIDHdelegation folder:
```
git clone https://github.com/gemeis/SIDHdelegation
cd SIDHdelegation
git clone https://github.com/microsoft/vOW4SIKE
```

To run, execute
``` 
magma run.m
```






### References
[1] R. Pedersen, O. Uzunkol: Delegating Supersingular Isogenies over $\mathbb{F}_{p^2}$ with Cryptographic Applications, 2021.

[2] https://github.com/microsoft/vOW4SIKE
