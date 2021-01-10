# MES

App written in C++ allows you to create 2D transient state thermal process simulations using finite element method.

Thermal phenomena occurring in the transient state are described by the Fourier equation in the
following form: 

![fourier equation](https://i.imgur.com/2swy97x.png)

The general solution can be descirbed by matrix equation:

![matrix form](https://i.imgur.com/LbMk0lB.png)

, where

![matrices](https://i.imgur.com/Uqq2Hra.png)

Example: (glass bottle, 6 degrees) -> -32 degrees Celsius
- cross-section

![nodes](https://i.imgur.com/dPOEOyt.png)

- output

![output](https://i.imgur.com/S3KPRJL.png)
