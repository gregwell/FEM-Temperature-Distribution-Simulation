# MES

This app allows you to create 2D transient state thermal process simulations using finite element method.

Thermal phenomena occurring in the transient state are described by the Fourier equation in the
following form: 

![fourier equation](https://i.imgur.com/2swy97x.png)

The general solution can be descirbed by matrix equation:

![matrix form](https://i.imgur.com/LbMk0lB.png)

, where

![matrices](https://i.imgur.com/Uqq2Hra.png)

Sample simulation: Glass bottle (6°C, filled with 40% vodka) was put in the freezer (-32°C). 
- part of the cross section in the center of a bottle with selected nodes

![nodes](https://i.imgur.com/dPOEOyt.png)

- output: temperatures of selected nodes after each 60 second step time

![output](https://i.imgur.com/S3KPRJL.png)
