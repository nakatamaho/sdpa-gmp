# sdpa-gmp
SDPA in arbitrary multiple precision.

# how to build
I verified build on Ubuntu 20.04
```
rm -rf sdpa-gmp
git clone https://github.com/nakatamaho/sdpa-gmp.git
cd sdpa-gmp
aclocal ; autoconf ; automake --add-missing
autoreconf --force --install
./configure --enable-openmp=yes --enable-shared=yes
make -j4
```

# Citation
```
@INPROCEEDINGS{SDPA-GMP,
author={Nakata, Maho},
booktitle={2010 IEEE International Symposium on Computer-Aided Control System Design},
title={A numerical evaluation of highly accurate multiple-precision arithmetic version of semidefinite programming solver: SDPA-GMP, -QD and -DD.},
year={2010},  volume={},  number={},
pages={29-34},
doi={10.1109/CACSD.2010.5612693}
}
```
# SDPLIB
binary64 optimals are taken from
https://github.com/vsdp/SDPLIB/blob/master/README.md
.

| Problem   |    m |    n | Optimal (binary64)      | Optimal (GMP) |
| --------- | ---: | ---: | ----------------------: | :---: |
| arch0     |  174 |  335 |  5.66517e-01            | 5.6651727321592959e-01 |
| arch2     |  174 |  335 |  6.71515e-01            | 6.7151540763990793e-01 |
| arch4     |  174 |  335 |  9.726274e-01           | 9.7262741740980893e-01 |
| arch8     |  174 |  335 |  7.05698e+00            | 7.0569800367002555e+00 |
| control1  |   21 |   15 |  1.778463e+01           | 1.7784626717523405e+01 |
| control2  |   66 |   30 |  8.300000e+00           | 8.2999999857902351e+00 |
| control3  |  136 |   45 |  1.363327e+01           | 1.3633266228377313e+01 |
| control4  |  231 |   60 |  1.979423e+01           | 1.9794230376537536e+01 |
| control5  |  351 |   75 |  1.68836e+01            | 1.6883599050955793e+01 |
| control6  |  496 |   90 |  3.73044e+01            | 3.7304412616333280e+01 |
| control7  |  666 |  105 |  2.06251e+01            | 2.0625072243801761e+01 |
| control8  |  861 |  120 |  2.0286e+01             | 2.0286360379531460e+01 |
| control9  | 1081 |  135 |  1.46754e+01            | 1.4675424692813939e+01 |
| control10 | 1326 |  150 |  3.8533e+01             | 3.8533032079581028e+01 |
| control11 | 1596 |  165 |  3.1959e+01             | 3.1958667450372498e+01 |
| eqaulG11  |  801 |  801 |  6.291553e+02           | 6.2915529282061428e+02 |
| equalG51  | 1001 | 1001 |  4.005601e+03           | 4.0056013154878550e+03 |
| gpp100    |  101 |  100 | -4.49435e+01            | -4.4943550775891146e+01 |
| gpp124-1  |  125 |  124 | -7.3431e+00             | -7.3430762652465377e+00 |
| gpp124-2  |  125 |  124 | -4.68623e+01            | -4.6862295072749908e+01 |
| gpp124-3  |  125 |  124 | -1.53014e+02            | -1.5301412730175306e+02 |
| gpp124-4  |  125 |  124 | -4.1899e+02             | -4.1898762587351130e+02 |
| gpp250-1  |  250 |  250 | -1.5445e+01             | -1.5444916882934067e+01 |
| gpp250-2  |  250 |  250 | -8.1869e+01             | -8.1868958840223643e+01 |
| gpp250-3  |  250 |  250 | -3.035e+02              | -3.0353932422884198e+02 |
| gpp250-4  |  250 |  250 | -7.473e+02              | -7.4732831101269269e+02 |
| gpp500-1  |  501 |  500 | -2.53e+01               | -2.5320543879075787e+01 |
| gpp500-2  |  501 |  500 | -1.5606e+02             | -1.5606038757941642e+02 |
| gpp500-3  |  501 |  500 | -5.1302e+02             | -5.1301760233182234e+02 |
| gpp500-4  |  501 |  500 | -1.56702e+03            | -1.5670187921561449e+03 |
| hinf1     |   13 |   14 |  2.0326e+00             |  2.6899806999311550e-05 |
| hinf2     |   13 |   16 |  1.0967e+01             |  1.0967055621049256e+01 |
| hinf3     |   13 |   16 |  5.69e+01               |  5.6940778009669388e+01 |
| hinf4     |   13 |   16 |  2.74764e+02            |  2.7149772617088246e+02 |
| hinf5     |   13 |   16 |  3.63e+02               |  3.3010372908768509e+02 |
| hinf6     |   13 |   16 |  4.490e+02              |  4.4892774532835125e+02 |
| hinf7     |   13 |   16 |  3.91e+02               |  1.5490470173390994e+02 |
| hinf8     |   13 |   16 |  1.16e+02               |     5 |
| hinf9     |   13 |   16 |  2.3625e+02             |     5 |
| hinf10    |   21 |   18 |  1.09e+02               |     5 |
| hinf11    |   31 |   22 |  6.59e+01               |     5 |
| hinf12    |   43 |   24 |  2e-1                   |     5 |
| hinf13    |   57 |   30 |  4.6e+01                |     5 |
| hinf14    |   73 |   34 |  1.30e+01               |     5 |
| hinf15    |   91 |   37 |  2.5e+01                |     5 |
| infd1     |   10 |   30 |  dual infeasible        |     6 |
| infd2     |   10 |   30 |  dual infeasible        |     6 |
| infp1     |   10 |   30 |  primal infeasible      |     6 |
| infp2     |   10 |   30 |  primal infeasible      |     6 |
| maxG11    |  800 |  800 |  6.291648e+02           |     7 |
| maxG32    | 2000 | 2000 |  1.567640e+03           |     7 |
| maxG51    | 1000 | 1000 |  4.003809e+03           |     7 |
| maxG55    | 5000 | 5000 |  9.999210e+03           |     7 |
| maxG60    | 7000 | 7000 |  1.522227e+04           |     7 |
| mcp100    |  100 |  100 |  2.261574e+02           |     8 |
| mcp124-1  |  124 |  124 |  1.419905e+02           |     8 |
| mcp124-2  |  124 |  124 |  2.698802e+02           |     8 |
| mcp124-3  |  124 |  124 |  4.677501e+02           |     8 |
| mcp124-4  |  124 |  124 |  8.644119e+02           |     8 |
| mcp250-1  |  250 |  250 |  3.172643e+02           |     8 |
| mcp250-2  |  250 |  250 |  5.319301e+02           |     8 |
| mcp250-3  |  250 |  250 |  9.811726e+02           |     8 |
| mcp250-4  |  250 |  250 |  1.681960e+03           |     8 |
| mcp500-1  |  500 |  500 |  5.981485e+02           |     8 |
| mcp500-2  |  500 |  500 |  1.070057e+03           |     8 |
| mcp500-3  |  500 |  500 |  1.847970e+03           |     8 |
| mcp500-4  |  500 |  500 |  3.566738e+03           |     8 |
| qap5      |  136 |   26 | -4.360e+02              |     9 |
| qap6      |  229 |   37 | -3.8144e+02             |     9 |
| qap7      |  358 |   50 | -4.25e+02               |     9 |
| qap8      |  529 |   65 | -7.57e+02               |     9 |
| qap9      |  748 |   82 | -1.410e+03              |     9 |
| qap10     | 1021 |  101 | -1.093e+01              |     9 |
| qpG11     |  800 | 1600 |  2.448659e+03           |    10 |
| qpG51     | 1000 | 2000 |  1.181000e+03           |    10 |
| ss30      |  132 |  426 |  2.02395e+01            |     1 |
| theta1    |  104 |   50 |  2.300000e+01           |    11 |
| theta2    |  498 |  100 |  3.287917e+01           |    11 |
| theta3    | 1106 |  150 |  4.216698e+01           |    11 |
| theta4    | 1949 |  200 |  5.032122e+01           |  5.0321221951837344e+01 |
| theta5    | 3028 |  250 |  5.723231e+01           |  5.7232307282180003e+01 |
| theta6    | 4375 |  300 |  6.347709e+01           |  6.3477087177964743e+01 |
| thetaG11  | 2401 |  801 |  4.000000e+02           |  4.0000000000000000e+02 |
| thetaG51  | 6910 | 1001 |  3.49000e+02            |  3.4900000000000000e+02 |
| truss1    |    6 |   13 | -8.999996e+00           | -8.9999963152868905e+00 |
| truss2    |   58 |  133 | -1.233804e+02           | -1.2338035636407390e+02 |
| truss3    |   27 |   31 | -9.109996e+00           | -9.1099962092020534e+00 |
| truss4    |   12 |   19 | -9.009996e+00           | -9.0099962910045294e+00 |
| truss5    |  208 |  331 | -1.326357e+02           | -1.3263567797250604e+02 |
| truss6    |  172 |  451 | -9.01001e+02            | -9.0100140477088096e+02 |
| truss7    |   86 |  301 | -9.00001e+02            | -9.0000140369343463e+02 |
| truss8    |  496 |  628 | -1.331146e+02           | -1.3311458915226341e+02 |



