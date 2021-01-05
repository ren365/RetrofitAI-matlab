# RetrofitAI-matlab
MATLAB version of safe control with supervised learning

requirement: OSQP (MATLAB version)

```matlab
install_osqp
```

The original repo: https://github.com/ren365/RetrofitAI

The structure of the codes. The file name remains almost the same as original python files.

| Name                                           | Type   | Decryption                                         | Status             |
| ---------------------------------------------- | ------ | -------------------------------------------------- | ------------------ |
| original python                                | folder | original python files                              | done               |
| osqp                                           | folder | installed lib for QP solver                        | done               |
| test codes                                     | folder | test whether matlab & python codes behave the same | under implementing |
| DynamicsAckermannZModified/ DynamicsAckermannZ | .m     | dynamics                                           | tested             |
| einsum                                         | .m     | helper function for numpy.einsum                   | tested             |
| LyapunovAckermannZ                             | .m     |                                                    | tested             |
| install_osqp                                   | .m     | automatically install the osqp on matlab           | done               |
| QPSolver                                       | .m     | CBF                                                | tested             |
| test_adaptive_clbf                             | .m     | main file to run                                   | untested           |
| gaussian                                       | .m     | gaussian instead of NN                             | undergoing         |
| adaptive_clbf                                  | .m     | combine NN and CBF                                 | untested           |
| BarrierAckermannPointZ                         | .m     | (cbf.py) update the barriers position              | tested             |
| BarrierAckermannVelocityZ                      | .m     | (cbf.py) update the barriers position              | tested             |
|                                                |        |                                                    |                    |


