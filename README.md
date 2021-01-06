# RetrofitAI-matlab
MATLAB version of safe control with supervised learning

requirement: OSQP (MATLAB version)

```matlab
install_osqp
```

The original repo: https://github.com/ren365/RetrofitAI

The structure of the codes. The file name remains almost the same as original python files.

| Name                                           | Type   | Decryption                                         | Status            |
| ---------------------------------------------- | ------ | -------------------------------------------------- | ----------------- |
| original python                                | folder | original python files                              | done              |
| osqp                                           | folder | installed lib for QP solver                        | done              |
| test codes                                     | folder | test whether matlab & python codes behave the same | undergoing        |
| install_osqp                                   | .m     | automatically install the osqp on matlab           | tested            |
| DynamicsAckermannZModified/ DynamicsAckermannZ | .m     | dynamics                                           | tested            |
| einsum                                         | .m     | helper function for numpy.einsum                   | tested            |
| LyapunovAckermannZ                             | .m     |                                                    | tested            |
| QPSolver                                       | .m     | CBF                                                | tested            |
| BarrierAckermannPointZ                         | .m     | (cbf.py) update the barriers position              | tested            |
| test_adaptive_clbf                             | .m     | main file to run                                   | untested, qp pass |
| adaptive_clbf                                  | .m     | combine NN and CBF                                 | untested, qp pass |
| ModelGP                                        | .m     | gaussian instead of NN                             | untested          |
| ScaledGP                                       | .m     | gaussian instead of NN                             | undergoing        |
|                                                |        |                                                    |                   |


