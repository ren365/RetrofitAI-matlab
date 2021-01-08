# RetrofitAI-matlab
MATLAB version of safe control with supervised learning

requirement: OSQP (MATLAB version)

Running instruction

```matlab
% go to the project path RetrofitAI-matlab
install_osqp % install the requirement
addpath("translated matlab")
run test_adaptive_clbf.m
```

The original repo: https://github.com/ren365/RetrofitAI

Reference paper: https://arxiv.org/abs/1910.02325 

The structure of the folder.

| Name              | Type   | Description                                               | Status |
| ----------------- | ------ | --------------------------------------------------------- | ------ |
| original python   | folder | [test] original python files                              | done   |
| osqp              | folder | installed lib for QP solver                               | done   |
| test codes        | folder | [test] test whether matlab & python codes behave the same | done   |
| translated matlab | folder | [main] matlab files                                       | done   |
| install_osqp      | .m     | automatically install the osqp on matlab                  | done   |
|                   |        |                                                           |        |

Inside the "translated matlab" folder, the descriptions of files are shown in the table below.

| Name                       | Type | Description                          |
| -------------------------- | ---- | ------------------------------------ |
| test_adaptive_clbf         | .m   | main file to run                     |
| adaptive_clbf              | .m   | combine GP and CBF                   |
| QPSolver                   | .m   | CBF                                  |
| ModelGP                    | .m   | Interface for gaussian process model |
| ScaledGP                   | .m   | Gaussian process model               |
| DynamicsAckermannZ         | .m   | Approximate dynamics                 |
| DynamicsAckermannZModified | .m   | True dynamics                        |
| BarrierAckermannPointZ     | .m   | barriers class                       |
| LyapunovAckermannZ         | .m   | CLF                                  |
| circles                    | .m   | outside helper function              |
| einsum                     | .m   | outside helper function              |
|                            |      |                                      |


