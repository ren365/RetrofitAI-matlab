#### Testing Files

**GIF explanation**: the first 4m shows the true dynamics is different from the dynamics and the model learns the error after 4m.

**QP status**(https://osqp.org/docs/interfaces/status_values.html#status-values):

| Status                     | Constant               | Value |
| -------------------------- | ---------------------- | ----- |
| solved                     | OSQP_SOLVED            | 1     |
| solved inaccurate          | OSQP_SOLVED_INACCURATE | 2     |
| maximum iterations reached | OSQP_MAX_ITER_REACHED  | -2    |

**File explanation**

| file        | status | description                                                  |
| ----------- | ------ | ------------------------------------------------------------ |
| movingTest1 | unsafe | Original settings (steering_limit=0.75,max_accel=1,min_accel=-1) |
| movingTest2 | unsafe | changing the vehicle limitation (steering_limit=1.25,max_accel=3,min_accel=-2) |
| movingTest3 | safe   | setting as movingTest1, use barrier_radius=1.2 as pseud-radius and 1.0 as actual radius |
| movingTest4 | unsafe | use predicted position to put into cbf                       |
|             |        |                                                              |

