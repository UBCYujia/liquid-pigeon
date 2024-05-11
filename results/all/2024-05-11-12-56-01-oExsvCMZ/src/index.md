## Local communication barrier 

When the global communication barrier is large, many chains may 
be required to obtain tempered restarts.

The local communication barrier can be used to visualize the cause 
of a high global communication barrier. For example, if there is a 
sharp peak close to a reference constructed from the prior, it may 
be useful to switch to a [variational approximation](https://pigeons.run/dev/variational/#variational-pt).

```@raw html
<iframe src="local_barrier.svg" style="height:500px;width:100%;"></iframe>
<a href="local_barrier.svg"> 🔍 Full page </a>  ⏐<a href="https://pigeons.run/dev/output-pt/#Local-communication-barrier">🔗 Info </a>
```


## GCB estimation progress 

Estimate of the Global Communication Barrier (GCB) 
as a function of 
the adaptation round. 

The global communication barrier can be used 
to set the number of chains. 
The theoretical framework of [Syed et al., 2021](https://academic.oup.com/jrsssb/article/84/2/321/7056147)
yields that under simplifying assumptions, it is optimal to set the number of chains 
(the argument `n_chains` in `pigeons()`) to roughly 2Λ.

Last round estimate: ``1.9984014443252818e-14``

```@raw html
<iframe src="global_barrier_progress.svg" style="height:500px;width:100%;"></iframe>
<a href="global_barrier_progress.svg"> 🔍 Full page </a>  ⏐<a href="https://pigeons.run/dev/output-pt/#Global-communication-barrier">🔗 Info </a>
```


## Evidence estimation progress 

Estimate of the log normalization (computed using 
the stepping stone estimator) as a function of 
the adaptation round. 

Last round estimate: ``0.0``

```@raw html
<iframe src="stepping_stone_progress.svg" style="height:500px;width:100%;"></iframe>
<a href="stepping_stone_progress.svg"> 🔍 Full page </a>  ⏐<a href="https://pigeons.run/dev/output-normalization/">🔗 Info </a>
```


## Pigeons summary 

| **round** | **n\_scans** | **n\_tempered\_restarts** | **global\_barrier** | **global\_barrier\_variational** | **last\_round\_max\_time** | **last\_round\_max\_allocation** | **stepping\_stone** |
|----------:|-------------:|--------------------------:|--------------------:|---------------------------------:|---------------------------:|---------------------------------:|--------------------:|
| 1         | 2            | missing                   | 7.10543e-15         | missing                          | 0.341824                   | 1.96259e7                        | 1.77636e-14         |
| 2         | 4            | missing                   | 1.06581e-14         | missing                          | 0.0385048                  | 48320.0                          | 3.55271e-15         |
| 3         | 8            | missing                   | 1.06581e-14         | missing                          | 0.0896436                  | 91568.0                          | 8.88178e-15         |
| 4         | 16           | missing                   | 1.42109e-14         | missing                          | 0.168872                   | 178000.0                         | 1.11022e-15         |
| 5         | 32           | missing                   | 1.19904e-14         | missing                          | 0.371182                   | 209344.0                         | -2.22045e-16        |
| 6         | 64           | missing                   | 1.54321e-14         | missing                          | 0.699563                   | 256608.0                         | -1.77636e-15        |
| 7         | 128          | missing                   | 1.31006e-14         | missing                          | 1.47585                    | 371296.0                         | 1.33227e-15         |
| 8         | 256          | missing                   | 1.60982e-14         | missing                          | 2.84649                    | 600672.0                         | -8.88178e-16        |
| 9         | 512          | missing                   | 1.16573e-14         | missing                          | 5.71326                    | 1.05942e6                        | 4.44089e-16         |
| 10        | 1024         | missing                   | 1.9984e-14          | missing                          | 11.3263                    | 1.97693e6                        | 0.0                 |
 

```@raw html
<a href="Pigeons_summary.csv">💾 CSV</a> ⏐<a href="https://pigeons.run/dev/output-reports/">🔗 Info </a>
```


## Pigeons inputs 

| **Keys**               | **Values**                                                                              |
|-----------------------:|:----------------------------------------------------------------------------------------|
| extended\_traces       | false                                                                                   |
| checked\_round         | 0                                                                                       |
| extractor              | nothing                                                                                 |
| record                 | Function[Pigeons.log\_sum\_ratio, Pigeons.timing\_extrema, Pigeons.allocation\_extrema] |
| multithreaded          | false                                                                                   |
| show\_report           | true                                                                                    |
| n\_chains              | 10                                                                                      |
| variational            | nothing                                                                                 |
| explorer               | nothing                                                                                 |
| n\_chains\_variational | 0                                                                                       |
| target                 | StanLogPotential(simple\_model\_model)                                                  |
| n\_rounds              | 10                                                                                      |
| exec\_folder           | nothing                                                                                 |
| reference              | nothing                                                                                 |
| checkpoint             | false                                                                                   |
| seed                   | 1                                                                                       |
 

```@raw html
<a href="Pigeons_inputs.csv">💾 CSV</a> ⏐<a href="https://pigeons.run/dev/reference/#Pigeons.Inputs">🔗 Info </a>
```

