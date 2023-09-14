crimCV
======

This repository contains the source code of a
[CRAN](https://cran.r-project.org/mirrors.html) package called
[crimCV](https://cran.r-project.org/web/packages/crimCV/index.html).  `crimCV`
fits trajectories using a finite mixture of zero-inflated Poisson models.  This
source code isn't of use to the majority of users and if you would just like to use it
you can can install it by starting [R](https://www.r-project.org/) and typing:

```R
> install.packages("crimCV")
```
Once installed you can use by loading the library:

```R
> library("crimCV")
```

Developing
----------

If you wish to help develop `crimCV` the best way to do so is with
[devtools](https://cran.r-project.org/web/packages/devtools/index.html).  With
`devtools` installed clone this repository to your current working directory,
start `R`, and type:
```R
> library("devtools")
> load_all("crimCV")
```
to build the local source and load the library into `R`.

License
-------

[GPL](https://github.com/drjdn/crimCV/blob/master/LICENSE.md)
