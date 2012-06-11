bhint
=====

bhint is a post-Newtonian, high-precision integrator for stellar systems surrounding a super-massive black hole.

 * Löckmann, U., &amp; Baumgardt, H. 2008 MNRAS, 384, 32 "Tracing intermediate-mass black holes in the Galactic Centre".

Original
---------

The **original** directory contains a clean version
of the ``bhint`` implementations by Löckmann and Baumgardt.

I fix some indentation issues,
to have a better looking code.

Slim
------

*On development*

The **slim** directory
aims to have a ``bhint`` version
with only the main idea of use
an Hermite Integrator Scheme
and the Kepler equation scenary,
it means, without the Post-newtonian corrections,
neither the SSE or GRAPE implementations.

GPU
-------

*On development*

The **gpu** directory
aims to contain an optimized ``bhint`` version.

The optimization was thinked to have
just a couple of modifications in the most
expensive functions (speaking of execution time).
