<div id="table-of-contents">
<h2>Table of Contents</h2>
<div id="text-table-of-contents">
<ul>
<li><a href="#orgheadline1">1. Files in this project</a></li>
<li><a href="#orgheadline2">2. Program files</a></li>
<li><a href="#orgheadline5">3. Compilation/installation</a>
<ul>
<li><a href="#orgheadline3">3.1. Prerequisite</a></li>
<li><a href="#orgheadline4">3.2. Compilation</a></li>
</ul>
</li>
</ul>
</div>
</div>


# Files in this project<a id="orgheadline1"></a>

-   **README.{org,md}:** this file.
-   **Doxyfile:** config file to build the automatic documentation.
-   **Makefile:** main Makefile for compilation.
-   **Makefile.inc:** macro common to all makefils: name and options of the compiler.
-   **doc:** directory for the Doxygen documentation. Directory create by doxygen.
-   **meshes:** exemple meshes in various format
-   **list.f90:** linked list, for integers and real numbers only.
-   **mat\_csc.f90:** matrices in the Compressed Sparse Column storage scheme.
-   **mat\_csr.f90:** matrices in the Compressed Sparse Row storage scheme.
-   **mat\_list.f90:** matrices as arrays of lists. Each list in the array represent
    a row.
-   **matrices.f90:** include file to use for all matrix operations.
-   **prec.f90:** declare the precision of numbers.
-   **sort.f90:** methods for sorting arrays of integers or reals.
-   **umfpack\_calls.f90:** module used to call the routines from the UMFPACK
    library.
-   **fe\_mesh.f90:** defines a mesh structure and methods on it, notably i/o and
    useful accessors to members.
-   **p1.f90:** assemble the mass and stiffness martrices with P1-Lagrange finite
    elements.
-   **vtkCellType.inc:** include file from the VTK library, used to write vtk files
    in fe\_mesh.f90

# Program files<a id="orgheadline2"></a>

-   **data.txt:** entry file that contains data for the laplace program
-   **laplace.f90:** program that solves the PDE alpha\*u-Delta(u) = f with
    homogeneous boundary conditions
-   **test\_csc.f90:** program used to test the CSC matrix format and resolution of
    linear system with UMFPACK.

# Compilation/installation<a id="orgheadline5"></a>

## Prerequisite<a id="orgheadline3"></a>

-   A fortran compiler
-   The gnu make utility
-   The UMFPACK library, which is part of the SuiteSparse package. You can find it
    on T. Davis web page: [<http://faculty.cse.tamu.edu/davis/suitesparse.html>](http://faculty.cse.tamu.edu/davis/suitesparse.html)

## Compilation<a id="orgheadline4"></a>

To compile, simply run 'make'. 

You can run 'make clean' to delete the local object files. The compiled exec files
and will remain. Type 'make distclean' to remove everything.
# TER_2A
