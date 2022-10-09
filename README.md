# PROJECT : Computing the Polar Decomposition
 
The aim of this project is to do some basic numerical linear algebra, programming in MATLAB and LaTeX.

---

## Project Outline

We first introduce 

1. Vector norms
2. Matrix norms
3. Singular value decomposition

as the preliminaries.

Then we introduce the main theory of our project, the **polar decomposition**.

Furthermore, we introduce two numerical methods and combine them together

1. Newton iteration
2. Newton Schulz iteration

Finally we test our code using several test matrices.

---

## File Outline 

LaTeX related files are all located at the `document` folder, and the MATLAB related files are located at the `code` folder.

* File `proj_description.pdf` is the description of this project.
* File `.\document\proj_polar.pdf` is the report for this project.
* File `.\code\poldec.m` is the MATLAB code for the polar decomposition.
* File `.\code\polsqrt.m` is the MATLAB code for computing matrix square root using the polar decomposition.

Notice that some of the figures are located at `.\code\` rather than `.\document\figs\`, and this is due to the `export_fig` command (for export figures in MATLAB) saves the figure to the current directory which is `.\code\`.

---
