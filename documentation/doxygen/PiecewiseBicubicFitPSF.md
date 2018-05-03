PiecewiseBicubicFitPSF {#PiecewiseBicubicFitPSF_main_page}
========================

  \brief Fitting for position dependent values and derivatives of piecewise
  bicubic PSF models.

  In the piecewise bicubic PSF model the PSF is split into a grid of cells
  and the intensity over each cell is assumed to be given by a bi-cubic
  function.

  Requiring the continuity of the first order derivatives with respect to x
  and y as well as the cross derivative of the PSF dramatically reduces the
  number of free parameters. In fact the PSF is fully specified by the
  values, the x and y derivatives and the xy cross derivatives at the grid
  points.

  Since we want the PSF to vary smoothly accross the image we will assume
  that these quantities are given by a polynomial expansion accross the
  image. Further, the grid points lying on the outside edge of the grid will
  be assumed to have a value and all derivatives equal to zero.

  The grid giving the splitting of the PSF will be defined by a set of
  vertical boundaries with coordinates \f$ \{x_i\} \f$ (i=0...M) relative to
  the source center and horizontal boundaries with coordinates
  \f$ \{y_j\} \f$ (j=0...N) relative to the source center. Let
  \f$ \{V_{i,j}\} \f$ be the values, \f$ \{D^x_{i,j}\} \f$ the horizontal
  derivatives, \f$ \{D^y_{i,j}\} \f$ the vertical derivatives, and
  \f$ \{D^{xy}_{i,j}\} \f$ the cross derivatives of the PSF at the grid
  points. Let \f$ \mathbf{Q} \f$ be a vector of these quantities excluding
  those on the outside grid edge.
  
  We need to impose that the integral of the PSF over the entire grid is one
  (scaling this by a factor is degenerate with scaling the amplitudes for all
  sources by the inverse factor). That can be imposed by noting that the
  integral of the PSF over the grid is given by \f$ \mathbf{I^TQ} \f$ for
  some column vector \f$ \mathbf{I} \f$. We then decompose
  \f$ \mathbf{Q}=\mathbf{I}/I^2+\sum q_i \mathbf{e}_i \f$ where \f$ I^2 \f$
  is the square norm of \f$ \mathbf{I} \f$ and \f$ \mathbf{e}_i \f$ are an
  orthonormal basis for the space perpendicular to \f$ \mathbf{I} \f$. Note
  that the SVD decomposition of \f$ \mathbf{I}^T \f$ has
  \f$ \mathbf{U}=(1) \f$, \f$ \mathbf{S}=(I, 0, 0, ...) \f$ and the first
  columns of \f$ \mathbf{V} \f$ is \f$ \mathbf{I}/I \f$ and the remaining
  columns can be used as the \f$ \mathbf{e}_i \f$. We will then fit for the
  spatially dependent \f$ q_i \f$:

  \f[
  	q_i=\sum_{k,l}q_{i,k,l}x^ky^l
  \f]

  where \f$ x \f$ and \f$ y \f$ are the image coordinates of the center of
  the source for which the PSF is being constructed.
  
  With these definitions, in order to fully specify PSF of each single source
  in an input image we need the quantities \f$ q_{i,k,l} \f$, where i is in
  the range 0 to (M-2)(N-2)-1, and the ranges for k and l depend on the order
  of the polynomial used to describe the PSF variability over the image, so
  we need to fit for those.

  Note that the coefficients of the bi-cubic polynomials specifying the PSF
  in each grid cell depend linearly on the values and derivatives at the grid
  nodes, which can be reconstructed from the quantities we fit for. In turn,
  the integrals of the PSF over CCD (sub-)pixels depends linearly on the
  bi-cubic coefficients.
  
  So if we somehow knew the fluxes of all the sources as well as their
  positions, we could constuct a matrix (\f$ \mathbf{M} \f$) such that:
  \f$ \mathbf{Mq}+\mathbf{\Delta} \f$ predicts the values of all image pixels
  assigned to any source scaled by the flux of that source. Where
  \f$ \mathbf{q} \f$ is the vector of the \f$ q_{i,k,l} \f$ unknowns and
  \f$ \Delta \f$ is the integral of a PSF specified entirely by
  \f$ \mathbf{I} \f$ over the corresponding pixel. So \f$ \mathbf{q} \f$ can
  be obtained by multiplying the rows of \f$ \mathbf{M} \f$ and
  \f$ \mathbf{\Delta} \f$ by the amplitude of the corresponding source and
  percforming a least squares fit. Since deriving the fluxes of the sources
  requires the PSF, we iterate. We begin by estimating the flux of each
  source by simply summing all the pixel values (in excess of the background)
  assigned to it or by aperture photometry. Then for each source we multiply
  the piece of the \f$ \mathbf{M} \f$ matrix that corresponds to that
  source's pixels by the current flux estimate and use the new matrix in a
  linear fit for the PSF.  Finally, we update the flux estimates by using the
  new PSF model.  This is repeated until the flux estimates stop changing too
  much.
  
  The Algorithm
  -------------

  Directly constructing \f$ \mathbf{M} \f$ is counter productive. For even
  modest grids (say 12x12 cells) and a reasonable number of sources used from
  a real image (say 3000 each with an average of 50 pixels assigned to it)
  and a reasonable expansion order (say 4) would result in unmanageable
  matrix sizes (\f$ 2\times10^9 \f$ elements in our example). Such matrices
  are hard to even allocate let alone solve. However, \f$ \mathbf{M} \f$ is
  of a special form which makes this unnecessary. In order to perform a least
  squares fit we need to multiply the RHS by \f$ \mathbf{M}^T \f$ (let the
  result be denoted by \f$ \mathbf{r'} \f$ and then solve 
  \f$ \mathbf{M}^T\,\mathbf{M}\,\mathbf{x}=\mathbf{r'} \f$. The matrix 
  \f$ \Lambda\equiv\mathbf{M}^T\,\mathbf{M} \f$ has only \f$ n^2K^2 \f$
  elements (\f$ 5.3 \times 10^7 \f$ in our example), which is reasonable to
  allocate and invert. Even more, calculating \f$ \mathbf{\Lambda} \f$ by
  first calculating \f$ \mathbf{M} \f$ is slower than the algorithm described
  below, and actually will dominate the computational time for non-trivial
  uses.

  Let: 
   * \f$ \mathcal{O} \f$ - the polynomial expansion order
   * \f$ G_x,G_y \f$ - number of grid boundaries in the x and y direction
                       respectively
   * S - number of sources
   * \f$ p_i \f$ - number of pixels assigned to the i-th source.
   * P - number of pixels in the input image assigned to all sources (i.e.
         \f$ P\equiv\sum_{i=0}^S p_i \f$)
   * K - number of terms in the polynomial expansion of a single PSF
         coefficient (i.e. K=\f$ \mathcal{O}+1)(\mathcal{O}+2)/2 \f$).
   * n - number of parameters defining a single psf (i.e.
         n=\f$ 4(G_x-2)(G_y-2)-1 \f$)
   * \f$ r \f$ - vector of the background excess values of all pixels
                 assigned to sources
   * \f$ \mathbf{I} \f$ - vector which when dotted with the PSF parameters
                          (\f$ \{V_{i,j}\} \f$, \f$ \{D^x_{i,j}\} \f$,
						   \f$ \{D^y_{i,j}\} \f$, and
						   \f$ \{D^{xy}_{i,j}\} \f$)
						  gives the integral of the PSF over the entire grid.
   * \f$ \Delta \f$ - Integrals of PSFs whose parameters are given by
                    \f$ \mathbf{I}/I^2 \f$ over the source pixels.
   * \f$ \mathbf{\tilde{M}}^i \f$ - a set of matrices, one for each source,
                                    with each row consisting of the integrals
									over the same source pixel but with only
									one \f$ q_i \f$ being non-zero. The
									dimensions of the i-th matrix
									are \f$ p_i\times n \f$.
   * \f$ \mathbf{\Lambda}^i\equiv(\mathbf{\tilde{M}}^i)^T\mathbf{\tilde{M}}^i \f$
   * \f$ \mathbf{\kappa^i} \f$ - a set of row-vectors, one for each source,
	                             containing the polynomial expansion terms
								 (i.e.\f$ \left(1, x, y, x^2, xy, \cdots,
								          y^\mathcal{O}\right) \f$)
   * \f$ f_i \f$ - current best estimate of the source flux

  With these definitions:
  \f[
    \mathbf{\Lambda}=\left(\begin{array}{ccccccc}
		\nwarrow & \uparrow & \nearrow &
		\nwarrow & \uparrow & \nearrow \\

		\longleftarrow & 
        \sum_{i=0}^S f_i^2 \Lambda^i_{1,1} 
				(\mathbf{\kappa^i})^T \mathbf{\kappa^i} &
        \longrightarrow &
        \longleftarrow &
        \sum_{i=0}^S f_i^2 \Lambda^i_{1,2} 
				(\mathbf{\kappa^i})^T\mathbf{\kappa^i} &
        \longrightarrow &
        \cdots \\
		\swarrow & \downarrow & \searrow &
		\swarrow & \downarrow & \searrow \\
        \\
		\nwarrow & \uparrow & \nearrow &
		\nwarrow & \uparrow & \nearrow \\
        \longleftarrow &
        \sum_{i=0}^S f_i^2 \Lambda^i_{2,1} 
				(\mathbf{\kappa^i})^T\mathbf{\kappa^i} &
        \longrightarrow &
        \longleftarrow &
        \sum_{i=0}^S f_i^2 \Lambda^i_{2,2} 
				(\mathbf{\kappa^i})^T\mathbf{\kappa^i} &
        \longrightarrow &
        \cdots \\
		\swarrow & \downarrow & \searrow &
		\swarrow & \downarrow & \searrow \\
        & \vdots & & & \vdots & & \ddots
    \end{array}\right)
  \f]

  So the algorithm we will use for PSF fitting is as follows. 
  
  ### One Time Initialization ###

  Overall Complexity: \f$ P\,n^2 \f$

  1. Obtain an initial estimate of the source fluxes (aperture photometry
  with flat PSF), the \f$ r \f$ and \f$ \Delta \f$ vectors.
  __Complexity__: large constant times S.

  2. For each source: 

     2.1 construct the corresponding \f$ \mathbf{\tilde{M}^i} \f$
         matrix and stack them together in \f$ \mathbf{\tilde{M}} \f$.
		 __Complexity__: large constant times \f$ n\,S \f$.

     2.2 For each source apply \f$ (\mathbf{\tilde{M}^i})^T \f$ to the
         corresponding piece of the \f$ r \f$ vector, storing the result as
		 \f$ \mathbf{\tilde{r}}^i \f$ a part of a permanent vector
		 \f$ \mathbf{\tilde{r}} \f$ (will have a total length of
		 \f$ S\,n \f$). __Complexity__: \f$ P\,n \f$.

     2.3 Repeat \f$ \Delta \f$ storing the result as
         \f$ \mathbf{\tilde{\Delta}} \f$.

     2.4 Calculate
		 \f$ \mathbf{\tilde{\Lambda}^i}
		 	 =
			 (\mathbf{\tilde{M}^i})^T\,\mathbf{\tilde{M}^i} \f$
		 and store the result as
		 part of the permament matrix \f$ \mathbf{\tilde{\Lambda}} \f$
		 (a stack of the individual \f$ \mathbf{\tilde{\Lambda}^i} \f$
		  matrices). __Complexity__: \f$ P\,n^2 \f$. 
  
  3. Calculate a stack of \f$ (\mathbf{\kappa^i})^T\mathbf{\kappa^i} \f$
  matrices (__complexity__: \f$ S\,K^2 \f$).

  ### Estimates for The Expansion Coefficients of the PSF Parameters ###

  Overall Complexity \f$ S\,n^2\,K^2 \f$

  1. Compute the \f$ \mathbf{\Lambda} \f$ matrix as it is written above
  (__complexity__: \f$ S\,n^2\,K^2 \f$).

  2. Compute:
  \f[
    \mathbf{r'}=\left(\begin{array}{c}
        \sum_{i=0}^S (\mathbf{\tilde{\kappa}^i})^T f_i
		(\tilde{r}^i_1-f_i\tilde{\Delta}^i_1)\\
        \sum_{i=0}^S (\mathbf{\tilde{\kappa}^i})^T f_i
		(\tilde{r}^i_2-f_i\tilde{\Delta}^i_2)\\
        \vdots\\
        \sum_{i=0}^S (\mathbf{\tilde{\kappa}^i})^T f_i
		(\tilde{r}^i_n-f_i\tilde{\Delta}^i_n\\
    \end{array}\right)
  \f]
  __complexity__: \f$ n\,K \f$
  3. Use LDLT decomposition to solve for the (updated) expansion
  coefficients (__complexity__: \f$ n^3\,K^3 \f$).

  ### Estimates for the Source Amplitudes. ###

  This is a simple linear regression for the background excesses of each
  source using its current best fit PSF.

  Smoothing.
  ----------

  In the case of fitting extremely sharp sources or very sparse images it is
  possible that there will not be enough stars or enough pixels within a star
  to constrain a grid, and yet there may be features of the PSF which require
  a fine grid to be used. In this case noise ends up being fitted and the
  resulting PSF has large amplitude high frequency components. This is
  clearly not physical and needs to be avoided in order to get stable
  photometry. In order to handle these cases piecewise bicubic PSF fitting
  support smoothing. This is implemented by adding the integral of the square
  of the second derivative in both x and y over the entire PSF to the
  function being minimized:

  \f[
	\tilde{\chi}^2 \equiv
  	\frac{10^\sigma}{S}\sum_{s=1}^S 
	\left(
		\frac{\partial^4 f_s(x,y)} {\partial x^2 \partial y^2}
	\right)^2
	+
	\sum_{s=1}^S \sum_{j=1}^{p_i} 
		\left(A_i\int_{pixel(i,j)} f_s(x,y) dx dy - r_i^j\right)^2
  \f]

  Adding the first term above, which is what is missing from the description
  above is actually quite simple. It can be expressed as the sum of terms
  containing constants and second order combinations of the PSF parameters.
  This in turn can be achieved by adding an extra
  \f$ \mathbf{\tilde{M}}^i \f$ matrix and a corresponding source with all
  zero fluxes. However, this extra matrix will not be multiplied by
  a \f$ \kappa^T\kappa \f$, but instead by the average of that over all
  sources.

  In order to add the first we need to write it as a sum of squares of the
  difference between a linear transformation of the \f$ q_i \f$ parameters
  and some right hand side vector.

  In the following let i and j be indices identifyin a grid cell (i is the
  grid column and j is the grid row) and let \f$ c^{i,j}_{m,n;s} \f$
  be the coefficient in front of \f$ x^m y^n \f$ term in the PSF funciton
  of source number \f$s\f$ over the \f$ (i,\ j) \f$ cell. In what
  follows we will drop the i, j and s indices for shorter notation. Then:

  \f{eqnarray*}{
    \int_{grid} dx dy
	\left(
		\frac{\partial^4 f_s(x,y)}{\partial x^2 \partial y^2}
	\right)^2
	=
	\sum_{i,j} \Bigg(&&					\\
		&& 16 c_{2,2}^2 w h 			\\
		&& + 48 c_{2,3}^2 w h^3		\\
		&& + 48 c_{3,2}^2 w^3 h		\\
		&& + 144 c_{3,3}^2 w^3 h^3		\\
		&& + 48 c_{2,2} c_{2,3} w h^2	\\
		&& + 48 c_{2,2} c_{3,2} h^2 w	\\
		&& + 72 c_{2,2} c_{3,3} w^2 h^2	\\
		&& + 72 c_{2,3} c_{3,2} w^2 h^2	\\
		&& + 144 c_{2,3} c_{3,3} w^2 h^3	\\
		&& + 144 c_{3,2} c_{3,3} w^3 h^2	\\
	\Bigg)								\\
	=
	\sum_{i,j} w h \Bigg(&&				\\
		&& 16 c_{2,2}^2					\\
		&& + 48 c_{2,3}^2 h^2			\\
		&& + 48 c_{3,2}^2 w^2			\\
		&& + 144 c_{3,3}^2 w^2 h^2		\\
		&& + 48 c_{2,2} c_{2,3} h		\\
		&& + 48 c_{2,2} c{3,2} w		\\
		&& + 72 c_{2,2} c_{3,3} w h		\\
		&& + 72 c_{2,3} c_{3,2} w h		\\
		&& + 144 c_{2,3} c_{3,3} w h^2	\\
		&& + 144 c_{3,2} c_{3,3} w^2 h	\\
	\Bigg)
  \f}
  where \f$ w,\ h \f$ are the width and height of the grid cell. Since this
  is a quadratic form in \f$ c_{m,n} \f$, it can be written as:
  \f[
  	\int_{grid} dx dy
	\left(
		\frac{\partial^4 f_{x_0, y_0}(x,y)}{\partial x^2 \partial y^2}
	\right)^2 = \mathbf{C}^T \mathbf{A} \mathbf{C}
  \f]
  where \f$ \mathbf{C}^T=(C^{0,0}_{2,2}, C^{0,0}_{2,3}, C^{0,0}_{3,2},
			C^{0,0}_{3,3}, C^{1,0}, \ldots) \f$ and A is a symmetric
  matrix which consists of matrices of the following form along its diagonal:
  \f[
	\mathbf{A_{cell}}=w h \left(\begin{array}{cccc}
		16		&	24 h		&	24 w		&	36 w h		\\
		24 h	&	48 h^2		&	36 w h		&	72 w h^2	\\
		24 w 	&	36 w h		&	48 w^2		&	72 w^2 h	\\
		36 w h 	&	72 w h^2	&	72 w^2 h	&	144 w^2 h^2
	\end{array}\right)
  \f]
  Since \f$ \mathbf{A} \f$ is symmetric it can be written as 
  \f$ \mathbf{E}^T \mathbf{\Lambda} \mathbf{E} \f$ where \f$ \mathbf{E} \f$
  is orthogonal and \f$ \mathbf{\Lambda} \f$ is diagonal containing the
  eigenvalues of \f$ \mathbf{A} \f$. Which can clearly be written as
  \f$ \mathbf{A} = \mathbf{B}^T\mathbf{B} \f$ with
  \f$ \mathbf{B}=\sqrt{\mathbf{\Lambda}} \mathbf{E} \f$.

  Further: 
  \f[
	\mathbf{C} = \mathbf{\Phi} \mathbf{q} = \mathbf{\Phi}\left(
		\mathbf{I}/I^2 + \sum_{i,k} q_{i,k} \kappa_k \mathbf{e}_i
	\right)
  \f]

  So we end up with:
  \f[
  	S\equiv \sum_{s=1}^S \int_{grid} dx dy
	\left(
		\frac{\partial^4 f_s(x,y)}{\partial x^2 \partial y^2}
	\right)^2 = \sum_{s=1}^S \sum_{c} \left[ \sum_{d,e}
		B_{c,d} \Phi_{d,e} \left(
			I_e/I^2 + \sum_{i,k} q_{i,k} \kappa^s_k \mathbf{e}_{i,e}
		\right)
	\right]^2
  \f]

  To find the minimum differentiate:
  \f{eqnarray*}{
  	\frac{\partial S}{\partial q_{\sigma,\mu}} &=&
	2 \sum_{s=1}^S \sum_{c} \left[
		\sum_{d,e} B_{c,d} \Phi_{d,e} \left(
			I_e/I^2 + \sum_{i,k} q_{i,k} \kappa^s_k \mathbf{e}_{i,e}
		\right)
	\right]
	\left[
		\sum_{\delta,\varepsilon} B_{c,\delta} \Phi_{\delta,\varepsilon}
		\left(
			\kappa^s_\mu e_{\sigma,\varepsilon}
		\right)
	\right]\\
	&=&
	2 \sum_{s=1}^S \kappa^s_{\mu}
		\mathbf{e}_{\sigma}^T\mathbf{\Phi}^T\mathbf{B}^T \mathbf{B}
		\mathbf{\Phi}
		\left(
			\mathbf{I}/I^2
			+
			\sum_{i,k} q_{i,k} \kappa^s_k \mathbf{e}_{i}
		\right)\\
	&=&
	\frac{2}{I^2} 
    \mathbf{e}_{\sigma}^T\mathbf{\Phi}^T\mathbf{A}\mathbf{\Phi} \mathbf{I}
    \sum_{s=1}^S \kappa^s_{\mu}
	+ 
	\frac{2}{I^2}
    \mathbf{e}_{\sigma}^T\mathbf{\Phi}^T\mathbf{A}\mathbf{\Phi} 
    \sum_{s=1}^S \sum_{i,k} q_{i,k} 
        \kappa^s_{\mu} \kappa^s_k \mathbf{e}_{i}

  \f}
  The first term needs to be added to the \f$ \mathbf{r}' \f$ vector and the
  matrix identified by the second term must be added to the
  \f$ \mathbf{\Lambda} \f$ matrix.

  Conveniently, neither correction depends on the current values of the
  parameters, so we can compute those once for a given grid and simply add
  them at each step, thus not resulting in any significant overhead to the
  fitting. Further the non-smoothing algorithm already computes the
  bicubic coefficients over each grid cell for the selected  basis vector so
  we have \f$ \mathbf{\Phi} \mathbf{I}/I^2 \f$ and
  \f$ \mathbf{\Phi} \mathbf{e}_i \f$ for all \f$ \mathbf{e}_i \f$.
