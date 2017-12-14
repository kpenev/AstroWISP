Integrals of the PSF {#PSF_integrals_page}
====================

In this package, integrals of the PSF are performed by using a local
polynomial approximation:

\f[
	\int dy \int dx PSF(x_0+x, y_0+y) dx dy \approx
	\int dy \int dx \sum_{m=0}^{m=M} \sum_{n=0}^{n=N} f_{m,n} x^m y^n
\f]
where \f$f_{m,n}\f$ could, but need not be the taylor expansion coefficients
around \f$(x_0, y_0)\f$.

As a result, all that is necessary is calculating integrals of \f$x^my^n\f$
for non-negative integer \f$m\f$ and \f$n\f$.

PSF fitting and aperture photometry, correcting for the subpixel structure
requires calculating the integral of the PSF over rectangles and over pieces
of a circle.

Integrating rectangles
----------------------

It is necessary to calculate integrals over rectangles defined as follows:
![Aperture photometry and PSF fitting need to calculate integrals over the shaded rectangle.](images/rectangle.png)

The integral over a rectangle is obviously given by:
\f{equation}{
	I^{rect}_{m,n}\equiv
	\int_{-\Delta y}^{\Delta y} dy \int_{-\Delta x}^{\Delta x} dx x^m y^n =
	\left\{
		\begin{array}{l@{\quad}l}
			0 & \mathrm{if}\ m\ \mathrm{is\ odd\ or}\ n\ \mathrm{is\ odd}\\
			\frac{4\Delta x^{m+1} \Delta y^{n+1}}{(m+1)(n+1)} &
				\mathrm{otherwise}
		\end{array}
	\right.
\f}

Integrating circle wedges
-------------------------

![We need to calculate integrals of the PSF over the shaded piece of the circle with x and y relative to the point (x_0, y_0).](images/circle_wedge.png)

Since we are using a local polynomial approximation only the following
integral needs to be computed:

\f{equation}{
	I^{circ}_{m,n}\equiv
	\int_{-\Delta y}^{\Delta y} dy 
		\int_{-\Delta x}^{\sqrt{r^2-y_0^2-y^2-2y_0y}-x_0} dx x^m y^n
\f}

### The solution:
\f{eqnarray}{
	I^{circ}_{m,n}
		&=& \frac{1}{m+1}\Bigg\{
				\sum_{i=0}^{m+1} \binom{m+1}{i} Q_{i,n} (-x_0)^{m+1-i}
				\nonumber\\&&{}+
				(-1)^m\frac{\Delta x^{m+1} \Delta y^{n+1}}{n+1}\left[1+(-1)^n\right]
			\Bigg\}\\
	Q_{m,n} &=& (r^2-y_0^2)Q_{m-2,n} - 2y_0Q_{m-2,n+1} - Q_{m-2,n+2}\\
	Q_{0,n}=P_n^{even}
			&=&\frac{2}{n+1}\left\{
				\begin{array}{ll}
					0 &, \quad n-\mathrm{odd}\\
					\Delta y^{n+1} &, \quad n-\mathrm{even}\\
				\end{array}
			\right.\\
	Q_{1,n}=P_n^{odd}&=&
		\frac{n-1}{n+2}(r^2-y_0^2)P_{n-2}^{odd}
		- \frac{2n+1}{n+2}y_0P_{n-1}^{odd}+\nonumber\\
		&&{}+\frac{2\Delta y^{n-1}}{n+2}\times
		\left\{\begin{array}{l@{,}l}
			-x_0^3-3x_0\Delta x^2 & \quad n-\mathrm{even}\\
			3x_0^2\Delta x + \Delta x^3 & \quad n-\mathrm{odd}\\
		\end{array}\right.\\
	P_0^{odd}
		&=& \frac{1}{2}\left[
				(y_0+\Delta y)(x_0-\Delta x) -
				(y_0-\Delta y)(x_0+\Delta x)
			\right]\nonumber\\
			&&{}+
			\frac{r^2}{2}\left[
				\tan^{-1} \frac{x_0+\Delta x}{y_0-\Delta y}
				-
				\tan^{-1} \frac{x_0-\Delta x}{y_0+\Delta y}
			\right]\\
	P_1^{odd}
		&=& \frac{y_0}{2}\left[
					(x_0+\Delta x)(y_0-\Delta y)
					-
					(x_0-\Delta x)(y_0+\Delta y)
			\right]\nonumber\\
			&&{}+
			\frac{(x_0+\Delta x)^3 - (x_0-\Delta x)^3}{3}\nonumber\\
			&&{}-
			 \frac{1}{2} y_0r^2 \left(
				\tan^{-1} \frac{x_0+\Delta x}{y_0 - \Delta y}
				-
				\tan^{-1} \frac{x_0-\Delta x}{y_0 + \Delta y}
			\right)
\f}

### Derivation of the solution:
\f{eqnarray}{
	I^{circ}_{m,n}
		&=& \frac{1}{m+1}\int_{-\Delta y}^{\Delta y} dy 
			y^n\left[(\sqrt{r^2-y_0^2-y^2-2y_0y}-x_0)^{m+1} -
					(-\Delta x)^{m+1}\right]\nonumber\\
		&=& \frac{1}{m+1}\Bigg\{
				\sum_{i=0}^{m+1} \binom{m+1}{i} Q_{i,n} (-x_0)^{m+1-i}
				\nonumber\\&&{}+
				(-1)^m\frac{\Delta x^{m+1} \Delta y^{n+1}}{n+1}\left[1+(-1)^n\right]
			\Bigg\}
\f}
Where:
\f{equation}{
	Q_{m,n}\equiv \int_{-\Delta y}^{\Delta y} y^n
					(r^2-y_0^2-y^2-2y_0y)^{m/2} dy
\f}

The \f$Q_{m,n}\f$ quantities clearly obey the following recursion relation:
\f{equation}{
	Q_{m,n}=(r^2-y_0^2)Q_{m-2,n} - 2y_0Q_{m-2,n+1} - Q_{m-2,n+2}
\f}

So a general solution only requires us to be able to calculate \f$Q_{0,n}\f$
and \f$Q_{1,n}\f$.

\f$Q_{0,n}\f$ is easy:
\f{eqnarray}{
	Q_{0,n}=P_n^{even}&=&
			\int_{-\Delta y}^{\Delta y} y^n dy\\
			&=&\frac{2}{n+1}\left\{
				\begin{array}{ll}
					0 &, \quad n-\mathrm{odd}\\
					\Delta y^{n+1} &, \quad n-\mathrm{even}\\
				\end{array}
			\right.
\f}

\f{equation}{
	Q_{1,n}=P_n^{odd}=
		\int_{-\Delta y}^{\Delta y} y^n \sqrt{r^2-y_0^2-y^2-2y_0y} dy\\
\f}

Consider the following combination:
\f{eqnarray*}{
	y_0 P_{n-1}^{odd} + P_n^{odd} 
		&=& \frac{1}{2} \int y^{n-1}2(y_0+y)\sqrt{r^2-y_0^2-y^2-2y_0y} dy\\
		&=& -\frac{1}{2} \int y^{n-1}\sqrt{r^2-y_0^2-y^2-2y_0y}
												d(r^2-y_0^2-y^2-2y_0y)\\
		&=& -\frac{1}{3} \int y^{n-1}d(r^2-y_0^2-y^2-2y_0y)^{3/2}\\
		&=& -\frac{1}{3} \left[y^{n-1}\Lambda(y)\right]_{-\Delta y}^{\Delta y} 
			+\frac{n-1}{3} \int y^{n-2}(r^2-y_0^2-y^2-2y_0y)^{3/2}dy\\
		&=& -\frac{1}{3} \Delta y^{n-1}\left[\Lambda(\Delta y) +
											(-1)^n \Lambda(-\Delta y)\right]
			+\frac{n-1}{3} \left[
				(r^2-y_0^2)P_{n-2}^{odd} - P_n^{odd} - 2y_0P_{n-1}^{odd}
			\right]\\
	\Rightarrow (n+2)P_{n}^{odd} &=& 
		-\Delta y^{n-1}\left[\Lambda^3(\Delta y) +
		(-1)^n \Lambda^3(-\Delta y)\right]
		+ (n-1)(r^2-y_0^2)P_{n-2}^{odd} - (2n+1)y_0P_{n-1}^{odd}
\f}

Where \f$\Lambda(y)\equiv \sqrt{r^2-y_0^2-y^2-2y_0y}\f$.

So finally all we need is \f$P_0^{odd}\f$ and \f$P_1^{odd}\f$.
From an online integrals table:

\f[
\int \sqrt{a x^2 + b x + c}\ dx = 
\frac{b+2ax}{4a}\sqrt{ax^2+bx+c}
+
\frac{4ac-b^2}{8a^{3/2}}\ln \left| 2ax + b + 2\sqrt{a(ax^2+bx+c)}\right |
\f]

Therefore:
\f{eqnarray*}{
	P_0^{odd}
		&=& \int_{-\Delta y}^{\Delta y} \sqrt{r^2-y_0^2-y^2-2y_0y} dy\\
		&=& \left\{\frac{y_0+y}{2}\sqrt{r^2-y_0^2-y^2-2y_0y} -
				\frac{i r^2}{2} \ln 2\left| y + y_0 - 
					i\sqrt{r^2-y_0^2-y^2-2y_0y}\right|
			\right\}_{-\Delta y}^{\Delta y}\\
		&=& \left\{\frac{y_0+y}{2}\sqrt{r^2-y_0^2-y^2-2y_0y} -
				\frac{i r^2}{2}\left[\ln 2 + \ln r -
					i \tan^{-1} \frac{\sqrt{r^2-y_0^2-y^2-2y_0y}}{y_0+y}
				\right]
			\right\}_{-\Delta y}^{\Delta y}\\
		&=& \left[\frac{y_0+y}{2}\Lambda(y) -
				\frac{r^2}{2}\tan^{-1} \frac{\Lambda(y)}{y_0+y}
			\right]_{-\Delta y}^{\Delta y}\\
		&=& \frac{y_0+\Delta y}{2}\Lambda(\Delta y) -
			\frac{y_0-\Delta y}{2}\Lambda(-\Delta y)
			+\frac{r^2}{2}\tan^{-1} \frac{\Lambda(-\Delta y)}{y_0-\Delta y}
			-\frac{r^2}{2}\tan^{-1} \frac{\Lambda(\Delta y)}{y_0+\Delta y}
\f}

From the same online integral table:
\f{equation}{
\begin{split}
\int &x \sqrt{a x^2 + bx + c}\ dx = \frac{1}{48a^{5/2}}\left ( 
2 \sqrt{a} \sqrt{ax^2+bx+c}
\right .  
  \left( - 3b^2 + 2 abx + 8 a(c+ax^2) \right)
\\ &  \left.
 + 3(b^3-4abc)\ln \left|b + 2ax + 2\sqrt{a}\sqrt{ax^2+bx+c} \right| \right)
 \end{split}
\f}

Therefore (dropping constants as we go):
\f{eqnarray*}{
	P_1^{odd}
		&=& \int_{-\Delta y}^{\Delta y} y\sqrt{r^2-y_0^2-y^2-2y_0y} dy\\
		&=& -\frac{i}{48}\Bigg\{
				2i\Lambda(y) \left[ -12y_0(y+y_0) - 8\Lambda^2(y)\right]
			 	- 24y_0r^2 \ln \left| y_0 + y - i\Lambda(y)\right|
			\Bigg\}_{-\Delta y}^{\Delta y}\\
		&=& -\frac{i}{48}\Bigg\{
				2i\Lambda(y) \left[ -12y_0(y+y_0) - 8\Lambda^2(y)\right]
			 	+ 24i y_0r^2 \tan^{-1} \frac{\Lambda(y)}{y_0 + y}
			\Bigg\}_{-\Delta y}^{\Delta y}\\
		&=& \left\{
				\frac{1}{2} y_0r^2 \tan^{-1} \frac{\Lambda(y)}{y_0 + y}
				-\frac{1}{6}\Lambda(y) \left[ 3y_0(y+y_0) + 2\Lambda^2(y)\right]
			\right\}_{-\Delta y}^{\Delta y}
\f}

Finally, note that from the drawing \f$\Lambda(\Delta y)=x-\Delta x\f$ and
\f$\Lambda(-\Delta y)=x+\Delta x\f$, so:

\f{eqnarray}{
	P_n^{odd}
		&=& \frac{n-1}{n+2}(r^2-y_0^2)P_{n-2}^{even}
		- \frac{2n+1}{n+2}y_0P_{n-1}^{even}+\nonumber\\
		&&{}+\frac{2\Delta y^{n-1}}{n+2}\times
		\left\{\begin{array}{l@{,}l}
			-x_0^3-3x_0\Delta x^2 & \quad n-\mathrm{even}\\
			3x_0^2\Delta x + \Delta x^3 & \quad n-\mathrm{odd}\\
		\end{array}\right.\\
	P_0^{odd}
		&=& \frac{1}{2}\left[
				(y_0+\Delta y)(x-\Delta x) -
				(y_0-\Delta y)(x+\Delta y)
			\right]\nonumber\\
			&&{}+
			\frac{r^2}{2}\left[
				\tan^{-1} \frac{x+\Delta x}{y_0-\Delta y}
				-
				\tan^{-1} \frac{x-\Delta x}{y_0+\Delta y}
			\right]\\
	P_1^{odd}
		&=& \frac{y_0}{2}\left[
					(x_0+\Delta x)(y_0-\Delta y)
					-
					(x_0-\Delta x)(y_0+\Delta y)
			\right]\nonumber\\
			&&{}+
			\frac{(x_0+\Delta x)^3 - (x_0-\Delta x)^3}{3}\nonumber\\
			&&{}-
			 \frac{1}{2} y_0r^2 \left(
				\tan^{-1} \frac{x_0+\Delta x}{y_0 - \Delta y}
				-
				\tan^{-1} \frac{x_0-\Delta x}{y_0 + \Delta y}
			\right)
\f}

### Checking the solution:

Checking this solution is much uglier, so we will only check the terms with
\f$m+n<=2\f$, since those we already have from the integrals local polynomial
PSFs, which have been extensively debugged.

#### m=n=0
Extracting from the code:
\f{eqnarray*}{
	I_{0,0}^{circ}
			&=&	\frac{1}{2}\left[
					r^2\gamma - 2(y_0-\Delta y)\Delta x
					-
					2(x_0-\Delta x)\Delta y
				\right]\\
			&=& \frac{r^2}{2}\gamma + 2 \Delta x \Delta y - y_0 \Delta x -
				x_0 \Delta y
\f}
Where \f$\gamma\equiv\tan^{-1} \frac{x_0+\Delta x}{y_0 - \Delta y} -
\tan^{-1} \frac{x_0-\Delta x}{y_0 + \Delta y}\f$

And from the solution:
\f{eqnarray*}{
	I_{0,0}^{circ}
		&=& -x_0 Q_{0,0} + Q_{1,0} + 2 \Delta x \Delta y\\
		&=& -2x_0\Delta y
			+
			\frac{1}{2}\left[
				(y_0+\Delta y)(x_0-\Delta x)
				-
				(y_0-\Delta y)(x_0+\Delta x)
			\right]
			+
			\frac{r^2}{2}\gamma + 2\Delta x \Delta y\\
		&=& -2x_0\Delta y
			+
			\frac{1}{2}\left[
				x_0 y_0 - y_0\Delta x + x_0\Delta y - \Delta x \Delta y
				- x_0 y_0 - y_0\Delta x + x_0 \Delta y + \Delta x\Delta y
			\right]
			+
			2\Delta x \Delta y + \frac{r^2}{2}\gamma\\
		&=& -2 x_0\Delta y - y_0 \Delta x + x_0\Delta y + 2\Delta x \Delta y
			+ \frac{r^2}{2}\gamma\\
		&=& \frac{r^2}{2}\gamma + 2 \Delta x \Delta y - y_0 \Delta x
			- x_0 \Delta y
\f}

#### m=1, n=0
Extracting from the code:
\f{eqnarray*}{
	I_{1,0}^{circ}
		&=& -x_0\left[
				\frac{r^2}{2}\gamma + 2\Delta x \Delta y - y_0\Delta x
				- x_0 \Delta y
			\right]
			+
			\frac{(y_0+\Delta y)^3-(y_0-\Delta y)^3}{3}
			-
			2x_0\Delta x(y_0-\Delta y)\\
		&=& -\frac{x_0 r^2}{2}\gamma -2 x_0 \Delta x\Delta y
			+ x_0 y_0\Delta x + x_0^2\Delta y + \frac{y_0^3}{3}
			+ y_0^2\Delta y + y_0 \Delta y^2 + \frac{\Delta y^3}{3}
			- \frac{y_0^3}{3} + y_0^2\Delta y - y_0\Delta y^2
			+ \frac{\Delta y^3}{3} - 2x_0y_0\Delta x + 2x_0\Delta x\Delta y\\
		&=& -\frac{x_0r^2}{2}\gamma + x_0^2\Delta y + 2y_0^2\Delta y +
			\frac{2}{3}\Delta y^3 - x_0y_0\Delta x\\
		&=& -\frac{x_0r^2}{2}\gamma + x_0^2\Delta y + x_0y_0\Delta x +
			\frac{2}{3}\Delta y^3\\
\f}
Where the last step uses \f$x_0\Delta x=y_0\Delta y\f$.

From the solution:
\f{eqnarray*}{
	I_{1,0}^{circ}
		&=& \frac{1}{2}\left\{
				x_0^2Q_{0,0} - 2x_0Q_{1,0} + Q_{2,0}
				- 2\Delta x^2 \Delta y
			\right\}\\
		&=& \frac{1}{2}\left\{
				2x_0^2\Delta y
				-
				x_0 \left[
					(y_0+\Delta y)(x_0-\Delta x)
					-
					(y_0-\Delta y)(x_0+\Delta x)
				\right]
				- x_0r^2\gamma + Q_{2,0} - 2\Delta x^2 \Delta y
			\right\}\\
		&=& \frac{1}{2}\left\{
				2x_0^2\Delta y
				-
				x_0(x_0 y_0 - y_0 \Delta x + x_0 \Delta y - \Delta x \Delta y
					-x_0 y_0 - y_0 \Delta x + x_0 \Delta y
					+ \Delta x \Delta y)
				- x_0 r^2\gamma + 2 r^2 \Delta y - 2y_0^2\Delta y
				- \frac{2\Delta y^3}{3} - 2 \Delta x^2 \Delta y
			\right\}\\
		&=& \frac{1}{2}\left\{
				2x_0^2\Delta y + 2 x_0 y_0 \Delta x - 2 x_0^2\Delta y
				- x_0r^2\gamma + 2r^2\Delta y - 2 y_0^2\Delta y
				- \frac{2\Delta y^3}{3} -2 \Delta x^2\Delta y
			\right\}\\
		&=& x_0 y_0 \Delta x + (r^2-y_0^2)\Delta y - \frac{\Delta y^3}{3}
			- \Delta x^2 \Delta y - \frac{x_0 r^2}{2} \gamma\\
		&=& x_0 y_0 \Delta x + (x_0^2 + \Delta x^2 + \Delta y^2)\Delta y
			- \frac{\Delta y^3}{3} - \Delta x^2 \Delta y
			- \frac{x_0 r^2}{2}\gamma\\
		&=& x_0 y_0 \Delta x + x_0^2 \Delta y + \frac{2\Delta y^3}{3}
			- \frac{x_0 r^2}{2}\gamma
\f}
Where the above uses:
\f{eqnarray*}{
	Q_{2,0} &=& (r^2-y_0^2) Q_{0,0} - 2 y_0 Q_{0,1} - Q_{0,2} =
		2(r^2-y_0^2)\Delta y - \frac{2\Delta y^3}{3}\\
	r^2 &=& x_0^2 + y_0^2 + \Delta x^2 + \Delta y^2
\f}

#### m=0, n=1
Extracting from the code:
\f[
	I_{1,0}^{circ}
		=	-y_0\left[
				\frac{r^2}{2}\gamma + 2\Delta x \Delta y - y_0\Delta x
				- x_0 \Delta y
			\right]
			+
			\frac{(x_0+\Delta x)^3-(x_0-\Delta x)^3}{3}
			-
			2y_0\Delta y(x_0-\Delta x)
\f]
Which is the same as \f$I_{1,0}^{circ}\f$, but with (\f$x_0\f$, \f$y_0\f$,
\f$\Delta x\f$, \f$\Delta y\f$) substituted with (\f$y_0\f$, \f$x_0\f$,
\f$\Delta y\f$, \f$\Delta x\f$), so the final expression must be:
\f[
	I_{1,0}^{circ} = -\frac{y_0r^2}{2}\gamma + y_0^2\Delta x
					+ x_0y_0\Delta y + \frac{2}{3}\Delta x^3
\f]

From the solution:
\f{eqnarray*}{
	I_{1,0}^{circ}
		&=& -x_0 Q_{0,1} + Q_{1,1} = Q_{1,1}\\
		&=& \frac{y_0}{2}\left[
				(x_0+\Delta x)(y_0-\Delta y)
				-
				(x_0-\Delta x)(y_0+\Delta y)
			\right]
			+
			\frac{(x_0+\Delta y)^3 - (x_0-\Delta x)^3}{3}
			-
			\frac{y_0 r^2}{2}\gamma\\
		&=& \frac{y_0}{2} \left[
				x_0 y_0 - x_0\Delta y + y_0 \Delta x - \Delta x\Delta y
				-x_0 y_0 - x_0\Delta y + y_0 \Delta x + \Delta x \Delta y
			\right]
			+ \frac{x_0^3}{3} + x_0^2\Delta x + x_0\Delta x^2
			+ \frac{\Delta x^3}{3} - \frac{x_0^3}{3} + x_0^2\Delta x
			- x_0\Delta x^2 + \frac{\Delta x^3}{3}
			- \frac{y_0 r^2}{2}\gamma\\
		&=& -x_0 y_0 \Delta y + y_0^2\Delta x + 2x_0^2\Delta x 
			+ \frac{2\Delta x^3}{3} - \frac{y_0 r^2}{2}\gamma\\
		&=& y_0^2\Delta x + x_0 y_0 \Delta y 
			+ \frac{2\Delta x^3}{3} - \frac{y_0 r^2}{2}\gamma
\f}

The (m, n)=(2, 0), (1, 1) and (0, 2) expressions contain way too many terms
to handle manually, so I fed the to Mathematica which confirmed that they
reproduce the expressions in the code.

Piecewise PSF pieces
--------------------------

The implementation of piecewise PSFs makes use of integrals over more general
areas than those defined above. The task can be split into rectangle and
wedge integrals, but with x and y defined relative to an arbitrary point.

![Integrals of the PSF need to be calculated over the shaded piece of the circle.](images/hcircle_piece_diagram.png)

The integral we wish to compute is:

\f{equation}{
	I^{circ}_{m,n}\equiv
	\int_{y_{min}}^{y_{max}} dy \int_{x_{min}}^{\sqrt{r^2-(y-y_c)^2}+x_c} dx
		x^m y^n
\f}

### The solution:
\f{eqnarray*}{
	I^{circ}_{m,n}
		&=& \frac{1}{m+1}\Bigg\{
				\sum_{i=0}^{m+1} \binom{m+1}{i} Q_{i,n} x_c^{m+1-i}
				\nonumber\\&&{}-
				\frac{x_{min}^{m+1}(y_{max}^{n+1}-y_{min}^{n+1})}{n+1}
			\Bigg\}\\
	Q_{m,n}&=&(r^2-y_c^2)Q_{m-2,n} + 2y_cQ_{m-2,n+1} - Q_{m-2,n+2}\\
	Q_{0,n}=P_n^{even}&=&
			\frac{y_{max}^{n+1}-y_{min}^{n+1}}{n+1}\\
	Q_{1,n}=P_n^{odd}&=&
		\frac{1}{n+2}\Big\{
			R_n(y_{min})-R_n(y_{max})
			+
			(2n+1)y_cP_{n-1}^{odd}
			+ 
			(n-1)(r^2-y_c^2)P_{n-2}^{odd}
		\Big\}\\
	P_0^{odd}&=&\frac{1}{2} \sqrt{r^2-(y_c-y_{max})^2}(y_{max}-y_c)
			-
			\frac{1}{2} \sqrt{r^2-(y_c-y_{min})^2} (y_{min}-y_c)
			-\\
			&&
			-\frac{1}{2} r^2 \left(
				\tan^{-1}\left[
					\frac{y_c-y_{max}}{\sqrt{r^2-(y_c-y_{max})^2}}
				\right]
				-\tan^{-1}\left[
					\frac{y_c-y_{min}}{\sqrt{r^2-(y_c-y_{min})^2}}
				\right]
			\right)\\
	P_1^{odd}&=&-\frac{1}{6} \sqrt{r^2-(y_c-y_{max})^2}
		\left(2 r^2+y_c^2+y_c y_{max}-2 y_{max}^2\right)+\\
		&&+\frac{1}{6} \sqrt{r^2-(y_c-y_{min})^2}
		\left(2 r^2+y_c^2+y_c y_{min}-2 y_{min}^2\right)-\\
		&&-
		\frac{1}{2} r^2 y_c \left(
			\tan^{-1}\left[
				\frac{y_c-y_{max}}{\sqrt{r^2-(y_c-y_{max})^2}}
			\right]
			-
			\tan^{-1}\left[
				\frac{y_c-y_{min}}{\sqrt{r^2-(y_c-y_{min})^2}}
			\right]
		\right)
\f}

### Derivation of the solution:
\f{eqnarray}{
	I^{circ}_{m,n}
		&=& \frac{1}{m+1}\int_{y_{min}}^{y_{max}} dy 
			y^n\left[(\sqrt{r^2-(y-y_c)^2}+x_c)^{m+1} -
					x_{min}^{m+1}\right]\nonumber\\
		&=& \frac{1}{m+1}\Bigg\{
				\sum_{i=0}^{m+1} \binom{m+1}{i} Q_{i,n} x_c^{m+1-i}
				\nonumber\\&&{}-
				\frac{x_{min}^{m+1}(y_{max}^{n+1}-y_{min}^{n+1})}{n+1}
			\Bigg\}
\f}
Where:
\f{equation}{
	Q_{m,n}\equiv \int_{y_{min}}^{y_{max}} y^n
					\left[r^2-(y-y_c)^2\right]^{m/2} dy
\f}

The \f$Q_{m,n}\f$ quantities clearly obey the following recursion relation:
\f{equation}{
	Q_{m,n}=(r^2-y_c^2)Q_{m-2,n} + 2y_cQ_{m-2,n+1} - Q_{m-2,n+2}
\f}

Just as above:
\f$Q_{0,n}\f$ is easy:
\f{eqnarray}{
	Q_{0,n}=P_n^{even}&=&
			\int_{y_{min}}^{y_{max}} y^n dy\\
			&=&\frac{1}{n+1}\left[y_{max}^{n+1}-y_{min}^{n+1}\right]
\f}

And:
\f{equation}{
	Q_{1,n}=P_n^{odd}=
		\int_{y_{min}}^{y_{max}} y^n \sqrt{r^2-(y-y_c)^2} dy\\
\f}

\f{eqnarray*}{
	P_n^{odd} - y_cP_{n-1}^{odd}
		&=& \int_{y_{min}}^{y_{max}} (y-y_c)y^{n-1}\sqrt{r^2-(y-y_c)^2}
				d(y-y_c)\\
		&=& \frac{1}{2}\int_{y_{min}}^{y_{max}} y^{n-1}\sqrt{r^2-(y-y_c)^2}
				d(y-y_c)^2\\
		&=& -\frac{1}{3}\int_{y_{min}}^{y_{max}} y^{n-1}
				d[r^2-(y-y_c)^2]^{3/2}\\
		&=& -\frac{1}{3}\left\{
				y^{n-1}\left[r^2-(y-y_c)^2\right]^{3/2}
			\right\}_{y_{min}}^{y_{max}}
			+
			\frac{n-1}{3}\int_{y_{min}}^{y_{max}} y^{n-2}
				\left[r^2-(y-y_c)^2\right]^{3/2} dy\\
		&=& -\frac{1}{3}\left\{
				y^{n-1}\left[r^2-(y-y_c)^2\right]^{3/2}
			\right\}_{y_{min}}^{y_{max}}
			+
			\frac{n-1}{3}\left\{
				(r^2-y_c^2)P_{n-2} + 2y_cP_{n-1} - P_n
			\right\}
\f}
Letting 
\f$R_n(y)\equiv\left\{
					y^{n-1}\left[r^2-(y-y_c)^2\right]^{3/2}
				\right\}_{y_{min}}^{y_{max}}\f$:
\f[
	(n+2)P_n^{odd} = -R_n(y) + (2n+1)y_cP_{n-1}^{odd} + 
					(n-1)(r^2-y_c^2)P_{n-2}^{odd}
\f]

Using Mathematica:
\f[
	P_0^{odd}=\frac{1}{2} \sqrt{r^2-(y_c-y_{max})^2}(-y_c+y_{max})
			+
			\frac{1}{2} \sqrt{r^2-(y_c-y_min)^2} (y_c-y_{min})
			-
			\frac{1}{2} r^2 \left(
				\tan^{-1}\left[
					\frac{y_c-y_{max}}{\sqrt{r^2-(y_c-y_{max})^2}}
				\right]
				-\text{ArcTan}\left[
					\frac{y_c-y_{min}}{\sqrt{r^2-(y_c-y_{min})^2}}
				\right]
			\right)
\f]

\f[
	P_1^{odd}=-\frac{1}{6} \sqrt{r^2-(y_c-y_{max})^2}
		\left(2 r^2+y_c^2+y_c y_{max}-2 y_{max}^2\right)
		+
		\frac{1}{6} \sqrt{r^2-(y_c-y_{min})^2}
		\left(2 r^2+y_c^2+y_c y_{min}-2 y_{min}^2\right)
		-
		\frac{1}{2} r^2 y_c \left(
			\tan^{-1}\left[
				\frac{y_c-y_{max}}{\sqrt{r^2-(y_c-y_{max})^2}}
			\right]
			-
			\tan^{-1}\left[
				\frac{y_c-y_{min}}{\sqrt{r^2-(y_c-y_{min})^2}}
			\right]
		\right)
\f]

The solution was checked by implementing it in Machematica and checknig
\f$P_n^{odd}\f$, \f$Q_{m,n}\f$ and \f$I_{m,n}^{circ}\f$ againts the actual
integral expressions for a few m, and n values.

Integrating __useless__ circle pieces
-------------------------------------
I derived the integral over the following area as well due to not thinking
through what is actually needed:

![We can calculate integrals of the PSF over the shaded piece of the circle.](images/useless_circle_wedge.png)

The integral we wish to compute is:

\f{equation}{
	I^{circ}_{m,n}\equiv
	\int_{y_0}^{y_{max}} dy \int_{x_0}^{\sqrt{r^2-y^2}} dx x^m y^n
\f}
where \f$y_{max}\equiv\sqrt{r^2-x_0^2}\f$. Similarly we define
\f$x_{max}\equiv\sqrt{r^2-y_0^2}\f$.

### The solution:
\f{eqnarray}{
	I^{circ}_{m,n} &=& \frac{1}{m+1}\left(Q_{m,n} - x_0^{m+1}P^{odd}_n
									\right)\\
	Q_{2k+1,n} &=& \sum_{i=0}^{k+1} {k+1 \choose i} r^{2(k+1-i)} (-1)^i
									P^{odd}_{n+2i}\\
	P^{odd}_n &=& \frac{y_{max}^{n+1} - y_0^{n+1}}{n+1}\\
	Q_{2k,n} &=& \sum_{i=0}^{k} {k \choose i} r^{2(k-i)} (-1)^i
								P^{even}_{n+2i}\\
	P^{even}_n &=& \frac{1}{n+2}
		\left(y_0^{n-1}x_{max}^3 - y_{max}^{n-1}x_0^3\right)
		+ \frac{n-1}{n+2}r^2P^{even}_{n-2}\\
	P^{even}_0 &=& \frac{1}{2}\left[y_{max}x_0 - y_0x_{max} +
				r^2\left(\tan^{-1}\frac{y_{max}}{x_0} -
				\tan^{-1}\frac{y_0}{x_{max}}\right)
			\right]\\
	P^{even}_1 &=& \frac{1}{3}\left(x_{max}^3 - x_0^3\right)
\f}

The \f$Q_{m,n}\f$ quantities satisfy the following identity:
\f[
	Q_{m+2,n}=r^2 Q_{m,n}-Q_{m,n+2}
\f]
See below for proof.

This allows the same solution to be written in a different form:
\f{eqnarray}{
	I^{circ}_{m,n} &=& \frac{1}{m+1}\left(Q_{m,n} - x_0^{m+1}P^{odd}_{n}
									\right)\\
	Q_{m,n} &=& r^2 Q_{m-2,n}-Q_{m-2,n+2}\\
	Q_{1,n} &=& r^2 P^{odd}_{n} - P^{odd}_{n+2}\\
	Q_{0,n} &=& P^{even}_n\\
\f}
Where \f$P^{even}_n\f$ and \f$P^{odd}_n\f$ are the same quantities as above.

### Derivation of the solution:

\f{eqnarray}{
	I^{circ}_{m,n}&=&\frac{1}{m+1}\int_{y_0}^{y_{max}} dy y^n \left[
		(r^2-y^2)^{(m+1)/2} - x_0^{m+1}\right]\\
	&=& \frac{1}{m+1}\left[Q_{m,n} -
		\frac{x_0^{m+1}}{n+1} \left(y_{max}^{n+1} - y_0^{n+1}\right)\right]
\f}
Where:
\f{equation}{
	Q_{m,n}\equiv\int_{y_0}^{y_{max}} dy y^n (r^2-y^2)^{(m+1)/2}
\f}

At this point we need to consider two separate cases.

#### Case 1: \f$m=2k+1\f$ with \f$k \in \mathbb{Z}^+\f$

\f{eqnarray}{
	Q_{2k+1,n}
		&=&\int_{y_0}^{y_{max}} y^n (r^2-y^2)^{k+1} dy\\
		&=&\int_{y_0}^{y_{max}} y^n \sum_{i=0}^{k+1} r^{2(k+1-i)}
			(-1)^i y^{2i} dy
\f}
Which gives:
\f{equation}{
	Q_{2k+1,n}=\sum_{i=0}^{k+1} {k+1 \choose i} r^{2(k+1-i)} (-1)^i
		\frac{y_{max}^{n+2i+1} - y_0^{n+2i+1}}{n+2i+1}
\f}

#### Case 2: \f$m=2k\f$ with \f$k \in \mathbb{Z}^+\f$
\f{eqnarray}{
	Q_{2k,n}
		&=&\int_{y_0}^{y_{max}} y^n (r^2-y^2)^k \sqrt{r^2-y^2} dy\\
		&=&\int_{y_0}^{y_{max}} y^n \sum_{i=0}^{k} {k \choose i} r^{2(k-i)}
			y^{2i}(-1)^i \sqrt{r^2-y^2}\\
		&=&\sum_{i=0}^{k} {k \choose i} r^{2(k-i)} (-1)^i P^{even}_{n+2i}
\f}
With:
\f{equation}{
	P^{even}_n\equiv\int_{y_0}^{y_{max}} y^n\sqrt{r^2-y^2} dy
\f}

\f{eqnarray}{
	P^{even}_n	&=&\frac{1}{2} \int_{y_0}^{y_{max}} y^{n-1} \sqrt{r^2-y^2}
													dy^2\\
		&=&-\frac{1}{3} \int_{y_0}^{y_{max}} y^{n-1} d(r^2-y^2)^{3/2}\\
		&=&\frac{1}{3}\bigg\{ y_0^{n-1}x_{max}^3 - y_{max}^{n-1} x_0^3 +
				\nonumber\\
		&&\quad\quad
			{}+(n-1)\int_{y_0}^{y_{max}} y^{n-2}(r^2-y^2)\sqrt{r^2-y^2}dy
			\bigg\}\\
		&=&\frac{1}{3}\left\{ y_0^{n-1} x_{max}^3 - y_max^{n-1}x_0^3 +
				(n-1)r^2P^{even}_{n-2} - (n-1)P^{even}_n
			\right\}
\f}

So we end up with:
\f{equation}{
	P^{even}_n=\frac{1}{n+2} \left(y_0^{n-1}x_{max}^3 -
									y_{max}^{n-1}x_0^3\right) +
		\frac{n-1}{n+2}r^2P^{even}_{n-2}
\f}

So in order to be able to calculate any \f$P^{even}_n\f$ we need
\f$P^{even}_0\f$ and \f$P^{even}_1\f$:
\f{eqnarray}{
	P^{even}_0 &=& \int_{y_0}^{y_{max}} \sqrt{r^2-y^2} dy \\
		&=& \frac{1}{2}\left[y_{max}x_0 - y_0x_{max} +
				r^2\left(\tan^{-1}\frac{y_{max}}{x_0} -
				\tan^{-1}\frac{y_0}{x_{max}}\right)
			\right]
\f}

\f{eqnarray}{
	P^{even}_1 &=& \int_{y_0}^{y_{max}} y\sqrt{r^2-y^2} dy\\
		&=& -\frac{1}{2}\int_{y_0}^{y_{max}}\sqrt{r^2-y^2} d(r^2-y^2)\\
		&=& \frac{1}{3}\left(x_{max}^3 - x_0^3\right)
\f}

### Checking the solution:

The best way to test if I made any errors in deriving the above is to
differentiate the result and see if we get what we expect.

In particular we expect:
\f[
	\frac{\partial I^{circ}_{m,n}}{\partial x_0} = \int_{y_0}^{y_{max}} dy
		x_0^m y^n = -x_0^m\frac{y_{max}^{n+1} - y_0^{n+1}}{n+1}
\f]
And
\f[
	\frac{\partial I^{circ}_{m,n}}{\partial y_0} =
		-y_0^n\frac{x_{max}^{m+1} - x_0^{m+1}}{m+1}
\f]

Two useful identities:
\f[
	\frac{\partial y_{max}}{\partial x_0}=-\frac{x_0}{y_max}
	\quad\mathrm{and}\quad
	\frac{\partial x_{max}}{\partial y_0}=-\frac{y_0}{x_max}
\f]

#### The derivative with respect to \f$x_0\f$:

From our solution:
\f{eqnarray*}{
	\frac{\partial I^{circ}_{m,n}}{\partial x_0} &=&
		\frac{\partial}{\partial x_0}\left\{
			\frac{1}{m+1}\int_{y_0}^{y_{max}} dy y^n \left[
				Q_{m,n} -
				\frac{x_0^{m+1}}{n+1} \left(y_{max}^{n+1} - y_0^{n+1}\right)
			\right]
		\right\}\\
	&=& \frac{1}{m+1}\frac{\partial Q_{m,n}}{\partial x_0} -
		x_0^m\frac{y_{max}^{n+1}-y_0^{n+1}}{n+1} +
		\frac{x_0^{m+2}y_{max}^{n-1}}{m+1}
\f}

Since the second term is exactly what we want, we need the following to be
satisfied:
\f[
	\frac{\partial Q_{m,n}}{\partial x_0} = -x_0^{m+2}y_{max}^{n-1}
\f]

From the solution:
\f{eqnarray*}{
	\frac{\partial Q_{2k+1,n}}{\partial x_0}
		&=& -\sum_{i=0}^{k+1} {k+1 \choose i} r^{2(k+1-i)} (-1)^i
			x_0 y_{max}^{n+2i-1}\\
		&=& -x_0 y_{max}^{n-1}
			\sum_{i=0}^{k+1} {k+1 \choose i} (r^2)^{k+1-i} (-y_{max}^2)^i\\
		&=& -x_0 y_{max}^{n-1}\left(r^2-y_{max}^2\right)^{k+1}\\
		&=& -x_0^{2k+3} y_{max}^{n-1}
\f}
Which is exactly what we needed.

Also from the solution:
\f{eqnarray*}{
	\frac{\partial Q_{2k,n}}{\partial x_0}
		&=& \sum_{i=0}^{k} {k \choose i} r^{2(k-i)} (-1)^i 
			\frac{\partial P^{even}_{n+2i}}{\partial x_0}\\
\f}

So we need:
\f{eqnarray*}{
	\sum_{i=0}^{k} {k \choose i} r^{2(k-i)} (-1)^i 
	\frac{\partial P^{even}_{n+2i}}{\partial x_0} &=&
		-x_0^{2k+2}y_{max}^{n-1}\\
	&=& -x_0^2 y_{max}^{n-1}
		\sum_{i=0}^{i=k} {k \choose i} r^{2(k-i)} (-1)^i y_{max}^{2i}
\f}
which translates to:
\f[
	\frac{\partial P^{even}_{n+2i}}{\partial x_0} = -x_0^2 y_{max}^{n+2i-1}
\f]
or
\f[
	\frac{\partial P^{even}_n}{\partial x_0} = -x_0^2 y_{max}^{n-1}
\f]
First:
\f{eqnarray*}{
	\frac{\partial P^{even}_0}{\partial x_0}
		&=& \frac{1}{2}\left\{ y_{max} - \frac{x_0^2}{y_{max}} -
				r^2 \frac{1}{1+\frac{y_{max}^2}{x_0^2}} \left[
					\frac{1}{y_{max}} + \frac{y_{max}}{x_0^2}
				\right]
			\right\}\\
		&=& \frac{1}{2y_{max}}\left\{ y_{max}^2 - x_0^2 -
				r^2 \frac{x_0^2}{x_0^2+y_{max}} \left[
					\frac{1} + \frac{y_{max}^2}{x_0^2}
				\right]
			\right\}\\
		&=& \frac{1}{2y_{max}}\left( y_{max}^2 - x_0^2 - r^2\right)\\
		&=& -\frac{x_0^2}{y_{max}}
\f}
which is exactly what we needed.

Next:
\f[
	\frac{\partial P^{even}_1}{\partial x_0}=-x_0^2
\f]
which is also what we needed.

Finally:
\f{eqnarray*}{
	\frac{\partial P^{even}_n}{\partial x_0}
		&=& \frac{1}{n+2}
			\left((n-1)y_{max}^{n-3}x_0^4 - 3y_{max}^{n-1}x_0^2\right)
			+ \frac{n-1}{n+2}r^2
				\frac{\partial P^{even}_{n-2}}{\partial x_0}\\
		&=& \frac{y_{max}^{n-3}x_0^2}{n+2}
			\left((n-1)x_0^2 - 3y_{max}^2\right)
			- \frac{n-1}{n+2}r^2 y_{max}^{n-3}x_0^2\\
		&=& \frac{y_{max}^{n-3}x_0^2}{n+2}
			\left((n-1)x_0^2 - 3y_{max}^2 - (n-1)r^2\right)\\
		&=& -y_{max}^{n-1}x_0^2
\f}
Which is the last piece we needed for the derivative with respect to
\f$x_0\f$.

#### The derivative with respect to \f$y_0\f$:

From our solution:
\f[
	\frac{\partial I^{circ}_{m,n}}{\partial y_0} = \frac{1}{m+1}\left(
		\frac{\partial Q_{m,n}}{\partial y_0} + x_0^{m+1} y_0^n \right)
\f]

So we need:
\f[
	\frac{\partial Q_{m,n}}{\partial y_0} = -x_{max}^{m+1} y_0^n
\f]

From the solution:
\f{eqnarray*}{
	\frac{\partial Q_{2k+1,n}}{\partial y_0}
		&=& -\sum_{i=0}^{k+1} {k+1 \choose i} r^{2(k+1-i)} (-1)^i
			y_0^{n+2i}\\
		&=& -y_0^n (r^2-y_0^2)^{k+1}\\
		&=& -x_max^{2k+2} y_0^n
\f}
Exactly what we need.

Also from the solution:
\f[
	\frac{\partial Q_{2k,n}}{\partial y_0} =
		\sum_{i=0}^{k} {k \choose i} r^{2(k-i)} (-1)^i
		\frac{\partial P^{even}_{n+2i}}{\partial y_0}
\f]

So we need:
\f{eqnarray*}{
	\sum_{i=0}^{k} {k \choose i} r^{2(k-i)} (-1)^i
			\frac{\partial P^{even}_{n+2i}}{\partial y_0}\\
		&=& -x_{max}^{2k+1}y_0^n\\
		&=& -x_{max}\sum_{i=0}^{k} {k \choose i} r^{2(k-i)} (-1)^i y_0^{n+2i}
\f}

In other words we must show that:
\f[
	\frac{\partial P^{even}_n}{\partial y_0}=-x_{max} y_0^n
\f]

First:
\f{eqnarray*}{
	\frac{\partial P^{even}_0}{\partial y_0}
		&=& \frac{1}{2}\left[\frac{y_0^2}{x_{max}} - x_{max} 
				-r^2 \frac{1}{1+\frac{y_0^2}{x_{max}^2}} \left(
					\frac{1}{x_{max}} + \frac{y_0^2}{x_{max}^3}
				\right)
			\right]\\
		&=& \frac{1}{2x_{max}}\left[ y_0^2 - x_{max}^2 -
				r^2\frac{x_{max}^2}{x_{max}^2+y_0^2}
				\frac{x_{max}^2 + y_0^2}{x_{max}^2}
			\right]\\
		&=& \frac{1}{2x_{max}}\left[ y_0^2 - x_{max}^2 - r^2 \right]\\
		&=& -x_{max}
\f}
Which matches our expectation.

Next:
\f[
	\frac{\partial P^{even}_1}{\partial y_0} = - x_{max} y_0
\f]
Which also matches.

Finally, assuming \f$\frac{\partial P^{even}_{n-2}}{\partial y_0} = -x_{max}
y_0^{n-2}\f$ :
\f{eqnarray*}{
	\frac{\partial P^{even}_n}{\partial y_0}
		&=& \frac{1}{n+2}
			\left((n-1) y_0^{n-2}x_{max}^3 - 3 y_0^n x_{max}\right)
			+ \frac{n-1}{n+2}r^2
				\frac{\partial P^{even}_{n-2}}{\partial y_0}\\
		&=&  \frac{1}{n+2}
			\left((n-1) y_0^{n-2}x_{max}^3 - 3 y_0^n x_{max}\right)
			- \frac{n-1}{n+2}r^2 x_{max} y_0^{n-2}\\
		&=& \frac{x_{max}}{n+2}\left((n-1) y_0^{n-2}x_{max}^2 - 3 y_0^n
				- (n-1)r^2 y_0^{n-2}
			\right)\\
		&=& -x_{max} y_0^n
\f}
Which was the last missing piece.

### Deriving the relation between \f$Q_{m,n}\f$ quantities.

For odd \f$m\f$:
\f{eqnarray*}{
	Q_{2k+3,n} &=& \sum_{i=0}^{k+2} {k+2 \choose i} r^{2(k+2-i)} (-1)^i
									P^{odd}_{n+2i}\\
			   &=& \sum_{i=0}^{k+2}
					\left[{k+1 \choose i} + {k+1 \choose i-1}\right]
					r^{2(k+2-i)} (-1)^i P^{odd}_{n+2i}\\
			   &=& r^2 Q_{2k+1,n} + \sum_{i=1}^{k+2} {k+1 \choose i-1}
					r^{2(k+2-i)} (-1)^i P^{odd}_{n+2i}\\
			   &=& r^2 Q_{2k+1,n} - \sum_{i=0}^{k+1} {k+1 \choose i}
					r^{2(k+1-i)} (-1)^i P^{odd}_{n+2+2i}\\
			   &=& r^2 Q_{2k+1,n} - Q_{2k+1,n+2}
\f}

For even \f$m\f$:
\f{eqnarray*}{
	Q_{2k+2,n} &=& \sum_{i=0}^{k+1} {k+1 \choose i} r^{2(k+1-i)} (-1)^i
								P^{even}_{n+2i}\\
			   &=& \sum_{i=0}^{k+1}
					\left[{k \choose i} + {k \choose i-1}\right]
					r^{2(k+1-i)} (-1)^i P^{even}_{n+2i}\\
			   &=& r^2Q_{2k,n} + \sum_{i=1}^{k+1} {k \choose i-1}
					r^{2(k+1-i)} (-1)^i P^{even}_{n+2i}\\
			   &=& r^2Q_{2k,n} - \sum_{i=0}^{k} {k \choose i}
					r^{2(k-i)} (-1)^i P^{even}_{n+2+2i}\\
			   &=& r^2Q_{2k,n} - Q_{2k,n+2}
\f}

Actually, this relation is obvious from the definition of \f$Q_{m,n}\f$.
