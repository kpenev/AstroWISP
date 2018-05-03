Implementing elliptical Gaussian integrals {#PSF_integrals_implementation_page}
==========================================

Calculating integrals:
-------------------------

We wish to calculate (up to some specified precision) integrals of elliptical
an Gaussian PSF (and their first and second order derivatives with respect to
S, D and K) over a rectangle specified as in the following diagram:

![Aperture photometry and PSF fitting need to calculate integrals over the shaded rectangle.](images/rectangle.png)

or over a piece of a circle like this:

![Aperture photometry and PSF fitting need to calculate integrals over the shaded rectangle.](images/rectangle.png)

Let us denote the integral with \f$\mathcal{I}\f$, and its derivatives with
subscripts.

For example:
\f{eqnarray*}{
	\mathcal{I}_S &\equiv& \frac{\partial \mathcal{I}}{\partial S}\\
	\mathcal{I}_{SK} &\equiv&
		\frac{\partial^2 \mathcal{I}}{\partial S\partial K}\\
\f}

Note that the current implementation never requires derivative information
for integrals over circle wedges, so even though those could be calculated,
the present implementation does not. If at any time in the future this is
required, it can be done in a fairly straightforward fashion.

The PSF model is:
\f{equation}{
	Psf(x,y)\equiv
	\exp\left\{-\frac{1}{2}\left[S(x^2+y^2)+D(x^2-y^2)+2Kxy\right]\right\}
\f}

We use the following quantities defined in <++>
\f{eqnarray}{
	C_{20}&\equiv&S+D\\
	C_{11}&\equiv&2K\\
	C_{02}&\equiv&S-D\\
	C_{10}&\equiv&2\left[(S+D)x+Ky\right]\\
	C_{01}&\equiv&2\left[(S-D)y+Kx\right]
\f}
And approximate the PSF as
\f{equation}{
	\frac{Psf(x+\delta x, y+\delta y)}{Psf(x,y)}=
		S_{20}S_{11}S_{02}S_{10}S_{01}
\f}
Where (as defined in <++>):
\f{eqnarray}{
	S_{20} & \equiv & \sum_{i=0}^I
		\frac{(-1)^i C_{20}^i \delta x^{2i}}{2^i i!}\\
	S_{11} & \equiv & \sum_{j=0}^J
		\frac{(-1)^j C_{11}^j \delta x^j \delta y^j}{2^j j!}\\
	S_{02} & \equiv & \sum_{k=0}^K
		\frac{(-1)^k C_{02}^k \delta y^{2k}}{2^k k!}\\
	S_{10} & \equiv & \sum_{l=0}^L
		\frac{(-1)^l C_{10}^l \delta x^l}{2^l l!}\\
	S_{01} & \equiv & \sum_{m=0}^M \frac{(-1)^m C_{01}^m \delta y^m}{2^m m!}
\f}

An upper limit to the error in the integral is:
\f{eqnarray}{
	\delta \mathcal{I} &<& Psf(x_0, y_0)A \big[\\
		&&(S_{20}+\Delta_{20}) (S_{11}+\Delta_{11}) (S_{02}+\Delta_{02})
			(S_{10}+\Delta_{10}) (S_{01}+\Delta_{01}) \\
		&&{} - S_{20}S_{11}S_{02}S_{10}S_{01}\big]
\f}

With:
\f{eqnarray}{
	\Delta_{20}(I) &\equiv& \frac{C_{20}^{I+1}\delta x^{2I+2}}{2^{I+1}(I+1)!}
		\exp\left(-\frac{C_{20}\delta x^2}{2}\right)\\
	\Delta_{11}(J) &\equiv& \left|\frac{C_{11}\delta x\delta y}{2}\right|^{J+1}
		\frac{1}{(J+1)!} \exp\left|\frac{C_{11}\delta x\delta y}{2}\right|\\
	\Delta_{02}(K) &\equiv& \frac{C_{02}^{K+1}\delta y^{2K+2}}{2^{K+1}(K+1)!}
		\exp\left(-\frac{C_{02}\delta y^2}{2}\right)\\
	\Delta_{10}(L) &\equiv& \left|\frac{C_{10}\delta x}{2}\right|^{L+1}
		\frac{1}{(L+1)!} \exp\left|\frac{C_{10}\delta x}{2}\right|\\
	\Delta_{01}(M) &\equiv& \left|\frac{C_{01}\delta y}{2}\right|^{M+1}
		\frac{1}{(M+1)!} \exp\left|\frac{C_{01}\delta y}{2}\right|
\f}

and where \f$A\f$ is the area over which the integral is being calculated:
\f{eqnarray}{
	A &=& 4\Delta x \Delta y \quad \mathrm{for\ a\ rectangle}\\
	A &=& P_0 - x_0(y_{max}-y_0) \quad \mathrm{for\ a\ circrle\ wedge}
\f}
See the [analytic expressions](@ref PSF_integrals_page) page for a
definition of \f$P_0\f$ and \f$y_{max}\f$.

At the start
------------

### For integrals over rectangles:

\f{eqnarray*}{
	\mathcal{I} &=& 4 Psf(x_0, y_0) \Delta x \Delta y\\
	\left.\begin{array}{l}
		\mathcal{I}_S,\mathcal{I}_D,\mathcal{I}_K,\mathcal{I}_{SS},
		\mathcal{I}_{SD},\\
		\mathcal{I}_{SK},\mathcal{I}_{DD},\mathcal{I}_{DK},\mathcal{I}_{KK}
		\end{array}\right\} &=& 0\\
	L_{2,0}=L_{1,1}=L_{0,2}=L_{1,0}=L_{0,1} &=& 1\\
	c_{2,0}&=&\frac{S+D}{2}\\
	c_{1,1}&=&K \\
	c_{0,2}&=&\frac{S-D}{2} \\
	c_{1,0}&=& \left((S+D) x_0 + K y_0\right)\\
	c_{0,1}&=& \left((S-D) y_0 + K x_0\right)\\
	f_{i,j}&=& c_{ij} \Delta x^i \Delta y^j\ ,\quad 
		(i,j) \in \left\{(2,0), (1,1), (0,2), (1,0), (0,1)\right\}\\
	\mathcal{O}_{i,j}&=&0\ ,\quad
		(i,j) \in \left\{(2,0), (1,1), (0,2), (1,0), (0,1)\right\}\\
	S_{i,j}&=&1\ ,\quad
		(i,j) \in \left\{(2,0), (1,1), (0,2), (1,0), (0,1)\right\}\\
	\Delta_{i,j}&=&\left\{
			\begin{array}{l@{\ ,\quad}l}
				|f_{ij}|\exp(-f_{ij}) & (i,j)\in\left\{(2,0), (0,2)\right\}\\
				|f_{ij}|\exp(|f_{ij}|) & i,j) \in
					\left\{(1,1), (1,0), (0,1)\right\}
			\end{array}
		\right.
\f}

\f$L_{i,j}\f$ refer to the last terms in the \f$S_{i,j}\f$ quantities,
\f$c_{i,j}\f$ and \f$f_{i,j}\f$ are just pre-stored shortcuts,
\f$\mathcal{O}_{i,j}\f$ are the expansion orders: \f$\mathcal{O}_{2,0}=I\f$,
\f$\mathcal{O}_{1,1}=J\f$, \f$\mathcal{O}_{0,2}=K\f$,
\f$\mathcal{O}_{1,0}=L\f$, \f$\mathcal{O}_{0,1}=M\f$, and \f$\Delta_{i,j}\f$
are the current limits on the error terms as given above.

### For integrals over circle wedges:

First, note that we start only requiring \f$\mathcal{I}_{0,0}\f$ and hence
\f$Q_{0,0}\f$. And at any step, we increase one of the \f$I\f$, \f$J\f$,
\f$K\f$, \f$L\f$ or \f$M\f$ indices by one, so for any new term we request,
we will have already requested all its prerequisites according to the
recursion relation outlined [here](@ref PSF_integrals_page), as long as the
summation is done with \f$M\f$ nested within \f$K\f$, nested within \f$L\f$,
nested within \f$J\f$ nested within \f$I\f$.


We store previously calculated \f$Q_{m,n}\f$ values in a pre-allocated array
(which can be expanded if higher order terms are necessary). In addition we
maintain two arrays (\f$\mathcal{N^odd}\f$ and \f$\mathcal{N^even}\f$) which
keeps track of the \f$Q_{m,n}\f$ terms calculated so far. In particular
\f$\mathcal{N}^{odd}_{2m+1+n}\f$ is the largest \f$m\f$ for which
\f$Q_{2m+1,n}\f$ has been calculated and \f$\mathcal{N}^{even}_{2m+n}\f$ is
the largest \f$m\f$ for which \f$Q_{2m,n}\f$ has been calculated.

The calculation of a new \f$Q_{m,n}\f$ value proceeds as follows:
 1. We make sure all \f$Q_{m\%1,(m+n)\%2+2*i}\f$ are calculated for
	\f$i<=(m+n)/2\f$.

 2. Use the recursion relation between \f$Q_{m,n}\f$ values to fill in the
    triangle of values which contribute to \f$Q_{m,n}\f$, that is
	\f$Q_{i,j}\f$ for which \f$i\%2=m\%2\f$, \f$j\%2=n\%2\f$,
	\f$(i+j)<=(m+n)\f$ and \f$j\geq n\f$ in orderof increasing \f$i+j\f$ and
	in order of increasing \f$i\f$ inside that.

\f{eqnarray*}{
	x_{max} &=& \sqrt{r^2-y_0^2}\\
	y_{max} &=& \sqrt{r^2-x_0^2}\\
	\Delta x &=& \frac{x_{max}-x_0}{2}\\
	\Delta y &=& \frac{y_{max}-y_0}{2}\\
	P_0 &=& \frac{1}{2}\left[y_{max}x_0 - y_0x_{max} +
				r^2\left(\tan^{-1}\frac{y_{max}}{x_0} -
				\tan^{-1}\frac{y_0}{x_{max}}\right)
			\right]\\
	P_1 &=& \frac{1}{3}\left(x_{max}^3-x_0^3\right)\\
	Q_{0,0} &=& Psf(x_0, y_0)\left[P_0 - x_0(y_{max}-y_0)\right]\\
	Q_{0,1} &=& \\
	Q_{1,0} &=& \\
	Q_{1,1} &=& \\
\f}
and the remaining variables are exactly the same as for rectangle integrals,
except for \f$Q_S\f$, \f$Q_D\f$, \f$Q_K\f$, \f$Q_{SS}\f$, \f$Q_{SD}\f$,
\f$Q_{SK}\f$, \f$Q_{DD}\f$, \f$Q_{DK}\f$ and \f$Q_{KK}\f$, which are
undefined.

The calculation
---------------

The first thing that is done is to ensure that \f$\Delta_{20}(I)\f$ and
\f$\Delta_{02}(J)\f$ are monotonically decreasing functions of the
corresponding expansion orders. That means increasing \f$I\f$ and \f$K\f$
until \f$I>f_{20}/2-1\f$ and \f$K>f_{02}/2-1\f$. This is done in order to
ensure that we do not stop prematurely refining the value of
\f$\mathcal{I}\f$. This increasing of expansion order is detailed below.

After that, the algorithm for estimating \f$\mathcal{I}\f$ is as follows:
 1. If the upper limit on the error is less than the tolerance: stop
 2. Pick the index (\f$(i,j)\f$) (one of
	\f$\{(2,0), (1,1), (0,2), (1,0), (0,1)\}\f$) to increment for which
	\f[
		\left|\Delta_{i,j}\right|\prod_{(k,l)\neq(i,j)} \left|S_{k,l}\right|
	\f]
	is the largest.
 3. Go to one order higher in the expansion in the corresponding index as
	outlined below, and go back to step 1.

Incrementing the expansion order of the \f$(p,q)\f$-th term 
-----------------------------------------------------------

\f{eqnarray*}{
	\mathcal{O}_{p,q} &=& \mathcal{O}_{p,q}+1\\
	L_{p,q} &=& -L_{p,q}\frac{f_{p,q}}{\mathcal{O}_{p,q}}
\f}

This way \f$L_{p,q}\f$ are indeed the last terms in the \f$S_{p,q}\f$
quantities.

Next, if \f$(p,q)\f$ are not (2,0) or (0,2): \f$S_{p,q}=S_{p,q}+L_{p,q}\f$,
and
\f$\Delta_{p,q}=\Delta_{p,q}\left|f_{p,q}\right|/(\mathcal{O}_{p,q}+1)\f$.
The reason why \f$S_{2,0}\f$ and \f$S_{0,2}\f$ are not updated is because for
high enough order \f$S_{2,0}\leq 1\f$ and \f$S_{0,2}\leq 1\f$, but due to the
alternating signs in their expansion they oscillate, so the limits on the
error will sometimes grow even though we have increased the order of the
expansion. So it is better to leave them to their initial value of 1.

Then for rectangles the \f$\mathcal{I}\f$ quantities are updated according
to:
\f[
	\mathcal{I}_*=\mathcal{I}_*+4 Psf(x_0, y_0) \Delta x \Delta y
		\sum^{\mathcal{O}_{2,0}}_{i=i_0} (-1)^i\frac{f_{2,0}^i}{i!}
		\sum^{\mathcal{O}_{1,1}}_{j=j_0} (-1)^j\frac{f_{1,1}^j}{j!}
		\sum^{\mathcal{O}_{0,2}}_{k=k_0} (-1)^k\frac{f_{0,2}^k}{k!}
		\sum^{\mathcal{O}_{1,0}}_
			{\begin{array}{c}l=l_0\\j+l:even\end{array}}
					(-1)^l\frac{f_{1,0}^l}{l!(2i+j+l+1)}
		\sum^{\mathcal{O}_{0,1}}_
			{\begin{array}{c}m=m_0\\j+m:even\end{array}}
					(-1)^m\frac{f_{0,1}^m}{m!(2k+j+m+1)}\Lambda_*
\f]

where:
\f{eqnarray*}{
	i_0 &=& \left\{\begin{array}{l@{\ ,\ }l}
							0 & (p,q)\neq(2,0)\\
							\mathcal{O}_{2,0} & (p,q)=(2,0)
				\end{array}\right.\\
	j_0 &=& \left\{\begin{array}{l@{\ ,\ }l}
							0 & (p,q)\neq(1,1)\\
							\mathcal{O}_{1,1} & (p,q)=(1,1)
				\end{array}\right.\\
	k_0 &=& \left\{\begin{array}{l@{\ ,\ }l}
							0 & (p,q)\neq(1,1)\\
							\mathcal{O}_{0,2} & (p,q)=(0,2)
				\end{array}\right.\\
	l_0 &=& \left\{\begin{array}{l@{\ ,\ }l}
					j\%2 & (p,q)\neq(1,0)\\
					\mathcal{O}_{1,0}+(\mathcal{O}_{1,0}+j)\%2 & (p,q)=(1,0)
			\end{array}\right.\\
	m_0 &=& \left\{\begin{array}{l@{\ ,\ }l}
					j\%2 & (p,q)\neq(0,1)\\
					\mathcal{O}_{0,1}+(\mathcal{O}_{0,1}+j)\%2 & (p,q)=(0,1)
			\end{array}\right.\\
	\Lambda &=& 1\\
	\Lambda_S &=& \frac{i}{2 c_{2,0}} + \frac{k}{2 c_{0,2}} +
					\frac{x_0 l}{c_{1,0}} + \frac{y_0 m}{c_{0,1}}\\
	\Lambda_D &=& \frac{i}{2 c_{2,0}} - \frac{k}{2 c_{0,2}} +
					\frac{x_0 l}{c_{1,0}} - \frac{y_0 m}{c_{0,1}}\\
	\Lambda_K &=& \frac{j}{c_{1,1}} + \frac{y_0 l}{c_{1,0}} +
					\frac{x_0 m}{c_{0,1}}\\
	\Lambda_{SS} &=& \Lambda_S^2 - \left(
		\frac{i}{4 c_{2,0}^2} + \frac{k}{4 c_{0,2}^2} +
		\frac{x_0^2 l}{c_{1,0}^2} + \frac{y_0^2 m}{c_{0,1}^2}\right)\\
	\Lambda_{SD} &=& \Lambda_S*\Lambda_D - \frac{i}{4 c_{2,0}^2} +
		\frac{k}{4 c_{0,2}^2} -\frac{x_0^2 l}{c_{1,0}^2} +
		\frac{y_0^2 m}{c_{0,1}^2}\\
	\Lambda_{SK} &=& \Lambda_S \Lambda_K - x_0 y_0\left(
			\frac{l}{c_{1,0}^2} + \frac{m}{c_{0,1}^2}\right)\\
	\Lambda_{DD} &=& \Lambda_D^2 - \left(
		\frac{i}{4 c_{2,0}^2} + \frac{k}{4 c_{0,2}^2} +
		\frac{x_0^2 l}{c_{1,0}^2} + \frac{y_0^2 m}{c_{0,1}^2}\right)\\
	\Lambda_{DK} &=& \Lambda_D \Lambda_K - x_0 y_0\left(
			\frac{l}{c_{1,0}^2} - \frac{m}{c_{0,1}^2}\right)\\
	\Lambda_{KK} &=& \Lambda_K^2 -\frac{x_0^2 m}{c_{0,1}^2} -
		\frac{y_0 l}{c_{1,0}} - \frac{j}{c_{1,1}^2}
\f}
