struct Amoeba {
	/*Multidimensional minimization by the downhill simplex methond of Nelder and Mead
	The concept of the method is written in Numerical Recipes 3rd Ed, Ch 10.5, p502
	The code is in C++, could be called in Matlab by matlab 'MEX' function using the
	include file 'nr3matlab.h' which is also included */
	const Doub ftol;  //tolerence of decrease in function value in terminating step
	Int nfunc;  //number of function evalutions
	Int mpts;  //number of points
	Int ndim;  //number of dimensions
	Doub fmin;  //function value at the minimum
	VecDoub y;  //function values at the vertices of the simplex
	MatDoub p;  //Current simplex
	Amoeba(const Doub ftoll) : ftol(ftoll) {}
	/*The constructor argument ftoll is the fractional convergence tolerance to be 
	achieved in the function value (n.b.!). */
	template <class T>
	VecDoub minimize(VecDoub_I &point, const Doub del, T &func)
	/*Multidimensional minimization of the function or functor func(x), where 
	x[0..ndim-1] is a vector in ndim dimensions, by the downhill simplex method of
	Nelder and Mead. The initial simplex is specified as in equation (10.5.1) by a 
	point[0..ndim-1] and a constant displacement del along each coordinate direction.
	Returned is the location of the minimum.*/
	{
		VecDoub dels(point.size(),del);
		return minimize(point,dels,func);
	}
	template <class T>
	VecDoub minimize(VecDoub_I &point, VecDoub_I &dels, T &func)
	/*Alternative interface that takes different displacements dels[0..ndim-1] in 
	different directions for the initial simplex.*/
	{
		Int ndim=point.size();
		MatDoub pp(ndim+1,ndim);
		for (Int i=0;i<ndim+1;i++) {
			for (Int j=0;j<ndim;j++)
				pp[i][j]=point[j];
			if (i !=0 ) pp[i][i-1] += dels[i-1];
		}
		return minimize(pp,func);
	}
	template <class T>
	VecDoub minimize(MatDoub_I &pp, T &func)
	/*Most general interface: initial simplex specified by the matrix 
	pp[0..ndim][0..ndim-1]. Its ndim+1 rows are ndim-dimensional vectors that are 
	the vertices of the starting simplex.*/
	{
		const Int NMAX=5000;  //Maximum allowed number of function evaluations
		const Doub TINY=1.0e-10;
		Int ihi,ilo,inhi;	//highest; lowset; next highest
		mpts=pp.nrows();
		ndim=pp.ncols();
		VecDoub psum(ndim),pmin(ndim),x(ndim);
		p=pp;
		y.resize(mpts);
		for (Int i=0;i<mpts;i++) {
			for (Int j=0;j<ndim;j++)
				x[j]=p[i][j];
			y[i]=func(x);
		}
		nfunc=0;
		get_psum(p,psum);
		for (;;) {
			ilo=0;
			/*First we must determine which point is the highest (worst), 
			next-highest, and lowest (best), by looping over the points in the 
			simplex.*/
			ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
			for (Int i=0;i<mpts;i++) {
				if (y[i] <= y[ilo]) ilo=i;
				if (y[i] > y[ihi]) {
					inhi=ihi;
					ihi=i;
				} else if (y[i] > y[inhi] && i != ihi) inhi=i;
			}
			Doub rtol=2.0*abs(y[ihi]-y[ilo])/(abs(y[ihi])+abs(y[ilo])+TINY);
			/*Compute the fractional range from highest to lowest and return 
			if satisfactory.*/
			if (rtol < ftol) {	//If returning, put best point and value in slot 0.
				SWAP(y[0],y[ilo]);
				for (Int j=0;j<ndim;j++) {
					SWAP(p[0][j],p[ilo][j]);
					pmin[j]=p[0][j];
				}
				fmin=y[0];
				return pmin;
			}
			if (nfunc >= NMAX) throw("NMAX exceeded");
			nfunc += 2;
			/*Begin a new iteration. First extrapolate by a factor -1 through the
			face of the simplex across from the high point, i.e., reflect the 
			simplex from the high point.*/
			Doub ytry=amotry(p,y,psum,ihi,-1.0,func);
			if (ytry <= y[ilo])
				/*Gives a result better than the best point, so try an additional
				extrapolation by a factor 2.*/
				ytry=amotry(p,y,psum,ihi,2.0,func);
			else if (ytry >= y[inhi]) {
				/*The reflected point is worse than the second-highest, so look 
				for an intermediate lower point, i.e., do a one-dimensional 
				contraction.*/
				Doub ysave=y[ihi];
				ytry=amotry(p,y,psum,ihi,0.5,func);
				if (ytry >= ysave) {  //worse than highest, Can’t seem to get rid of that high point.
					for (Int i=0;i<mpts;i++) {	//Better contract around lowest (best) pt
						if (i != ilo) {
							for (Int j=0;j<ndim;j++)
								p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
							y[i]=func(psum);
						}
					}
					nfunc += ndim;	//Keep track of function evaluations.
					get_psum(p,psum);	//Recompute psum.
				}
			} else --nfunc;	//Correct the evaluation count.
		}	//Go back for the test of doneness and the next iteration
	}
	inline void get_psum(MatDoub_I &p, VecDoub_O &psum)
	/*Utility function, sum of rows for each col of p*/
	{
		for (Int j=0;j<ndim;j++) {
			Doub sum=0.0;
			for (Int i=0;i<mpts;i++)
				sum += p[i][j];
			psum[j]=sum;
		}
	}
	template <class T>
	Doub amotry(MatDoub_IO &p, VecDoub_O &y, VecDoub_IO &psum,
		const Int ihi, const Doub fac, T &func)
	/*Helper function: Extrapolates by a factor fac through the face of the 
	simplex across from the high point, tries it, and replaces the high point
	if the new point is better.*/	
	{
		VecDoub ptry(ndim);
		Doub fac1=(1.0-fac)/ndim;
		Doub fac2=fac1-fac;
		for (Int j=0;j<ndim;j++)
			ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
		Doub ytry=func(ptry); //Evaluate the function at the trial point.
		if (ytry < y[ihi]) {  //If it’s better than the highest, then replace the highest
			y[ihi]=ytry;
			for (Int j=0;j<ndim;j++) {
				psum[j] += ptry[j]-p[ihi][j];
				p[ihi][j]=ptry[j];
			}
		}
		return ytry;
	}
};
