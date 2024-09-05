#include <RcppArmadillo.h>
#include<omp.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]

using namespace arma;

#define tolr 1e-5
#define crossprod(x) symmatu((x).t() * (x))
#define tcrossprod(x) symmatu((x) * (x).t())

#define check_symmetric(A) cout << A.is_symmetric() << endl;

struct theta{
  arma::mat Phi;
  arma::field<arma::mat> A;
  arma::vec ps;
};

struct hyperparams{
  arma::vec ps,A;
  arma::mat Phi;
};


void copy_theta(const theta &src,theta &dest){
  dest.Phi=src.Phi;
  dest.A=src.A;
  dest.ps=src.ps;
}

double theta_ssq(const theta &x){
  double s= dot( x.Phi, x.Phi ) +dot( x.ps, x.ps );
  for(const auto &it:x.A)
    s+=dot(it, it);
  return s;
}

// [[Rcpp::export]]
double rig(double mu){
  double y = randn<double>();
  y *= y;
  double mu2 = pow(mu,2); //gsl_pow_2(mu);
  double quad = 4 * mu * y + mu2 * pow(y,2);
  // double x = mu + y * mu2 / 2 - mu / 2  * sqrt(quad);
  double  x_div_mu=(2+mu*y - sqrt(quad) )/2;
  double  x = (mu* x_div_mu);
  if(x<=0)
    return(1e-5);
  
  double u = log (randu<double>());
  // if(u <= (mu / (x + mu))) return x;
  if(u <=  -log1p(x_div_mu) ) return x;
  else return mu / x_div_mu;
}

double    _unur_bessel_k_nuasympt (double x, double nu, int islog, int expon_scaled)    {
#define M_LNPI     1.14472988584940017414342735135      /* ln(pi) */
  
  double z;                   /* rescaled argument for K_nu() */
  double sz, t, t2, eta;      /* auxiliary variables */
  double d, u1t,u2t,u3t,u4t;  /* (auxiliary) results for Debye polynomials */
  double res;                 /* value of log(K_nu(x)) [= result] */
  
  /* rescale: we comute K_nu(z * nu) */
  z = x / nu;
  
  /* auxiliary variables */
  sz = hypot(1,z);   /* = sqrt(1+z^2) */
  t = 1. / sz;
  t2 = t*t;
  
  eta = (expon_scaled) ? (1./(z + sz)) : sz;
  eta += log(z) - log1p(sz);                  /* = log(z/(1+sz)) */
  
  /* evaluate Debye polynomials u_j(t) */
  u1t = (t * (3. - 5.*t2))/24.;
  u2t = t2 * (81. + t2*(-462. + t2 * 385.))/1152.;
  u3t = t*t2 * (30375. + t2 * (-369603. + t2 * (765765. - t2 * 425425.)))/414720.;
  u4t = t2*t2 * (4465125.
                   + t2 * (-94121676.
                   + t2 * (349922430.
                   + t2 * (-446185740.
                   + t2 * 185910725.)))) / 39813120.;
                   d = (-u1t + (u2t + (-u3t + u4t/nu)/nu)/nu)/nu;
                   
                   /* log(K_nu(x)) */
                   res = log(1.+d) - nu*eta - 0.5*(log(2.*nu*sz) - M_LNPI);
                   
                   return (islog ? res : exp(res));
}

double _gig_mode(double lambda, double omega)
  /*---------------------------------------------------------------------------*/
  /* Compute mode of GIG distribution.                                         */
  /*                                                                           */
  /* Parameters:                                                               */
  /*   lambda .. parameter for distribution                                    */
  /*   omega ... parameter for distribution                                    */
  /*                                                                           */
  /* Return:                                                                   */
  /*   mode                                                                    */
  /*---------------------------------------------------------------------------*/
{
  if (lambda >= 1.)
    /* mode of fgig(x) */
    return (sqrt((lambda-1.)*(lambda-1.) + omega*omega)+(lambda-1.))/omega;
  else
    /* 0 <= lambda < 1: use mode of f(1/x) */
    return omega / (sqrt((1.-lambda)*(1.-lambda) + omega*omega)+(1.-lambda));
} /* end of _gig_mode() */
    
    void
    _rgig_ROU_shift_alt (double *res, int n, double lambda, double lambda_old, double omega, double alpha)
  /*---------------------------------------------------------------------------*/
  /* Type 8:                                                                   */
  /* Ratio-of-uniforms with shift by 'mode', alternative implementation.       */
  /*   Dagpunar (1989)                                                         */
  /*   Lehner (1989)                                                           */
  /*---------------------------------------------------------------------------*/
    {
      double xm, nc;     /* location of mode; c=log(f(xm)) normalization constant */
  double s, t;       /* auxiliary variables */
  double U, V, X;    /* random variables */
  
  int i;             /* loop variable (number of generated random variables) */
  int count = 0;     /* counter for total number of iterations */
  
  double a, b, c;    /* coefficent of cubic */
  double p, q;       /* coefficents of depressed cubic */
  double fi, fak;    /* auxiliary results for Cardano's rule */
  
  double y1, y2;     /* roots of (1/x)*sqrt(f((1/x)+m)) */
  
  double uplus, uminus;  /* maximum and minimum of x*sqrt(f(x+m)) */
  
  /* -- Setup -------------------------------------------------------------- */
  
  /* shortcuts */
  t = 0.5 * (lambda-1.);
  s = 0.25 * omega;
  
  /* mode = location of maximum of sqrt(f(x)) */
  xm = _gig_mode(lambda, omega);
  
  /* normalization constant: c = log(sqrt(f(xm))) */
  nc = t*log(xm) - s*(xm + 1./xm);
  
  /* location of minimum and maximum of (1/x)*sqrt(f(1/x+m)):  */
  
  /* compute coeffients of cubic equation y^3+a*y^2+b*y+c=0 */
  a = -(2.*(lambda+1.)/omega + xm);       /* < 0 */
  b = (2.*(lambda-1.)*xm/omega - 1.);
  c = xm;
  
  /* we need the roots in (0,xm) and (xm,inf) */
  
  /* substitute y=z-a/3 for depressed cubic equation z^3+p*z+q=0 */
  p = b - a*a/3.;
  q = (2.*a*a*a)/27. - (a*b)/3. + c;
  
  /* use Cardano's rule */
  fi = acos(-q/(2.*sqrt(-(p*p*p)/27.)));
  fak = 2.*sqrt(-p/3.);
  y1 = fak * cos(fi/3.) - a/3.;
  y2 = fak * cos(fi/3. + 4./3.*M_PI) - a/3.;
  
  /* boundaries of minmal bounding rectangle:                  */
  /* we us the "normalized" density f(x) / f(xm). hence        */
  /* upper boundary: vmax = 1.                                 */
  /* left hand boundary: uminus = (y2-xm) * sqrt(f(y2)) / sqrt(f(xm)) */
  /* right hand boundary: uplus = (y1-xm) * sqrt(f(y1)) / sqrt(f(xm)) */
  uplus  = (y1-xm) * exp(t*log(y1) - s*(y1 + 1./y1) - nc);
  uminus = (y2-xm) * exp(t*log(y2) - s*(y2 + 1./y2) - nc);
  
  /* -- Generate sample ---------------------------------------------------- */
  
  for (i=0; i<n; i++) {
    do {
      ++count;
      U = uminus + unif_rand() * (uplus - uminus);    /* U(u-,u+)  */
  V = unif_rand();                                /* U(0,vmax) */
  X = U/V + xm;
    }                                         /* Acceptance/Rejection */
  while ((X <= 0.) || ((log(V)) > (t*log(X) - s*(X + 1./X) - nc)));
    
    /* store random point */
    res[i] = (lambda_old < 0.) ? (alpha / X) : (alpha * X);
  }
  
  /* -- End ---------------------------------------------------------------- */
  
  return;
    } /* end of _rgig_ROU_shift_alt() */
  
  void _rgig_newapproach1 (double *res, int n, double lambda, double lambda_old, double omega, double alpha)
  /*---------------------------------------------------------------------------*/
  /* Type 4:                                                                   */
  /* New approach, constant hat in log-concave part.                           */
  /* Draw sample from GIG distribution.                                        */
  /*                                                                           */
  /* Case: 0 < lambda < 1, 0 < omega < 1                                       */
  /*                                                                           */
  /* Parameters:                                                               */
  /*   n ....... sample size (positive integer)                                */
  /*   lambda .. parameter for distribution                                    */
  /*   omega ... parameter for distribution                                    */
  /*                                                                           */
  /* Return:                                                                   */
  /*   random sample of size 'n'                                               */
  /*---------------------------------------------------------------------------*/
  {
    /* parameters for hat function */
    double A[3], Atot;  /* area below hat */
    double k0;          /* maximum of PDF */
    double k1, k2;      /* multiplicative constant */
    
    double xm;          /* location of mode */
    double x0;          /* splitting point T-concave / T-convex */
    double a;           /* auxiliary variable */
    
    double U, V, X;     /* random numbers */
    double hx;          /* hat at X */
    
    int i;              /* loop variable (number of generated random variables) */
    int count = 0;      /* counter for total number of iterations */
    
    /* -- Check arguments ---------------------------------------------------- */
    
    if (lambda >= 1. || omega >1.)
      Rcpp::stop ("invalid parameters");
    
    /* -- Setup -------------------------------------------------------------- */
    
    /* mode = location of maximum of sqrt(f(x)) */
    xm = _gig_mode(lambda, omega);
    
    /* splitting point */
    x0 = omega/(1.-lambda);
    
    /* domain [0, x_0] */
    k0 = exp((lambda-1.)*log(xm) - 0.5*omega*(xm + 1./xm));     /* = f(xm) */
    A[0] = k0 * x0;
    
    /* domain [x_0, Infinity] */
    if (x0 >= 2./omega) {
      k1 = 0.;
      A[1] = 0.;
      k2 = pow(x0, lambda-1.);
      A[2] = k2 * 2. * exp(-omega*x0/2.)/omega;
    }
    
    else {
      /* domain [x_0, 2/omega] */
      k1 = exp(-omega);
      A[1] = (lambda == 0.)
        ? k1 * log(2./(omega*omega))
          : k1 / lambda * ( pow(2./omega, lambda) - pow(x0, lambda) );
      
      /* domain [2/omega, Infinity] */
      k2 = pow(2/omega, lambda-1.);
      A[2] = k2 * 2 * exp(-1.)/omega;
    }
    
    /* total area */
    Atot = A[0] + A[1] + A[2];
    
    /* -- Generate sample ---------------------------------------------------- */
    
    for (i=0; i<n; i++) {
      do {
        ++count;
        
        /* get uniform random number */
        V = Atot * unif_rand();
        
        do {
          
          /* domain [0, x_0] */
          if (V <= A[0]) {
            X = x0 * V / A[0];
            hx = k0;
            break;
          }
          
          /* domain [x_0, 2/omega] */
          V -= A[0];
          if (V <= A[1]) {
            if (lambda == 0.) {
              X = omega * exp(exp(omega)*V);
              hx = k1 / X;
            }
            else {
              X = pow(pow(x0, lambda) + (lambda / k1 * V), 1./lambda);
              hx = k1 * pow(X, lambda-1.);
            }
            break;
          }
          
          /* domain [max(x0,2/omega), Infinity] */
          V -= A[1];
          a = (x0 > 2./omega) ? x0 : 2./omega;
          X = -2./omega * log(exp(-omega/2. * a) - omega/(2.*k2) * V);
          hx = k2 * exp(-omega/2. * X);
          break;
          
        } while(0);
        
        /* accept or reject */
        U = unif_rand() * hx;
        
        if (log(U) <= (lambda-1.) * log(X) - omega/2. * (X+1./X)) {
          /* store random point */
          res[i] = (lambda_old < 0.) ? (alpha / X) : (alpha * X);
          break;
        }
      } while(1);
      
    }
    
    /* -- End ---------------------------------------------------------------- */
    
    return;
  } /* end of _rgig_newapproach1() */
    
    void _rgig_ROU_noshift (double *res, int n, double lambda, double lambda_old, double omega, double alpha)
  /*---------------------------------------------------------------------------*/
  /* Tpye 1:                                                                   */
  /* Ratio-of-uniforms without shift.                                          */
  /*   Dagpunar (1988), Sect.~4.6.2                                            */
  /*   Lehner (1989)                                                           */
  /*---------------------------------------------------------------------------*/
    {
      double xm, nc;     /* location of mode; c=log(f(xm)) normalization constant */
  double ym, um;     /* location of maximum of x*sqrt(f(x)); umax of MBR */
  double s, t;       /* auxiliary variables */
  double U, V, X;    /* random variables */
  
  int i;             /* loop variable (number of generated random variables) */
  int count = 0;     /* counter for total number of iterations */
  
  /* -- Setup -------------------------------------------------------------- */
  
  /* shortcuts */
  t = 0.5 * (lambda-1.);
  s = 0.25 * omega;
  
  /* mode = location of maximum of sqrt(f(x)) */
  xm = _gig_mode(lambda, omega);
  
  /* normalization constant: c = log(sqrt(f(xm))) */
  nc = t*log(xm) - s*(xm + 1./xm);
  
  /* location of maximum of x*sqrt(f(x)):           */
  /* we need the positive root of                   */
  /*    omega/2*y^2 - (lambda+1)*y - omega/2 = 0    */
  ym = ((lambda+1.) + sqrt((lambda+1.)*(lambda+1.) + omega*omega))/omega;
  
  /* boundaries of minmal bounding rectangle:                   */
  /* we us the "normalized" density f(x) / f(xm). hence         */
  /* upper boundary: vmax = 1.                                  */
  /* left hand boundary: umin = 0.                              */
  /* right hand boundary: umax = ym * sqrt(f(ym)) / sqrt(f(xm)) */
  um = exp(0.5*(lambda+1.)*log(ym) - s*(ym + 1./ym) - nc);
  
  /* -- Generate sample ---------------------------------------------------- */
  
  for (i=0; i<n; i++) {
    do {
      ++count;
      U = um * unif_rand();        /* U(0,umax) */
  V = unif_rand();             /* U(0,vmax) */
  X = U/V;
    }                              /* Acceptance/Rejection */
  while (((log(V)) > (t*log(X) - s*(X + 1./X) - nc)));
    
    /* store random point */
    res[i] = (lambda_old < 0.) ? (alpha / X) : (alpha * X);
  }
  
  /* -- End ---------------------------------------------------------------- */
  
  return;
    } /* end of _rgig_ROU_noshift() */
  
  
  
#define ZTOL (datum::eps*10.0)
  /*---------------------------------------------------------------------------*/
  double rgig( double lambda, double chi, double psi)
  {
    double omega, alpha;     /* parameters of standard distribution */
  // SEXP sexp_res;           /* results */
  double *res;
  int i;
  int n=1;
  
  /* check sample size */
  if (n<=0) {
    Rcpp::stop("sample size 'n' must be positive integer.");
  }
  
  /* check GIG parameters: */
  if ( !(R_FINITE(lambda) && R_FINITE(chi) && R_FINITE(psi)) ||
  (chi <  0. || psi < 0)      ||
  (chi == 0. && lambda <= 0.) ||
  (psi == 0. && lambda >= 0.) ) {
    // cout<<"lambda="<<lambda<<", chi="<<chi<<", psi"=psi<<endl;
    printf("lambda= %lf, chi=%lf, psi= %lf\n",lambda, chi,psi);
    Rcpp::stop("invalid parameters for GIG distribution!!");
  }
  
  /* allocate array for random sample */
  // PROTECT(sexp_res = NEW_NUMERIC(n));
  res = (double *)malloc(n*sizeof(double)) ;
  
  if (chi < ZTOL) {
    /* special cases which are basically Gamma and Inverse Gamma distribution */
    if (lambda > 0.0) {
      for (i=0; i<n; i++) res[i] = R::rgamma(lambda, 2.0/psi);
    }
    else {
      for (i=0; i<n; i++) res[i] = 1.0/R::rgamma(-lambda, 2.0/psi);
    }
  }
  
  else if (psi < ZTOL) {
    /* special cases which are basically Gamma and Inverse Gamma distribution */
    if (lambda > 0.0) {
      for (i=0; i<n; i++) res[i] = 1.0/R::rgamma(lambda, 2.0/chi);
    }
    else {
      for (i=0; i<n; i++) res[i] = R::rgamma(-lambda, 2.0/chi);
    }
    
  }
  
  else {
    double lambda_old = lambda;
    if (lambda < 0.) lambda = -lambda;
    alpha = sqrt(chi/psi);
    omega = sqrt(psi*chi);
    
    /* run generator */
    do {
      if (lambda > 2. || omega > 3.) {
        /* Ratio-of-uniforms with shift by 'mode', alternative implementation */
        _rgig_ROU_shift_alt(res, n, lambda, lambda_old, omega, alpha);
        break;
      }
      
      if (lambda >= 1.-2.25*omega*omega || omega > 0.2) {
        /* Ratio-of-uniforms without shift */
        _rgig_ROU_noshift(res, n, lambda, lambda_old, omega, alpha);
        break;
      }
      
      if (lambda >= 0. && omega > 0.) {
        /* New approach, constant hat in log-concave part. */
        _rgig_newapproach1(res, n, lambda, lambda_old, omega, alpha);
        break;
      }
      
      /* else */
      Rcpp::stop("parameters must satisfy lambda>=0 and omega>0.");
      
    } while (0);
  }
  
  /* return result */
  // UNPROTECT(1);
  double ret=res[0];
  free(res);
  return ret;
  } /* end of do_rgig() */
  
  double quad_form(arma::vec x, arma::mat S){
    arma::mat ch=chol(symmatu(S),"lower");
    arma::vec tmp=solve(trimatl(ch),x);
    cout<<"ch="<<ch<<"tmp="<<tmp<<endl;
    
    return as_scalar (crossprod(tmp) );
  }
  
  void update_loading(const arma::mat &eta, arma::mat &Lambda,arma::mat &T_phi,double &T_sum,  const arma::mat &Y,double &tau, const arma::vec &ps,const double a){
    const unsigned p=Lambda.n_rows, k=Lambda.n_cols;
    // --- UPDATE PSI --- //
    arma::mat abs_theta=abs(Lambda);
    arma::mat psi_inv=(tau/T_sum)* (T_phi/   abs_theta) ;
    psi_inv.transform( [](double val) { return ( rig(val) ); } );
    
    // --- UPDATE TAU --- //
    double temp= T_sum *accu(abs_theta/T_phi );
    // cout<<"T_Phi.max= "<<T_phi.max()<<" T_Phi.min= "<<T_phi.min()<<" chi ="<<temp<<endl;
    tau = rgig((double) ((a-1)*k*p), 2* temp, 1);
    // tau = rgig((double) (1- k*p), 2* temp, 1);
    Rcpp::Rcout<<"tau updated; tau= "<<tau<<endl;
    
    
    // --- UPDATE PHI --- //
    T_phi=abs_theta;
    T_phi.transform( [&a](double val) { double tmp=rgig(a-1, 2* val, 1); return(  std::max(tmp,9e-6)) ; } );
    T_sum=accu(T_phi);
    if(T_sum<1){
      T_phi/=T_sum;
      T_sum=1;
    }
    
    Rcpp::Rcout<<"phi updated T_sum="<<T_sum<<endl;
    if(T_phi.has_nan()){
      T_phi.print("T_phi :");
      // abs_theta.print("abs_theta :");
      // arma::uvec tmppp=abs_theta.elem( find_nonfinite( T_phi ) );
      // cout<<"corresponding abs_theta"<<  tmppp <<"at that rgig ";
      Rcpp::Rcout<<"corresponding abs_theta"<<  abs_theta.elem( find_nonfinite( T_phi ) ) <<endl;
      // for(int h=0; h<tmppp.size();++h )
      //   cout<<rgig(a-1, 2* tmppp(h), 1)<<endl;
      
      Rcpp::Rcout<<"has inf in phi "<<T_phi.has_inf()<<endl;
      Rcpp::Rcout<<"T_sum"<<T_sum<<endl;
      Rcpp::stop("phi has nan!!");
    }
    
    
    // --- UPDATE PLAM --- //
    temp=pow(T_sum/tau,2);//gsl_pow_2(T_sum/tau);
    arma::mat Plam= ( temp*(psi_inv/square(T_phi)));
    Plam.transform( [](double val) { return ( val/(1+tolr*val) ); } );
    // cout<<"Plam updated; accu_plam"<<accu(Plam) <<endl;
    if(Plam.has_nan()){
      cout<<"corresponding psi_inv"<<  psi_inv.elem( find_nonfinite( Plam ) ) <<endl;
      // cout<<"corresponding T_phi"<<  T_phi.elem( find_nonfinite( Plam ) ) <<endl;
      Rcpp::stop("plam has nan!!");
    }
    
    // --- UPDATE LAMBDA --- //
    arma::cube Llamt(k,k,p);
    arma::mat eta2 = symmatu(eta.t() * eta);    // prepare eta crossproduct before the loop
#pragma omp parallel for
    for(unsigned j = 0; j < p; ++j) {
      Llamt.slice(j) = chol(diagmat(Plam.row(j)) + ps(j)*eta2 );
      Lambda.row(j) = (solve(trimatu(Llamt.slice(j)), randn<arma::vec>(k) +ps(j) *solve(trimatl((Llamt.slice(j)).t()),  (eta.t() * Y.col(j)),solve_opts::no_approx ),solve_opts::no_approx ) ).t();
    }
    // Rcpp::Rcout<<"Lambda updated"<<endl;
    if(Lambda.has_nan())
      Rcpp::stop("has nan!!");
    
    // Plam.reset(); psi_inv.reset(); abs_theta.reset();Llamt.reset();eta2.reset();
  }
  
  void update_loading_params( arma::mat &Lambda,arma::mat &T_phi, arma::mat &Plam,const double a){
    const unsigned p=Lambda.n_rows, k=Lambda.n_cols;
    // --- UPDATE PSI --- //
    arma::mat abs_theta=abs(Lambda);
    arma::mat psi_inv= (T_phi/   abs_theta) ;
    psi_inv.transform( [](double val) { return ( rig(val) ); } );
    
    
    // --- UPDATE PHI --- //
    T_phi=abs_theta;
    T_phi.transform( [&a](double val) { double tmp=rgig(a-1, 2* val, 1); return(  std::max(tmp,9e-6)) ; } );
    
    if(T_phi.has_nan()){
      T_phi.print("T_phi :");
      // abs_theta.print("abs_theta :");
      // arma::uvec tmppp=abs_theta.elem( find_nonfinite( T_phi ) );
      // cout<<"corresponding abs_theta"<<  tmppp <<"at that rgig ";
      Rcpp::Rcout<<"corresponding abs_theta"<<  abs_theta.elem( find_nonfinite( T_phi ) ) <<endl;
      // for(int h=0; h<tmppp.size();++h )
      //   cout<<rgig(a-1, 2* tmppp(h), 1)<<endl;
      
      Rcpp::Rcout<<"has inf in phi "<<T_phi.has_inf()<<endl;
      Rcpp::stop("phi has nan!!");
    }
    
    // --- UPDATE PLAM --- //
    Plam= ( psi_inv/square(T_phi));
    // Plam.transform( [](double val) { return ( val/(1+tolr*val) ); } );
    Plam=( Plam/(1+tolr*Plam) );
      // Rcpp::Rcout<<"Plam updated; accu_plam"<<accu(Plam) <<endl;
      if(Plam.has_nan()){
        Rcpp::Rcout<<"corresponding psi_inv"<<  psi_inv.elem( find_nonfinite( Plam ) ) <<endl;
        // cout<<"corresponding T_phi"<<  T_phi.elem( find_nonfinite( Plam ) ) <<endl;
        Rcpp::stop("plam has nan!!");
      }
  }
  
  void update_eta(const arma::mat &Lambda,const arma::vec &ps, const arma::mat &Y, arma::mat &eta){
    // --- UPDATE eta --- //
    unsigned k=eta.n_cols, n=eta.n_rows;
    
    // arma::mat Lmsg = Lambda.each_col() % ps;
    // arma::mat Veta1 = eye<arma::mat>(k,k) + Lmsg.t() * lambda;
    // arma::mat S = inv(trimatu(chol(Veta1)));
    // arma::mat Veta = S * S.t();
    // arma::mat Meta = Y * Lmsg * Veta;
    // arma::mat noise = randn<arma::mat>(n, k);
    // arma::mat eta = Meta + noise * S.t();
    
    arma::mat Lmsg = Lambda.each_col() % ps;
    arma::mat Lmsg2 =symmatu( Lmsg.t() * Lambda);
    // arma::vec z=Lmsg2.diag();
    // (z.t()).print("diagvec:");
    Lmsg2.diag()+=1.0;
    
    arma::mat noise(k, n,fill::randn);
    arma::mat Meta=(Y * Lmsg).t(); //of order k x n
    
    arma::mat chol_precmat_eta = chol( Lmsg2);//Cholesky of posterior prec mt of eta
    eta = (solve(trimatu( chol_precmat_eta), noise +solve(trimatl(chol_precmat_eta.t()),  Meta,solve_opts::no_approx) ,solve_opts::no_approx) ).t();
    
    // Lmsg.reset(); Lmsg2.reset(); noise.reset(); Meta.reset(); chol_precmat_eta.reset();
  }
  
//' Matrix inversion using low-rank plus diagonal decomposition
//'
//' Fast inversion of \eqn{\Sigma=\Lambda \Lambda^T+\mathrm{diag}(\texttt{ps})} using the \href{https://en.wikipedia.org/wiki/Woodbury_matrix_identity}{Woodbury matrix identity}
//'
//' @param Lambda A \eqn{d \times q} order matrix \eqn{\Lambda} typically \eqn{q\ll d} 
//' @param ps A vector of length \eqn{d} comprising positive elements
//'
//' @return The matrix \eqn{(\Lambda \Lambda^T+\mathrm{diag}(\texttt{ps}) )^{-1}}
//' @export
//'
//' @examples
//' d=1e3; q=50
//' lam=matrix(rnorm(d*q,sd=.25),nrow=d,ncol=q ); ps=rgamma(d,1,1) 
//' sig=tcrossprod(lam)+diag(ps)
//' rbenchmark::benchmark( chol2inv(chol(sig )),##Cholesky based inverse
//'  fast_fact_inv(lam,ps), ##fast_fact_inv
//'  replications = 10  ##benchmark using 10 replicates
//'  ) 
// [[Rcpp::export]]
arma::mat fast_fact_inv(const arma::mat &Lambda,const arma::vec &ps){
   arma::mat Lmsg = Lambda.each_col() / ps;
   arma::mat Lmsg2 =symmatu( Lmsg.t() * Lambda);
   Lmsg2.diag()+=1;
   
   // Lmsg2.diag().t().print("Lmsg2.diag() :");
   
   /****** VERSION 1 ******/
   arma::mat S, exp_mat2;
   bool err=chol(S,Lmsg2,"lower");
   /*if(err)
    Rcpp::stop ("Cholesky failed in fast_fact_inv!");*/
   err=solve(exp_mat2,trimatl( S) ,Lmsg.t()) ; //order of exp_mat2= d x p
   /*if(err)
    Rcpp::stop("Solve failed in fast_fact_inv!");*/
   return diagmat(1/ps)- symmatu(crossprod(exp_mat2)) ;
   /************************/
   
   /****** VERSION 2 ******/
   /*arma::mat exp_mat2=Lmsg * solve(Lmsg2 ,Lmsg.t(),solve_opts::likely_sympd) ;
    return diagmat(1/ps)- exp_mat2 ;*/
   /************************/
}
  
  arma::mat fast_fact_inv_exten(const arma::mat &Lambda_orig,const arma::mat &A,const arma::vec &sq_ps, arma::vec &s){
    arma::mat U,V;
    
    bool b=svd(U,s,V,A);
    if(!b)
      Rcpp::stop("First SVD failed in fast_fact_inv_exten!");
    
    s=sqrt(1+square(s));
    U.cols(0,s.n_elem-1).each_row()%=s.t();
    arma::mat Lambda=Lambda_orig* U;
    
    // cout<<"2nd SVD 0"<<endl;
    // Lambda.save("Lambda.csv", csv_ascii);
    b=svd(U,s,V,Lambda);
    if(!b)
      Rcpp::stop("Second SVD failed in fast_fact_inv_exten!");
    // cout<<"2nd SVD 1"<<endl;
    s%=s;
    
    arma::vec s1=sqrt(1+s);
    // d.subvec(0,s.n_elem-1)=sqrt(1+square(s));
    U.each_col() /=sq_ps;
    U.cols(0,s1.n_elem-1).each_row()/=s1.t();
    
    /*arma::vec d(k,fill::ones);
     d.subvec(0,s.n_elem-1)=sqrt(1+square(s));
     U.each_row()*=d.t();
     arma::mat Lambda=Lambda_orig* U;
     
     svd(U,s,V,Lambda);
     
     d.set_size(p);d.fill(1.0);
     d.subvec(0,s.n_elem-1)=sqrt(1+square(s));
     U.each_col() /=sq_ps;
     U.each_row()/=d.t();*/
    
    return tcrossprod(U);
  }
  
//' Determinant of a matrix using low-rank plus diagonal decomposition
//'
//' Fast determinant of \eqn{\Sigma=\Lambda \Lambda^T+\mathrm{diag}(\texttt{ps})} using the \href{https://en.wikipedia.org/wiki/Matrix_determinant_lemma}{Matrix determinant lemma}
//'
//' @param Lambda_orig A \eqn{d \times q} order matrix \eqn{\Lambda} typically \eqn{q\ll d} 
//' @param ps A vector of length \eqn{d} comprising positive elements
//' @param lg A logical variable; if \code{TRUE} (default), \code{log} of the determinant will be returned.
//'
//' @return If \code{lg=TRUE} then \code{log}-determinant of \eqn{\Lambda \Lambda^T+\mathrm{diag}(\texttt{ps}) }; else determinant of \eqn{\Lambda \Lambda^T+\mathrm{diag}(\texttt{ps}) }.
//' @export
//'
//' @examples
//' d=1e3; q=50
//' lam=matrix(rnorm(d*q,sd=.25),nrow=d,ncol=q ); ps=rgamma(d,1,1) 
//' sig=tcrossprod(lam)+diag(ps)
//' rbenchmark::benchmark( log(det(sig)), ##naive determinant
//' 2*sum(log(diag(chol(sig ))) ),##Cholesky based determinant
//'  fast_fact_det(lam,ps), ##fast_fact_det
//'  replications = 10  ##benchmark using 10 replicates
//'  ) 
// [[Rcpp::export]]
double fast_fact_det(const arma::mat &Lambda_orig,const arma::vec &ps, const bool lg=1){
  arma::vec s, sq_ps=sqrt(ps);
  arma::mat Lambda=Lambda_orig.each_col() / sq_ps;
  
  bool err=svd(s,Lambda);
  /*if(err)
  Rcpp::stop("svd failed in fast_fact_det");*/
  
  double log_det= sum(log1p(square(s)))+ sum(log(ps)) ;
  
  return (lg ? log_det : exp(log_det));
}
  
  inline double fast_fact_det_exten(const arma::vec &s, const arma::vec &ps){
    return sum(log1p(s))+ sum(ps);
  }
  
  arma::mat solve_precond(arma::mat X, arma::mat L){
    arma::vec diag_X=sqrt(X.diag());
    L.each_col()/=diag_X;
    
    X.each_col()/=diag_X;
    X.each_row()/=diag_X.t();
    X=symmatu(X);
    
    arma::mat res=solve( X, L,solve_opts::likely_sympd);
    res.each_col()/=diag_X;
    
    return res;
  }
  
  arma::mat del_v(const arma::mat &Lambda, const arma::vec &ps, const arma::mat &D, const arma::mat &S, arma::mat &del_inv_S, const int n){
    arma::mat del_inv=fast_fact_inv(Lambda, ps);
    arma::mat del_inv_L=del_inv* Lambda; del_inv_S= del_inv* S;
    
    arma::mat tmpmat= del_inv_S *del_inv_L -n*del_inv_L;
    return (Lambda % D - tmpmat  );
  }
  
  void leapfrog(const unsigned nstep,const double delta,arma::mat &Lambda, arma::mat &p_lam,const arma::vec &ps,const arma::mat &S,const unsigned n,
                arma::mat &del_inv_S,const arma::mat &D, arma::mat &v_old){
    for(unsigned i=0;i<nstep;++i){
      Lambda+= delta*(p_lam-(delta/2)*v_old);
      arma::mat v_new=del_v( Lambda, ps, D, S, del_inv_S, n);
      p_lam-=(delta/2)* ( v_old+v_new);
      v_old=v_new;
    }
  }
  
  void del_v_sig(const arma::mat &Lambda, const arma::vec &ps, const arma::mat &D, const arma::mat &S, arma::mat &del_inv_S, const unsigned n, arma::mat &del_lam,arma::vec &del_sig,const arma::vec &sig_param){
    arma::vec sig=exp(ps);
    arma::mat del_inv=fast_fact_inv(Lambda, sig);
    del_inv_S= del_inv* S;
    
    arma::mat tmpmat= n*del_inv-del_inv_S *del_inv ;
    del_sig= (tmpmat.diag() % sig) /2+ (ps-sig_param(0))/sig_param(1);
    del_lam=Lambda % D+tmpmat * Lambda;
  }
  
  void del_v_msfa(const theta &param, theta &grads, const hyperparams &hypers, const arma::cube &S,
                  arma::cube &del_inv_S, const arma::uvec &ns, arma::field<arma::mat> &Phi_sq_As, const int nthreads ){
    // Rcpp::Rcout<<"del_v_msfa flag start "<<endl;
    const unsigned nstudy=S.n_slices, p=param.Phi.n_rows, k=param.Phi.n_cols;//, k_eff=col_inds.n_elem;
    grads.A.set_size(nstudy);
    Phi_sq_As.set_size(nstudy);
    // Rcpp::Rcout<<"del_v_msfa flag mid "<<endl;
    
    arma::vec sig=exp(param.ps);
    arma::cube tmpmat(p,k,nstudy);
    arma::mat del_ps(p,nstudy);
    del_inv_S.set_size(p,p,nstudy);
    // Rcpp::Rcout<<"del_v_msfa flag loop "<<endl;
    
    omp_set_num_threads(1);
#pragma omp parallel for num_threads(nstudy)
    for(unsigned s=0;s<nstudy;++s){
      // cout<<"del_v_msfa flag s= "<<s<<endl;
      mat tmpmat2=symmatu(tcrossprod(param.A[s])); tmpmat2.diag()+=1; //C_s
      // mat sqrt_A= chol( tmpmat2, "lower" ); //chol(C_s)
      mat sqrt_A; chol(sqrt_A, tmpmat2, "lower" ); //chol(C_s)
      // cout<<"del_v_msfa sqrt_A done s= "<<s<<endl;
      Phi_sq_As(s)=param.Phi * trimatl(sqrt_A); // Lambda \times chol(I_q + A_s A_s^T, "lower")
      mat del_inv=fast_fact_inv(Phi_sq_As(s), sig);//\Sigma_s^{-1}
      // cout<<"del_v_msfa del_inv done s= "<<s<<endl;
      del_inv_S.slice(s)=del_inv * S.slice(s); //\Sigma_s^{-1} \times W_s
      mat Gs = ns(s) *del_inv -del_inv_S.slice(s) *del_inv;
      
      del_ps.col(s)=Gs.diag();
      
      tmpmat.slice(s)=Gs * param.Phi;
      grads.A[s]=param.Phi.t() * tmpmat.slice(s) * param.A[s]+ (param.A[s]-hypers.A(0))/hypers.A(1);
      
      tmpmat.slice(s)*=tmpmat2;
    }
    /* arma::cube del_inv(p,p,nstudy), Gs(p,p,nstudy),  tmpmat2(k,k,nstudy), sqrt_A(k,k,nstudy); //mGs is negative of Gs
#pragma omp parallel for num_threads(nthreads)
     for(unsigned s=0;s<nstudy;++s){
     // cout<<"del_v_msfa flag s= "<<s<<endl;
     tmpmat2.slice(s)=symmatu(tcrossprod(param.A[s])); tmpmat2.slice(s).diag()+=1.; //C_s
     // sqrt_A.slice(s)= chol( tmpmat2.slice(s), "lower" ); //chol(C_s)
     chol(sqrt_A.slice(s), tmpmat2.slice(s), "lower" ); //chol(C_s)
     Phi_sq_As(s)=param.Phi * trimatl(sqrt_A.slice(s)); // Lambda \times chol(I_q + A_s A_s^T, "lower")
     del_inv.slice(s)=fast_fact_inv(Phi_sq_As(s), sig);//\Sigma_s^{-1}
     // cout<<"del_v_msfa flag1 s= "<<s<<endl;
     del_inv_S.slice(s)=del_inv.slice(s) * S.slice(s); //\Sigma_s^{-1} \times W_s
     Gs.slice(s)= ns(s) *del_inv.slice(s)-del_inv_S.slice(s) *del_inv.slice(s) ;
     del_ps.col(s)=Gs.slice(s).diag();
     
     tmpmat.slice(s)=Gs.slice(s) * param.Phi;
     grads.A[s]=param.Phi.t() * tmpmat.slice(s) * param.A[s]+ (param.A[s]-hypers.A(0))/hypers.A(1);
     
     tmpmat.slice(s)*=tmpmat2.slice(s);
     }*/
    omp_set_num_threads(nthreads);
    // Rcpp::Rcout<<"del_v_msfa flag 1 "<<endl;
    
    grads.Phi=param.Phi % hypers.Phi;
    grads.Phi+= sum(tmpmat,2) ;
    
    grads.ps= (param.ps-hypers.ps(0))/hypers.ps(1)+(sum(del_ps,1) % sig) /2.0 ;
    // Rcpp::Rcout<<"del_v_msfa flag end "<<endl;
  }
  
  void leapfrog_msfa(const unsigned nstep, const double delta,theta &param,theta &v_old, theta &p_lam, const arma::cube &S, const arma::uvec &ns,
                     arma::cube &del_inv_S, const hyperparams &hyper, arma::field<arma::mat> &Phi_sq_As, const unsigned k_eff, const int nthreads){
    const unsigned nstudy=ns.n_elem, k=param.Phi.n_cols;
    const arma::uvec col_inds=randperm(k , k_eff );
    // col_inds.print("col_inds");
    /*theta v_old;
     del_v_msfa(param, v_old, hyper, S, del_inv_S, ns, Phi_sq_As );*/
    
    for(unsigned i=0;i<nstep;++i){
      // cout<<"flag -1 leap i="<<i<<endl;
      if(k_eff==k)
        param.Phi += delta*(p_lam.Phi-(delta/2)*v_old.Phi);
      else param.Phi.cols(col_inds) += delta*(p_lam.Phi-(delta/2)*v_old.Phi.cols(col_inds));
      param.ps+=delta*(p_lam.ps-(delta/2)*v_old.ps);
      // Rcpp::Rcout<<"flag 0 leap i="<<i<<endl;
#pragma omp parallel for num_threads(nthreads)
      for(unsigned s=0;s<nstudy;++s)
        param.A[s]+=delta*(p_lam.A[s]-(delta/2)*v_old.A[s]);
      // Rcpp::Rcout<<"flag 1 leap i="<<i<<endl;
      
      theta v_new;
      del_v_msfa(param, v_new, hyper, S, del_inv_S, ns, Phi_sq_As  , nthreads);
      // Rcpp::Rcout<<"flag 1.5 leap i="<<i<<endl;
      if(k_eff==k)
        p_lam.Phi -=(delta/2)* ( v_old.Phi+v_new.Phi);
      else{
        mat v_phi=v_old.Phi+v_new.Phi;
        p_lam.Phi -=(delta/2)* v_phi.cols(col_inds);
      }
      p_lam.ps -=(delta/2)* ( v_old.ps+v_new.ps);
      // cout<<"flag 2 leap i="<<i<<endl;
#pragma omp parallel for num_threads(nthreads)
      for(unsigned s=0;s<nstudy;++s)
        p_lam.A[s]-=(delta/2)*(v_old.A[s]+v_new.A[s]);
      // cout<<"flag 3 leap i="<<i<<endl;
      copy_theta(v_new,v_old);
    }
  }
  
  void leapfrog_sig(const unsigned nstep,const double delta,arma::mat &Lambda, arma::mat &p_lam,arma::vec &p_sig,arma::vec &ps,const arma::mat &S,const unsigned n,
                    arma::mat &del_inv_S,const arma::mat &D, arma::mat &v_old,arma::vec &v_old_sig,const arma::vec &sig_param){
    for(unsigned i=0;i<nstep;++i){
      Lambda+= delta*(p_lam-(delta/2)*v_old);
      ps+=delta*(p_sig-(delta/2)*v_old_sig);
      arma::mat v_new; arma::vec v_new_sig;
      del_v_sig( Lambda, ps, D, S, del_inv_S, n,v_new,v_new_sig,sig_param);
      p_lam-=(delta/2)* ( v_old+v_new);
      p_sig-=(delta/2)* ( v_old_sig+v_new_sig);
      v_old=v_new;  v_old_sig=v_new_sig;
    }
  }
  
//' Cholesky factorization based fast determinant of positive defnite matrices
//'
//' @param X A positive definite matrix
//' @param lg A logical variable; if \code{TRUE} (default), \code{log} of the determinant will be returned.
//'
//' @return If \code{lg=TRUE} then \code{log}-determinant of \code{X}; else determinant of \code{X}
//' @export
// [[Rcpp::export]]
double log_det_pd(const arma::mat &X, const bool lg=1){
  // (X.diag().t()).print("X.diag()");
  arma::mat ch=chol(X);
  double ret=2*sum(log(ch.diag() ) );
  return (lg ? ret : exp(ret));
}
  
  // [[Rcpp::export]]
  Rcpp::List cov_est_HMC(const double a, const arma::vec ps_hyper,
                         const int nrun,const int thin, unsigned nleapfrog,arma::vec del_range,
                         arma::mat phimat,const arma::mat Y,const unsigned leapmin=5,const unsigned leapmax=15){
    const unsigned p=phimat.n_rows,k=phimat.n_cols, n=Y.n_rows;
    unsigned acceptance=0; double acpt_prob;
    
    
    arma::vec ps(p,fill::zeros), sig_param(2);
    sig_param(1)=log1p(ps_hyper(1) /pow(ps_hyper(0),2));sig_param(0)=log(ps_hyper(0))-sig_param(1)/2;
    
    arma::mat T_phi = randg<arma::mat>( p, k, distr_param(a,0.5) );
    T_phi.transform( [](double val) { return(  std::max(val,9e-17)) ; } );
    
    //--- Initiate PLAM --- //
    arma::mat Plam= (1/square(T_phi));
    
    // --- initialise loop objects --- //
    arma::mat p_lam,del_inv_S,S= crossprod(Y);
    double H_ll_new, H_ll_old=0, H_prior_new, H_prior_old ;
    
    //get H_ll_old
    arma::mat covmat=diagmat(exp(ps))+tcrossprod(phimat);
    
    H_ll_old = (n*fast_fact_det(phimat,exp(ps))+ trace(solve(covmat,S,solve_opts::likely_sympd+solve_opts::equilibrate+solve_opts::refine	) )
                  + accu(square(ps-sig_param(0)))/sig_param(1))/2 ;
    ///////////////
    
    unsigned n_mc=std::floor((nrun)/thin);
    arma::mat ps_mat(n_mc,p);
    arma::cube  Phi_mc(p,k, n_mc,fill::zeros);
    
    Rcpp::Rcout<<"loop starts"<<endl;
    // --- loop --- //
    for(int i=0; i<nrun; ++i){
      // cout<<"Flag start!!"<<endl;
      // --- UPDATE \PHI prior params--- // //This is not the DL prior phi...this is the common loading matrix
      update_loading_params( phimat,T_phi,Plam,a);
      // Plam.print("Plam ");
      // cout<<"Phi params updated"<<endl;
      
      
      // --- UPDATE \PHI --- // //This is not the DL prior phi...this is the common loading matrix
      arma::mat p_lam= randn<arma::mat>(p,k); arma::vec p_sig=randn<arma::vec>(p);
      double kin_energy= (dot(p_lam,p_lam)+dot(p_sig,p_sig))/2;
      
      arma::mat phimat_new=phimat;arma::vec ps_new=ps;
      
      arma::mat del_v_lam; arma::vec del_v_ps;
      del_v_sig(phimat_new, ps_new, Plam, S, del_inv_S, n, del_v_lam, del_v_ps, sig_param);
      
      int nstep=R::rpois(nleapfrog);
      nstep=std::clamp(nstep,1,(int) leapmax);
      double delta= randu( distr_param(del_range(0),del_range(1)) );
      leapfrog_sig(nstep, delta,phimat_new, p_lam,p_sig, ps_new, S, n, del_inv_S,Plam, del_v_lam, del_v_ps, sig_param);
      
      H_prior_old=  dot(square(phimat) ,Plam)/2;
      H_prior_new=  dot(square(phimat_new) ,Plam)/2;
      
      //get H_ll_new
      // covmat=diagmat(ps)+tcrossprod(phimat_new);
      H_ll_new = (n*fast_fact_det(phimat_new,exp(ps_new))+ trace(del_inv_S ) + accu(square(ps_new-sig_param(0)))/sig_param(1) )/2 ;
      
      double H_new= H_ll_new+ H_prior_new+ (dot(p_lam,p_lam)+ dot(p_sig,p_sig )) /2, H_old=H_ll_old+H_prior_old+kin_energy;
      
      /*cout<<"H_ll_new= "<<H_ll_new+H_prior_new<<" H_ll_old= "<<H_ll_old+H_prior_old<<endl;
       cout<<"H_new= "<<H_new<<" H_old= "<<H_old<<endl;*/
      
      if(log(randu())< -(H_new-H_old) ){
        phimat=phimat_new;
        ps=ps_new;
        H_ll_old=H_ll_new;
        ++acceptance;
      }
      acpt_prob= ((double) acceptance  )/ ((double)(i+1) );
      
      int remainder= (i+1 );
      int quotient= (int) std::floor(remainder/thin);
      remainder-= (quotient*thin) ;
      
      if(remainder==0 ){
        Rcpp::Rcout<<"Acceptance probability="<<acpt_prob<<" delta="<<delta<<" nleapfrog= "<<nleapfrog<<endl;
        if(acpt_prob<.5){
          --nleapfrog;
          nleapfrog=std::max(leapmin,nleapfrog);
        }
        if(acpt_prob>.85){
          ++nleapfrog;
          nleapfrog=std::min(leapmax,nleapfrog);
        }
        
        Rcpp::Rcout<<"Iteration: "<<i<<endl;
        
        ps_mat.row(quotient-1)=ps.t();
        Phi_mc.slice(quotient-1) =phimat;
      }
    }
    
    return Rcpp::List::create(Rcpp::Named("Lambda") =Phi_mc,
                              Rcpp::Named("residuals") = exp(ps_mat),
                              Rcpp::Named("acceptance_probability")=acpt_prob);
  }
  
  
  /*Rcpp::List DL_MSFA(unsigned nmix, double dir_prec, double diag_psi_iw, double niw_kap,
   double temper,int approx,arma::uvec del,
   double as, double bs, double a,
   int nrun, int burn, int thin,
   Rcpp::List Lambda_list, Rcpp::List f_list,
   Rcpp::List Y_list, Rcpp::List l_list,
   arma::mat phimat, unsigned S,
   std::string current_dir
   ,bool dofactor,bool docluster){*/
  
  // // [[Rcpp::export]]
  // Rcpp::List DL_MSFA(double as, double bs, double a,
  //                    int nrun, int burn, int thin,
  //                    Rcpp::List Lambda_list, Rcpp::List f_list,
  //                    Rcpp::List Y_list, Rcpp::List l_list,
  //                    arma::mat phimat, unsigned S){
  //   unsigned j;
  //   bool dofactor=1;
  //   arma::uvec ns(S),ks(S);
  //   arma::vec taus(S),T_sums(S);
  //   // arma::field<arma::mat> Lambda(S), Y(S), l(S), T_phis(S);
  //   
  //   std::ofstream sigma_diag_file;
  //   
  //   // sigma_diag_file.open("sigma_diags.txt");
  //   
  //   arma::mat *Lambda, *Y, *l, *T_phis;
  //   T_phis=new arma::mat[S];
  //   Lambda=new arma::mat[S ];
  //   // f= new f[S];
  //   Y= new arma::mat[S];
  //   l= new arma::mat[S];
  //   
  //   const unsigned p=phimat.n_rows,k=phimat.n_cols;
  //   
  //   for(j=0;j<S;++j){
  //     Lambda[j]=Rcpp::as<arma::mat>(Lambda_list[j]);
  //     Y[j]=Rcpp::as<arma::mat>(Y_list[j]);
  //     // f[j]=as<arma::mat>(f_list[j]);
  //     l[j]=Rcpp::as<arma::mat>(l_list[j]);
  //     
  //     // const int n=Y.n_rows,p=Y.n_cols,k=eta.n_cols;
  //     
  //     ns(j)=Y[j].n_rows; ks(j)=l[j].n_cols;
  //     
  //     T_phis[j]=randg<arma::mat>( p, ks(j), distr_param(a,1.0) );
  //     T_sums(j) =accu(T_phis[j]);
  //     
  //     taus(j)=randg<double>(distr_param(p*ks(j) *a ,0.5)   );
  //   }
  //   arma::uvec cs=cumsum(ns);
  //   unsigned n=sum(ns), count=0;
  //   
  //   arma::mat E(n,p), eta(n,k), J_MC(n,p,fill::zeros),J;
  //   
  //   for(j=0;j<S;++j){
  //     E.rows(cs(j)-ns(j),cs(j)-1)= Y[j]- l[j] * Lambda[j] .t();
  //     eta.rows(cs(j)-ns(j),cs(j)-1)=Rcpp::as<arma::mat>( f_list[j]);
  //   }
  //   
  //   arma::vec ps(p,fill::ones), ps_MC(p,fill::zeros); //Diagonal variance parameters
  //   
  //   arma::mat T_phi = randg<arma::mat>( p, k, distr_param(a,0.5) );
  //   T_phi.transform( [](double val) { return(  std::max(val,9e-17)) ; } );
  //   double T_sum=accu(T_phi);
  //   if(T_sum<1){
  //     T_phi/=T_sum;
  //     T_sum=1;
  //   }
  //   double tau=randg<double>(distr_param(p*k*a ,0.5)   );
  //   
  //   // --- Initiate PLAM --- //
  //   // Plam=1/( 1+psi % square (phi*tau) );
  //   
  //   // arma::mat pi_mat(std::floor((nrun+burn)/thin),nmix);
  //   // umat alloc_var_mat(std::floor((nrun+burn)/thin), n);
  //   
  //   // --- initialise loop objects --- //
  //   arma::mat Ytil;
  //   
  //   
  //   unsigned n_mc=std::floor((nrun+burn)/thin);
  //   arma::mat ps_mat(n_mc,p);
  //   
  //   arma::cube *lambda_s_mc= new arma::cube[S], *l_s_mc= new arma::cube[S],  Phi_mc(p,k, n_mc,fill::zeros), eta_mc(n,k, n_mc,fill::zeros);
  //   for(j=0;j<S;++j){
  //     (lambda_s_mc[j]).zeros(p,ks(j), n_mc);
  //     (l_s_mc[j]).zeros(ns(j),ks(j),n_mc);
  //   }
  //   
  //   
  //   ///////////////////////////////////////////////////////////////////
  //   
  //   Rcpp::Rcout<<"loop starts"<<endl;
  //   // --- loop --- //
  //   for(int i=0; i<nrun+burn; ++i){
  //     if(dofactor){
  //       cout<<"Flag start!!"<<endl;
  //       // --- UPDATE \PHI --- // //This is not the DL prior phi...this is the common loading matrix
  //       update_loading(eta, phimat,T_phi,T_sum,  E,tau, ps,a);
  //       cout<<"Phi updated"<<endl;
  //       
  //       // --- UPDATE ETA (f in cd msMSFA notation) --- //
  //       update_eta(phimat,ps, E, eta);
  //       cout<<"eta updated"<<endl;
  //       J=eta * phimat.t()  ;
  //       
  //       for(j=0;j<S;++j){
  //         arma::mat tmpmat=Y[j]- J.rows(cs(j)-ns(j),cs(j)-1);
  //         // --- UPDATE Lambda_s  --- //
  //         update_loading(l[j], Lambda[j],T_phis[j],T_sums(j),  tmpmat,taus(j), ps,a);
  //         cout<<"Lambda["<<j<< "] updated"<<endl;
  //         // -----------------------//
  //         
  //         // ----- UPDATE L_s  ---- //
  //         update_eta(Lambda[j],ps, tmpmat, l[j]);
  //         cout<<"l["<<j<< "] updated"<<endl;
  //         // -----------------------//
  //         
  //         E.rows(cs(j)-ns(j),cs(j)-1)= Y[j]- l[j] * ((Lambda[j]) .t());
  //         cout<<"Flag end!!"<<endl;
  //       }
  //       
  //       // --- UPDATE SIGMA --- //
  //       Ytil = E-J;
  //       ps =  bs + 0.5 * sum(square(Ytil), 0).t();
  //       ps.transform( [&as,&n](double val) { return  randg(distr_param(as + 0.5*n, 1 / val ) ); ; } );
  //       cout<<"Sigma updated"<<endl;
  //     }
  //     
  //     int remainder= (i+1 );
  //     int quotient= (int) std::floor(remainder/thin);
  //     remainder-= (quotient*thin) ;
  //     
  //     if(remainder==0 ){
  //       // sigma_diag_file<<1/ps.t();
  //       Rcpp::Rcout<<"Iteration: "<<i<<endl;
  //       // pi_mat.row(quotient-1)=pi.t();
  //       // alloc_var_mat.row(quotient-1)=del.t();
  //       
  //       ps_mat.row(quotient-1)=(1/ps).t();
  //       eta_mc.slice(quotient-1) =eta;
  //       Phi_mc.slice(quotient-1) =phimat;
  //       
  //       for(j=0;j<S;++j){
  //         (l_s_mc[j]).slice(quotient-1)=l[j];
  //         (lambda_s_mc[j]).slice(quotient-1)=Lambda[j];
  //       }
  //       
  //       if( i>burn){
  //         J_MC+=J;
  //         ++count;
  //       }
  //     }
  //     
  //     // printcheck = (start+1) % 1000;
  //     // if(!printcheck){
  //     //   cout << (start+1) << endl;
  //     // }
  //   }
  //   
  //   Rcpp::List list_lambdas(S),list_ls(S);
  //   arma::field<arma::mat> J_avg(S);
  //   for(j=0;j<S;++j){
  //     // cout<<"save s= "<<j<<endl;
  //     J_avg(j) = J_MC.rows(cs(j)-ns(j),cs(j)-1)/count;
  //     list_lambdas[j]=lambda_s_mc[j];
  //     list_ls[j]=l_s_mc[j];
  //   }
  //   
  //   
  //   // delete[] inds_eq_j;
  //   
  //   delete[] T_phis;
  //   delete[] Lambda;
  //   delete[] Y;
  //   delete[] l;
  //   delete[] lambda_s_mc;
  //   delete[] l_s_mc;
  //   
  //   return Rcpp::List::create(Rcpp::Named("Denoised_data") = J_avg,
  //                             Rcpp::Named("Phi") =Phi_mc,
  //                             Rcpp::Named("eta") =eta_mc,
  //                             Rcpp::Named("Lambda_s") =list_lambdas,
  //                             Rcpp::Named("l_s") =list_ls,
  //                             Rcpp::Named("residuals") = ps_mat);
  // }
  
  // [[Rcpp::export]]
  Rcpp::List cov_est_Gibbs(double as, double bs, double a,
                     int nrun, int burn, int thin,
                     arma::mat phimat, arma::mat eta, arma::mat Y){
    const unsigned p=phimat.n_rows,k=phimat.n_cols, n=Y.n_rows;
    arma::vec ps(p,fill::ones);
    
    arma::mat T_phi = randg<arma::mat>( p, k, distr_param(a,0.5) );
    T_phi.transform( [](double val) { return(  std::max(val,9e-17)) ; } );
    double T_sum=accu(T_phi);
    if(T_sum<1){
      T_phi/=T_sum;
      T_sum=1;
    }
    double tau=randg<double>(distr_param(p*k*a ,0.5)   );
    
    // --- Initiate PLAM --- //
    // Plam=1/( 1+psi % square (phi*tau) );
    
    // arma::mat pi_mat(std::floor((nrun+burn)/thin),nmix);
    // umat alloc_var_mat(std::floor((nrun+burn)/thin), n);
    
    // --- initialise loop objects --- //
    arma::mat Ytil;
    
    unsigned n_mc=std::floor((nrun+burn)/thin);
    arma::mat ps_mat(n_mc,p);
    arma::cube  Phi_mc(p,k, n_mc,fill::zeros); //, eta_mc(n,k, n_mc,fill::zeros);
    
    Rcpp::Rcout<<"loop starts"<<endl;
    // --- loop --- //
    for(int i=0; i<nrun+burn; ++i){
      if(1){
        Rcpp::Rcout<<"Flag start!!"<<endl;
        // --- UPDATE \PHI --- // //This is not the DL prior phi...this is the common loading matrix
        update_loading(eta, phimat,T_phi,T_sum,  Y,tau, ps,a);
        Rcpp::Rcout<<"Phi updated"<<endl;
        
        // --- UPDATE ETA (f in MSFA notation) --- //
        update_eta(phimat,ps, Y, eta);
        Rcpp::Rcout<<"eta updated"<<endl;
        arma::mat J=eta * phimat.t()  ;
        
        // --- UPDATE SIGMA --- //
        Ytil = Y-J;
        ps =  bs + 0.5 * sum(square(Ytil), 0).t();
        ps.transform( [&as,&n](double val) { return  randg(distr_param(as + 0.5*n, 1 / val ) ); ; } );
        Rcpp::Rcout<<"Sigma updated"<<endl;
      }
      
      int remainder= (i+1 );
      int quotient= (int) std::floor(remainder/thin);
      remainder-= (quotient*thin) ;
      
      if(remainder==0 ){
        Rcpp::Rcout<<"Iteration: "<<i<<endl;
        // pi_mat.row(quotient-1)=pi.t();
        // alloc_var_mat.row(quotient-1)=del.t();
        
        ps_mat.row(quotient-1)=1/ps.t();
        
        Phi_mc.slice(quotient-1) =phimat;
        // eta_mc.slice(quotient-1) =eta;
      }
    }
    
    return Rcpp::List::create(Rcpp::Named("Phi") =Phi_mc,
                              Rcpp::Named("residuals") = ps_mat);
  }
  
  // [[Rcpp::export]]
  Rcpp::List SUFA_HMC(const unsigned nrun,const unsigned thin, unsigned nleapfrog,const arma::vec del_range,
                      const arma::vec ps_hyper,const arma::vec A_hyper,const double a,
                      Rcpp::List Y_list, const arma::uvec ks,  arma::mat phi_init,
                      unsigned leapmax=18,unsigned leapmin=5, const double col_prob=1.0, const int nthreads=4){
    
    const unsigned p=phi_init.n_rows, nstudy=ks.n_elem, k=phi_init.n_cols;
    arma::uvec ns(nstudy);
    theta p_lam, params, v_old;
    params.A.set_size(nstudy); p_lam.A.set_size(nstudy);
    
    arma::field<arma::mat> Phi_sq_As;
    
    // const unsigned p=Rcpp::ncol(Y_list[ 0]);
    
    cube S(p,p,nstudy), del_inv_S;
    
    
    if(col_prob<=0 || col_prob>1)
      Rcpp::stop("Enter valid proportion of columns to update!");
    // unsigned k_eff= GSL_MAX_INT (ceil(k*col_prob),1);
    const unsigned k_eff= (col_prob==1.0) ? k: (std::max ((int) ceil(k*col_prob),1));
    Rcpp::Rcout<< "k_eff="<<k_eff<<endl;
    
    params.Phi=phi_init;
    // params.Phi.randn(p,k);
    params.ps.set_size(p); params.ps.fill(0.0);
    
    hyperparams hyper;
    hyper.ps.set_size(2);
    
    hyper.ps(1)=log1p(ps_hyper(1) /pow(ps_hyper(0),2));hyper.ps(0)=log(ps_hyper(0))-hyper.ps(1)/2;
    hyper.A=A_hyper;
    
    double H_ll_new, H_ll_old=accu(square(params.ps-hyper.ps(0)))/(2*hyper.ps(1))  , H_prior_new, H_prior_old ;
    for(unsigned j=0;j<nstudy;++j){
      arma::mat tmpmat=Rcpp::as<arma::mat>(Y_list[j]);
      ns(j)=tmpmat.n_rows;
      S.slice(j)=crossprod(tmpmat);
      params.A[j].randn(k,ks(j));
      
      //get H_ll_old
      arma::mat covmat=diagmat(exp(params.ps))+tcrossprod(params.Phi*(tcrossprod(params.A[j]) ));
      
      H_ll_old+= (ns(j) * log_det_pd(covmat)+ trace(solve(covmat,S.slice(j),solve_opts::likely_sympd+solve_opts::equilibrate+solve_opts::refine	) )
                    + accu(square(params.A[j]-hyper.A(0)))/hyper.A(1))/2 ;
    }
    
    unsigned acceptance=0; double acpt_prob;
    
    arma::mat T_phi = randg<arma::mat>( p, k, distr_param(a,0.5) );
    T_phi.transform( [](double val) { return(  std::max(val,9e-17)) ; } );
    
    //--- Initiate PLAM --- //
    hyper.Phi= (1/square(T_phi));
    
    // --- initialise loop objects --- //
    // omp_set_num_threads(1);
    del_v_msfa(params, v_old, hyper, S, del_inv_S, ns , Phi_sq_As, nthreads);
    v_old.Phi-= (params.Phi %hyper.Phi); //subtract the DL prior part since they are updated separately
    
    ///////////////
    
    unsigned n_mc=std::floor( nrun/thin);
    arma::mat ps_mat(n_mc,p);
    arma::cube *A= new arma::cube[nstudy], Phi_mc(p,k, n_mc,fill::zeros) ;
    for(unsigned j=0;j<nstudy;++j)
      A[j].set_size(k,ks(j),n_mc);
    
    Rcpp::Rcout<<"loop starts"<<endl;
    // --- loop --- //
    for(unsigned i=0; i<nrun; ++i){
      // --- UPDATE \PHI prior params--- // //This is not the DL prior phi...this is the common loading matrix
      update_loading_params( params.Phi,T_phi,hyper.Phi,a);
      // hyper.Phi.print("hyper.Phi ");
      // cout<<"Phi params updated"<<endl;
      
      // --- UPDATE \PHI --- // //This is not the DL prior phi...this is the common loading matrix
      p_lam.Phi.randn(p,k_eff ); p_lam.ps.randn(size(params.ps));
      double kin_energy= dot(p_lam.Phi,p_lam.Phi)+dot(p_lam.ps,p_lam.ps);
      for(unsigned j=0;j<nstudy;++j){
        p_lam.A[j].randn(size(params.A[j]));
        kin_energy+= dot(p_lam.A[j], p_lam.A[j] );
      }
      kin_energy/=2;
      
      theta params_new, v_new;
      copy_theta(params,params_new);
      copy_theta(v_old,v_new);
      v_new.Phi+= (params.Phi %hyper.Phi);
      
      int nstep=R::rpois(nleapfrog);
      nstep=std::clamp(nstep,1,(int) leapmax);
      double delta=randu( distr_param(del_range(0),del_range(1)) );
      // Rcpp::Rcout<<"leapfrog starts"<<endl;
      leapfrog_msfa(nstep, delta, params_new, v_new, p_lam, S, ns, del_inv_S, hyper, Phi_sq_As,k_eff,nthreads);
      // Rcpp::Rcout<<"leapfrog ends"<<endl;
      
      H_prior_old= dot(square(params.Phi) , hyper.Phi)/2;
      H_prior_new= dot(square(params_new.Phi) , hyper.Phi)/2;
      
      //get H_ll_new
      H_ll_new =  accu(square(params_new.ps-hyper.ps(0) ))/(2*hyper.ps(1)) ;
      double Hvec=0;
      vec tmp_exp=exp(params_new.ps);
      omp_set_num_threads(1);
#pragma omp parallel for num_threads(nthreads) reduction(+:Hvec)
      for(unsigned j=0;j<nstudy;++j)
        Hvec+= ( ns(j) *fast_fact_det(Phi_sq_As[j],tmp_exp)+ trace(del_inv_S.slice(j) )
                   + accu(square(params_new.A[j]-hyper.A(0) ) )/hyper.A(1) ) ;
      omp_set_num_threads(nthreads);
      H_ll_new+=( Hvec/2.0 );
      
      double H_new= H_ll_new+ H_prior_new+ theta_ssq(p_lam) /2, H_old=H_ll_old+H_prior_old+kin_energy;
      /*Rcpp::Rcout<<"H_ll_new= "<<H_ll_new<<" H_ll_old= "<<H_ll_old<<endl;
       Rcpp::Rcout<<"H_new= "<<H_new<<" H_old= "<<H_old<<endl;*/
      
      if(log(randu())< -(H_new-H_old) ){
        copy_theta(params_new,params);
        H_ll_old=H_ll_new;
        copy_theta(params_new,params);
        copy_theta(v_new,v_old);
        v_old.Phi-= (params.Phi %hyper.Phi);
        ++acceptance;
      }
      acpt_prob= ((double) acceptance  )/ ((double)(i+1) );
      
      int remainder= (i+1 );
      int quotient= (int) std::floor(remainder/thin);
      remainder-= (quotient*thin) ;
      
      if(remainder==0 ){
        if(acpt_prob<.5){
          --nleapfrog;
          nleapfrog=std::max(leapmin,nleapfrog);
        }
        if(acpt_prob>.85){
          ++nleapfrog;
          nleapfrog=std::min(leapmax,nleapfrog);
        }
        
        /*if(acpt_prob<.3){
         delta-= 1e-4;
         delta= GSL_MAX_DBL(delta,.006);
        }
         if(acpt_prob>.67){
         delta+= 1e-4;
         }*/
        
        Rcpp::Rcout<<"Iteration: "<<i<<" Acceptance probability="<<acpt_prob<<" delta="<<delta<<" nleapfrog= "<<nleapfrog<<endl;
        
        ps_mat.row(quotient-1)=(params.ps).t();
        Phi_mc.slice(quotient-1) =params.Phi;
        
        for(unsigned j=0;j<nstudy;++j)
          A[j].slice(quotient-1)=params.A[j];
      }
    }
    
    Rcpp::List list_A(nstudy);
    for(unsigned j=0;j<nstudy;++j){
      // cout<<"save s= "<<j<<endl;
      list_A[j]=A[j];
    }
    
    delete [] A;
    
    return Rcpp::List::create(Rcpp::Named("Lambda") =Phi_mc,
                              Rcpp::Named("A") = list_A,
                              Rcpp::Named("residuals") = exp(ps_mat),
                              Rcpp::Named("acceptance_probability")=acpt_prob);
  }
  
  // [[Rcpp::export]]
  arma::mat matchalign_permmat(arma::mat lambda, arma::mat pivot) {
    arma::mat refr = join_rows(lambda, -lambda);
    int k = lambda.n_cols;
    uvec ind = regspace<uvec> (0, k-1);
    uvec perm(k);
    arma::vec signs(k);
    rowvec norms(2*k);
    unsigned int w, c, wc;
    arma::mat diff, diffsq;
    
    for(int i=0; i<k; i++){
      diff = refr.each_col() - pivot.col(i);
      diffsq = square(diff);
      norms = sum(diffsq);
      w = index_min(norms);
      c = refr.n_cols / 2;
      if(w>=c){
        wc = w-c;
        signs(i) = -1;
        perm(i) = ind(wc);
        refr.shed_col(w);
        refr.shed_col(wc);
        ind.shed_row(wc);}
      else {
        wc = w+c;
        signs(i) = 1;
        perm(i) = ind(w);
        refr.shed_col(wc);
        refr.shed_col(w);
        ind.shed_row(w);}
    }
    
    arma::mat permmat = zeros<arma::mat>(k,k);
    for(int i=0; i<k; i++){
      permmat(perm(i), i) = signs(i);
    }
    
    /*lambda *= permmat;
     return Rcpp::wrap(lambda);*/
    return permmat;
  }
  