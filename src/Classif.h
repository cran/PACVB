class Hinge{
  public:
	double Obj(mat m,double s, const mat &X, const mat &Y,double v,double l,int B=0)
	{
        boost::math::normal_distribution<> d(0,1);
		int n=X.n_rows;
		int p=X.n_cols;
		double sum=0;
        int k;
        boost::random::uniform_int_distribution<> U1(0,n-1);
        RandomG::Random<boost::random::uniform_int_distribution<> > runif1(U1);
		//double s=exp(ls);
        if(B==0){B=n;}        
		for(int i=0;i<B;i++)
		{
            if(n!=B)
            {
                k=runif1();
            }else{
                k=i;
            }
			mat G=Y(k,0)*X(k,span::all);
			double nG=norm(G);
			double mG=dot(G,m);
			double ratio=(double)(1-mG)/(s*nG);
			sum+=(1-mG)*cdf(d,ratio)+s*nG*pdf(d,ratio);
		}
		sum=(l/B)*sum+0.5*dot(m,m)/v-0.5*p*(log(s*s)-s*s/v);
		return sum;		
	}	
	
	mat Deriv(mat m,double s, const mat &X, const mat &Y,double v,double l,int B=0)
	{
		int n=X.n_rows;
		int p=X.n_cols;
        int k;
		mat sum1(1,p,fill::zeros);
		double sum2=0;
		boost::math::normal_distribution<> d(0,1);
        boost::random::uniform_int_distribution<> U1(0,n-1);
        RandomG::Random<boost::random::uniform_int_distribution<> > runif1(U1);
        if(B==0){ B=n;}        
		for(int i=0;i<B;i++)
		{
            if(n!=B)
            {
                k=runif1();
            }else{
                k=i;
            }
			mat G=Y(k,0)*X(k,span::all);
			double nG=norm(G);
			double mG=dot(G,m);
			double ratio=(double)(1-mG)/(s*nG);
			sum1=sum1-G*cdf(d,ratio);
			sum2+=pdf(d,ratio)*nG;
		}
		sum1=(l/B)*sum1+m.t()/v;
		sum2=(l/B)*sum2-p*(1.0/(s)-s/v);
		mat deriv(1,p+1);
		deriv(0,span(0,p-1))=sum1;
		deriv(0,p)=sum2;
		return deriv;
	}
};

class HingeAUC{
  public:
    HingeAUC(const mat &Y)
    {
        _n1=sum(sum(0.5*(Y+1)));
        _n=Y.n_rows;
        _n0=_n-_n1;
        _Y0=mat(_n0,1);
        _Y1=mat(_n1,1);
        int k1=0,k0=0;
        for(int i=0;i<_n;i++)
        {
            if(Y(i,0)==1)
            {
                _Y1(k1,0)=i;
                k1++;
            }else{
                _Y0(k0,0)=i;
                k0++;
            }
        }
    }
	double Obj(mat m,double s, const mat &X, const mat &Y,double v,double l,int B=0)
	{
		boost::math::normal_distribution<> d(0,1);
		int n=X.n_rows;
		int p=X.n_cols;
		double sum=0;
        int n0=n-_n1;
        boost::random::uniform_int_distribution<> U1(0,_n1-1);
        RandomG::Random<boost::random::uniform_int_distribution<> > runif1(U1);
        boost::random::uniform_int_distribution<> U2(0,_n0-1);
        RandomG::Random<boost::random::uniform_int_distribution<> > runif0(U2);
		//double s=exp(ls);
        if(B==0) B=n;        
        int n01=0;
		for(int i=0;i<B;i++)
		{
            if(n!=B)
            {
               int k=runif0();
               int l=runif1();
               mat G=(_Y1(l,0)-_Y0(k,0))*(X(k,span::all)-X(l,span::all));
               double nG=sqrt(dot(G,G));
               double mG=dot(G,m);
               double ratio=(1-mG)/(s*nG);
               sum+=(1-mG)*cdf(d,ratio)+s*nG*pdf(d,ratio);
               n01++;
            }else{
               if(Y(i,0)==1)
               {
                for(int j=0;j<n;++j)
                {
                    if(Y(j,0)==-1)
                    {             
                        mat G=(Y(i,0)-Y(j,0))*(X(i,span::all)-X(j,span::all));
                        double nG=sqrt(dot(G,G));
                        double mG=dot(G,m);
                        double ratio=(1-mG)/(s*nG);
                        sum+=(1-mG)*cdf(d,ratio)+s*nG*pdf(d,ratio);
                        n01++;
                    }
                }
               }
            }
        }
		sum=(l/n01)*sum+dot(m,m)/v-0.5*p*(log(s*s)-s*s/v);
		return sum;		
	}	
	
	mat Deriv(mat m,double s, const mat &X, const mat &Y,double v,double l, int B=0)
	{
		int n=X.n_rows;
		int p=X.n_cols;
		mat sum1(1,p,fill::zeros);
		double sum2=0;
		boost::math::normal_distribution<> d(0,1);
        int n0=n-_n1;
        boost::random::uniform_int_distribution<> U1(0,_n1-1);
        RandomG::Random<boost::random::uniform_int_distribution<> > runif1(U1);
        boost::random::uniform_int_distribution<> U2(0,_n0-1);
        RandomG::Random<boost::random::uniform_int_distribution<> > runif0(U2);

        int n01=0;
		if(B==0){ B=n;}
        for(int i=0;i<B;i++)
		{
            if(n!=B)
            {
                    int k=runif1();
                    int l=runif0();
                    mat G=(_Y1(k,0)-_Y0(l,0))*(X(k,span::all)-X(l,span::all));
                    double nG=sqrt(dot(G,G));
                    double mG=dot(G,m);
                    double ratio=(1-mG)/(s*nG);
                    sum1=sum1-G*cdf(d,ratio);
                    sum2+=pdf(d,ratio)*nG;
                    n01++;
            }else{
                if(Y(i,0)==1)
                {
                  for(int j=0;j<n;++j)
                  {
                    if(Y(j,0)==-1)
                    {
                        mat G=(Y(i,0)-Y(j,0))*(X(i,span::all)-X(j,span::all));
                        double nG=sqrt(dot(G,G));
                        double mG=dot(G,m);
                        double ratio=(1-mG)/(s*nG);
                        sum1=sum1-G*cdf(d,ratio);
                        sum2+=pdf(d,ratio)*nG;
                        n01++;
                    }
                  }
                }
            }
		}
		sum1=(l/n01)*sum1+m.t()/v;
		sum2=(l/n01)*sum2-p*(1.0/(s)-s/v);
//		sum2=sum2*s;
		mat deriv(1,p+1);
		deriv(0,span(0,p-1))=sum1;
		deriv(0,p)=sum2;
		return deriv;
	}
  private:
    int _n;
    int _n1;
    int _n0;
    mat _Y0;
    mat _Y1; 
};

template<class Funct>
class GradientDescent
{
  public:
    GradientDescent(int K, Funct *funct)
    {
      _K=K;
      _funct=funct;
    }
    mat operator()(mat m, double ls,const mat &X,const mat &Y,double v,double l,bool linesearch, int B)
    {
      int k=0;
      mat mp=m;
      int p=m.n_rows;
      double eta=sqrt(p+1)/sqrt(_K+1);
      double si=exp(ls);
      double sip=si;
      double err=1;
      double alpha;
      while(k<_K)
      {
        if(B==0){eta=1.0/pow(k+100,0.7);}
        
        mat grad=(*_funct).Deriv(m,si,X,Y,v,l,B);
        mat grad1=grad(0,span(0,p-1));
        double grad2=grad(0,p);
	    double foo=norm(grad,2);
        mat pkm=-eta*grad1.t()/foo;
        double pks=-eta*grad2/foo;
	    if(linesearch){
            alpha=findalpha(m,si,X,Y,v,l, pks,pkm);
        }else{
            alpha=1;
        }
        m=mp+alpha*pkm;
        si=sip+alpha*pks;
        if(si<0.000001) si=0.000001;
	    //cout << "obj" << (*_funct).Obj(m,si,X,Y,v,l) << "\n";
        err=norm(m-mp);
        mp=m;
        sip=si;
        k++;
      }
      mat theta(1,p+1);
      theta(0,span(0,p-1))=m.t();
      theta(0,p)=log(si);
      _F=(*_funct).Obj(m,si,X,Y,v,l,B);
      return theta;
    }
	double findalpha(mat m,double si,const mat &X,const mat &Y,double v,double l,double pks, mat pkm)
	{
		double alpha=1,Tau=0.99,c1=0.0001,c2=0.9,fap,f,pf,pfap;
		bool Wolfe=false;
		f=(*_funct).Obj(m,si,X,Y,v,l);
		mat df=(*_funct).Deriv(m,si,X,Y,v,l);
        int p=pkm.n_elem;
		mat pk(1,p+1);
        pk(0,span(0,p-1))=pkm.t();
        pk(0,p)=pks;
        pf=dot(df,pk);
		int count=0;
		do{
			mat nm=m+alpha*pkm;
			double ns=si+alpha*pks;
			fap=(*_funct).Obj(nm,ns,X,Y,v,l);
			mat dfap=(*_funct).Deriv(nm,ns,X,Y,v,l);
			pfap=dot(dfap,pk);
			if((fap<=f+c1*alpha*pf)&(pfap>=c2*pf))
			{
				Wolfe=true;
			}else{
				alpha*=Tau;
			}
			count++;
		}while(!Wolfe & count <200);
		if(count>=200){ alpha=0.01;}
		return alpha;
	
	}
    double Get_F(){return _F;}
  private:
    int _K;
    Funct *_funct;
    double _F;
};
template<class Funct>
class GradientDescentF0
{
  public:
    GradientDescentF0(int K, Funct *funct)
    {
      _K=K;
      _funct=funct;
    }
    mat operator()(mat m,const mat &X,const mat &Y,double v,double l,bool linesearch,int B)
    {
      int k=0;
      mat mp=m;
      int n=X.n_rows;
      int p=m.n_rows;
      double eta=sqrt(p+1)/sqrt(_K+1);
      double si=1.0/n;
      double sip=si;
      double err=1;
      double alpha;
      while(k<_K)
      {
        if(B==0){eta=1.0/pow(k+100,0.7);}
        mat grad=(*_funct).Deriv(m,si,X,Y,v,l,B);
        mat grad1=grad(0,span(0,p-1));
	    double foo=norm(grad1,2);
        mat pkm=-eta*grad1.t()/foo;
	    if(linesearch){
            alpha=findalpha(m,X,Y,v,l,pkm);
        }else{
            alpha=1;
        }
        m=mp+alpha*pkm;
	    //cout << "obj" << (*_funct).Obj(m,si,X,Y,v,l) << "\n";
        err=norm(m-mp);
        mp=m;
        sip=si;
        k++;
      }
      mat theta(1,p);
      theta(0,span(0,p-1))=m.t();
      _F=(*_funct).Obj(m,si,X,Y,v,l,B);
      return theta;
    }
	double findalpha(mat m,const mat &X,const mat &Y,double v,double l, mat pkm)
	{
		double alpha=1,Tau=0.99,c1=0.0001,c2=0.9,fap,f,pf,pfap;
		bool Wolfe=false;
        int n=X.n_rows;
        double si=1.0/n;
		f=(*_funct).Obj(m,si,X,Y,v,l);
		mat df=(*_funct).Deriv(m,si,X,Y,v,l);
        int p=pkm.n_elem;
        pf=dot(df,pkm);
		int count=0;
		do{
			mat nm=m+alpha*pkm;
			fap=(*_funct).Obj(nm,si,X,Y,v,l);
			mat dfap=(*_funct).Deriv(nm,si,X,Y,v,l);
			pfap=dot(dfap,pkm);
			if((fap<=f+c1*alpha*pf)&(pfap>=c2*pf))
			{
				Wolfe=true;
			}else{
				alpha*=Tau;
			}
			count++;
		}while(!Wolfe & count <200);
		if(count>=200){alpha=0.01;}
		return alpha;
	
	}
    double Get_F(){return _F;}
  private:
    int _K;
    Funct *_funct;
    double _F;
};
