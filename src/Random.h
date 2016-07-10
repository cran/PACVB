#ifndef RANDOM_HPP
#define RANDOM_HPP
namespace RandomG{
	
using namespace std;
using namespace arma;
//using namespace boost::random;

inline int getSeed()
{
   // ifstream rand("/dev/urandom");
   /* GetRNGstate
    char tmp[sizeof(int)];
    rand.read(tmp,sizeof(int));
    rand.close();
    int* number = reinterpret_cast<int*>(tmp);
    return (*number);*/
    int output=R::runif(0,INT_MAX);
    return output;
}

template<class distribution>
class Random 
{
	public:
	Random(distribution nd)
	{
		boost::random::mt19937 randgen(getSeed());
		rng=new variate_generator<boost::random::mt19937,distribution >(randgen, nd);	
	}
	Random<distribution>& operator=(const Random<distribution>& X)
	{
		rng=X.rng;
	}
	virtual ~Random()
	{
		delete rng;
	}
	double operator()(void)
	{
		double x;
//		#pragma omp critical
//		{
			 x=(*rng)();
//		}
		return x;
	}
	variate_generator<boost::random::mt19937,distribution > *rng;
};



}
#endif

