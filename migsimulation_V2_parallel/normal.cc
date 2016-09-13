#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions/normal.hpp>
using namespace std;
using namespace boost;
/* Copied form http://en.wikibooks.org/wiki/C++_Programming/Libraries/Boost */
/* Returns a random number from normal distribution with mean of "mean" and
   standard deviation sigma.  Uses boost random number generator that 
   is doubly random because of use of current time to generate seed to
   random number generator */

int SampleNormalAddN(double* randN, double mean, double sigma, double N)
{
    /* Generate N random number from the current time */
    // the memory should be preallocated for *randN!!!.
    static mt19937 rng(static_cast<unsigned> (std::time(0)));
    /* This defines a normal distribution to generate numbers */
    normal_distribution<double> norm_dist(mean,sigma);
    // bind random number generator to distribution to forma function object
    variate_generator<mt19937&, normal_distribution<double> >
        normal_sampler(rng, norm_dist);
    for(int i=0; i<N;i++) { 
        randN[i]+=normal_sampler();
    }
   // return(normal_sampler());
    return 0;
}

int SampleNormalN(double* randN, double mean, double sigma, double N)
{
    /* Generate N random number from the current time */
    // the memory should be preallocated for *randN!!!.
    static mt19937 rng(static_cast<unsigned> (std::time(0)));
    /* This defines a normal distribution to generate numbers */
    normal_distribution<double> norm_dist(mean,sigma);
    // bind random number generator to distribution to forma function object
    variate_generator<mt19937&, normal_distribution<double> >
        normal_sampler(rng, norm_dist);
    for(int i=0; i<N;i++) {
	randN[i]=normal_sampler();
    }
   // return(normal_sampler());
    return 0;
}

double SampleNormal(double mean, double sigma)
{
    /* Generate a random number from the current time */
    static mt19937 rng(static_cast<unsigned> (std::time(0)));
    /* This defines a normal distribution to generate numbers */
    normal_distribution<double> norm_dist(mean,sigma);
    // bind random number generator to distribution to forma function object
    variate_generator<mt19937&, normal_distribution<double> >
        normal_sampler(rng, norm_dist);
    return(normal_sampler());
}

