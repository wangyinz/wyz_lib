#ifndef _SYNTHETICSEISMOGRAM_H_
#define _SYNTHETICSEISMOGRAM_H_
#include "TimeSeries.h"
#include "ThreeComponentSeismogram.h"
namespace SEISPP {
using namespace std;
using namespace SEISPP;
/*! Abstract base class for a synthetic seismogram generator of any kinds.
Design is aimed to be generic enough that any synthetic generator can
be constructed as an implementation of this abstraction.  Most 
synthetic calculators will have a large amount of overhead in the constructors.
As a minimum a constructor for a generator will need to load a model.  It may
or may not need to do some initial calculations to allow one to call the
methods that return seismogram objects.  For example, a finite difference
synthetic program will do a ton of calculations in tne constructor of 
have to load a previously computed grid to do it's thing.  In contrast
a simple layered calculator will only need to load a 1d model and be
ready to go.  
*/
class SyntheticSeismogram
{
public:
    /*! Compute a single component synthetic.  

      Sometimes synthetics are naturally only single component things.
      e.g. if you want the pressure field for an acoustic approximation
      the output is naturally a scalar field.  This method should be used
      to return that type of synthetic.  

      \param hypo is an object that contains the source coordinates
      \param rlat receiver latitude (radians)
      \param rlon receiver longitude (radians)
      \param relev receiver elevation (km)
      \param type is a string used to sort out possible options for scalar
        output.  e.g. an implementation could use P for pressure, Z for
        vertical, R for radial, and T for transverse.  
      \param units can be optionally used to specify units of output.
      */
    virtual TimeSeries ComputeScalar(Hypocenter& hypo,
            double rlat, double rlon, double relev,string type)=0;
    /*! Similar to simpler method above, but clone parent. */
    virtual TimeSeries ComputeScalar(const TimeSeries& parent,
            Hypocenter& hypo,
            double rlat, double rlon, double relev,string type)=0;
    /*! Compute a three component synthetic.

      Sometimes we want a synthetic that computes a full three component
      seismogram.  This method should be implemented to do that.  

      \param hypo is an object that contains the source coordinates
      \param rlat receiver latitude (radians)
      \param rlon receiver longitude (radians)
      \param relev elevation (km)
      \param units can be optionally used to specify units of output.
      */
    virtual ThreeComponentSeismogram Compute3C(Hypocenter& hypo,
            double rlat, double rlon, double relev,string units)=0;
    /*! Similar to simpler method above, but clone parent. */
    virtual ThreeComponentSeismogram Compute3C(const ThreeComponentSeismogram& parent,
            Hypocenter& hypo,
            double rlat, double rlon, double relev,string units)=0;
};

}  // End SEISPP namespace encapsulation
#endif
