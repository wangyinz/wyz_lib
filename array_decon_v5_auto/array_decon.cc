#include <vector>
#include <list>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <time.h>
#include <cmath>
#include <complex>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <omp.h>
#include "perf.h"
#include "DeconOperator.h"
#include "SeisppKeywords.h"
#include "TimeSeries.h"
#include "ThreeComponentSeismogram.h"
#include "ensemble.h"
#include "Metadata.h"
#include "Hypocenter.h"
#include "dbpp.h"
#include "filter++.h"
#include "resample.h"
#include "seispp.h"
#include "stack.h"
#include "TraceEditPlot.h"
#include "SignalToNoise.h"

/*
This the deconvolution program takes the stack of all vertical components
as the wavelet to deconvolve with each three components. Input and output
are all based on antelope database. 
All parameters are defined in the trace_decon.pf.
*/

using namespace std;
using namespace SEISPP;

bool SEISPP::SEISPP_verbose(true);
/*! \brief Generic algorithm to load arrival times from a database.

In passive array processing a very common need is to extract time windows
around marked phase pick or a theoretical arrival time.  For measured
arrival times the css database has an awkward link to waveforms that 
causes problems when dealing with continuous data.  This is one element of
a set of functions aimed at dealing with this problem.  

This particular procedure aims to take an input ensemble of data and 
find matches for arrival times in an external database.  For each member
of the ensemble with a matching arrival it posts the arrival time to the
generalized header (Metadata) of the parent object.  To avoid testing 
for valid arrivals in the ensemble after this procedure is run a index
of data with valid arrivals is returned as an stl list container of ints.

\param dat is the input ensemble passed as a reference to the vector 
	container of data objects.  The type of the objects in the container
	is generic except class T MUST be an object that inherits Metadata.
	At time of this writing T could be TimeSeries, ThreeComponentSeismogram,
	or ComplexTimeSeries. 
\param dbh match handle to a Datascope database containing arrival. 
	Normally this should reference the standard catalog view.  That
	is:  event->origin->assoc->arrival subsetted to orid==prefor.
	The algorithm actually only cares that a find using the Metadata
	area of the vector data components will provide a unique match
	into the table the handle references (normally arrival).  This 
	tacitly assumes the handle has been properly constructed and the
	proper attributes have been loaded with the data to provide a unique
	match into arrival.  With css this REQUIRES that each component of
	dat MUST have had evid and/or orid posted before calling this function.  There
	is no other unambiguous way to match waveforms in this context to
	a unique arrival. 
\param keyword is the attribute name used to posted the arrival time data.

\return list of members of input vector with valid arrival times posted.
	The contents of the list can be traversed and used with operator[]
	to extract valid data from the ensemble.  
*/
list<long> LoadArrivalTimes(vector<ThreeComponentSeismogram>& dat,
                DatascopeMatchHandle& dbh,
		    const string keyword)
{
        std::vector<ThreeComponentSeismogram>::iterator d;
	long i;
	list<long> data_with_arrivals;
	const string base_error("Warning (LoadArrivalTimes): ");
	for(d=dat.begin(),i=0;d!=dat.end();++d,++i)
	{
		double atime;  //  arrival time.  
		if(d->live)
		{
		// First see if there is an arrival for this
		// station.  If not, skip it. 
			list<long> records
				=dbh.find(dynamic_cast<Metadata&>(*d),false);
			// if no arrival silently skip data for this station
			if(records.size()<=0) continue;
			if(records.size()>1)
			{
				string sta=d->get_string("sta");
				cerr << base_error 
					<< "found "
					<< records.size()
					<< " arrivals for station "
					<< sta <<endl
					<< "Using first found "
					<< "in database view"<<endl;
			}
			Dbptr db=dbh.db;
			// tricky usage here.  begin() returns
			// an iterator so the * operator gets the
			// value = record number of the match
			db.record=*(records.begin());
			char csta[10];
			if(dbgetv(db,0,"arrival.time",&atime,"sta",csta,0)
				== dbINVALID) 
			{
				string sta=d->get_string("sta");
				cerr << base_error
					<< "dbgetv failed"
					<< " in attempt to obtain"
					<< " arrival time for station"
					<< sta << endl
					<< "Data from this station"
					<< " will be dropped"<<endl;	
			}
			d->put(keyword,atime);
			data_with_arrivals.push_back(i);
		}
	}
	return(data_with_arrivals);
} 
ThreeComponentEnsemble *BuildRegularGather(ThreeComponentEnsemble& raw,
	DatascopeMatchHandle& dbh, 
	ResamplingDefinitions& rdef,
	double target_dt,
	TimeWindow processing_window)
{
	const string arrival_keyword("arrival.time");
	const double samprate_tolerance(0.01);  // fractional sample rate tolerance
	int nmembers=raw.member.size();
	auto_ptr<TimeSeries> x1,x2,x3;
	ThreeComponentEnsemble *result;
	result = new ThreeComponentEnsemble(raw);
	// An inefficiency here, but this allow us to discard dead
	// traces and problem data from ensemble as we assemble the
	// new one.
	result->member.clear();
	result->member.reserve(raw.member.size());
	// Load arrivals from database.  List returned is index into raw of
	// data with valid arrivals loaded
	list<long> data_with_arrivals;
        data_with_arrivals=LoadArrivalTimes(raw.member,dbh,arrival_keyword);
	list<long>::iterator index;
	for(index=data_with_arrivals.begin();index!=data_with_arrivals.end();++index)
	{
		ThreeComponentSeismogram d=raw.member[*index];
		if(d.live)
		{
		try {
cout << d.get_string("sta")<<" has arrival time ="
	<<strtime(d.get_double("arrival.time"))<<endl;
			d.rotate_to_standard();	
			// partial clone used to hold result
			ThreeComponentSeismogram d3c(d);  
			x1=auto_ptr<TimeSeries>(ExtractComponent(d,0));
			x2=auto_ptr<TimeSeries>(ExtractComponent(d,1));
			x3=auto_ptr<TimeSeries>(ExtractComponent(d,2));
			// resample if necessary.  Using auto_ptr to avoid temporary pointer
			// and as good practice to avoid memory leaks
			if( (abs( (d.dt)-target_dt)/target_dt) > samprate_tolerance)
			{
				*x1=ResampleTimeSeries(*x1,rdef,target_dt,false);
				*x2=ResampleTimeSeries(*x2,rdef,target_dt,false);
				*x3=ResampleTimeSeries(*x3,rdef,target_dt,false);
			}
			// This procedure returns an auto_ptr.  An inconsistency in
			// SEISPP due to evolutionary development
			x1=ArrivalTimeReference(*x1,arrival_keyword,processing_window);
			x2=ArrivalTimeReference(*x2,arrival_keyword,processing_window);
			x3=ArrivalTimeReference(*x3,arrival_keyword,processing_window);
			// safer than using attribute ns in x1
			// assumes all three components are equal length,
			// which is pretty much guaranteed here
			int ns=x1->s.size();
			d3c.ns=ns;
			d3c.dt=x1->dt;
			d3c.t0=x1->t0;
			d3c.tref=x1->tref;
			d3c.u=dmatrix(3,ns);
			// Using blas here for speed
			dcopy(ns,&(x1->s[0]),1,d3c.u.get_address(0,0),3);
			dcopy(ns,&(x2->s[0]),1,d3c.u.get_address(1,0),3);
			dcopy(ns,&(x3->s[0]),1,d3c.u.get_address(2,0),3);
			result->member.push_back(d3c);
		} catch (SeisppError serr)
		{
			// Minor maintenance issue here.  Frozen name is
			// assumed to be in Metadata area.  Avoiding second
			// handler for possible MetadataError intentionally
			string sta=d.get_string("sta");
			cerr << "Problem assembling 3C seismogram for station "
				<< sta <<endl;
			serr.log_error();
			raw.member[*index].live=false;
			cerr << "Data for this station dropped"<<endl;
		}
		}
			
	}
	return(result);
}
void PostEvid(ThreeComponentEnsemble *d,int evid)
{
	vector<ThreeComponentSeismogram>::iterator dptr;
	for(dptr=d->member.begin();dptr!=d->member.end();++dptr)
		dptr->put("evid",evid);
}
/*! \brief Builds a standard catalog view from a CSS3.0 database.

Much passive array processing is built around the css3.0 schema.
The css3.0 schema has a standard view that defines the definitive
catalog for a network.  This is formed by the join of:
	event->origin->assoc->arrival
and is usually (as here) subsetted to only keep rows with
orid equal to the "preferred origin" (prefor).  This procedure
builds a handle to this database view.

\param dbh is a handle to the parent Datascope database 
	from which the view is to be constructed.  It need
	only reference the correct database. 
\return new handle that refers to the "standard" catalog view
	described above.
*/
DatascopeHandle StandardCatalogView(DatascopeHandle& dbh)
{
	DatascopeHandle result(dbh);
	dbh.lookup("event");
	dbh.natural_join("origin");
	string ss_to_prefor("orid==prefor");
	dbh.subset(ss_to_prefor);
	dbh.natural_join("assoc");
	dbh.natural_join("arrival");
	return(dbh);
}

void SaveResult(DatascopeHandle& dbh,
	ThreeComponentEnsemble* gather,
		AttributeMap& amo,
			MetadataList& mdlo,
				bool use_wfdisc,
					string dir,
						string dfile,
							vector<double> &all_meta,
								string dbout)
{
	const string wfdtable("wfdisc");
	const string wfptable("wfprocess");
	const string x0_name("IE");
	const string x1_name("IN");
	const string x2_name("IZ");
	DatascopeHandle dbsclink(dbh);
	DatascopeHandle dbevlink(dbh);
	DatascopeHandle dbdecon(dbh);
	
	if(use_wfdisc)
	{
		dbdecon.lookup("decon");
		vector<ThreeComponentSeismogram>::iterator d;
		int i;
		for(d=gather->member.begin(), i=0;d!=gather->member.end();++d,++i)
		{
			try {
				for(int j=0;j<3;j++)
				{
					TimeSeries *data=ExtractComponent(*d,j);
					switch(j)
					{
						case 0:
							data->put(string("chan"),x0_name);
						break;
						case 1:
							data->put(string("chan"),x1_name);
						break;
						case 2:
							data->put(string("chan"),x2_name);
						break;
					}
					data->put("dir",dir);
					data->put("dfile",dfile);
					int rnum=dbsave(*data,dbh.db,wfdtable,mdlo,amo);
					if(rnum!=-1)
					{
						rnum++;
                        dbdecon.append();
                        dbdecon.put("pwfid",rnum);
                        dbdecon.put("sta",d->get_string("sta"));
                        //cout<<"sta = "<<d->get_string("sta")<<endl;
                        switch(j)
			{
				case 0:
					dbdecon.put("chan",x0_name);
					//cout<<"rawsnr0 = "<<d->get_double("rawsnr0")<<endl;
					dbdecon.put("rawsnr",d->get_double("rawsnr0"));
				break;
				case 1:
					dbdecon.put("chan",x1_name);
					//cout<<"rawsnr1 = "<<d->get_double("rawsnr1")<<endl;
					dbdecon.put("rawsnr",d->get_double("rawsnr1"));
				break;
				case 2:
					dbdecon.put("chan",x2_name);
					//cout<<"rawsnr2 = "<<d->get_double("rawsnr2")<<endl;
					dbdecon.put("rawsnr",d->get_double("rawsnr2"));
				break;
			}
                        dbdecon.put("niteration",   int(all_meta[i*15+j*5]));
                        dbdecon.put("nspike",       int(all_meta[i*15+j*5+1]));
                        dbdecon.put("epsilon",      all_meta[i*15+j*5+2]);
                        dbdecon.put("peakamp",      all_meta[i*15+j*5+3]);
                        dbdecon.put("averamp",      all_meta[i*15+j*5+4]);
						//dbdecon.put("snr", 			1.0);
					}
					delete data;
				}
			} catch (SeisppError& serr) {
				string sta=d->get_string("sta");
				cerr << "Error saving station "<<sta<<endl;
				serr.log_error();
			}
		}
	}
	else
	{
		dbsclink.lookup("sclink");
		dbevlink.lookup("evlink");
		dbdecon.lookup("decon");
		vector<ThreeComponentSeismogram>::iterator d;
		int i;
		for(d=gather->member.begin(),i=0;d!=gather->member.end();++d,++i)
		{
			try {
				d->put("dir",dir);
				d->put("dfile",dfile);
				d->put("timetype",string("a"));
				d->put("wfprocess.algorithm",string("trace_decon"));
				int rnum;
				if((rnum=dbsave(*d,dbh.db,wfptable,mdlo,amo))!=-1)
				{
                    rnum++;
					for(int j=0;j<3;j++)
					{
                        dbdecon.append();
                        dbdecon.put("pwfid",rnum);
                        dbdecon.put("sta",d->get_string("sta"));
                        switch(j)
						{
							case 0:
                                dbdecon.put("chan",x0_name);
                                //cout<<"rawsnr0 = "<<d->get_double("rawsnr0")<<endl;
								dbdecon.put("rawsnr",d->get_double("rawsnr0"));
                                break;
							case 1:
                                dbdecon.put("chan",x1_name);
                                //cout<<"rawsnr1 = "<<d->get_double("rawsnr1")<<endl;
								dbdecon.put("rawsnr",d->get_double("rawsnr1"));
                                break;
							case 2:
                                dbdecon.put("chan",x2_name);
                                //cout<<"rawsnr2 = "<<d->get_double("rawsnr2")<<endl;
								dbdecon.put("rawsnr",d->get_double("rawsnr2"));
                                break;
						}
                        dbdecon.put("niteration",   int(all_meta[i*15+j*5]));
                        dbdecon.put("nspike",       int(all_meta[i*15+j*5+1]));
                        dbdecon.put("epsilon",      all_meta[i*15+j*5+2]);
                        dbdecon.put("peakamp",      all_meta[i*15+j*5+3]);
                        dbdecon.put("averamp",      all_meta[i*15+j*5+4]);
					}
					dbevlink.append();
					dbevlink.put("evid",d->get_int("evid"));
					dbevlink.put("pwfid",rnum);
					dbsclink.append();
					dbsclink.put("sta",d->get_string("sta"));
					dbsclink.put("chan","3C");
					dbsclink.put("pwfid",rnum);
				}
			} catch (SeisppError& serr) {
				string sta=d->get_string("sta");
				cerr << "Error saving station "<<sta<<endl;
				serr.log_error();
			}
		}
	}
}
void TestNaN(ThreeComponentSeismogram& dat)
{
	for(int i=0;i<dat.u.rows();i++)
		for(int j=0;j<dat.u.columns();j++)
		{
			if(std::isnan(dat.u(i,j)))
				{
					dat.live=false;
					cout<<"NaN found in "
						<<dat.get_string("sta")
						<<endl;
					return;
				}
		}
}
/*! \brief Multiple stage TimeInvariantFilter operator.

Sometimes it is necessary to apply a series of stages of
filtering to equalize different data sets or just to 
simply the process of defining a multiple filter chain.
This object simplifies that process by allowing the definition
of a chain of filters and a set of supplied methods to apply
these filters to data.
*/
class MultiStageFilter
{
public:
	/*! \brief Default constructors.  

	Loads a null filter definition.  That is the default is
	a do nothing operator. */
	MultiStageFilter();
	/*! \brief Construct the filter definitions from a string.

	This is currently the prime constructor.  A string is parsed
	into tokens that describe a series of filters that will be
	chained together.  The strings that define each individual
	filter type are parsed into blocks determined by the separator
	argument.  The parsed strings are currently used to construct
	a set of TimeInvariantFilter processing objects.
	An example helps explain how this would be used.  If we
	passed "DEMEAN; BW 0.5 5 2 5" and define ";" as the separator
	this would yield a two stage filter:  demean followed by a 
	0.5 to 2 Hz bandpass filter.

	\param filterlist contains the set of filter recipes to use.
	\param separator is the string used as a separator between
		the recipes for each filter description.
	*/
	MultiStageFilter(string filterlist,string separator);
	template <class T> void  apply(T& d);
private:
	list<TimeInvariantFilter> stages;
};
MultiStageFilter::MultiStageFilter()
{
	TimeInvariantFilter f(string("none"));
	stages.push_back(f);
}
MultiStageFilter::MultiStageFilter(string filterlist,string separator)
{
    try {
	const string white(" \t\n");
	int current=0,end_current;
	string stmp;
	string filterparam;
	// Strip any leading white space
	stmp=filterlist;
	if((current=stmp.find_first_not_of(white,0)) != 0)
	{
		stmp.erase(0,current);
		current=0;
	}
	int endstmp=stmp.size();
	do {
		end_current=stmp.find_first_of(separator,current);
		if(end_current<0)
		{
			filterparam.assign(stmp,current,endstmp);
			stages.push_back(TimeInvariantFilter(filterparam));
			break;
		}			
		filterparam.assign(stmp,current,end_current-current);
		stages.push_back(TimeInvariantFilter(filterparam));
		current=stmp.find_first_not_of(separator,end_current);
		current=stmp.find_first_not_of(white,current);
	}
	while(current<endstmp && current>=0);
    } catch (...) {throw;};
}
		
	
	
/* For now this only works on ensembles.  */
template <class T> void MultiStageFilter::apply(T& d)
{
    try {
	list<TimeInvariantFilter>::iterator filt;
	for(filt=stages.begin();filt!=stages.end();++filt)
		FilterEnsemble(d,*filt);
	/* Note this template could be made to work with TimeSeries
	 or ThreeComponentSeismogram components individually if 
	we used a typeid check on T and then used this line
	for that type of beast:  filt->apply(d);
	*/
    } catch (...) {throw;};
}
void ApplyFST(ThreeComponentEnsemble& e,Hypocenter& hypo)
{
	double vp0(6.0),vs0(3.5);
	for(int i=0;i<e.member.size();++i)
	{
		double lat,lon,elev;
		lat=e.member[i].get_double("lat");
		lon=e.member[i].get_double("lon");
		elev=e.member[i].get_double("elev");
		SlownessVector u=hypo.pslow(lat,lon,elev);
		e.member[i].free_surface_transformation(u,vp0,vs0);
	}
}
void ApplyFST(ThreeComponentSeismogram& e,Hypocenter& hypo)
{
	double vp0(6.0),vs0(3.5);
	double lat,lon,elev;
	lat=e.get_double("lat");
	lon=e.get_double("lon");
	elev=e.get_double("elev");
	SlownessVector u=hypo.pslow(lat,lon,elev);
	e.free_surface_transformation(u,vp0,vs0);
	
}
void ApplyKills(ThreeComponentEnsemble *d, set<int>& kills)
{
    set<int>::iterator kptr;
    int nmembers=d->member.size();
    for(kptr=kills.begin();kptr!=kills.end();++kptr)
    {
        int i=(*kptr);
        if((i<0) || (i>=nmembers) ) throw i;
        d->member[i].live=false;
    }
}
unsigned int nextPowerOf2(unsigned int n)
{
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n++;
    return n;
}
void usage()
{
	cerr << "array_decon db [-b beamdb -o outputdb -pf pfname]" << endl;
	exit(-1);
}

int main(int argc, char **argv)
{
	// This is a C++ thing.  Not required on all platforms, but
	// recommended by all the books.  
	ios::sync_with_stdio();
	// This is a standard C construct to crack the command line
	if(argc<2) usage();
	string pfin(argv[0]);
	string dbin(argv[1]);
	string beamdb=dbin;
	string dbout=dbin;
	for(int i=2;i<argc;++i)
	{
		if(!strcmp(argv[i],"-V"))
			usage();
		else if(!strcmp(argv[i],"-pf"))
		{
			++i;
			if(i>=argc) usage();
			pfin = string(argv[i]);
		}
		else if(!strcmp(argv[i],"-b"))
		{
			++i;
			if(i>=argc) usage();
			beamdb=string(argv[i]);
		}
		else if(!strcmp(argv[i],"-o"))
		{
			++i;
			if(i>=argc) usage();
			dbout=string(argv[i]);
		}
		else
			usage();
	}
	// Standard Antelope C construct to read a parameter file
	// Well almost standard.  const_cast is a C++ idiom required
	// due to a collision of constantness with Antelope libraries.
	Pf *pf;
        if(pfread(const_cast<char *>(pfin.c_str()),&pf)) 
	{
		cerr << "pfread error for pf file="<<pfin<<".pf"<<endl;
		exit(-1);
	}
	const string arrival_keyword("arrival.time");
	Metadata control(pf);
	MetadataList mdbeam=pfget_mdlist(pf,"Beam_mdlist");
	MetadataList mdens=pfget_mdlist(pf,"Ensemble_mdlist");
	MetadataList mdtrace=pfget_mdlist(pf,"station_mdlist");
	MetadataList mdlo=pfget_mdlist(pf,"output_mdlist");
	try {
		// Get set of control parameters
		string netname=control.get_string("netname");
		string phase=control.get_string("phase");
		double ts,te;
		double tpad=control.get_double("data_time_pad");
		ts=control.get_double("noise_window_start");
		te=control.get_double("noise_window_end");
		TimeWindow noise_twin(ts,te);
		StationChannelMap stachanmap(pf);
		string schemain=control.get_string("InputAttributeMap");
		string schemaout=control.get_string("OutputAttributeMap");
		double target_dt=control.get_double("target_sample_interval");
		control.put("shaping_wavelet_dt",target_dt);
		double proclen=control.get_double("processing_length");
		control.put("processing_sample_length",(int)(proclen/target_dt));
		double snr_signal_window_start=control.get_double("snr_signal_window_start");
		double snr_signal_window_end=control.get_double("snr_signal_window_end");
		double snr_noise_window_start=control.get_double("snr_noise_window_start");
		double snr_noise_window_end=control.get_double("snr_noise_window_end");
		TimeWindow noise(snr_noise_window_start,snr_noise_window_end),
			signal(snr_signal_window_start,snr_signal_window_end);
		ResamplingDefinitions rdef(pf);
		
		string output_dir=control.get_string("output_dir");
		bool use_wfdisc=control.get_bool("use_wfdisc");
		
		// First get all the database components assembled.
		// input data, output data
		DatascopeHandle dbh(dbin,false);
		/* Build a match handle for arrivals using the standard catalog
		view of the join of event->origin->assoc->arrival */
		AttributeMap am(schemain);  
		AttributeMap amo(schemaout);  
		DatascopeHandle dbcatalog=StandardCatalogView(dbh);
		list<string> matchkeys;
		matchkeys.push_back(string("sta"));
		matchkeys.push_back(string("evid"));
		DatascopeMatchHandle dbhm(dbcatalog,string(""),matchkeys,am);
		
		/*Now construct a handle to get beam traces */
		DatascopeHandle dbhbeam(dbh);
		if(beamdb!=dbin)
			dbhbeam=DatascopeHandle(beamdb,true);
		DatascopeHandle dbhv(dbhbeam,pf,
			string("dbprocess_commands_for_beams"));
		
		/* Define output database*/
		DatascopeHandle dbho(dbh);
		if(dbout!=dbin)
			dbho=DatascopeHandle(dbout,false);
		dbhv.rewind();  // Not essential, but good practice to initialize
		/* New stuff using seismic plot abstraction.   Three 
		   seperate windows will be generated.  One of beam, one to 
		   display deconvolution waveforms (expected to be inverse
		   filter and actual output along with a copy of the beam). 
		   May evolve */
		// This may need to be a TimeWindowPicker.  For now 
		// define it as a simple static plot.
		//Metadata beamdef(pf,string("beam_plot_properties"));
		//SeismicPlot beamplot(beamdef); 
		/* Need an editor window for data plotting.   Deconvolved
		   data will be replace original data. */
		//Metadata dpdef(pf,string("data_plot_properties"));
		//TraceEditPlot dataplot(dpdef);
		/* these raw pointers are used to hold dynamic ensembles
		   below.  Usual warning about pointers. */
		ThreeComponentEnsemble *rawdata=NULL, 
				*regular_gather=NULL,
				*decondata=NULL;
		SeismicArray *stations=NULL;
		/* Assorted variables needed in the loop below */
		double lat,lon,depth,otime;
		int evid;
		int record;
		string filter_param;
		ArrayDecon *decon_op;
		for(record=0;record<dbhv.number_tuples();++record,++dbhv)
		{
			// Read the beam trace and load it as beam
			TimeSeries beam(dynamic_cast<DatabaseHandle&>(dbhv),
				mdbeam,am);
			try{
				if( (abs( (beam.dt)-target_dt)/target_dt) > 0.01)
					beam=ResampleTimeSeries(beam,rdef,target_dt,false);
			}catch(SeisppError& serr)
			{
				cout<<"problem with resampling beam"<<endl;
				serr.log_error();
				continue;
			}
			double beamnorm=fabs(beam.s[0]);
			for(int i=1;i<beam.s.size();i++)
				if(fabs(beam.s[i])>beamnorm)
					beamnorm=fabs(beam.s[i]);
			//sqrt(ddot(beam.s.size(),&(beam.s[0]),1,&(beam.s[0]),1));
			dscal(beam.s.size(), 1.0/beamnorm, &(beam.s[0]), 1);
			
			int time_shift=-dbhv.get_double(string("wfprocess.time"))/target_dt;
			control.put("sample_shift",time_shift);
			int num_gather;
			lat=beam.get_double(string("origin.lat"));
			lon=beam.get_double(string("origin.lon"));
			depth=beam.get_double(string("origin.depth"));
			otime=beam.get_double(string("origin.time"));
			evid=beam.get_int(string("evid"));
			// origin coordinates are degrees in the db,
			// but must be radians for hypocenter object
			lat=rad(lat);
			lon=rad(lon);
			
			Hypocenter hypo(lat,lon,depth,otime,
				string("tttaup"),string("iasp91"));
			// Define the filter
			try {
				filter_param=beam.get_string("xcorbeam.filter");
			} catch (MetadataGetError mdge) {
				cerr	<<"filter attribute not defined"
					<<" using default of DEMEAN"<<endl;
				filter_param=string("DEMEAN");
			}
			MultiStageFilter filt(filter_param,string(";"));
			// On the first record we need to load the station
			// geometry object
			if(record==0)
			{
				stations = new SeismicArray(
					dynamic_cast<DatabaseHandle&>(dbh),
					hypo.time,netname);
			}
			else
			{
				TimeWindow test(hypo.time,hypo.time+2000.0);
				if(!stations->GeometryIsValid(test))
				{
					delete stations;
					stations = new SeismicArray(
					dynamic_cast<DatabaseHandle&>(dbh),
					hypo.time,netname);
				}
			}
			ts=dbhv.get_double(string("wfprocess.time"));
			te=round((dbhv.get_int(string("wfprocess.nsamp"))-1)/dbhv.get_double(string("wfprocess.samprate")))+ts;
            if(te<proclen)
                te=proclen;
			TimeWindow datatwin(ts,te);
			TimeWindow processingtwin(ts,proclen);
			// Read the raw data using the time window based constructor
			try{
			rawdata=array_get_data(*stations,hypo,
				phase,datatwin,tpad,dynamic_cast<DatabaseHandle&>(dbh),
				stachanmap,mdens,mdtrace,am);
			}catch (SeisppError& serr)
			{
				cout<<"bad event"<<endl;
				serr.log_error();
				continue;
			}
			if(rawdata->member.size()==0)
			{
				cout<<"bad event with problem wf data"<<endl;
				continue;
			}
			try{
				// Filter the raw data.
				filt.apply(*rawdata);
				PostEvid(rawdata,evid);
		
				num_gather=rawdata->member.size();
				for(int i=0;i<num_gather;i++)
				{
					double lat=stations->array[rawdata->member[i].get_string("sta")].lat;
					double lon=stations->array[rawdata->member[i].get_string("sta")].lon;
					double elev=stations->array[rawdata->member[i].get_string("sta")].elev;
					rawdata->member[i].put("lat",lat);
					rawdata->member[i].put("lon",lon);
					rawdata->member[i].put("elev",elev);
					//cout<<rawdata->member[i].get_string("sta")<<":"<<lat<<":"<<lon<<":"<<elev<<endl;
				}

				regular_gather=BuildRegularGather(*rawdata, dbhm,rdef,target_dt,
						datatwin);
						
				//this may not be a good choice, but works fine for most of the times.
				for(int i=0;i<num_gather;i++)
				{
					if(regular_gather->member[i].ns<time_shift)
					{
						cout<<"bad event with problem waveform that's too short"<<endl;
						throw i;
					}
				}
			            
				if(control.get_bool("apply_free_surface_transform"))
					ApplyFST(*regular_gather,hypo);
				
				// regular gather is now assumed to contain data with
				// a common start time.  
				regular_gather->put(string("moveout"),0.0);
				//dataplot.plot(*regular_gather,true);
			
				//string filt_str;
				//cout<<"change filter parameters?"<<endl;
				//std::getline(cin,filt_str);
				//if(filt_str=="no")
				//	break;
				//else
				//{
				//	delete filt;
				//	filt=new MultiStageFilter(filt_str,string(";"));
				//	delete rawdata_filt;
				//}
			}catch (...)
			{
				cout<<"bad event with problem waveform"<<endl;
				//serr.log_error();
				continue;
			}
            
			num_gather=regular_gather->member.size();
			int maxns=0;
			for(int i=0;i<num_gather;i++)
			{
				if(regular_gather->member[i].ns>maxns)
					maxns=regular_gather->member[i].ns;
			}
			if(maxns<=0)
			{
				cout<<"bad event that has a number of sample eq or lt 0"<<endl;
				continue;
			}
			
			control.put("operator_nfft",static_cast<int>(nextPowerOf2(maxns)));
			control.put("taper_length",maxns);
			string decon_type=control.get_string("deconvolution_type");
			
			if(decon_type=="least_square")
				decon_op=new ArrayLeastSquareDecon(control);
			else if(decon_type=="water_level")
				decon_op=new ArrayWaterLevelDecon(control);
			else if(decon_type=="multi_taper")
				decon_op=new ArrayMultiTaperDecon(control);
			else if(decon_type=="iterative")
				decon_op=new ArrayGeneralIterDecon(control);
			else
			{
				cout<<"WARNING: deconvolution_type is unavailable, set to least_square by default."<<endl;
				decon_type="least_square";
				decon_op=new ArrayLeastSquareDecon(control);
			}
			
			bool multi_taper_on=false;
			if(decon_type=="multi_taper" ||
				(decon_type=="iterative" && control.get_string("iterative_decon_type")=="multi_taper"))
				multi_taper_on=true;
			
			vector<double> noise_vector;
			if(multi_taper_on)
			{
				auto_ptr<TimeSeries> noise;
				TimeSeries beamtmp(beam);
				beamtmp.tref=absolute;
				noise=ArrivalTimeReference(beamtmp,arrival_keyword,noise_twin);
				noise_vector=noise->s;
				if(decon_type=="multi_taper")
					dynamic_cast<ArrayMultiTaperDecon *>(decon_op)->loadnoise(noise_vector);
				else
					dynamic_cast<ArrayGeneralIterDecon *>(decon_op)->loadnoise(noise_vector);
			}
			
			//auto_ptr<TimeSeries> beam_tmp;
			//beam_tmp=ArrivalTimeReference(beam_stack.stack,arrival_keyword,beam_twin);
			//beamplot.plot(*beam_tmp,false);
			//beam_tmp.release();
			
			// Deconvolution operator called here
			// set decondata
			decon_op->load(beam.s);
			
			/*// Show the actual output
			decon_op->loaddata(beam.s);
			TimeSeries actual_output(beam);
			actual_output.s=decon_op->getresult();
			beamplot.plot(actual_output,false);*/
			
			decondata=new ThreeComponentEnsemble(*regular_gather);
			//omp_set_num_threads(4);
			num_gather=regular_gather->member.size();
			vector<vector<double> > all_data;
			vector<double> all_scale;
			all_data.reserve(num_gather*3);
			for(int i=0;i<num_gather;i++)
			{
				for(int j=0;j<3;j++)
				{
					TimeSeries *data=ExtractComponent(regular_gather->member.at(i),j);
					all_data.push_back(data->s);
					if(j==2)
					{
						int len=beam.s.size()<=data->s.size()?beam.s.size():data->s.size();
						double datamp=ddot(len,&(beam.s[0]),1,&(data->s[0]),1);
						double denom=sqrt(ddot(len,&(beam.s[0]),1,&(beam.s[0]),1));
						all_scale.push_back((datamp/denom));
						all_scale.push_back((datamp/denom));
						all_scale.push_back((datamp/denom));
					}
					delete data;
				}
			}
			for(int i=0;i<num_gather*3;i++)
				for(int j=0;j<all_data[i].size();j++)
					all_data[i][j]/=all_scale[i];
			vector<double> all_meta;
			if(decon_type=="iterative")
			{
				all_meta.reserve(num_gather*3*5);
				all_meta.resize(num_gather*3*5);
			}
			
			vector<vector<double> > all_result;
			all_result.reserve(num_gather*3);
			all_result.resize(num_gather*3);
			int paralleltime=num_gather*3;
			#pragma omp parallel
			{
			#pragma omp for schedule(guided, 1)
				for(int i=0;i<paralleltime;i++)
				{
					ArrayDecon *para_decon_op;
					if(decon_type=="least_square")
						para_decon_op=new ArrayLeastSquareDecon(*dynamic_cast<ArrayLeastSquareDecon *>(decon_op));
					else if(decon_type=="water_level")
						para_decon_op=new ArrayWaterLevelDecon(*dynamic_cast<ArrayWaterLevelDecon *>(decon_op));
					else if(decon_type=="multi_taper")
						para_decon_op=new ArrayMultiTaperDecon(*dynamic_cast<ArrayMultiTaperDecon *>(decon_op));
					else if(decon_type=="iterative")
						para_decon_op=new ArrayGeneralIterDecon(*dynamic_cast<ArrayGeneralIterDecon *>(decon_op));
					else
						para_decon_op=new ArrayLeastSquareDecon(*dynamic_cast<ArrayLeastSquareDecon *>(decon_op));
					para_decon_op->loaddata(all_data[i]);
					all_result[i]=para_decon_op->getresult();
					if(decon_type=="iterative")
					{
						int spikecount=0;
						double peakamp=0;
						double rms=0;
						for(int j=0;j<all_result[i].size();j++)
						{
							if(abs(all_result[i][j])>1e-15)
								spikecount++;
							if(abs(all_result[i][j])>peakamp)
								peakamp=abs(all_result[i][j]);
							rms+=all_result[i][j]*all_result[i][j];
						}
						all_meta[i*5]=dynamic_cast<ArrayGeneralIterDecon *>(para_decon_op)->numberiter()+1;
						all_meta[i*5+1]=spikecount;
						all_meta[i*5+2]=dynamic_cast<ArrayGeneralIterDecon *>(para_decon_op)->epsilon();
						all_meta[i*5+3]=peakamp;
						all_meta[i*5+4]=sqrt(rms/all_result[i].size());
					}
					cout<<"THREAD:"<<omp_get_thread_num()<<" processed "
						<<regular_gather->member.at(i/3).get_string("sta")
						<<":"<<i%3
						<<endl;
					delete para_decon_op;
				}
			}
			
			delete rawdata;
			
			for(int i=0;i<num_gather;i++)
			{
				vector<TimeSeries> decon_result_vector;
				for(int j=0;j<3;j++)
				{
					double snr_rms;
					TimeSeries *data=ExtractComponent(regular_gather->member.at(i),j);
					//get SNR before replacing the raw data with the deconvolved data. XIAOTAO YANG
					//cout<<"t0 in timeseries = "<<data->t0<<endl;
					snr_rms=SNR_rms(*data,signal,noise);
					//data->put(SEISPP::snr_keyword,snr_rms);

					data->s=all_result[i*3+j];
					double arrtime=regular_gather->member.at(i).get_double(arrival_keyword);
					data->rtoa(arrtime);
					data->put(arrival_keyword,arrtime);
					auto_ptr<TimeSeries> datacut=ArrivalTimeReference(*data,arrival_keyword,processingtwin);
					decon_result_vector.push_back(*datacut);
					//CHECK SNR VALUES. xiaotao yang
					//DEBUG
					/*
					cout<<"sta = "<<data->get_string("sta")<<", comp = "
						<<j<<", snr = "
						<<decon_result_vector[j].get_double(SEISPP::snr_keyword)<<endl;
					*/
					//put snr to decondata. Xiaotao Yang
					switch(j)
					{
						case 0:
							decondata->member[i].put("rawsnr0",snr_rms);
							break;
						case 1:
							decondata->member[i].put("rawsnr1",snr_rms);
							break;
						case 2:
							decondata->member[i].put("rawsnr2",snr_rms);
							break;
					}
					
					delete data;
				}
				ThreeComponentSeismogram decon_result_seismogram(decon_result_vector,2);
				decondata->member[i].u=decon_result_seismogram.u;
				decondata->member[i].ns=decon_result_seismogram.ns;
				TestNaN(decondata->member[i]);
				decondata->member[i].rotate_to_standard();
			}
			/*cout << "Edit deconvolved data to remove bad traces"
			  <<endl
			  <<"Type X in the data window or use the menu to continue"
			  <<endl;
			dataplot.plot(*decondata,true);
			set<int> kill_list=dataplot.report_kills();
			try {
				ApplyKills(decondata,kill_list);
			}catch(int ierr)
			{
				cerr << "Fatal problem with TraceEditPlot return"
					<<endl
					<< "Editor returned invalid kill index of "
					<< ierr<<endl
					<< "Must be in range 0 to "
					<< decondata->member.size()-1<<endl;
				exit(-1);
			}*/
			
			//convert time reference to absolute time
			for(int i=0;i<num_gather;i++)
			{
				decondata->member.at(i).rtoa(regular_gather->member.at(i).get_double(arrival_keyword));
			}
			
			stringstream ss;
			ss << evid;
			string dfile = ss.str();
			SaveResult(dbho,decondata,amo,mdlo,use_wfdisc,output_dir,decon_type+"_"+dfile,all_meta,dbout);
			delete regular_gather;
			delete decondata;
			delete decon_op;
		}
	}
	catch (SeisppError& serr)
	{
		serr.log_error();
	}
	catch (...)
	{
		cerr << "Something threw an unexpected exception"<<endl;
	}
}

