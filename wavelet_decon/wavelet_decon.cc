#include <vector>
#include <list>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <time.h>
#include <cmath>
#include <complex>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <omp.h>

#include "DeconOperator.h"
#include "SeisppKeywords.h"
#include "TimeSeries.h"
#include "ThreeComponentSeismogram.h"
#include "ensemble.h"
#include "Metadata.h"
#include "dbpp.h"
#include "seispp.h"
#include "resample.h"

using namespace std;
using namespace SEISPP;

bool SEISPP::SEISPP_verbose(true);

void SaveResult(DatascopeHandle& dbh,
	TimeSeries& dcdata,
		AttributeMap& amo,
			MetadataList& mdlo,
				string dir,
					string dfile)
{
	const string wfdtable("wfdisc");
	dcdata.put("dir",dir);
	dcdata.put("dfile",dfile);
	dbsave(dcdata,dbh.db,wfdtable,mdlo,amo);
}
void TestNaN(TimeSeries& dat)
{
	for(int i=0;i<dat.s.size();i++)
		if(std::isnan(dat.s[i]))
			{
				dat.s[i]=0;
			}
	return;
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
	cerr << "wavelet_decon db -o outputdb [-n start_number -pf pfname]" << endl;
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
	string dbout=dbin;
    string stnumstring("NULL");
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
		else if(!strcmp(argv[i],"-o"))
		{
			++i;
			if(i>=argc) usage();
			dbout=string(argv[i]);
		}
		else if(!strcmp(argv[i],"-n"))
		{
			++i;
			if(i>=argc) usage();
			stnumstring=string(argv[i]);
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
	Metadata control(pf);
	MetadataList mdl=pfget_mdlist(pf,"input_mdlist");
	MetadataList mdlo=pfget_mdlist(pf,"output_mdlist");
	try {
		// Get set of control parameters
		string schemain=control.get_string("InputAttributeMap");
		string schemaout=control.get_string("OutputAttributeMap");
		double gau_width=control.get_double("gaussian_width");
		long grouprecord=control.get_long("records_in_group");
		int time_shift=0;
		control.put("sample_shift",time_shift);
		
		string output_dir=control.get_string("output_dir");
		
        bool resample_flag=control.get_bool("resample_flag");
        ResamplingDefinitions rdef;
        double rdt;
        if(resample_flag)
        {
            ResamplingDefinitions rdeftemp(pf);
            rdef=rdeftemp;
            rdt=control.get_double("target_sample_interval");
        }
        
		// First get all the database components assembled.
		// input data, output data
		DatascopeHandle dbh(dbin,false);
		AttributeMap am(schemain);  
		AttributeMap amo(schemaout);
        
		/* Build a view to drive the process*/
		DatascopeHandle dbhv(dbh);
		dbhv.lookup(string("wfdisc"));
		
		/* Define output database*/
		DatascopeHandle dbho(dbh);
		if(dbout!=dbin)
			dbho=DatascopeHandle(dbout,false);
		dbho.lookup(string("wfdisc"));
		
		control.put("deconvolution_type",string("iterative"));
		control.put("iterative_decon_type",string("least_square"));
		control.put("shaping_wavelet_frequency_for_inverse",double(0.0));
		 
		
		/* Assorted variables needed in the loop below */
		dbhv.rewind();
		if(stnumstring.compare("NULL")!=0)
		{
			long idst=atol(stnumstring.c_str());
			for(int i=0;i<dbhv.number_tuples();++i,++dbhv)
			{
				if(i==idst)
					break;
			}
		}
		long recordnum=dbhv.number_tuples()-dbhv.current_record();
		for(long k=0;k<recordnum;k+=grouprecord)
		{
			vector<TimeSeries> datagroup;
			vector<TimeSeries> resultgroup;
			vector<string> dfilegroup;
			datagroup.reserve(grouprecord);
			for (long record=k;record<((k+grouprecord<recordnum)?k+grouprecord:recordnum);++record,++dbhv)
			{
				try{
					TimeSeries d(dbhv,mdl,am);
                    if(resample_flag)
                    {
                        try{
                            if( (abs( (d.dt)-rdt)/rdt) > 0.01)
                            {
                                d=ResampleTimeSeries(d,rdef,rdt,false);
                            }
                        }catch(...)
                        {
                            cerr<<"Failed resample record No."<<record<<endl;
                            continue;
                        }
                    }
					datagroup.push_back(d);
				}catch(SeisppError& serr)
				{
					cerr<<serr.message<<endl;
					cerr<<"Can not read record No."<<record<<endl;
					continue;
				}
			}
			resultgroup.reserve(datagroup.size());
			resultgroup.resize(datagroup.size());
			dfilegroup.reserve(datagroup.size());
			dfilegroup.resize(datagroup.size());
			#pragma omp parallel
			{
				#pragma omp for schedule(guided, 1)
				for(long record=0;record<datagroup.size();record++)
				{
					Metadata controlpara(control);
					try{
						SimpleDecon *decon_op;
			
						int nfft=static_cast<int>(nextPowerOf2(datagroup[record].get_int(string("nsamp"))));
						controlpara.put("operator_nfft",nfft);
						double dt=1.0/datagroup[record].get_double(string("samprate"));
						controlpara.put("shaping_wavelet_dt",dt);
			
						decon_op=new SimpleGeneralIterDecon(controlpara);
			
						double *w=gaussian(gau_width,dt,nfft);
						vector<double> wavelet;
						wavelet.reserve(nfft);
						for(int i=0;i<nfft;i++)
							wavelet.push_back(w[i]);
						TimeSeries decondata(datagroup[record]);
					
						decon_op->load(wavelet,datagroup[record].s);
						decondata.s=decon_op->getresult();
			
			
						delete [] w;
						delete decon_op;
			
						TestNaN(decondata);
						resultgroup[record]=decondata;
						dfilegroup[record]=datagroup[record].get_string(string("dfile"));
						
					}catch (SeisppError& serr)
					{
						cerr<<"bad trace at No."<<record+k<<endl;
						continue;
					}
				}
			}
			for(long record=0;record<resultgroup.size();record++)
			{
				SaveResult(dbho,resultgroup[record],amo,mdlo,output_dir,dfilegroup[record]+"c");
			}
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

