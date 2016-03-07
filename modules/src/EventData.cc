#include "EventData.hh"
#include "Message.hh"

void EventData::Print(int verbosity)
{
    if (verbosity >= 1)
    {
	Message m(INFO);
	m<<std::endl;
	m<<"************************************************************************"<<std::endl;
	m<<"************************* EVENT INFORMATION ****************************"<<std::endl;
	m<<"************************************************************************"<<std::endl;
	m<<"Run ID: "<<run_id<<std::endl;
	m<<"Event ID: "<<event_id<<std::endl;
	m<<"Status: "<<status<<std::endl;
	m<<"Trigger Count: "<<trigger_count<<std::endl;
	m<<"Time Stamp: "<<timestamp<<std::endl;
	m<<"Dt: "<<dt<<std::endl;
	m<<"Event Time: "<<event_time<<std::endl;
	m<<"N. Chans: "<<nchans<<std::endl;
	m<<"Saturated: "<<saturated<<std::endl;
	m<<"************************************************************************"<<std::endl;
    }
    if (verbosity == 2)
    {
	{
	    Message m(INFO);
	    m<<std::endl;
	    m<<"************************************************************************"<<std::endl;
	    m<<"************************** CHANNEL SUMMARY *****************************"<<std::endl;
	    m<<"************************************************************************"<<std::endl;
	    m<<"\t\tCh. Baseline\t\t\t\t\t"<<std::endl;
	    m<<"Ch.\tPulses\tFound\tMean\tInt.Min"<<std::endl;
	    m<<"************************************************************************"<<std::endl;
	    for (size_t ch = 0; ch < channels.size(); ch++)
	    {
		m<<std::setw(2)<<ch<<"\t"<<std::setw(2)<<channels[ch].npulses<<"\t"
		 <<std::setw(2)<<channels[ch].baseline.found_baseline<<"\t"<<std::setw(6)<< std::setprecision(4)<<channels[ch].baseline.mean<<"\t"	
		 <<std::setw(7)<< std::setprecision(5) <<channels[ch].integral_min<<"\t"
		 <<std::endl;
	    }
	    m<<"***********************************************************************"<<std::endl;

	    //Print details for sum channel
	    if (GetChannelByID(ChannelData::CH_SUM))
	    {
		m<<std::endl;
		m<<"************************************************************************"<<std::endl;
		m<<"*********************** SUM CHANNEL INFORMATION ************************"<<std::endl;
		m<<"************************************************************************"<<std::endl;
	    }
	}
	if (GetChannelByID(ChannelData::CH_SUM))
	    GetChannelByID(ChannelData::CH_SUM)->Print(verbosity);
    }
    if (verbosity == 3)
    {
	{
	    Message m(INFO);
	    m<<std::endl;
	    m<<"************************************************************************"<<std::endl;
	    m<<"************************* CHANNEL INFORMATION **************************"<<std::endl;
	    m<<"************************************************************************"<<std::endl;
	}
	for (size_t ch = 0; ch < channels.size(); ch++)
	{
	    channels[ch].Print(verbosity);
	}
    }
}
