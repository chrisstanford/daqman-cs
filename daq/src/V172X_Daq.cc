/* --------------------------------------------------------
CAEN_V172XDAQ.cc
This is the implementation for the CAEN_V172XDAQ class, 
which inherits from the WARP_VetoDAQ class.
----------------------------------------------------------*/

#include "V172X_Daq.hh"
//#include "CAEN_V172XEvent.hh"
#include "CAENVMElib.h"
#include "RawEvent.hh"
#include "Message.hh"
#include "ConfigHandler.hh"
#include "EventHandler.hh"
#include <string>
#include <time.h>
#include <bitset>
#include <algorithm>
#include "boost/ref.hpp"
#include "boost/timer.hpp"
#include "boost/date_time/posix_time/posix_time_duration.hpp"
#include <sstream>

//declare some useful constants
const int event_size_padding = 8;

enum VME_REGISTERS{
  VME_ChZSThresh =         0x1024,
  VME_ChZSNsamples =       0x1028,
  VME_ChTrigThresh =       0x1080,
  VME_ChTrigSamples =      0x1084,
  VME_ChStatus =           0x1088,
  VME_ChBuffersFull =      0x1094,
  VME_ChDAC =              0x1098,
  VME_ChannelsConfig =     0x8000,
  VME_BufferCode =         0x800C,
  VME_BufferFree =         0x8010,
  VME_CustomSize =         0x8020,
  VME_AcquisitionControl = 0x8100,
  VME_AcquisitionStatus =  0x8104,
  VME_SWTrigger =          0x8108,
  VME_TrigSourceMask =     0x810C,
  VME_TrigOutMask =        0x8110,
  VME_PostTriggerSetting = 0x8114,
  VME_FrontPanelIO =       0x811C,
  VME_ChannelMask =        0x8120,
  VME_DownsampleFactor =   0x8128,
  VME_EventsStored =       0x812C,
  VME_BoardInfo =          0x8140,
  VME_EventSize =          0x814C,
  VME_VMEControl =         0xEF00,
  VME_VMEStatus =          0xEF04,
  VME_BoardID =            0xEF08,
  VME_RelocationAddress =  0xEF10,
  VME_InterruptID =        0xEF14,
  VME_InterruptOnEvent =   0xEF18,
  VME_BLTEvents =          0xEF1C,
  VME_SWReset =            0xEF24,
  VME_SWClear =            0xEF28,
};

V172X_Daq::V172X_Daq() : BaseDaq(), _initialized(false), 
			 _params(), _triggers(0), _vme_mutex()
{
  ConfigHandler::GetInstance()->RegisterParameter(_params.GetDefaultKey(),
						  _params);
  _handle_vme_bridge = 0;
  std::fill (_handle_board, _handle_board + _params.nboards, 0);
}

V172X_Daq::~V172X_Daq()
{
  if (_params.vme_bridge_link >= 0) CAENVME_End(_handle_vme_bridge);

  for(int i=0; i < _params.nboards; i++)
    if(_params.board[i].enabled && _params.board[i].link > 0 && _params.board[i].link != _params.vme_bridge_link) CAENVME_End(_handle_board[i]);
}

int init_link (int link, int board, bool usb, int32_t *handle) {
  CVErrorCodes err = CAENVME_Init(usb ? cvV1718 : cvV2718, link, board, handle);
  if(err != cvSuccess){
    Message m(ERROR);
    m<<"Unable to initialize CAEN VME bridge for link " << link;
    if(usb) m<<" on USB ";
    m<< ": "<<std::endl;
    m<<"\t"<<CAENVME_DecodeError(err)<<std::endl;
    return -1;
  }
  char message[100];
  Message(DEBUG)<<"CAEN VME bridge successfully initialized for link " 
		<< link << "!"<<std::endl;
  CAENVME_BoardFWRelease(*handle,message);
  Message(DEBUG)<<"\tFirmware Release: "<<message<<std::endl;
  CAENVME_DriverRelease(*handle,message);
  Message(DEBUG)<<"\tDriver Release: "<<message<<std::endl;
  Message(INFO)<< "Link " << link << " initialized on handle " 
	       << *handle <<std::endl;
  return 0;
}


int V172X_Daq::Initialize()
{
  if(_initialized){
    Message(WARNING)<<"Reinitializing V172X_Daq..."<<std::endl;
    if (_params.vme_bridge_link >= 0) CAENVME_End(_handle_vme_bridge);

    for(int i=0; i < _params.nboards; i++)
      if(_params.board[i].enabled && _params.board[i].link > 0 && _params.board[i].link != _params.vme_bridge_link) CAENVME_End(_handle_board[i]);
  }
      
  if (_params.vme_bridge_link >= 0) {
    if (init_link (_params.vme_bridge_link, 0, false, &_handle_vme_bridge) < 0){
      _status=INIT_FAILURE;
      return -1;
    }
    //reset everything
    //CAENVME_SystemReset(_handle_vme_bridge);
    //CAENVME_DeviceReset(_handle_vme_bridge);
    if(_params.send_start_pulse) {
      CAENVME_ClearOutputRegister(_handle_vme_bridge, 0xFFFF);
      for(int line = 0; line<5; line++){
	CVErrorCodes err = CAENVME_SetOutputConf(_handle_vme_bridge, (CVOutputSelect)line, cvDirect, 
	    cvActiveHigh, cvManualSW);
	if(err != cvSuccess){
	  Message(ERROR)<<"Unable to configure the V2718 output registers.\n";
	  return -2;
	}
      }
    }
  } 
  else if (_params.send_start_pulse) {
    Message(WARNING)<<"send_start_pulse enabled but no bridge link present, disabling pulses" << std::endl;
    _params.send_start_pulse = 0;
  }

  for(int i=0; i < _params.nboards; i++) {
    if(_params.board[i].enabled && 
       ( ( _params.board[i].link >= 0 && 
	 _params.board[i].link != _params.vme_bridge_link )
	 || _params.board[i].usb ) ){
      if (init_link (_params.board[i].link, _params.board[i].chainindex, 
		     _params.board[i].usb, _handle_board + i) < 0) {
	_status=INIT_FAILURE;
	return -3;
      }
    } 
    else _handle_board[i] = _handle_vme_bridge;
  }

  for(int i=0; i < _params.nboards; i++){
    //Initialize each enabled board
    if(!_params.board[i].enabled) 
      continue;
    try{
      if(InitializeBoard(i)){
	_status = INIT_FAILURE;
	return -5;
      }
      
    }
    catch(std::exception &error)
      {
	Message(EXCEPTION)<<error.what()<<std::endl;
	return -6;
      }
  }
  _initialized = true;
  return Update();
}

int V172X_Daq::InitializeBoard(int boardnum)
{
  V172X_BoardParams& board = _params.board[boardnum];
  WriteVMERegister(board.address + VME_SWReset, 1, _handle_board[boardnum]);
  uint32_t data = ReadVMERegister(board.address+VME_BoardInfo, _handle_board[boardnum]);
  board.board_type = (BOARD_TYPE)(data&0xFF);
  board.nchans = (data>>16)&0xFF;
  board.mem_size = (data>>8)&0xFF;
  if(board.UpdateBoardSpecificVariables()){ //returns -1 on error
    Message(CRITICAL)<<"Board "<<boardnum<<" with address "
		     <<std::hex<<std::showbase
		     <<board.address<<std::dec<<std::noshowbase
		     <<" is not a V172X digitizer!"<<std::endl;
    _status = INIT_FAILURE;
    return -2;
  }
  WriteVMERegister(board.address+VME_BoardID,
		   board.id, _handle_board[boardnum]);
  WriteVMERegister(board.address+VME_InterruptID,
		   board.id, _handle_board[boardnum]);
  return 0;
}

int V172X_Daq::Update()
{
  if(!_initialized){
    Message(ERROR)<<"Attempt to update parameters before initializations!"
		  <<std::endl;
    return -1;
  }
  if(_is_running){
    Message(ERROR)<<"Attempt to update parameters while in run."
		  <<std::endl;
    return -2;
  }
  try{
    for(int iboard=0; iboard < _params.nboards; iboard++){
      V172X_BoardParams& board = _params.board[iboard];
      //downsample factor is not implemented in modern firmware
      board.downsample_factor=1; 
      if(!board.enabled) 
	continue;
      //determine the trigger acquisition window for the database
      //WARNING: Assumes it is the same for all boards!!!
      runinfo* info = EventHandler::GetInstance()->GetRunInfo();
      if(info){
	std::stringstream ss;
	ss<<"board"<<iboard<<".pre_trigger_time_us";
	info->SetMetadata(ss.str(), board.pre_trigger_time_us);
	ss.str("");
	ss<<"board"<<iboard<<".post_trigger_time_us";
	info->SetMetadata(ss.str(), board.post_trigger_time_us);
	
      }
       
      uint32_t channel_mask = 0;
      uint32_t trigger_mask = 0;
      uint32_t trigger_out_mask = 0;
      //need to know total_nsamps to estimate event size
      //do the per-channel stuff
      for(int i=0; i<board.nchans; i++){
	V172X_ChannelParams& channel = board.channel[i];
	channel_mask += (1<<i) * channel.enabled;
	trigger_mask += (1<<i) * channel.enable_trigger_source;
	trigger_out_mask += (1<<i) * channel.enable_trigger_out;
	//write the per-channel stuff
	//Zero suppression threshold
	uint32_t zs_thresh = (1<<31) * channel.zs_polarity +
	  channel.zs_threshold;
	WriteVMERegister(board.address+VME_ChZSThresh+i*0x100,zs_thresh, _handle_board[iboard]);
	//zero suppression time over threshold
	uint32_t nsamp = channel.zs_thresh_time_us * board.GetSampleRate();
	if(nsamp >= (1<<20)) nsamp = (1<<20) -1;
	if(board.zs_type == ZLE){
	  //nsamp contains the pre and post samples
	  if(channel.zs_pre_samps>=(1<<16)) channel.zs_pre_samps = (1<<16)-1;
	  if(channel.zs_post_samps>=(1<<16)) channel.zs_post_samps = (1<<16)-1;
	  uint32_t npre = 
	    std::ceil(channel.zs_pre_samps/board.stupid_size_factor);
	  uint32_t npost = 
	    std::ceil(channel.zs_post_samps/board.stupid_size_factor);
	  nsamp = (npre<<16) + npost;
	}
	//warning: what shall we do with the ZS n samples? 
	if(board.board_type == V1725)//we just use the x4 gain, or 0.5V input range
	  WriteVMERegister(board.address+ VME_ChZSNsamples+i*0x100, 1, _handle_board[iboard]);
	else 
	  WriteVMERegister(board.address+ VME_ChZSNsamples+i*0x100, nsamp, _handle_board[iboard]);

	//trigger threshold
	WriteVMERegister(board.address+ VME_ChTrigThresh+i*0x100, 
			 channel.threshold, _handle_board[iboard]);
	//time over trigger threhsold
	nsamp = std::ceil(channel.thresh_time_us * board.GetSampleRate()) 
	  / board.stupid_size_factor;
	
	if(nsamp >= (1<<12)) nsamp = (1<<12) - 1;
	if(board.board_type == V1725) {//warning: 1n84 is trigger logic now, maybe write to 1n70?
	  //we may need to write to 1n84 to make this meaningful
	  WriteVMERegister(board.address+ 0x1070+i*0x100, nsamp, _handle_board[iboard]);
	}
	else
	  WriteVMERegister(board.address+ VME_ChTrigSamples+i*0x100, nsamp, _handle_board[iboard]);
	//dc offset
	WriteVMERegister(board.address+ VME_ChDAC+i*0x100, channel.dc_offset, _handle_board[iboard]);
	//wait until the dac has updated
	uint32_t status = 0x4;
	while( channel.enabled &&( (status&0x4) || !(status&0x2)) ){
	  //std::cerr<<"Waiting for channel "<<i<<" of board "<<iboard<<" to update DAC...";
	  status = ReadVMERegister(board.address+VME_ChStatus +i*0x100, _handle_board[iboard]);
	  //std::cerr<<" status is now "<<status<<std::endl;
	}
      }
      //finish up with the board parameters
      uint32_t channel_config = (1<<16) * board.zs_type + 
	(1<<6) * board.trigger_polarity + 
	(1<<4) + //Memory Sequential access
	(1<<3) * board.enable_test_pattern + 
	(1<<1) * board.enable_trigger_overlap;
      WriteVMERegister(board.address+ VME_ChannelsConfig, channel_config, _handle_board[iboard]);
      //Buffer code (determines total trigger time
      WriteVMERegister(board.address+VME_BufferCode, board.GetBufferCode(), _handle_board[iboard]);
      //Custom size of register
      
      WriteVMERegister(board.address+ VME_CustomSize, 
		       board.GetCustomSizeSetting(), _handle_board[iboard]);
      
      //Acquisition Control
      uint32_t acq_control =  
	(1<<3) * board.count_all_triggers +
	_params.send_start_pulse;
//      acq_control = (1<<3) * board.count_all_triggers + 4;
      WriteVMERegister(board.address+VME_AcquisitionControl,acq_control, _handle_board[iboard]);
      board.acq_control_val = acq_control;
      //trigger mask
      if(board.local_trigger_coincidence >7) 
	board.local_trigger_coincidence = 7;
      if(board.coincidence_window_nclocks>15)
	board.coincidence_window_nclocks = 15;
      trigger_mask += (1<<31) * board.enable_software_trigger 
	+ (1<<30) * board.enable_external_trigger
	+ (1<<24) * board.local_trigger_coincidence
	+ (1<<20) * board.coincidence_window_nclocks;
      WriteVMERegister(board.address+ VME_TrigSourceMask, trigger_mask, _handle_board[iboard]);
      //trigger out mask
      trigger_out_mask += (1<<31) * board.enable_software_trigger_out +
	(1<<30) * board.enable_external_trigger_out;
      WriteVMERegister(board.address+ VME_TrigOutMask, trigger_out_mask, _handle_board[iboard]);
      //post trigger setting
      
      WriteVMERegister(board.address+VME_PostTriggerSetting, 
		       board.GetPostTriggerSetting(), _handle_board[iboard]);
      //signal logic and front panel programming
      WriteVMERegister(board.address+ VME_FrontPanelIO,
		       board.signal_logic + (1<<6), _handle_board[iboard]);
      //channel mask
      WriteVMERegister(board.address+ VME_ChannelMask, channel_mask, 
		       _handle_board[iboard]);
      //VME control
      uint32_t vme_control = (1<<5) * _params.align64 + 
	(1<<4) + //enable bus error
	(1<<3) + //enable optical link error
	1; //interrupt level
      WriteVMERegister(board.address+ VME_VMEControl, vme_control, _handle_board[iboard]);
      //Interrupt num, BLT event num
      WriteVMERegister(board.address+VME_InterruptOnEvent, 0, _handle_board[iboard]);
      WriteVMERegister(board.address+VME_BLTEvents, 1, _handle_board[iboard]);
      //wait until the board is ready to take data
      uint32_t status = 0;
      int count = 0;
      while( !((status&0x100) && (status&0xc0)) ){
	status = ReadVMERegister(board.address+VME_AcquisitionStatus, _handle_board[iboard]);
	if(count++ > 500){
	  Message(ERROR)<<"Unable to initialize board "<<iboard<<" at address "
			<<std::hex<<board.address<<std::dec<<"\n";
	  return 1;
	}
      }
      
    }
    /*Find the max expected event size in bytes
    The header size is 16 bytes per board
    The data size is 2 bytes per sample per channel
    */
    Message(DEBUG)<<"The expected event size is "<<_params.GetEventSize()
		  <<" bytes."<<std::endl;   
    //wait 2 seconds for DC offset levels to adjust
    time_t now = time(0);
    while(time(0) - now < 2) {}
  }
  catch(...){ 
    return -3;
  }
  return 0;
}

//predicate for find_if function used to test if any board has data
bool DataAvailable(uint32_t status)
{
  return (status & 0x8);
}

void V172X_Daq::DataAcquisitionLoop()
{
  if(!_initialized) Initialize();
  //prepare some variables
  _triggers = 0;
  std::vector<uint32_t> acq_status(_params.enabled_boards);
  //Arm the boards
  try{
    std::vector<uint32_t> acq_write;
    for(int i=0; i<_params.nboards; i++){
      if(_params.board[i].enabled){
	acq_write.push_back(_params.board[i].acq_control_val+0x4);
      }
    }
    WriteVMERegisters(VME_SWClear,1);
    WriteVMERegisters(VME_AcquisitionControl,&(acq_write[0]));
    if(_params.send_start_pulse)
      CAENVME_SetOutputRegister(_handle_vme_bridge, 
				cvOut0Bit | cvOut1Bit | cvOut2Bit | 
				cvOut3Bit | cvOut4Bit );
  }
  catch(std::exception& e){
    Message(ERROR)<<"Unable to arm the board for run!\n";
    _initialized=false;
    _status=INIT_FAILURE;
    return;
  }
  int32_t irq_handle = _handle_vme_bridge;
  if (_params.vme_bridge_link) 
    for(int i=0; i < _params.nboards; i++) { 
      irq_handle = _handle_board[i]; 
      if (_params.board[i].enabled && _params.board[i].link >= 0) break;
    }
  
  //figure out whether we can use interrupts
  //usb interface can't do interrupts, so use polling only
  
  bool use_interrupt = true;
  for(int i=0; i<_params.nboards; i++){
    if(_params.board[i].enabled && _params.board[i].usb) 
      use_interrupt = false;
  }

  
  while(_is_running ){

    CVErrorCodes err = cvTimeoutError;
    
    if(use_interrupt){
      //Enable IRQ lines and wait for an interrupt
      CAENVME_IRQEnable(irq_handle,0xFF);
      err = CAENVME_IRQWait(irq_handle,0xFF,
					 _params.trigger_timeout_ms);
      CAENVME_IRQDisable(irq_handle,0xFF);
    }
    else{
      for(int i=0; i<_params.nboards; i++){
	if(!_params.board[i].enabled) continue;
	//see if there is an event ready on the board
	if(DataAvailable(ReadVMERegister(_params.board[i].address+
					 VME_AcquisitionStatus, 
					 _handle_board[i])) ){
	  err = cvSuccess;
	  break;
	}
      }
    }

    //Check to see if there was an interrupt or just the timeout
    switch(err){
    case cvSuccess:
      break;
    case cvGenericError:
      Message(DEBUG)<<"Generic error occurred. Probably just a timeout...\n";
      //notice: no break here!
    case cvTimeoutError:
      if(!use_interrupt){
	//sleep for the timeout interval
	boost::this_thread::sleep(
	    boost::posix_time::millisec(_params.trigger_timeout_ms));
      }
      if(_params.auto_trigger){
	Message(DEBUG)<<"Triggering...\n";
	WriteVMERegisters(VME_SWTrigger,1);
      }
      else
	Message(DEBUG)<<"Waiting for trigger..."<<std::endl;
      continue;
      break;
    default:
      Message(ERROR)<<"Unknown error waiting for trigger interrupt\n";
      break;
    }
    
    //if we get here, there is an event ready for download
    //get a new event ready 
    RawEventPtr next_event(new RawEvent);
    size_t blocknum = 
      next_event->AddDataBlock(RawEvent::CAEN_V172X,
			       _params.event_size_bytes+event_size_padding);
    unsigned char* buffer = next_event->GetRawDataBlock(blocknum);
    const uint32_t UNSET_EVENT_COUNTER = 0xFFFFFFFF;
    uint32_t event_counter = UNSET_EVENT_COUNTER;
    //get the data
    int data_transferred = 0;
    for(int i=0; i<_params.nboards; i++){
      if(!_params.board[i].enabled) continue;
      //wait until the event is ready on the board
      
      int tries = 0;
      const int maxtries = 50;
      while(!DataAvailable(ReadVMERegister(_params.board[i].address+
					   VME_AcquisitionStatus, _handle_board[i])) &&
	    tries++ < maxtries) {}
      if(tries >= maxtries){
	Message(DEBUG)<<"No trigger received on board "<<i<<"\n";
	continue;
      }
      
      int this_dl_size = 0;
      tries=0;
      CVErrorCodes err = cvSuccess;
      while(this_dl_size == 0 && ++tries<maxtries ){
	err = 
	  CAENVME_FIFOBLTReadCycle(_handle_board[i], _params.board[i].address,
				   buffer + data_transferred,
				   _params.board[i].event_size_bytes,
				   cvA32_U_MBLT, cvD64, &this_dl_size);
	if(err != cvSuccess && err != cvBusError){
	  Message(ERROR)<<"Error generated while downloading event from board"
			<<i<<": "<<CAENVME_DecodeError(err)<<"\n";
	  continue;
	}
      }
      
      if(this_dl_size<=0){
	Message(ERROR)<<"0 bytes downloaded for board "<<i<<std::endl;
	Message(DEBUG)<<"Events stored on this board: "
                      <<ReadVMERegister(_params.board[i].address+
                                        VME_EventsStored, _handle_board[i])<<"\n";
        uint32_t out;
        out = ReadVMERegister(_params.board[i].address+VME_VMEStatus, 
			      _handle_board[i]);
        Message(DEBUG)<<"VME status:"
                      <<"\n\tBERR flag: "<< (out&4)
                      <<"\n\tOutput buffer full: "<< (out&2)
                      <<"\n\tData ready: "<< (out&1)<<"\n";
        out = ReadVMERegister(_params.board[i].address+VME_AcquisitionStatus, 
			      _handle_board[i]);
        Message(DEBUG)<<"Acquisition status:"
                      <<"\n\tReady for acquisition: "<< (out&256)
                      <<"\n\tPLL Status: "<< (out&128)
                      <<"\n\tPLL Bypass: "<< (out&64)
                      <<"\n\tClock source: "<< (out&32)
                      <<"\n\tEvents full: "<< (out&16)
                      <<"\n\tEvent ready: "<< (out&8)
                      <<"\n\t Run on: "<< (out&4) <<"\n";
        out = ReadVMERegister(_params.board[i].address+VME_EventSize, 
			      _handle_board[i]);
	Message(DEBUG)<<"Next event size: "<<out<<"\n";
        Message(DEBUG)<<"Expected event size "
                      <<_params.board[i].event_size_bytes/4<<"\n";
	
	
	//free the buffer so we don't re-trigger spuriously
	WriteVMERegister(_params.board[i].address+VME_BufferFree,1, 
			 _handle_board[i]);
	Message(DEBUG)<<"Events stored on this board: "
		      <<ReadVMERegister(_params.board[i].address+
					VME_EventsStored, _handle_board[i])
		      <<"\n";
	WriteVMERegister(_params.board[i].address+VME_BufferFree,2, 
			 _handle_board[i]);
	Message(DEBUG)<<"Events stored on this board: "
		      <<ReadVMERegister(_params.board[i].address+
					VME_EventsStored, _handle_board[i])
		      <<"\n";
	Message(ERROR)<<"Boards don't usually recover from this error! "
		      <<"Aborting!\n";
	_status = COMM_ERROR;
	_is_running = false;
	break;
      }
      else{
	
	//this could ignore potential align64 extra bits:
	//data_transferred += this_dl_size;
	//instead, we check the actual event
	long ev_size = (*((uint32_t*)(buffer+data_transferred)) & 
			     0x0FFFFFFF) * sizeof(uint32_t);
	if(std::abs(ev_size - this_dl_size) > 5){
	  Message(WARNING)<<"Event size does not match download count!\n\t"
			  <<"Event size: "<<ev_size<<"; download size: "
			  <<this_dl_size<<"; requested download "
			  <<_params.board[i].event_size_bytes<<std::endl;
	  Message(ERROR)<<"Boards don't usually recover from this error! "
		      <<"Aborting!\n";
	  _status = COMM_ERROR;
	  _is_running = false;
	  break;
	}
	//check for event ID misalignment
	uint32_t evct = (((uint32_t*)(buffer+data_transferred))[2])&0xFFFFFF;
	if(event_counter == UNSET_EVENT_COUNTER)
	  event_counter = evct;
	else if(evct != event_counter){
	  Message(CRITICAL)<<"Mismatched event ID on board "<<i
			   <<"; received "<<evct<<", expected "<<event_counter
			   <<"; Aboring run\n";
	  _status = GENERIC_ERROR;
	  _is_running=false;
	  break;
	}
      
	data_transferred += ev_size;
	
      }
    }
    if(GetStatus() != NORMAL){
      _is_running = false;
      break;
    }
    if(data_transferred <= 0){
      Message(ERROR)<<"No boards transferred usable data this event!\n";
      continue;
    }
    else{
      _triggers++;
      next_event->SetDataBlockSize(blocknum, data_transferred);
      PostEvent(next_event);
    }    
  }//end while(_is_running)
   //We only reach here once the event has stopped, so clean up any mess
  //First, end the run, and clear the buffers
  if(_params.send_start_pulse)
    CAENVME_ClearOutputRegister(_handle_vme_bridge, 0xFFFF);
  for(int i=0; i<_params.nboards; i++){
    if(_params.board[i].enabled){
      WriteVMERegister(_params.board[i].address+ VME_AcquisitionControl,
		       _params.board[i].acq_control_val, _handle_board[i]);
      WriteVMERegister(_params.board[i].address+VME_SWClear,0x1, _handle_board[i]);
    }
  }
  WriteVMERegisters(VME_SWReset,1);
  Message(DEBUG)<<_triggers<<" total triggers downloaded."<<std::endl;
}  
