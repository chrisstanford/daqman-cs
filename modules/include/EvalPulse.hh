#ifndef EVAL_PULSE_h
#define EVAL_PULSE_h

//this files provides a function to evaluate pulse objects
//it is used by Rois, Pulses, and TOFs -- Jingke Xu

class Pulse;
class ChannelData;

//function to calculate pulse parameters, if max_peak_time>0
//the peak of the pulse should be within (start time, start time + max_peak_time)
int EvaluatePulse(Pulse& pulse, ChannelData* chdata, double max_peak_time=-1);
				
#endif
