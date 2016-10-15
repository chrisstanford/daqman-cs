#include "EventHandler.hh"
#include "BaselineFinder.hh"
#include "ConvertData.hh"
#include "PulseFinder.hh"
#include "TOF.hh"
#include "PSD.hh"
#include "SpeFinder.hh"
#include "SumChannels.hh"
#include "EvalRois.hh"
#include "Integrator.hh"
#include "AverageWaveforms.hh"

int EventHandler::AddCommonModules()
{
  AddModule<ConvertData>();
  AddModule<SumChannels>();
  AddModule<BaselineFinder>();
  AddModule<Integrator>();
  AddModule<SpeFinder>();
  AddModule<EvalRois>();
  AddModule<PulseFinder>();
  AddModule<TOF>();
  AddModule<PSD>();
  AddModule<AverageWaveforms>();
  return 9;
}
