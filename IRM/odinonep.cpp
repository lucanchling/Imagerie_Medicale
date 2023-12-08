#include <odinseq/seqall.h>


class METHOD_CLASS : public SeqMethod {

 public:
  METHOD_CLASS(const STD_string& label);

  void method_pars_init();
  void method_seq_init();
  void method_rels();
  void method_pars_set();

 private:
  LDRenum   OnepulseMode;
  LDRint    NumOfSamples;
  LDRdouble PulseDuration;
  LDRdouble FreqOffset;


  SeqPulsarBP exc;
  SeqPulsarBP refoc;
  SeqAcq acq;
  SeqDelay relaxdelay;
  SeqDelay exc2refoc;
  SeqDelay refoc2acq;
  SeqObjLoop acculoop;
  SeqObjList scanpart;
  SeqGradConstPulse spoiler;
  SeqVector averagevec; // index vector for averages
  SeqVecIter phaseiter;


};


////////////////////////////////////////////////////////////////////////////


METHOD_CLASS::METHOD_CLASS (const STD_string& label)
                                   : SeqMethod(label) {
  set_description("A simple onepulse experiment (non-selective RF pulse + acquisition).");
}

////////////////////////////////////////////////////////////////////////////


void METHOD_CLASS::method_pars_init() {

  OnepulseMode.add_item("FID");
  OnepulseMode.add_item("SpinEcho");
  OnepulseMode.set_actual("FID");

  commonPars->set_EchoTime(0.0);

  NumOfSamples=1024;
  PulseDuration=1.0;
  FreqOffset=0.0;

  append_parameter(OnepulseMode,"OnepulseMode");
  append_parameter(NumOfSamples,"NumOfSamples");
  append_parameter(PulseDuration,"PulseDuration");
  append_parameter(FreqOffset,"FreqOffset");
}

////////////////////////////////////////////////////////////////////////////

void METHOD_CLASS::method_seq_init() {


  exc=SeqPulsarBP("exc",PulseDuration,commonPars->get_FlipAngle(),systemInfo->get_main_nucleus());
  exc.set_pulse_type(excitation);

  refoc=SeqPulsarBP("refoc",PulseDuration,180.0,systemInfo->get_main_nucleus());
  refoc.set_pulse_type(refocusing);

  acq=SeqAcq("acq",(unsigned int)NumOfSamples,commonPars->get_AcqSweepWidth(),1.0,systemInfo->get_main_nucleus());

  exc.set_freqoffset(FreqOffset);
  refoc.set_freqoffset(FreqOffset);
  acq.set_freqoffset(FreqOffset);


  if(commonPars->get_RFSpoiling()) {

    exc.set_phasespoiling(4, 90.0);
    acq.set_phasespoiling(4, 90.0);
    refoc.set_phasespoiling(4, 90.0, 90.0);

    phaseiter=SeqVecIter("phaseiter");
    phaseiter.add_vector(exc.get_phaselist_vector());
    phaseiter.add_vector(acq.get_phaselist_vector());
    phaseiter.add_vector(refoc.get_phaselist_vector());

  } else {

    refoc.set_phase(90.0);
  }



  spoiler=SeqGradConstPulse("spoiler",sliceDirection,0.5*systemInfo->get_max_grad(),2.0);

  relaxdelay=SeqDelay("relaxdelay");
  exc2refoc=SeqDelay("exc2refoc");
  refoc2acq=SeqDelay("refoc2acq");

  acculoop=SeqObjLoop("acculoop");

  if(OnepulseMode=="SpinEcho") scanpart= exc + exc2refoc + spoiler + refoc + spoiler + refoc2acq + acq  + relaxdelay;
  else scanpart= exc  + acq  + spoiler + relaxdelay;

  if(commonPars->get_RFSpoiling()) scanpart += phaseiter;


  averagevec=SeqVector("averagevec",commonPars->get_NumOfRepetitions());

  set_sequence( acculoop(scanpart)[averagevec] );
}

////////////////////////////////////////////////////////////////////////////

void METHOD_CLASS::method_rels() {

   float min_echo_time1=0.0;
   float min_echo_time2=0.0;

   if(OnepulseMode=="SpinEcho") {
     min_echo_time1=(exc.get_duration()-exc.get_magnetic_center())+spoiler.get_duration()+refoc.get_magnetic_center();
     min_echo_time2=(refoc.get_duration()-refoc.get_magnetic_center())+spoiler.get_duration()+acq.get_acquisition_center();
   } else {
     min_echo_time1=(exc.get_duration()-exc.get_magnetic_center())+acq.get_acquisition_center();
   }

   if(commonPars->get_EchoTime()<min_echo_time1) commonPars->set_EchoTime(min_echo_time1);
   if(commonPars->get_EchoTime()<min_echo_time2) commonPars->set_EchoTime(min_echo_time2);


   exc2refoc=commonPars->get_EchoTime()-min_echo_time1;
   refoc2acq=commonPars->get_EchoTime()-min_echo_time2;

  if(scanpart.get_duration()>commonPars->get_RepetitionTime()) commonPars->set_RepetitionTime(scanpart.get_duration());
  relaxdelay.set_duration(commonPars->get_RepetitionTime()-scanpart.get_duration());
}

////////////////////////////////////////////////////////////////////////////

void METHOD_CLASS::method_pars_set() {

  acq.set_reco_vector(average,averagevec);

  recoInfo->set_Recipe("averagecoll | averagesum | kspace | fft | slicecoll | image | store"); // Quick hack to display a single spectrum from multiple averages

}


////////////////////////////////////////////////////////////////////////////

// entry point for the sequence module
ODINMETHOD_ENTRY_POINT
