#include <odinseq/seqall.h>


class METHOD_CLASS : public SeqMethod {

 private:
  SeqPulsar exc;
  SeqPulsarReph exc_reph;

  SeqPulsar refoc;

  SeqAcqRead acqread;
  SeqAcqRead acqread_templ;
  SeqAcqDeph readdeph;

  SeqGradPhaseEnc pe1,pe2;
  SeqGradPhaseEnc pe3d1,pe3d2;

  SeqGradConstPulse spoiler;

  SeqDelay echopad;
  SeqDelay excpad;
  SeqDelay relaxdelay;

  SeqObjList echopart;
  SeqObjList echopart_templ;

  SeqGradTrapezParallel postexc;

  SeqObjList scan;
  SeqObjList readout;
  SeqObjList readout_templ;

  SeqObjList pe1part;
  SeqObjList pe2part;
  SeqDelay pe1delay;
  SeqDelay pe2delay;

  SeqObjLoop echoloop;
  SeqObjLoop segloop;
  SeqObjLoop sliceloop;
  SeqObjLoop reploop;
  SeqObjLoop loop3d;

  SeqVector echoindex; // index vector for echoes



  LDRenum   PhaseEncoding;
  LDRint    NumOfSegments;
  LDRdouble PulseDur;
  LDRdouble refocFlipAngle;
  LDRdouble SpoilerStrength;
  LDRdouble SpoilerDuration;
  LDRbool   TemplScan;


 public:
  METHOD_CLASS(const STD_string& label) : SeqMethod(label) {
    set_description("Segmented 2D/3D RARE (Turbo SE) Sequence. "
                    "Phase correction is possible by using a prescan with phase-encoding gradients switched off. ");
  }

  unsigned int numof_testcases() const {return 3;}

  void method_pars_init() {

    PhaseEncoding.add_item("Linear",linearEncoding);
    PhaseEncoding.add_item("Reverse",reverseEncoding);
    PhaseEncoding.add_item("CenterOut",centerOutEncoding);
    PhaseEncoding.add_item("CenterIn",centerInEncoding);
    PhaseEncoding=centerOutEncoding;
    PhaseEncoding.set_description("Phase encoding order");

    NumOfSegments=4;
    NumOfSegments.set_description("Number of excitation to acquire one slice/partition");

    PulseDur=5.0;
    PulseDur.set_unit(ODIN_TIME_UNIT).set_description("Duration of excitation and refocusing pulse");

    refocFlipAngle=180.0;
    refocFlipAngle.set_unit(ODIN_ANGLE_UNIT).set_description("Refocusing flip angle");

    SpoilerStrength=50.0;
    SpoilerStrength.set_unit("%").set_description("Spoiler-gradient strength in percent of max gradient amplitude");

    SpoilerDuration=1.0;
    SpoilerDuration.set_unit(ODIN_TIME_UNIT).set_description("Spoiler-gradient duration");

    TemplScan=false;
    TemplScan.set_description("Acquire and use template scan for phase correction");


    // alternative default settings for sequence test
    if(get_current_testcase()==1) {
      geometryInfo->set_Mode(voxel_3d);
      PulseDur=6.0;
      refocFlipAngle=120.0;
      SpoilerStrength=30.0;
      SpoilerDuration=2.0;
    }
    if(get_current_testcase()==2) {
      commonPars->set_ReductionFactor(3);
      TemplScan=false;
      NumOfSegments=1;
    }

    append_parameter(PhaseEncoding,"PhaseEncoding");
    append_parameter(NumOfSegments,"NumOfSegments");
    append_parameter(PulseDur,"PulseDur");
    append_parameter(TemplScan,"TemplScan");


    append_parameter(SpoilerStrength,"SpoilerStrength");
    append_parameter(SpoilerDuration,"SpoilerDuration");

  }


/////////////////////////////////////////////////////////////////////////////////////////////////


void method_seq_init(){
  Log<Seq> odinlog(this,"method_seq_init");

  float gamma=systemInfo->get_gamma("");

  ///////////////// Pulses: /////////////////////

  float slicethick=geometryInfo->get_sliceThickness();

  if(geometryInfo->get_Mode()==voxel_3d) {
    slicethick=geometryInfo->get_FOV(sliceDirection);
  }

  float spatres=slicethick/4.0;

  // Excitation Pulse
  exc=SeqPulsarSinc("exc",slicethick,false,PulseDur,commonPars->get_FlipAngle(),spatres,256);
  if(geometryInfo->get_Mode()==slicepack)  exc.set_freqlist  (gamma * exc.get_strength() / (2.0*PII) *  geometryInfo->get_sliceOffsetVector()  );
  if(geometryInfo->get_Mode()==voxel_3d)   exc.set_freqoffset(gamma * exc.get_strength() / (2.0*PII) *  geometryInfo->get_offset(sliceDirection)  );
  exc.set_pulse_type(excitation);


  // rephasing lobe for excitation pulse
  exc_reph=SeqPulsarReph("exc_reph",exc);


  // Refocusing Pulse
  refoc=SeqPulsarSinc("refoc",slicethick,false,PulseDur,refocFlipAngle,spatres,256);
  refoc.set_phase(90.0);
  if(geometryInfo->get_Mode()==slicepack) refoc.set_freqlist  (gamma * refoc.get_strength() / (2.0*PII) *  geometryInfo->get_sliceOffsetVector()  );
  if(geometryInfo->get_Mode()==voxel_3d)  refoc.set_freqoffset(gamma * refoc.get_strength() / (2.0*PII) *  geometryInfo->get_offset(sliceDirection)  );
  refoc.set_pulse_type(refocusing);


  ////////////////// Geometry: /////////////////////////////////

  // calculate the resolution in the read Channel and set the number of phase encoding
  // steps so that we will obtain a uniform resolution in read and phase Channel:
  float resolution=secureDivision(geometryInfo->get_FOV(readDirection),commonPars->get_MatrixSize(readDirection));

  commonPars->set_MatrixSize(phaseDirection,int(secureDivision(geometryInfo->get_FOV(phaseDirection),resolution)+0.5),noedit);

  commonPars->set_MatrixSize(sliceDirection,int(secureDivision(geometryInfo->get_FOV(sliceDirection),resolution)+0.5),noedit);





  //////////////// Phase Encoding: //////////////////////////

  pe1=SeqGradPhaseEnc("pe1",commonPars->get_MatrixSize(phaseDirection),geometryInfo->get_FOV(phaseDirection),
                  phaseDirection,0.25*systemInfo->get_max_grad(),
                  encodingScheme(int(PhaseEncoding)),interleavedSegmented,NumOfSegments,
                  commonPars->get_ReductionFactor(), DEFAULT_ACL_BANDS, commonPars->get_PartialFourier());

  pe2=pe1;
  pe2.set_label("pe2");
  pe2.invert_strength();

  // Echo indices for template correction
  unsigned int nechoes=((const SeqVector&)pe1).get_numof_iterations();
  ODINLOG(odinlog,significantDebug) << "nechoes=" << nechoes << STD_endl;
  echoindex=SeqVector("echoindex",nechoes);

  //////////////// Phase Encoding (3D): //////////////////////////

  if(geometryInfo->get_Mode()==voxel_3d) {

    pe3d1=SeqGradPhaseEnc("pe3d1",commonPars->get_MatrixSize(sliceDirection),geometryInfo->get_FOV(sliceDirection),pe1.get_constduration(),
                          sliceDirection, linearEncoding, noReorder, 1, commonPars->get_ReductionFactor());

    pe3d2=pe3d1;
    pe3d2.set_label("pe3d2");
    pe3d2.invert_strength();
  }


  //////////////// Readout: //////////////////////////////

  acqread=SeqAcqRead("acqread",commonPars->get_AcqSweepWidth(),commonPars->get_MatrixSize(readDirection),
                               geometryInfo->get_FOV(readDirection),readDirection);

  acqread_templ=SeqAcqRead(acqread);
  acqread_templ.set_label("acqread_templ");
  acqread_templ.set_template_type(phasecorr_template);

  readdeph=SeqAcqDeph("readdeph",acqread,spinEcho);



  //////////////// Spoiler: //////////////////////////////

  float spoiler_strength_phys=(double)SpoilerStrength/100.0*systemInfo->get_max_grad();

  spoiler=SeqGradConstPulse("spoiler",phaseDirection,spoiler_strength_phys,SpoilerDuration);


  //////////////// several padding delays ////////////////////

  excpad=SeqDelay("excpad");
  echopad=SeqDelay("echopad");

  relaxdelay=SeqDelay("relaxdelay");

  pe1delay=SeqDelay("pe1delay");
  pe2delay=SeqDelay("pe2delay");


  //////////////// total sequence: //////////////////////////////

  echopart=SeqObjList("echopart");

  echopart_templ=SeqObjList("echopart_templ");

  scan=SeqObjList("scan");

  echoloop=SeqObjLoop("echoloop");

  segloop=SeqObjLoop("segloop");

  sliceloop=SeqObjLoop("sliceloop");

  reploop=SeqObjLoop("reploop");

  loop3d=SeqObjLoop("loop3d");

  readout=SeqObjList("readout");


  fvector gradint(3); gradint=0.0;
  gradint+=exc_reph.get_gradintegral();
  gradint+=readdeph.get_gradintegral();
  gradint+=spoiler.get_gradintegral();
  postexc=SeqGradTrapezParallel("postexc",gradint[0],gradint[1],gradint[2],0.5*systemInfo->get_max_grad());
  ODINLOG(odinlog,significantDebug) << "postexc. get_gradintegral()=" << postexc.get_gradintegral().printbody() << STD_endl;

  float spoiler_k_diff=systemInfo->get_gamma()*(spoiler.get_gradintegral()[1]-postexc.get_gradintegral()[1]);
  ODINLOG(odinlog,significantDebug) << "spoiler_k_diff=" << spoiler_k_diff << STD_endl;


  if(geometryInfo->get_Mode()==voxel_3d) {
    pe1part+=pe1/pe3d1;
    pe2part+=pe2/pe3d2;
  } else {
    pe1part+=pe1;
    pe2part+=pe2;
 }

  pe1delay=SeqDelay("pe1delay",pe1part.get_duration());
  pe2delay=SeqDelay("pe2delay",pe2part.get_duration());

  echopart_templ=      refoc + spoiler + echopad + pe1delay      + acqread_templ  +  pe2delay      + echopad + spoiler;
  echopart      =      refoc + spoiler + echopad + pe1part       + acqread        +  pe2part       + echopad + spoiler;

  readout_templ = exc + postexc + excpad + echoloop(echopart_templ)[echoindex];
  readout       = exc + postexc + excpad + echoloop(echopart)      [pe1][pe2][echoindex];


  if(geometryInfo->get_Mode()==voxel_3d) {

    if(TemplScan) {
      scan+= readout_templ + relaxdelay;
    }

    scan+=  reploop(
              segloop(
                loop3d(
                  readout + relaxdelay
                )[pe3d1][pe3d2]
              )[pe1.get_reorder_vector()][pe2.get_reorder_vector()]
            )[commonPars->get_NumOfRepetitions()];

  } else {

    if(TemplScan) {
      scan+=  sliceloop(
                readout_templ + relaxdelay
             )[exc][refoc];
    }

    scan+=  reploop(
              segloop(
                sliceloop(
                  readout + relaxdelay
                )[exc][refoc]
              )[pe1.get_reorder_vector()][pe2.get_reorder_vector()]
            )[commonPars->get_NumOfRepetitions()];

  }


  set_sequence( scan );

  }


/////////////////////////////////////////////////////////////////////////////////////////////////


  void method_rels(){
    Log<Seq> odinlog(this,"method_rels");

    // CPMG timings
    double TE1=exc.get_duration()-exc.get_magnetic_center()+
               postexc.get_duration()+
               refoc.get_magnetic_center();

    double TE2=(echopart.get_duration())/2.;

    ODINLOG(odinlog,significantDebug) << "TE1/TE2=" << TE1 << "/" << TE2 << STD_endl;

    if(TE2-TE1 < 0) echopad.set_duration(echopad.get_duration()-(TE2-TE1));
    else excpad.set_duration(excpad.get_duration()+(TE2-TE1));
    ODINLOG(odinlog,significantDebug) << "echopad=" << echopad.get_duration() << STD_endl;
    ODINLOG(odinlog,significantDebug) << "excpad=" << excpad.get_duration() << STD_endl;

    commonPars->set_EchoTime(echopart.get_duration(),noedit);

    // calculate relaxdelay to get the desired repetition time
    float scandur=readout.get_duration()*float(geometryInfo->get_nSlices());
    if(scandur>commonPars->get_RepetitionTime()) commonPars->set_RepetitionTime(scandur);
    relaxdelay=(commonPars->get_RepetitionTime()-scandur)/float(geometryInfo->get_nSlices());


  }


/////////////////////////////////////////////////////////////////////////////////////////////////

  void method_pars_set(){

    // inform the readout about the used phase encoding and slice vector (for automatic reconstruction)

    acqread.set_reco_vector(line,pe1);
    acqread_templ.set_reco_vector(echo,echoindex);
    acqread.set_reco_vector(echo,echoindex);

    if(geometryInfo->get_Mode()==voxel_3d) {
      acqread.set_reco_vector(line3d,pe3d1);
    }

    acqread.set_reco_vector(slice,exc);
    acqread_templ.set_reco_vector(slice,exc);
  }

};



/////////////////////////////////////////////////////////////////////////////////////////////////

// entry point for the sequence module
ODINMETHOD_ENTRY_POINT


