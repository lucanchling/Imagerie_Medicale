//////////////////////////////////////////////////////////////////////////////
//
// A simple gradient echo sequence for the ODIN framework
//
//////////////////////////////////////////////////////////////////////////////



// includes for the sequence objects
#include <odinseq/seqall.h>



// The method itself is a class that is derived from
// an abstract base class 'SeqMethod' that contains all
// routines common to all methods:

class METHOD_CLASS : public SeqMethod {

 public:

  // Constructor that takes the methods identifier (unique string) as its argument
  METHOD_CLASS(const STD_string& label);

  // virtual functions that are overwritten in this method to build the sequence,
  // calculate parameter relations, etc.:
  void method_pars_init();
  void method_seq_init();
  void method_rels();
  void method_pars_set();

 private:

  // Parameters for this method:
  LDRfloat T1Ernst;
  LDRint   DummyCycles;
  LDRfloat PartialFourierRead;

  // Sequence objects for this method:
  SeqPulsar exc;

  SeqDelay relaxdelay;

  SeqObjLoop peloop;
  SeqObjLoop peloop3d;
  SeqObjLoop sliceloop;
  SeqObjLoop reploop;
  SeqObjLoop dummyloop;

  SeqObjList scanpart;
  SeqObjList slicepart;
  SeqObjList dummypart;

  SeqDelay dummypad;

  SeqGradEcho grech;

  SeqVecIter phaseiter;

  SeqGradTrapezParallel crusher;
  SeqDelay   crusherdelay;

};

////////////////////////////////////////////////////////////////////////////////////

METHOD_CLASS::METHOD_CLASS (const STD_string& label)
                                   : SeqMethod(label) {
  // Put in here all stuff that will once be initialised and
  // never changed

  // Specify a short description of the method
  set_description("This is a simple spoiled gradient echo sequence where "
                "each excitation is used to scan one line "
                "in k-space via a gradient echo with the "
                "apropriate phase encoding. The flip angle is "
                "the Ernst angle for the given T1 (to obtain maximum SNR). " );
}


////////////////////////////////////////////////////////////////////////////////////

void METHOD_CLASS::method_pars_init() {
  // In this function the methods parameters will be initialised


  // Assign default values:
  commonPars->set_RepetitionTime(1000.0);
  commonPars->set_AcqSweepWidth(25.0);
  commonPars->set_MatrixSize(readDirection,128);


  // set default values and a short a description (which
  // will appear as a tooltip in the UI) for each parameter
  T1Ernst=1300.0;
  T1Ernst.set_minmaxval(0.0,5000.0).set_description("For optimum SNR, the flip angle will be set to the Ernst angle using this T1");

  DummyCycles=3;
  DummyCycles.set_description("Number of dummy shots before actual acquisition");

  PartialFourierRead=0.0;
  PartialFourierRead.set_description("The amount of partial Fourier undersampling in readout direction (0=no undersampling, 1=half fourier)");

  // Make method specific parameters visible in the user interface
  append_parameter(T1Ernst,"T1Ernst");
  append_parameter(DummyCycles,"DummyCycles");
  append_parameter(PartialFourierRead,"PartialFourierRead");

}


////////////////////////////////////////////////////////////////////////////////////


void METHOD_CLASS::method_seq_init() {
  // Put in here all stuff to create the layout of the sequence


  ///////////////// Excitation Pulse: /////////////////////

  // Get the slice thickness from the global geometry handler 'geometryInfo':
  float slicethick=geometryInfo->get_sliceThickness();
  if(geometryInfo->get_Mode()==voxel_3d) slicethick=0.9*geometryInfo->get_FOV(sliceDirection);

  float spatres=slicethick/4.0;


  // Create the sinc shaped excitation pulse.
  //exc=SeqPulsarSinc("exc",slicethick,false,2.0,commonPars->get_FlipAngle(),spatres);
	exc=SeqPulsarGauss("exc",slicethick,false,4.0,commonPars->get_FlipAngle(), 256);
  // Set the frequency list of the excitation pulse so that we will excite all slices in the slicepack
  exc.set_freqlist( systemInfo->get_gamma("") * exc.get_strength() / (2.0*PII) *  geometryInfo->get_sliceOffsetVector() );

  // This is useful for simulating/visualization of the sequence
  exc.set_pulse_type(excitation);

  // This loop object is used to loop over the slices
  sliceloop=SeqObjLoop("sliceloop");



  ////////////////// Geometry: /////////////////////////////////


  // calculate the resolution in the read Channel and set the number of phase encoding
  // steps so that we will obtain a isotropic resolution in read and phase Channel:
  float resolution=secureDivision(geometryInfo->get_FOV(readDirection),commonPars->get_MatrixSize(readDirection));
//  commonPars->set_MatrixSize(phaseDirection,int(secureDivision(geometryInfo->get_FOV(phaseDirection),resolution)+0.5),noedit);

  int size3d=1;
  if(geometryInfo->get_Mode()==voxel_3d) {
    size3d=int(secureDivision(geometryInfo->get_FOV(sliceDirection),resolution)+0.5); // isotropic resolution
  }
  commonPars->set_MatrixSize(sliceDirection,size3d,noedit);


  //////////////// Delays: //////////////////////////////

  // relaxation delay after each readout
  relaxdelay=SeqDelay("relaxdelay");

  //////////////// Gradient Echo Module: //////////////////////////////

  // This is the gradient-recalled echo kernel of the sequence.
  // Calculations of read/phase gradient strengts and durations are
  // performed by this module itself so we just have to pass in the size and FOV
  // in each direction. This module has a tight timing scheme so that short TE is possible.
  if(geometryInfo->get_Mode()==voxel_3d) {

    // 3D mode
    grech=SeqGradEcho("grech",
              commonPars->get_MatrixSize(readDirection), geometryInfo->get_FOV(readDirection),
              commonPars->get_MatrixSize(phaseDirection), geometryInfo->get_FOV(phaseDirection),
              commonPars->get_MatrixSize(sliceDirection), geometryInfo->get_FOV(sliceDirection),
              exc, commonPars->get_AcqSweepWidth(),
              commonPars->get_ReductionFactor(), DEFAULT_ACL_BANDS, false,
              commonPars->get_PartialFourier(),PartialFourierRead);

  } else {

    // 2D mode
    grech=SeqGradEcho("grech", exc, commonPars->get_AcqSweepWidth(),
              commonPars->get_MatrixSize(readDirection), geometryInfo->get_FOV(readDirection),
              commonPars->get_MatrixSize(phaseDirection), geometryInfo->get_FOV(phaseDirection),
              linearEncoding, noReorder, 1, commonPars->get_ReductionFactor(), DEFAULT_ACL_BANDS, false,
              commonPars->get_PartialFourier(),PartialFourierRead);
  }


  if(commonPars->get_RFSpoiling()) {

    exc.set_phasespoiling();
    grech.set_phasespoiling();

    phaseiter=SeqVecIter("phaseiter");
    phaseiter.add_vector(exc.get_phaselist_vector());
    phaseiter.add_vector(grech.get_phaselist_vector());
  }


  //////////////// Crusher Gradient: //////////////////////////////

  float crusher_strength=0*0.5*systemInfo->get_max_grad();
  float crusher_integral=0*4.0*fabs(grech.get_gradintegral().sum());

  crusher=SeqGradTrapezParallel("crusher",crusher_integral,crusher_integral,crusher_integral, crusher_strength);
  crusherdelay=SeqDelay("crusherdelay",0.1); // Small delay to avoid gradient-induced stimulation


  //////////////// Loops: //////////////////////////////

  // Construct a loop object to iterate through the phase encoding steps
  peloop=SeqObjLoop("peloop");

  // Construct a loop object to iterate through the 2nd phase encoding steps
  peloop3d=SeqObjLoop("peloop3d");

  // Construct a loop object to perform repetitions of the experiment
  reploop=SeqObjLoop("reploop");

  // Construct a loop object to perform dummy cylces
  dummyloop=SeqObjLoop("dummyloop");


  //////////////// Constructing the Sequence: ///////////////


  // This sequence container will hold the objects for one readout:
  scanpart=SeqObjList("scanpart");
  scanpart = grech + crusherdelay + crusher;
  if(commonPars->get_RFSpoiling()) scanpart+=phaseiter;


  // Contains kernel of one dummy cycle
  dummypart=SeqObjList("dummypart");
  dummypad=SeqDelay("dummypad",scanpart.get_duration()-exc.get_duration()); // Padding delay for dummy cycle
  dummypart = exc + dummypad + relaxdelay;


  // This sequence container will hold the objects for one slicepack:
  slicepart=SeqObjList("slicepart");


  if(geometryInfo->get_Mode()==voxel_3d) {

    slicepart = peloop3d(
                  peloop(
                    sliceloop(
                      scanpart + relaxdelay   // the sequence kernel
                    )[grech.get_exc_vector()]
                  )[grech.get_pe_vector()]
                )[grech.get_pe3d_vector()];
  } else {

    slicepart = peloop(
                  sliceloop(
                    scanpart + relaxdelay   // the sequence kernel
                  )[grech.get_exc_vector()]
                )[grech.get_pe_vector()];

  }


  // Finally, build the whole sequence
  set_sequence(
    dummyloop(
      sliceloop(
        dummypart
      )[grech.get_exc_vector()]
    )[DummyCycles]
  + reploop(
      slicepart
     )[commonPars->get_NumOfRepetitions()]
  );

}


////////////////////////////////////////////////////////////////////////////////////


void METHOD_CLASS::method_rels() {
  // Put in here all stuff that has to be performed whenever one of the sequence parameters
  // has been changed by the user


  // calculate relaxdelay to get the desired repetition time
  float scandur=scanpart.get_duration()*float(geometryInfo->get_nSlices());
  if(scandur>commonPars->get_RepetitionTime()) commonPars->set_RepetitionTime(scandur);
  relaxdelay.set_duration( (commonPars->get_RepetitionTime()-scandur)/float(geometryInfo->get_nSlices()) );


  // calculate Ernst angle accordng to TR
  float flipangle=180.0/PII * acos( exp ( -secureDivision ( commonPars->get_RepetitionTime(), T1Ernst) ) );
  commonPars->set_FlipAngle( flipangle, noedit );
  exc.set_flipangle( flipangle );


  // retrieve echo time from gradient echo module to display it in the user interface
  commonPars->set_EchoTime(grech.get_echo_time(), noedit );
}


////////////////////////////////////////////////////////////////////////////////////

void METHOD_CLASS::method_pars_set() {
  // Put in here all stuff that has to be performed after the parameters have been edited by the user
  // and before the sequence is played out

}


////////////////////////////////////////////////////////////////////////////////////

// entry point for the sequence module
ODINMETHOD_ENTRY_POINT


