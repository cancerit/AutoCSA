package uk.ac.sanger.cgp.autocsa.util;

import java.io.IOException;
import java.util.ArrayList;
import java.beans.*;
import java.lang.reflect.*;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import uk.ac.sanger.cgp.autocsa.analysis.*;
import uk.ac.sanger.cgp.autocsa.beans.*;
import uk.ac.sanger.cgp.autocsa.exceptions.*;
import org.biojava.bio.chromatogram.UnsupportedChromatogramFormatException;

/**
 *<p> Class provides the interface between autoCSA and hotCSA.</p>
 *<p> Public static methods analyseTrace() & performComparison().</p>
 *
 *Original author:  emd
 *@author $Author$
 *@version $Revision$
 */
public class CSALoader {
	
	protected static Log log = LogFactory.getLog(CSALoader.class.getName());

  private static final boolean debug=false;
// definition of wait time for trace file locking problem
  private static final long WAIT_TIME=5000L;

/**
 * null constructor.
 */
  public CSALoader() {
// constructor not needed as all methods are static
  }
/**
  * Performs a CSA analysis step on the specified trace
  * @param csatrace CSATrace bean
  * @return traceDet TraceDetails bean containing trace analysis parameters, Note traceDet is an updated version of TraceDetails beans obtained from input CSATrace bean.
  * @throws BadTraceException if any trace data Channel has zero peaks
  * @throws CSAException if any (miscelleaneous) IOExceptions occur
  **/

  public static TraceDetails analyseTrace(CSATrace csatrace ) throws BadTraceException, CSAException, BadCommentException {

// get input parameters from csatrace
// n.b. traceDet is updated with trace analysis info then output
  TraceDetails traceDet=csatrace.getTraceDetails();

  String refseq=csatrace.getRefSeq();
  String sts = traceDet.getStsName();
  char forRev=traceDet.getForwardOrReverseStrand();
  boolean reverse=false;
  if( forRev == 'r' ) {reverse=true; }

// n.b. minAmp,minDist,refStart,refEnd hardired until we're
// ready to set them directly from the database
  int minAmp=csatrace.getMinPeakIntensity();
  int minDist=csatrace.getMinPeakSpacing();
  minAmp=500;  // dummy value - not used
  int refStart=csatrace.getRefSearchStart();
  int refEnd=refseq.length()+csatrace.getRefSearchEnd();

// set Mobility-correction type
  Adjustment.setAdjustmentType(csatrace.getMobilityCorrection());

// read in trace file name in traceDet - data is set in Trace object
  Trace myTrace;
  try {
    myTrace=getTraceObj(traceDet);
  } catch ( UnsupportedChromatogramFormatException e ) {
    if(log.isWarnEnabled()) log.warn("Caught UnsupportedChromatogramFormatException: "+e.getMessage());
    throw new CSAException("CSAOutput.analyseTrace() ",e);
  } catch ( IOException e ) {
    if(log.isWarnEnabled()) log.warn("Caught IOException: "+e.getMessage());
    throw new CSAException("CSAOutput.analyseTrace() ",e);
  }
  if( csatrace.getMobilityCorrection() == 0 ) {
    baseline(myTrace);
  }

// ready for AutoCSA analysis - construct AutoCSA object
  AutoCSA autoCSA = new AutoCSA(sts,refseq,reverse,refStart,refEnd);
// n.b. we're now ready to reset default CSA parameters
  CSAParameters csaParams = getCSAParameters(csatrace);
//reset refEnd in csaParams as absolute value needed rather than offset
  csaParams.setRefSearchEnd(refEnd);
  autoCSA.setCSAParameters(csaParams);  // pass input params bean

// perform CSA trace analysis step and set all analysis results
  SeqTraceAnalysis traceAnal;
  try {
    traceAnal=autoCSA.doCSAAnalysis(myTrace,"G",minAmp,minDist,traceDet.isNormal());
    traceDet.setCsaVersion(autoCSA.getCSAVersion());
    setTraceAnalInfo(traceAnal,traceDet);
    setTraceParams(traceAnal,traceDet);
  } catch ( BadTraceException e ) {
    if(log.isWarnEnabled()) log.warn("Caught BadTraceException: "+e.getMessage());
    throw new BadTraceException("CSAOutput.load() ",e);
  }

// results are returned via updated TraceDetails bean
  return traceDet;
  }
/**
  * Performs a CSA Comparison step on the specified pair of traces
  * @param csainput CSAInput bean
  * @return csaout CSAOutput bean containing trace comparison parameters
  * @throws BadTraceException if any trace data Channel has zero peaks
  * @throws CSAException if any (miscelleaneous) IOExceptions occur
  **/
  public static CSAOutput performComparison(CSAInput csainput ) throws BadTraceException, CSAException, BadCommentException {

// get input parameters from csainput
  TraceDetails traceWt=csainput.getWtTraceDetails();
  TraceDetails traceMut=csainput.getMutTraceDetails();

  String refseq=csainput.getRefSeq();
  String sts = traceMut.getStsName();
  char forRev=traceWt.getForwardOrReverseStrand();
  boolean reverse=false;
  if( forRev == 'r' ) {reverse=true; }

// n.b. minAmp,minDist,refStart,refEnd hardired until we're
// ready to set them directly from the database
  int minAmp=csainput.getMinPeakIntensity();
  int minDist=csainput.getMinPeakSpacing();
  minAmp=500; // dummy value - not used
  int refStart=csainput.getRefSearchStart();
  int refEnd=refseq.length()+csainput.getRefSearchEnd();
// set Mobility-correction type
  Adjustment.setAdjustmentType(csainput.getMobilityCorrection());

// read in trace file names in traceDet objs - data is set in Trace object
  Trace myTrace1,myTrace2;
  try {
    myTrace1=getTraceObj(traceWt);
    myTrace2=getTraceObj(traceMut);
  } catch ( UnsupportedChromatogramFormatException e ) {
    if(log.isWarnEnabled()) log.warn("Caught UnsupportedChromatogramFormatException: "+e.getMessage());
    throw new CSAException("CSAOutput.performComparison() ",e);
  } catch ( IOException e ) {
    if(log.isWarnEnabled()) log.warn("Caught IOException: "+e.getMessage());
    throw new CSAException("CSAOutput.performComparison() ",e);
  }
  if( csainput.getMobilityCorrection() == 0 ) {
    baseline(myTrace1);
    baseline(myTrace2);
  }

// ready for AutoCSA comparison - construct AutoCSA object
  AutoCSA autoCSA = new AutoCSA(sts,refseq,reverse,refStart,refEnd);
  CSAParameters csaParams = getCSAParameters(csainput);
//reset refEnd in csaParams as absolute value needed rather than offset
  csaParams.setRefSearchEnd(refEnd);
  autoCSA.setCSAParameters(csaParams);  // pass input params bean

// perform CSA trace analyses steps, then comparison step (if possible)
// and set all comparison results (in CSAOutput obj)
  SeqTraceComparison comparison=null;
  SeqTraceAnalysis cTraceAnal,sTraceAnal;
  CSAOutput csaout = null;
  try {
    cTraceAnal=autoCSA.doCSAAnalysis(myTrace1,"N",minAmp,minDist,traceWt.isNormal());
    sTraceAnal=autoCSA.doCSAAnalysis(myTrace2,"S",minAmp,minDist,traceMut.isNormal());
//  n.b. QC is not applied on the analyses here so we are passing
//  the trace based on the result from it's original analysis.
    comparison= autoCSA.setupAutoCSA(cTraceAnal,sTraceAnal);
    if( comparison == null ) {
      csaout=new CSAOutput();
    } else {
// n.b. firstly set q array and search for dye blobs
// these must be set as are referenced by comparison methods
      setTraceAnalQualityAndDyeBlobs(cTraceAnal,traceWt);
      setTraceAnalQualityAndDyeBlobs(sTraceAnal,traceMut);
      int[] roi=getROICoordsForAnalysis(traceWt);
// set roi coords for ALL objects WT, TT & comparison
      cTraceAnal.setROICoords(roi[0],roi[1]);
      sTraceAnal.setROICoords(roi[0],roi[1]);
      autoCSA.setROICoords(roi[0],roi[1]);
// set comparison parameters (drop,bin) dependent on dna type
      if( traceMut.getDnaType() == Constants.DNA_MRX ) {
        comparison.setComparisonParams(0.5f,6.0f);
      }
      autoCSA.runAutoCSA(comparison);
// set association between CSAOutput and comparison
      comparison.setCSAOutput(new CSAOutput());
      autoCSA.mutationScan(comparison);
// hand back updated CSAOutput object
      csaout = comparison.getCSAOutput();
    }
// set all trace info in output object csaout
    csaout.setMutTraceDetails(traceMut);
    csaout.setWtTraceDetails(traceWt);
    //csaout.setCSAMutTraceAnalInfo(sTraceAnal);
    //csaout.setCSAWtTraceAnalInfo(cTraceAnal);
    csaout.setComparisonInfo(comparison,autoCSA.getCSAVersion());
//  n.b. these calls now done by analyseTrace
//  if ( sTraceAnal != null ) {
//    csaout.setMutTraceHoles(sTraceAnal.getTraceHoles());
//  }
//  if ( cTraceAnal != null ) {
//    csaout.setWtTraceHoles(cTraceAnal.getTraceHoles());
//  }
//  setTraceParams(sTraceAnal,traceMut);
//  setTraceParams(cTraceAnal,traceWt);
//  extract required fields from CSAInput bean and copy into CSAOutput
    csaout.setNormalNormalComp(csainput.getNormalNormalComp());
  } catch ( BadTraceException e ) {
    if(log.isWarnEnabled()) log.warn("Caught BadTraceException: "+e.getMessage());
    throw new BadTraceException("CSAOutput.performComparison() ",e);
  }

// results are returned via CSAOutput object
  return csaout;
  }
/**
  * Instantiates a Trace object from a specified trace file
  * @param td input TraceDetails object (bean)
  * @return trace Trace object formed from specified trace in input
  * @throws UnsupportedChromatogramFormatException if trace file is wrong type
  * @throws IOException for any other miscelleaneous exceptions
  **/
  private static Trace getTraceObj(TraceDetails td) throws IOException, UnsupportedChromatogramFormatException, BadTraceException, BadCommentException {

  int nChans=4;
  Trace trace = null;

  if(td.useChromatogram()) {
  	trace = new Trace(td.getChromatogram(), nChans);
  }
  else {
	  trace = new Trace(td.getFilePositionFromRoot(),nChans);
  }
//Trace trace=new Trace(file.toString(),nChans);
// use generic trace loader in Canutil package
  try {
    trace.loadSCFFile(1);    // uses generic loader
//  trace.loadTraceFile();    // uses ab1 loader
  } catch ( UnsupportedChromatogramFormatException e ) {
    if(log.isWarnEnabled()) log.warn("Caught UnsupportedChromatogramFormatException: "+e.getMessage());
    throw new IOException("CSALoader.getTraceObj() "+e.getMessage());
  } catch ( IOException e ) {
// take account of case where trace file might be temporarily locked
// by another user - sleep and then try again
    try {
      if(debug){if(log.isInfoEnabled()) log.info("Trace access problem waiting and retrying");}
      Thread.sleep(WAIT_TIME); // wait time is in milliseconds
    } catch (InterruptedException ie) {
      throw new IOException("CSALoader.getTraceObj() (InterruptedException)"+ie.getMessage());
    }
    try {
      trace.loadSCFFile(1);    // uses generic loader
    } catch ( IOException e1 ) {
      if(log.isWarnEnabled()) log.warn("Caught IOException: "+e1.getMessage());
      throw new IOException("CSALoader.getTraceObj() "+e1.getMessage());
    }
  }
  
// return Trace object complete with data Channels loaded
  return trace;
  }
/**
  * Performs a baselining of a trace
  * @param trace input Trace object specifying trace file
  **/
  private static void baseline(Trace trace) {

// baseline channels if necessary
  int i,traceLimit1,traceLimit2;
  if ( trace.AbiType.equals("3730") ) {
    traceLimit1=0;
    traceLimit2=0;
    for( i=0; i < trace.nChannels; ++i ) {
      trace.baseline(i,traceLimit1,traceLimit2);
    }
  }

  }
/**
  * Derives trace analysis information from a SeqTraceAnalysis object
  * and creates a set of TraceBaseParam beans for each data point, these
  * are then added to an ArrayList and set in the TraceDetails object
  * @param traceAnal input SeqTraceAnalysis object
  * @param traceDet input TraceDetails object
  **/
  private static void setTraceParams(SeqTraceAnalysis traceAnal,
                                        TraceDetails traceDet) {
    if ( traceAnal == null ) { return; }
    ArrayList traceParams = new ArrayList();
    float[] quality=traceAnal.getTraceQuality();
    int[] bpos=traceAnal.getbpPos();
    int[] scan=traceAnal.getScan();
    int[] hgts=traceAnal.getIntensity();
// set up ArrayList of position, q value, scan index in bean
    for( int i=0; i < bpos.length ; ++i ) {
      if ( bpos[i] == -1 ) {continue; }
      traceParams.add(new TraceBaseParam(bpos[i],quality[i],scan[i],hgts[i]));
    }
// update traceDet with traceParams ArrayList
    traceDet.setTraceBaseParams(traceParams);
  }
/**
  * Instantiates a Trace object from a specified trace file
  * @param td input TraceDetails object (bean)
  * @return trace Trace object formed from specified trace in input
  **/
  private static int[] getROICoordsForAnalysis(TraceDetails traceDet) {
   int[] temp=traceDet.getROICoords();
    int[] roi=new int[2];
    roi[0]=temp[0];
    roi[1]=temp[1];
    if( traceDet.getForwardOrReverseStrand() == 'r' ) {
      int stsLen=traceDet.getRefSeq().length();
      if ( roi[0] > 0 && roi[1] > 0 ) {
        roi[0]=stsLen-temp[1]+1;
        roi[1]=stsLen-temp[0]+1;
      }
    }
// return ROI coordinates appropriate for required Seq orientation
    return roi;
  }
  public static void setTraceAnalQualityAndDyeBlobs(SeqTraceAnalysis traceAnal, TraceDetails traceDet) {
    if ( traceAnal != null ) {
      traceAnal.setTumourSample( ! traceDet.isNormal());
      traceAnal.setTraceQuality(Constants.QUALITY_WINDOW_SIZE);
      int ndb=0;
      if( Constants.DYE_BLOB_WINDOW_START1 > 0 ) {
        ndb+=traceAnal.searchForDyeBlobs(Constants.DYE_BLOB_WINDOW_START1,
                Constants.DYE_BLOB_WINDOW_WIDTH,1.75f);
      }
      if( Constants.DYE_BLOB_WINDOW_START2 > 0 ) {
        ndb+=traceAnal.searchForDyeBlobs(Constants.DYE_BLOB_WINDOW_START2,
                Constants.DYE_BLOB_WINDOW_WIDTH,1.75f);
      }
    }
  }
  private static void setTraceAnalInfo(SeqTraceAnalysis traceAnal,
                                        TraceDetails traceDet) {
    int[] limits={0, 0};
    float quality=0.0f, qProfile=0.0f;
    int statusCode;
    if ( traceAnal != null ) {
      traceAnal.setTumourSample( ! traceDet.isNormal());
      traceAnal.setTraceQuality(Constants.QUALITY_WINDOW_SIZE);
      int[] roi=getROICoordsForAnalysis(traceDet);
      traceAnal.setROICoords(roi[0],roi[1]);
// reset any q values appropriately for dye blobs, under & over loading
      setTraceAnalQualityAndDyeBlobs(traceAnal,traceDet);
//    traceAnal.setModifiedQwhereUnderloaded();
//    traceAnal.setModifiedQwhereOverloaded();
// extract the results from the analysis
      quality=traceAnal.getAverageTraceQuality();
      qProfile=traceAnal.getQualityProfile(3.0f);
      limits=traceAnal.getAnalysisLimits();
      statusCode=traceAnal.getTraceStatus();
//  if trace is a Normal and has status of 1 (passed analysis) then
//  test to check if normal has signature of an Indel affecting the ROR
      if( traceDet.isNormal() && statusCode == Constants.ORA_ID_RERUN_ANALYSIS ) {
        int indelPos=traceAnal.checkQualityProfileForInDel(false);
        if( indelPos > 0 && ( roi[0] == 0 ||
            (roi[0] > 0 && indelPos < roi[1]) ) ) {
          statusCode=Constants.ORA_ID_SUSPECT_INDEL_IN_ROR;
        }
      }
      traceDet.setTraceHoles(traceAnal.getTraceHoles());
      traceDet.setMaxScanIndex(traceAnal.getTraceObj().getChannel(0).getChannelLength());
    } else {
// the analysis has failed only if traceAnal == null
      statusCode=Constants.ORA_ID_RERUN_MATCH_FAILED;     // code for failure to locate parent seq
    }
// update traceDet bean with all analysis parameters
    traceDet.setOverallTraceQuality(new Float(quality));
    traceDet.setTraceQualityPercentage(new Float(qProfile));
    traceDet.setStsCoverageStart(new Integer(limits[0]));
    traceDet.setStsCoverageStop(new Integer(limits[1]));
    traceDet.setTraceStatus(statusCode);
    if ( debug ) {
      if(log.isInfoEnabled()) log.info("AutoCSA Found: Status: "+statusCode+" Q: "+quality);
      if(log.isInfoEnabled()) log.info("AutoCSA Found: Analysis Limits: "+limits[0]+" "+limits[1]);
    }
  }
  private static CSAParameters getCSAParameters(CSAExtractableParameters
params) throws CSAException {
// extract the parameters we know are definately in the input bean
// and set these in the output bean
    CSAParameters csaParams = new CSAParameters();
    csaParams.setPeakSearchBin(params.getPeakSearchBin());
    csaParams.setRefSearchStartInc(params.getRefSearchStartInc());
    csaParams.setMaxBasesMissed(params.getMaxBasesMissed());
    csaParams.setMobilityCorrection(params.getMobilityCorrection());
    csaParams.setRefSearchEnd(params.getRefSearchEnd());
    csaParams.setRefSearchStart(params.getRefSearchStart());
    csaParams.setMinPeakSpacing(params.getMinPeakSpacing());
    csaParams.setMinPeakIntensity(params.getMinPeakIntensity());

// now test input bean for presence of a getCritMutRatio() method
// n.b this bean parameter will only be needed for CSA comparisons
    BeanInfo info;
    try {
      info = Introspector.getBeanInfo(params.getClass());
    } catch (IntrospectionException e) {
      if(log.isWarnEnabled()) log.warn("Caught IntrospectionException: "+e.getMessage());
      throw new CSAException("CSALoader.getCSAParameters() ",e);
    }
    PropertyDescriptor[] pd = info.getPropertyDescriptors();

// test CSAExtractableParameters for presence of getCritMutRatio method
// if present then set the CritMutRatio
    String value;
    for(int i = 0 ; i< pd.length ; i++) {
      Method meth = pd[i].getReadMethod();
      if( meth.getName().equals("getCritMutRatio") ) {
        try {
          value = meth.invoke(params,null).toString();
        } catch (IllegalAccessException e) {
          if(log.isWarnEnabled()) log.warn("Caught IllegalAccessException: "+e.getMessage());
          throw new CSAException("CSALoader.getCSAParameters() ",e);
        } catch (InvocationTargetException e) {
          if(log.isWarnEnabled()) log.warn("Caught InvocationTargetException: "+e.getMessage());
          throw new CSAException("CSALoader.getCSAParameters() ",e);
        }
        csaParams.setCritMutRatio(Float.valueOf(value).floatValue());
        break;
      }
    }

// return the set CSAParameters bean
    return csaParams;
  }

}
