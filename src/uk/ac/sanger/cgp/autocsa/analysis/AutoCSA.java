package uk.ac.sanger.cgp.autocsa.analysis ;

import java.util.ArrayList;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import uk.ac.sanger.cgp.autocsa.beans.*;
import uk.ac.sanger.cgp.autocsa.exceptions.*;
import uk.ac.sanger.cgp.autocsa.util.*;

/**
 *<p> Top level class for AutoCSA package.</p>
 *
 *Original author:  emd
 *@author $Author$
 *@version $Revision$
 */
public class AutoCSA {
	
	protected static Log log = LogFactory.getLog(AutoCSA.class.getName());

/**
 * version string of AutoCSA package hardwired into class.
 */
  public static final String VERSION="1.81";
/**
 * sts name.
 */
  public String analysedSts;
/**
 * sts amplimer sequence (reverse complemented if necessary).
 */
  public String refSeq;
  private int refSearchStart;
  private int refSearchEnd;
  private SeqTraceAnalysis controlAnalysis=null;
  private SeqTraceAnalysis sampleAnalysis=null;
  private CSAParameters csaParams=null;
  private int[][] missingBases=null;
  private int[] minPeakAmps;
  private int[] minShoulderAmps;
  private int[] minShoulderAmpsRelaxed;
  private int ROIStartCoord=0;
  private int ROIEndCoord=0;
// constants very specific to AutoCSA
// constants indicating a description of trace holes
  public static final int MISSING_IN_NORMAL=1;
  public static final int MISSING_IN_SAMPLE=2;
  public static final int MISSING_IN_BOTH=3;
  public static final int ASSIGNED_DELETION=4;
  public static final int ASSIGNED_INSERTION=5;
  public static final int SPECULATIVE_INDEL=6;
  public static final int MATCHING_ERROR=7;
// constants indicating peak classification
  public static final int NOISE_PEAK=0;
  public static final int CALLED_PEAK=1;
  public static final int HET_SNP_PEAK=2;
  public static final int HOM_SNP_PEAK=3;
  public static final int HET_INS_PEAK=4;
  public static final int HOM_INS_PEAK=5;
  public static final int DYE_BLOB_PEAK=6;

/**
 * Constructor for AutoCSA.
 * @param sts sts name
 * @param refseq sts amplimer sequence (supply in sense direction)
 * @param reverse boolean indicating whether an anti-sense trace(s)
 * @param refStart base number to start analysis
 * @param refEnd base number to finish analysis
 */
  public AutoCSA(String sts, String refseq, boolean reverse, int refStart, int refEnd) {
    analysedSts=sts;
    refseq=refseq.toUpperCase();
    if ( reverse ) {
      refSeq=reverseComp(refseq);
    } else {
      refSeq=refseq;
    }
    refSearchStart=refStart;
    refSearchEnd=refEnd;
  }
/**
 * Returns the version of AutoCSA.
 */
  public String getCSAVersion() {
    return VERSION;
  }
/**
 * Set the CSA analysis and comparison parameters.
 * @param params input CSAParameters bean to override the defaults
 */
  public void setCSAParameters(CSAParameters params) {
    csaParams=params;
  }
/**
 * Perform a CSA analysis step on the suplied trace.
 * @param trace Trace object specifying input trace
 * @param minAmp NOT relevant (unused)
 * @param minDist minimum allowed distance between peaks (in scans)
 * @param traceID identifier for the trace
 * @param isNormal boolean indicating whether Trace is a Normal
 * @return traceAnalysis SeqTraceAnalysis object
 * @throws BadTraceException If 1 or more Channels contains zero peaks
 */
  public SeqTraceAnalysis analyseTrace(Trace trace, int minAmp, int minDist,
      String traceID, boolean isNormal) throws BadTraceException {

    SeqTraceAnalysis traceAnalysis = null;
    int i,iChan;

    int[] apeaks = null;
    int[] cpeaks = null;
    int[] gpeaks = null;
    int[] tpeaks = null;
    int[] apeaksPos = null;
    int[] cpeaksPos = null;
    int[] gpeaksPos = null;
    int[] tpeaksPos = null;

// derive the min intensity for peak detection for each channel
    setMinIntensityForPeaks(trace,isNormal);
    int AIndex=Constants.INPUT_BASE_ORDERING.indexOf("A");
    int CIndex=Constants.INPUT_BASE_ORDERING.indexOf("C");
    int GIndex=Constants.INPUT_BASE_ORDERING.indexOf("G");
    int TIndex=Constants.INPUT_BASE_ORDERING.indexOf("T");
    int[] numPeaks = new int[] { 0, 0, 0, 0 };

// find peaks then set intensities and locations to individual arrays
    trace.getChannel(GIndex).findPeaks(minPeakAmps[GIndex],minShoulderAmps[GIndex],minShoulderAmpsRelaxed[GIndex],minDist);
    gpeaks=trace.getChannel(GIndex).getPeaks();
    gpeaksPos=trace.getChannel(GIndex).getPeaksPos();

    trace.getChannel(AIndex).findPeaks(minPeakAmps[AIndex],minShoulderAmps[AIndex],minShoulderAmpsRelaxed[AIndex],minDist);
    apeaks=trace.getChannel(AIndex).getPeaks();
    apeaksPos=trace.getChannel(AIndex).getPeaksPos();

    trace.getChannel(TIndex).findPeaks(minPeakAmps[TIndex],minShoulderAmps[TIndex],minShoulderAmpsRelaxed[TIndex],minDist);
    tpeaks=trace.getChannel(TIndex).getPeaks();
    tpeaksPos=trace.getChannel(TIndex).getPeaksPos();

    trace.getChannel(CIndex).findPeaks(minPeakAmps[CIndex],minShoulderAmps[CIndex],minShoulderAmpsRelaxed[CIndex],minDist);
    cpeaks=trace.getChannel(CIndex).getPeaks();
    cpeaksPos=trace.getChannel(CIndex).getPeaksPos();

    if ( gpeaks != null ) {
      numPeaks[GIndex]=gpeaks.length;
    }
    if ( apeaks != null ) {
      numPeaks[AIndex]=apeaks.length;
    }
    if ( tpeaks != null ) {
      numPeaks[TIndex]=tpeaks.length;
    }
    if ( cpeaks != null ) {
      numPeaks[CIndex]=cpeaks.length;
    }
    int npeaks=numPeaks[0]+numPeaks[1]+numPeaks[2]+numPeaks[3];
    if(log.isInfoEnabled()) log.info("Total Peaks: "+npeaks);
    if(log.isInfoEnabled()) log.info("No A,C,G,T Peaks: "+numPeaks[AIndex]+" "+
            numPeaks[CIndex]+" "+numPeaks[GIndex]+" "+numPeaks[TIndex]);

//  assume that if we have zero of any base type then have a bad trace
//  in this case throw a BadTraceException
    if( numPeaks[0] == 0 || numPeaks[1] == 0 || numPeaks[2] == 0 
        || numPeaks[3] == 0 ) {
      throw new BadTraceException("Bad Trace File: Zero peaks in a Channel");
    }

// construct the SeqTraceAnalysis object (workhorse trace analysis obj)
    traceAnalysis = new SeqTraceAnalysis(npeaks,refSearchStart,refSearchEnd,minPeakAmps,minDist);
// associate the input Trace obj with the output SeqTraceAnalysis obj
    traceAnalysis.setTraceObj(trace);

//  set CSA analysis parameters
    if( csaParams != null ) {
      traceAnalysis.setAnalysisParams(csaParams.getPeakSearchBin(),
                                      csaParams.getRefSearchStart(),
                                      csaParams.getRefSearchEnd(),
                                      csaParams.getRefSearchStartInc(),
                                      csaParams.getMaxBasesMissed(),
                                      csaParams.getMinPeakSpacing());
    }
// merge individual base scan & intensity arrays
    int[] scan = new int[npeaks+1];
    scan[0]=0;
    System.arraycopy(gpeaksPos,0,scan,1,numPeaks[GIndex]);
    System.arraycopy(apeaksPos,0,scan,1+numPeaks[GIndex],numPeaks[AIndex]);
    System.arraycopy(tpeaksPos,0,scan,1+numPeaks[GIndex]+numPeaks[AIndex],numPeaks[TIndex]);
    System.arraycopy(cpeaksPos,0,scan,1+numPeaks[GIndex]+numPeaks[AIndex]+numPeaks[TIndex],numPeaks[CIndex]);

    int[] intensity = new int[npeaks+1];
    intensity[0]=0;
    System.arraycopy(gpeaks,0,intensity,1,numPeaks[GIndex]);
    System.arraycopy(apeaks,0,intensity,1+numPeaks[GIndex],numPeaks[AIndex]);
    System.arraycopy(tpeaks,0,intensity,1+numPeaks[GIndex]+numPeaks[AIndex],numPeaks[TIndex]);
    System.arraycopy(cpeaks,0,intensity,1+numPeaks[GIndex]+numPeaks[AIndex]+numPeaks[TIndex],numPeaks[CIndex]);

    // set ref seq
    traceAnalysis.setRefSeq(refSeq);
    // set array (individual bases as elements)
    traceAnalysis.setPeakToBaseData(numPeaks[GIndex],numPeaks[AIndex],numPeaks[TIndex],numPeaks[CIndex]);
    traceAnalysis.setAnalysisData(scan, intensity);

    // sort the arrays based on Scan index
    int[] sortIndices=traceAnalysis.sortDataset(traceAnalysis.getScan());

    // now apply sort indices to rest of data
    // n.b only curIntensity & curBase are defined at this stage
    traceAnalysis.setScan(traceAnalysis.sortAnalysisData(sortIndices,traceAnalysis.getScan()));
    traceAnalysis.setIntensity(traceAnalysis.sortAnalysisData(sortIndices,traceAnalysis.getIntensity()));
    traceAnalysis.setBase(traceAnalysis.sortBaseData(sortIndices,traceAnalysis.getBase()));

    // apply method to identify ref seq in the trace
    traceAnalysis.matchPeaksToRefSeq();

    // n.b all data must be sorted wrt bpPos before SeqTraceComparison
    // can be instantiated
    sortIndices=traceAnalysis.sortDataset(traceAnalysis.getbpPos());
    traceAnalysis.setbpPos(traceAnalysis.sortAnalysisData(sortIndices,traceAnalysis.getbpPos()));
    traceAnalysis.setScan(traceAnalysis.sortAnalysisData(sortIndices,traceAnalysis.getScan()));
    traceAnalysis.setIntensity(traceAnalysis.sortAnalysisData(sortIndices,traceAnalysis.getIntensity()));
    traceAnalysis.setTracePeakNo(traceAnalysis.sortAnalysisData(sortIndices,traceAnalysis.getTracePeakNo()));
    traceAnalysis.setRealPeak(traceAnalysis.sortAnalysisData(sortIndices,traceAnalysis.getRealPeak()));
    traceAnalysis.setFilter(traceAnalysis.sortAnalysisData(sortIndices,traceAnalysis.getFilter()));
    traceAnalysis.setBase(traceAnalysis.sortBaseData(sortIndices,traceAnalysis.getBase()));

    int[] bpPos=traceAnalysis.getbpPos();
// n.b. reset npeaks as analysis adds extra peak slots to arrays
    npeaks=bpPos.length-1;      // as arrays are len npeaks+1
    // find extent of ref seq analysed
    for( i=0; i < npeaks+1; ++i ) {
      if( bpPos[i] > -1 ) { break; }
    }
    if( i < npeaks ) {

    // set individual trace holes (i.e. points where analysis failed)

      traceAnalysis.setTraceHoles(findMissingBases(bpPos,bpPos[i],bpPos[npeaks],traceAnalysis.getAnalysisIndex(bpPos[i])));
    }

// return the SeqTraceAnalysis object
    return traceAnalysis;
  }
/**
 * Outputs plain text files containing columns of trace analysis parameters.
 * @param idNormal Identifier for Normal trace
 * @param idSample Identifier for Sample trace
 */
  public void writeAnalysisTables(String idNormal, String idSample) {

    if( controlAnalysis != null ) {
      String outTableN=analysedSts+"_csa.dat"+idNormal;
      controlAnalysis.createAnalysisTable(outTableN);
    }
    if( sampleAnalysis != null ) {
      String outTableS=analysedSts+"_csa.dat"+idSample;
      sampleAnalysis.createAnalysisTable(outTableS);
    }

  }
/**
 * Top level caller for CSA trace analysis.
 * @param trace Trace object specifying input trace
 * @param traceID identifier for the trace
 * @param minAmp NOT relevant (unused)
 * @param minDist minimum allowed distance between peaks (in scans)
 * @param isNormal boolean indicating whether Trace is a Normal
 * @return traceAnalysis SeqTraceAnalysis object
 * @throws BadTraceException If 1 or more Channels contains zero peaks
 */
  public SeqTraceAnalysis doCSAAnalysis (Trace trace, String traceID,
                                        int minAmp, int minDist,
                                        boolean isNormal)
                                        throws BadTraceException {

    SeqTraceAnalysis traceAnalysis=null;
    try {
// analyse the trace i.e. match peaks to amplimer seq
      traceAnalysis=analyseTrace(trace,minAmp,minDist,traceID,isNormal);
    } catch ( BadTraceException e ) {
      throw new BadTraceException(e.getMessage(),e);
    }
    return traceAnalysis;
  }

/**
 * Sets the CSA analysis objects for the autoCSA class.
 * @param cTraceAnal the Normal analysis object
 * @param sTraceAnal the Sample analysis object
 */
  public void setAnalyses(SeqTraceAnalysis cTraceAnal,
                          SeqTraceAnalysis sTraceAnal) {
    controlAnalysis=cTraceAnal;
    sampleAnalysis=sTraceAnal;
  }
/**
 * Setup method defines a comparison object from individual analyses.
 * @param cTraceAnal the Normal analysis object
 * @param sTraceAnal the Sample analysis object
 * @return comparison the SeqTraceComparison object
 */
  public SeqTraceComparison setupAutoCSA(SeqTraceAnalysis cTraceAnal,
                                         SeqTraceAnalysis sTraceAnal) {

    SeqTraceAnalysis[] traceAnalysis = new SeqTraceAnalysis[2];
    SeqTraceComparison comparison = null;

    int[] bpLimit1 = new int[2];
    int[] bpLimit2 = new int[2];
    int[] bpIndex = new int[2];
    int i,j;

    traceAnalysis[0]=cTraceAnal;
    traceAnalysis[1]=sTraceAnal;
    controlAnalysis=cTraceAnal;
    sampleAnalysis=sTraceAnal;

// loop over both SeqTraceAnalysis objects for convenience
    for( j=0; j < 2; ++j ) {

      if( traceAnalysis[j] == null ) {
        if(log.isInfoEnabled()) log.info("Trace Analysis Failed");
        return comparison;
      }
      int[] bpPos=traceAnalysis[j].getbpPos();
      int npeaks=bpPos.length-1;      // as arrays are len npeaks+1
    // find extent of ref seq analysed
      for( i=0; i < npeaks+1; ++i ) {
        if( bpPos[i] > -1 ) { break; }
      }
      // n.b. sort on bpPos moves all analysed data to end
      if( i > npeaks ) {
        if(log.isInfoEnabled()) log.info("Failed to Locate Parent Seq");
        return comparison;
      }
      bpIndex[j]=i;
      bpLimit1[j]=bpPos[i];
      bpLimit2[j]=bpPos[npeaks];

    }

    // set individual trace holes (i.e. points where analysis failed)
    // n.b this is not wrt the comparison extent!
//   comment as done in analyseTrace

//  traceAnalysis[0].setTraceHoles(findMissingBases(traceAnalysis[0].getbpPos(),bpLimit1[0],bpLimit2[0],traceAnalysis[0].getAnalysisIndex(bpLimit1[0])));
//  traceAnalysis[1].setTraceHoles(findMissingBases(traceAnalysis[1].getbpPos(),bpLimit1[1],bpLimit2[1],traceAnalysis[1].getAnalysisIndex(bpLimit1[1])));

    // Trace Comparison section

    //int bpLim1=Math.max(bpLimit1[0],bpLimit1[1]);
    //int bpLim2=Math.min(bpLimit2[0],bpLimit2[1]);
      int[] lims = new int[2];
      lims=findComparisonLimits(traceAnalysis[0].getbpPos(),
                                traceAnalysis[1].getbpPos());
      int bpLim1=lims[0];
      int bpLim2=lims[1];
      int npoints=bpLim2-bpLim1+1;
    //int ctrlStart=bpIndex[0]+bpLim1-bpLimit1[0];
    //int sampStart=bpIndex[1]+bpLim1-bpLimit1[1];
      int ctrlStart=traceAnalysis[0].getAnalysisIndex(bpLim1);
      int sampStart=traceAnalysis[1].getAnalysisIndex(bpLim1);
// n.b. normally ctrlStart & sampStart will not be -1 but if there
// is no overlap between analysed seq then these will be -1
    if ( ctrlStart == -1 || sampStart == -1 ) {
      if(log.isInfoEnabled()) log.info("RefSeq Match failure (No Coverage Overlap between traces)");
      return comparison;
    }

    // set trace holes common to portion of seq used for comparison
    ArrayList ctrlList = findMissingBases(traceAnalysis[0].getbpPos(),bpLim1,bpLim2,ctrlStart);
    ArrayList sampList = findMissingBases(traceAnalysis[1].getbpPos(),bpLim1,bpLim2,sampStart);

    int count1=0;
    int count2=0;
    int count3=0;
    ArrayList missing=null;
    if ( ctrlList.size() + sampList.size() == 0 ) {
      missing=new ArrayList(1);
    if(log.isInfoEnabled()) log.info("NO Bases Missing");
    } else {
      missingBases = findMissingInfo(ctrlList,sampList);
      missing=new ArrayList();
    // clear ctrlList & sampList so can reset with bp pos not in both
      ctrlList.clear();
      sampList.clear();
      for( i=0; i < missingBases[0].length; ++i ) {
        if ( missingBases[1][i] == MISSING_IN_NORMAL ) {
          ctrlList.add(count1++,new Integer(missingBases[0][i]));
        }
        if ( missingBases[1][i] == MISSING_IN_SAMPLE ) {
          sampList.add(count2++,new Integer(missingBases[0][i]));
        }
        if ( missingBases[1][i] == MISSING_IN_BOTH ) {
          ++count3;
        }
        missing.add(i,new Integer(missingBases[0][i]));
    if(log.isInfoEnabled()) log.info("Missing: "+missingBases[0][i]+" "+missingBases[1][i]);
      }
    }
    // take account of missing bases in the overall analysed range
    npoints-=(count1+count2+count3);
    int nmiss=count1+count2+count3;

// n.b. must have a minimum coverage i.e. comparison region
    if ( npoints < Constants.MIN_ALLOWABLE_BASES ) {
      if(log.isInfoEnabled()) log.info("RefSeq Match failure (< Constants.MIN_ALLOWABLE_BASES Bases to Analyse)");
      return comparison;
    }
    if(log.isInfoEnabled()) log.info("Seq Range: "+bpLim1+" "+bpLim2+" Npoints: "+npoints+
                       " WT: "+bpLimit1[0]+" "+bpLimit2[0]+
                       " MT: "+bpLimit1[1]+" "+bpLimit2[1]+
                       " Missed: "+nmiss);
//  if(log.isInfoEnabled()) log.info("Index: "+ctrlStart+" "+sampStart);

//  construct SeqTraceComparison object (workhorse trace comparison obj)
    comparison = new SeqTraceComparison(npoints);

//  set CSA comparison parameters
    if( csaParams != null ) {
      comparison.setComparisonParams(csaParams.getCritMutRatio(),
                                     csaParams.getPeakSearchBin());
    }
    // set missing base info array - assumed as HOM mutations
    comparison.setMissingBases(missingBases);

    int[] tempScan = null;
    int[] tempHt = null;

    tempScan=extractData(traceAnalysis[0].getScan(),ctrlStart,npoints,bpLim1,1,sampList.size());
    tempHt=extractData(traceAnalysis[0].getIntensity(),ctrlStart,npoints,bpLim1,1,sampList.size());

    // set scan and intensity arrays for control data
    comparison.setControlData(tempScan,tempHt);

    // n.b. need to reset these as arraycopy only copies into empty arrays
    tempScan = null;
    tempHt = null;

    tempScan=extractData(traceAnalysis[1].getScan(),sampStart,npoints,bpLim1,2,ctrlList.size());
    tempHt=extractData(traceAnalysis[1].getIntensity(),sampStart,npoints,bpLim1,2,ctrlList.size());

    // set scan and intensity arrays for sample data
    comparison.setSampleData(tempScan,tempHt);

    // set arrays for other required data
    comparison.initialiseData(bpLim1,missing);

    // set portion of ref seq analysed from original refSeq string
    comparison.setRefSeq(refSeq);
    comparison.setAnalysedRefSeqData(bpLim1,missing);

// return the SeqTraceComparison object
    return comparison;
  }
/**
 * Runs the autoCSA comparison step, normalises and calculates peak drops.
 * @param comparison the SeqTraceComparison object
 */
  public void runAutoCSA(SeqTraceComparison comparison ) {
    // set peak ratios (factors)
    comparison.setPeakFactors();

// if ROI is set then set it in the SeqTraceComparison object
    if (ROIStartCoord > 0 ) {
      comparison.setROICoords(ROIStartCoord,ROIEndCoord);
    }
    // apply normalisation to peak factors (i.e. derives scale)
    // n.b. first test tumour trace for signature of a het Indel and
    // if found set location parameter in SeqTraceComparison obj
    int indelPos=sampleAnalysis.checkQualityProfileForInDel(true);
    comparison.setIndelBasePosition(indelPos);
    comparison.normalisePeaks(controlAnalysis,sampleAnalysis);

  }
/**
 * Sets a Region of Interest (ROI) for the amplimer in which to optimize the data analysis.
 * @param start the start base number for the ROI
 * @param end the finish base number for the ROI
 */
  public void setROICoords(int start, int end) {
      ROIStartCoord=start;
      ROIEndCoord=end;
  }
/**
 * Performs a search for all mutation types on the SeqTraceComparison object.
 * @param comparison the SeqTraceComparison to analyse
 */
  public void mutationScan(SeqTraceComparison comparison) {

    // search array scale for Heterozygous point mutations
    // n.b. hetMutationScan needs control & sample trace analysis object
    // to search over and also to reset certain fields
    comparison.hetMutationScan(controlAnalysis,sampleAnalysis);

    comparison.finalHetMutationScan(controlAnalysis,sampleAnalysis);

//  Leave out finding Hets in Normals 
//  comparison.controlHetScan(controlAnalysis,sampleAnalysis);

    // search for Heterozygous INS/DEL type mutations
    comparison.testHetIndel();

    // search for Homozygous point mutations & InDels
    // n.b hom snp's and del's are dependent on missingBases
    if ( missingBases != null) {
      comparison.homDeletionScan(controlAnalysis,sampleAnalysis,
                                 missingBases);
      comparison.homMutationScan(controlAnalysis,sampleAnalysis,
                                 missingBases);
    }
    // n.b hom ins's are independent on missingBases
    comparison.homInsertionScan(controlAnalysis,sampleAnalysis,
                                missingBases);
  }
  private int[] findComparisonLimits(int[] bpPos1,int[] bpPos2) {
    int[] basePos = new int[refSeq.length()+1];
    int[] limits = new int[2];
    int i;
    for( i=0; i < basePos.length; ++i ) {
      basePos[i]=0;
    }
// increment basePos[i] for each trace where we have analysed the base
    for( i=0; i < bpPos1.length; ++i ) {
      if( bpPos1[i] == -1 ) { continue; }
      ++basePos[bpPos1[i]];
    }
    for( i=0; i < bpPos2.length; ++i ) {
      if( bpPos2[i] == -1 ) { continue; }
      ++basePos[bpPos2[i]];
    }
// basePos[i] == 2 defines which bases we have analysed in both traces
    for( i=0; i < basePos.length; ++i ) {
      if( basePos[i] == 2 ) { break; }
    }
    limits[0]=i;
    for( i=basePos.length-1; i > 0; --i ) {
      if( basePos[i] == 2 ) { break; }
    }
    limits[1]=i;
// return the limits of the comparison range
    return limits;
  }
  private ArrayList findMissingBases(int[] bpPos,int bpLim1, int bpLim2, int startIndex) {
// returns an ArrayList of any missing base No.s
// n.b. will always return in ascending order
    int nbase=bpLim1;
    int count=0;
    int nmiss;

    ArrayList missing = new ArrayList();
    for( int i=startIndex; i < bpPos.length; ++i ) {
      if( bpPos[i] > bpLim2 ) { break; }
      if( bpPos[i] != nbase++ ) { 
// check how many bases are missing
        nmiss=bpPos[i]-bpPos[i-1]-1;
        for (int j=0; j< nmiss ; ++j) {
          missing.add(count++,new Integer(nbase-1));
          ++nbase;
        }
      }
    }
// return ArrayList of missing base No.s
    return missing;
  }
  private int[] extractData(int[] data, int startIndex, int npoints, int startBp, int thisTrace, int nmissOther) {
  // returns an extracted array to take into account any missing base No.s
  // i.e. due to missing base No.s in the other trace

    int[] newData = new int[npoints];
    if ( nmissOther == 0 ) {
      System.arraycopy(data,startIndex,newData,0,npoints);
      return newData;
    }
    int count=0;
    int index=0;
    int tempLen=npoints+nmissOther;
    int[] tempData = new int[tempLen];
    System.arraycopy(data,startIndex,tempData,0,tempLen);
    int save,otherTrace;
    int iarray=0;
    if( thisTrace == 1 ) {
      otherTrace=2;
    } else {
      otherTrace=1;
    }
    int nmissTot=0;
    if( missingBases != null ) {
      nmissTot=missingBases[0].length;
    }
    for( int bp=startBp; bp < startBp+npoints+nmissTot ; ++bp ) {
      save=1;
      if ( index < nmissTot && missingBases[0][index] == bp ) {
        save=0;        // don't save if data missing
// don't save this trace's data point if not in the other trace
        if ( missingBases[1][index] == otherTrace ) {
          ++iarray;
        }
        ++index;
      }

      if( save == 0 ) { continue; }
      newData[count++]=tempData[iarray++];
    }
    return newData;
  }
  private int[][] findMissingInfo(ArrayList ctrlList,ArrayList sampList) {
  // returns an array of amalgamated information about missing base No.s
  
    int nctrl=ctrlList.size();
    int nsamp=sampList.size();
    int i,j;
  
    if ( nctrl == 0 ) {
      int[][] missingInfo = new int[2][nsamp];
      for( i=0; i < nsamp; ++i ) {
        missingInfo[0][i]=((Integer) sampList.get(i)).intValue();
        missingInfo[1][i]=MISSING_IN_SAMPLE;    // 2 indicates Sample only
      }
      return missingInfo;
    } else if ( nsamp == 0 ) {
      int[][] missingInfo = new int[2][nctrl];
      for( i=0; i < nctrl; ++i ) {
        missingInfo[0][i]=((Integer) ctrlList.get(i)).intValue();
        missingInfo[1][i]=MISSING_IN_NORMAL;    // 1 indicates Control only
      }
      return missingInfo;
    } else {
  // remove duplicates
      ArrayList ctype = new ArrayList(nctrl);
      ArrayList stype = new ArrayList(nsamp);
      int ictrl,isamp,styp;
      for( j=0; j < nsamp; ++j ) {
        stype.add(j,new Integer(MISSING_IN_SAMPLE));
      }
      for( i=0; i < nctrl; ++i ) {
        ictrl=((Integer) ctrlList.get(i)).intValue();
        ctype.add(i,new Integer(MISSING_IN_NORMAL));
        j=0;
        while (sampList.size() > 0 && j < sampList.size()) {
          isamp=((Integer) sampList.get(j)).intValue();
          if ( isamp == ictrl ) {
            sampList.remove(j);
            stype.remove(j);
            --nsamp;
  // reset code type of duplicates
            ctype.set(i,new Integer(MISSING_IN_BOTH));
            continue;
          }
          ++j;
        }
      }
  // now insert samp info into control
      nsamp=sampList.size();
      for( j=0; j < nsamp; ++j ) {
        isamp=((Integer) sampList.get(j)).intValue();
        styp=((Integer) stype.get(j)).intValue();
        for( i=0; i < nctrl; ++i ) {
          ictrl=((Integer) ctrlList.get(i)).intValue();
          if ( isamp < ictrl ) { break; }
        }
        i=Math.min(i,nctrl);
        ++nctrl;
        ctype.add(i,new Integer(styp));
        ctrlList.add(i,new Integer(isamp));
      }
  // dimension with ctrl as sample copied into ctrl
      int[][] missingInfo = new int[2][nctrl];
      for( i=0; i < nctrl; ++i ) {
        missingInfo[0][i]=((Integer) ctrlList.get(i)).intValue();
        missingInfo[1][i]=((Integer) ctype.get(i)).intValue();
      }
      return missingInfo;
    }
  }
  private String reverseComp(String seq) {
  //  reverse compliment sequence
    String base;
    StringBuffer buf= new StringBuffer();
        for( int i = 0 ; i< seq.length() ; ++i) {
          base = seq.substring(i, i+1);
          if ( base.equals("C") ) {
            buf.append("G");
          }
          if ( base.equals("A") ) {
            buf.append("T");
          }
          if ( base.equals("T") ) {
            buf.append("A");
          }
          if ( base.equals("G") ) {
            buf.append("C");
          }
          if ( base.equals("N") ) {
            buf.append("N");
          }
        }
    buf.reverse();
    return buf.toString();
  }
  private void setMinIntensityForPeaks(Trace trace,boolean isNormal) {
    minPeakAmps = new int[] { 0, 0, 0, 0 };
    minShoulderAmps = new int[] { 0, 0, 0, 0 };
    minShoulderAmpsRelaxed = new int[] { 0, 0, 0, 0 };
//  set multiplication factors for definition of min peak intensities
//  n.b. G channel can contain smaller peaks so factor is set smaller
    float[] factors = new float[] { Constants.A_CHANNEL_PEAK_RATIO,
                                    Constants.C_CHANNEL_PEAK_RATIO,
                                    Constants.G_CHANNEL_PEAK_RATIO,
                                    Constants.T_CHANNEL_PEAK_RATIO };
    float[] sfactors = new float[] { Constants.A_CHANNEL_SHOULDER_RATIO,
                                     Constants.C_CHANNEL_SHOULDER_RATIO,
                                     Constants.G_CHANNEL_SHOULDER_RATIO,
                                     Constants.T_CHANNEL_SHOULDER_RATIO };
// get base ordering in input file
    int[] indices =getBaseIndices();
    int len=trace.getChannel(0).getChannelLength();
    int i;
    int strictMinPeak;
    float decrease=Constants.UL_DECREASE_FACTOR;
    float shFac=Constants.SHOULDER_RELAXATION;
    int cut=Constants.SCAN_SUBTRACTOR_FOR_PEAKS;
// cut off 1000 scans from each end of the trace to avoid dye blobs
// and larger spikes at end of trace - to get a reasonable max intensity
    for( i=0; i < 4; ++i ) {
      minPeakAmps[i]=trace.getChannel(i).getMaxIntensity(cut,len-cut);
    }
// multiply base factors and max I's to get min I's for peak detection
    for( i=0; i < 4; ++i ) {
//if(log.isInfoEnabled()) log.info("Channel Max: "+minPeakAmps[i]);
      minShoulderAmps[i]=(int) (sfactors[indices[i]]*(float) minPeakAmps[i]);
      minShoulderAmpsRelaxed[i]=(int) (shFac * sfactors[indices[i]] *(float) minPeakAmps[i]);
      strictMinPeak=Constants.STRICT_MIN_PEAK;
// n.b. Do NOT make allowances for Underloaded data in Normals
      if ( ! isNormal ) {
        if( minPeakAmps[i] < Constants.CUTOFF_UNDERLOAD_DECREASE ) strictMinPeak= (int) (decrease* (float) Constants.STRICT_MIN_PEAK);
      }
      minPeakAmps[i]=(int) (factors[indices[i]]*(float) minPeakAmps[i]);
      minPeakAmps[i]=Math.max(strictMinPeak,minPeakAmps[i]);
    }
  }
  private int[] getBaseIndices() {
    int[] indices = new int[4];
    int AIndex=Constants.INPUT_BASE_ORDERING.indexOf("A");
    int CIndex=Constants.INPUT_BASE_ORDERING.indexOf("C");
    int GIndex=Constants.INPUT_BASE_ORDERING.indexOf("G");
    int TIndex=Constants.INPUT_BASE_ORDERING.indexOf("T");
    indices[AIndex]=0;
    indices[CIndex]=1;
    indices[GIndex]=2;
    indices[TIndex]=3;
// return array of indices defining base ordering in input trace file
    return indices;
  }
}
