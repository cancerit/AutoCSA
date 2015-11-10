package uk.ac.sanger.cgp.autocsa.analysis ;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.Collections;
import java.io.IOException;
import java.io.FileWriter;
import java.io.File;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import uk.ac.sanger.cgp.autocsa.exceptions.*;
import uk.ac.sanger.cgp.autocsa.util.*;

/**
 *<p> Workhorse Class for Sequence Trace Analysis operations.</p>
 *
 *Original author:  Ed Dicks
 *@author $Author$
 *@version $Revision$
 */
public class SeqTraceAnalysis {
	
	protected static Log log = LogFactory.getLog(SeqTraceAnalysis.class.getName());

  public int nPeaks;
  public String stsName;
  public String well;
  public String refSeq;
  private boolean isTumourSample=false;
  private int[] fork;
  private int[] bpPos;
  private char[] curBase;
  private int[] curScan;
  private int[] curIntensity;
  private int[] filter;
  private int[] tracePeakNo;
  private int[] realPeak;
  private float[] quality;
  private int refSearchStart;
  private int refSearchEnd;
  private int refSearchStartInc=40;
  private int maxBasesMissed=50;
  private int[] minPeakIntensity;
  private int minPeakSpacing;
  private float peakSearchBin=6.0f;
  private ArrayList traceHoles=null;
  private Trace traceObj;
  private boolean matchedSomeRefSeq=false;
  private int ROIStartCoord=0;
  private int ROIEndCoord=0;
  private float maxAllowableQuality=25.0f;
  private int reservePeakIndex=0;
  private int totalReservedPeaks=10;
  private String nonSearchableBases="NRYKMBDHV";

/**
 * Allocates a SeqTraceAnalysis object.
 * @param np No of data points in analysis, i.e. No of peaks in trace.
 * @param refStart base No to start Amplimer matching
 * @param refEnd base No to stop Amplimer matching
 * @param minAmp int[4] array of minimum intensities for each channel that were used for peak detection.
 * @param minDist minimum distance between peaks used in peak detection.
 */
  public SeqTraceAnalysis(int np, int refStart, int refEnd, int[] minAmp, int minDist) {
       nPeaks=np+totalReservedPeaks;
       refSearchStart=refStart;
       refSearchEnd=refEnd;
       minPeakIntensity=minAmp;
       minPeakSpacing=minDist;
  }
/**
 * Setter for isTumour boolean.
 * @param isTumour flag whether trace is a Normal or Tumour
 */
  public void setTumourSample(boolean isTumour) {
       isTumourSample=isTumour;
  }
/**
 * Getter for array of Min peak intensities used in peak detection
 * @return array (4) of min peak intensities (order is GATC)
 */
  public int[] getMinPeakIntensity() {
       return minPeakIntensity;
  }
/**
 * Getter for Min peak intensity used in peak detection for a channel
 * @param base channel to select i.e. one of A,C,G,T
 * @return min peak intensity for the specified channel
 */
  public int getMinPeakIntensity(String base) {
       String order=Constants.INPUT_BASE_ORDERING;   // order of channels
       return minPeakIntensity[order.indexOf(base)];
  }
/**
 * Setter for Amplimer sequence
 * @param rSeq amplimer sequence for current trace
 */
  public void setRefSeq(String rSeq) {
       refSeq=rSeq.toUpperCase();
  }
/**
 * <p>Setter for trace analysis data.</p>
 * Initialises all other class analysis arrays.
 * @param scan int array of scan indices of all peaks detected
 * @param intensity int array of intensities of all peaks detected
 */
  public void setAnalysisData(int[] scan, int[] intensity){
       if( totalReservedPeaks > 0 ) {
         curScan=new int[nPeaks+1];
         curIntensity=new int[nPeaks+1];
         curScan[0]=scan[0];
         curIntensity[0]=intensity[0];
         for(int n=0 ; n < totalReservedPeaks ; ++n ) {
           curScan[n+1]=n+1;
           curIntensity[n+1]=1;
         }
         int i=totalReservedPeaks+1;
         for(int n=1 ; n < scan.length ; ++n ) {
           curScan[i]=scan[n];
           curIntensity[i++]=intensity[n];
         }
       } else {
         curScan = scan;
         curIntensity = intensity;
       }
       fork= new int[nPeaks+1];
       bpPos= new int[nPeaks+1];
       filter= new int[nPeaks+1];
       tracePeakNo= new int[nPeaks+1];
       realPeak= new int[nPeaks+1];
       quality= new float[nPeaks+1];
  }
/**
 * Setter for base types of detected peaks
 * @param nG No of G peaks detected
 * @param nA No of A peaks detected
 * @param nT No of T peaks detected
 * @param nC No of C peaks detected
 */
  public void setPeakToBaseData(int nG, int nA, int nT, int nC){
      curBase = new char[nPeaks+1];
      curBase[0]='n';
      int nN=totalReservedPeaks; 
      if( totalReservedPeaks > 0 ) Arrays.fill(curBase,1,1+nN,'N');
      Arrays.fill(curBase,1+nN,1+nN+nG,'G');
      Arrays.fill(curBase,1+nN+nG,1+nN+nG+nA,'A');
      Arrays.fill(curBase,1+nN+nG+nA,1+nN+nG+nA+nT,'T');
      Arrays.fill(curBase,1+nN+nG+nA+nT,nPeaks+1,'C');
  }
/**
 * Sorts array of data and returns sort indices
 * @param data array to derive sort indices
 * @return array of sort indices
 */
  public int[] sortDataset(int[] data){

    List dpList = new ArrayList();

    for(int n=1 ; n < nPeaks+1 ; ++n ) {
      dpList.add(new DataPoint(data[n],n));
    }

    Collections.sort(dpList);
    int[] sortIndices = new int[nPeaks+1];

//  curScan[0]=0;
    sortIndices[0]=0;
    for(int n=1 ; n < nPeaks+1 ; ++n ) {
      sortIndices[n]=((DataPoint) dpList.get(n-1)).getIndex();
    }
    return sortIndices;
  }
/**
 * Sorts an int array of data based on a previous sort
 * @param indices sort indices to define the sorting order
 * @param data array to sort
 * @return array of sorted data
 */
  public int[] sortAnalysisData(int[] indices, int[] data ){

    int[] newData=new int[nPeaks+1];
    for(int n=0 ; n < nPeaks+1 ; ++n ) {
	newData[n] = data[indices[n]];
    }
    return newData;
  }
/**
 * Sorts an bases array (chars) based on a previous sort
 * @param indices sort indices to define the sorting order
 * @param data bases array (char) to sort
 * @return array of sorted bases
 */
  public char[] sortBaseData(int[] indices, char[] data ){

    char[] newData=new char[nPeaks+1];
    for(int n=0 ; n < nPeaks+1 ; ++n ) {
       newData[n] = data[indices[n]];
    }
    return newData;
  }

/**
 * Setter for an individual array element of realPeak.
 * @param index index of the realPeak array
 * @param value new value to use for the array element
 */
  public void setRealPeak(int index, int value ){
    realPeak[index] = value;
  }

/**
 * Setter for the curScan array.
 * @param scan the new value for the curScan array
 */
  public void setScan(int[] scan){
    curScan=scan;
  }
/**
 * Setter for the intensity array.
 * @param intensity the new value for the intensity array
 */
  public void setIntensity(int[] intensity){
    curIntensity=intensity;
  }
/**
 * Setter for the bpPos array.
 * @param pos the new value for the bpPos array
 */
  public void setbpPos(int[] pos){
    bpPos=pos;
  }
/**
 * Setter for the curBase array.
 * @param base the new value for the base array
 */
  public void setBase(char[] base){
    curBase=base;
  }
/**
 * Setter for the realPeak array.
 * @param rpeak the new value for the realPeak array
 */
  public void setRealPeak(int[] rpeak){
    realPeak=rpeak;
  }
/**
 * Setter for the filter array.
 * @param filt the new value for the filter array
 */
  public void setFilter(int[] filt){
    filter=filt;
  }
/**
 * Setter for the tracePeakNo array.
 * @param peakno the new value for the tracePeakNo array
 */
  public void setTracePeakNo(int[] peakno){
    tracePeakNo=peakno;
  }

/**
 * Getter for the curScan array.
 * @return the current value of the curScan array.
 */
  public int[] getScan(){
    return curScan;
  }
/**
 * Getter for the curIntensity array.
 * @return the current value of the curIntensity  array.
 */
  public int[] getIntensity(){
    return curIntensity;
  }
/**
 * Getter for the bpPos array.
 * @return the current value of the bpPos array.
 */
  public int[] getbpPos(){
    return bpPos;
  }
/**
 * Getter for the curBase array.
 * @return the current value of the curBase array.
 */
  public char[] getBase(){
    return curBase;
  }
/**
 * Getter for the realPeak array.
 * @return the current value of the realPeak array.
 */
  public int[] getRealPeak(){
    return realPeak;
  }
/**
 * Getter for the tracePeakNo array.
 * @return the current value of the tracePeakNo array.
 */
  public int[] getTracePeakNo(){
    return tracePeakNo;
  }
/**
 * Getter for the filter array.
 * @return the current value of the filter array.
 */
  public int[] getFilter(){
    return filter;
  }

/**
 * Top level method to organise matching of trace peaks with amplimer
 */
  public void matchPeaksToRefSeq() {
    
    //  set variables
    int refLen = refSeq.length();
    int prev_peak[][]=new int[refLen+1][ 6];
    
    String baseList[]={"A","C","G","T"};
    int i,j,k,l,m,n,x,y,z,x1;
    int peak,ht,spacing,scan;
    float searchLimit1,searchLimit2,baseoffset;
    int searchType,searchFilter;
    String base,base1,base2,searchBase;
    float b1_no,b2_no,max_filter;
    float b1_curve,b2_curve,b1_slope,b2_slope,b1_int,b2_int;
    float b1_loc,b2_loc,b1_real,b2_real;
    int position=0,startBase=0,index=0,scanInc=0,prevScan=0;
    boolean startedMatching=false;

//  define an upper search limit which may encompass any HOM INS of
//  the set maximum length, n.b. use a reduced average base spacing as
//  we are using a majority of sts's with len < 500 i.e lower av spacing
    int maxSearchLen=Math.max(200,(Constants.AVERAGE_BASE_SPACING-1)*Constants.MAX_POSSIBLE_HOM_INS);
//  n.b the search limit saul below limits the max size of any
//  possible HOM Insertion (n.b. after matching started)
    int sall = 150;           //  search area lower limit (scans)
    int saul = Math.max(275,maxSearchLen);  //  search area upper limit (scans)
    int matchSeqLenCrit=15;   //  success if matched this No. of bases
   
// initialise arrays & initially set all peaks to noise peaks
    for( j = 0 ; j<= nPeaks ; ++j) {
        fork[j] = -1;
        bpPos[j] = -1;
        tracePeakNo[j] = j;
        realPeak[j] = AutoCSA.NOISE_PEAK;
    }
    
//  filter out low peaks by looking at local intensity ratios
    x = 1;
    for( i = 1 ; i<= nPeaks ; ++i) {
        scan = curScan[i];
        for( j = x ; j<= nPeaks ; ++j) {
          if ( curScan[j] > scan - 5 ) { break; }
        }
        x = j;
        x1 = j;
        for( j = x1+1 ; j<= nPeaks ; ++j) {
          if ( curScan[j] > scan + 5 ) { break; }
        }
        x1 = Math.min(j-1,nPeaks);
        max_filter=0.0f;
        for( j = x ; j<= x1 ; ++j) {
          max_filter = Math.max((float) curIntensity[j],max_filter);
        }
        if ( curIntensity[i] > max_filter * 0.2f ) {
            filter[i] = 1;
        }
    }
    
    for( j = 1 ; j<= refLen ; ++j) {
        for( k = 1 ; k<= 5 ; ++k) {
            prev_peak[j][ k] = 0;
        }
    }
    
// n.b. prev_peak[0][ 0] - current max matched base No found
// n.b. prev_peak[0][ 1] - current length of portion of matched bases found
// prev_peak[j][ 1] - 'fork' value
// prev_peak[j][ 2] - index to map base No to arays
// prev_peak[j][ 3] - scan of base
// prev_peak[j][ 4] - intensity of base
// prev_peak[j][ 5] - fork switch
    prev_peak[0][ 1] = 0;       //  max length found
    
// define model parameters for base No Vs scan relationship
    int startOffset=50;
    int startInc=10;
    int startBaseNo=50;

    for ( int mainLoop=0; mainLoop < 2; ++mainLoop ) {

// if we have not matched any amplimer after first pass of this
// maion loop then assume we have a truncation of bases at the
// start of the trace and alter parameterisation so that we use
// the very start of the trace as the search window
    if( mainLoop == 1 ) {
      startOffset=1;
      startInc=0;
    }
    prevScan=0;   // this setting causes trace param model to be used
    int adjustmentType=Adjustment.getAdjustmentType();
    //  loop matchSequence to find first batch of bases in trace
    int refSearchStartLimit=refSearchStart+refSearchStartInc;

    for( i = refSearchStart ; i<= refSearchStartLimit ; ++i) {

// can use estimated start scan if data is mobility corrected
// n.b use fact that base 50 in seq is never before scan 50
// also use a conservative base spacing of 10 from base 50 onwards
// to make sure search bin does not miss required base
        if( adjustmentType == 1 ) {
          prevScan=startOffset+startInc*(i-startBaseNo);
        }
        prev_peak[0][ 0] = i;    //  current max peak position
        for( j = 0 ; j<= nPeaks ; ++j) {
            fork[j] = -1;
            bpPos[j] = -1;
            realPeak[j] = AutoCSA.NOISE_PEAK;
        }
        position = i;
//if(log.isInfoEnabled()) log.info("Search (L1) base no "+position+" of "+refLen);
        matchSequence(position, sall, saul, prev_peak, prevScan, matchSeqLenCrit,true);
        
// assume matching ok if >= matchSeqLenCrit bases are found
        if ( prev_peak[0][ 0] - position >= matchSeqLenCrit ) {
//if(log.isInfoEnabled()) log.info("Breaking "+prev_peak[0][ 0]+" i: "+i);
            startedMatching=true;
            break;   // jump out of i loop
        }
    }
//  write ref seq info once first batch of contiguous bases identified
    if( startedMatching ) {
      startBase=position;
      for( j = position ; j<= prev_peak[0][ 0] ; ++j) {
        if ( prev_peak[j][ 2] > 0 ) {
            index = prev_peak[j][ 2];
            fork[index] = prev_peak[j][ 1];
            realPeak[index] = AutoCSA.CALLED_PEAK;
            bpPos[index] = j;
        }
      }
      sall = 0;             // redefine search limits
      saul = maxSearchLen;  // for efficiency
    }

// assign next base No to search for n.b always assign next
// as this becomes new start base which is always recorded
    position = prev_peak[0][ 0] + 1;  // try next base

  whileloop:
    while ( position < refSearchEnd ) {
        for( i = 1 ; i<= maxBasesMissed  ; ++i) {
            for( j = index + 1 ; j<= nPeaks ; ++j) {
                fork[j] = -1;
                bpPos[j] = -1;
                realPeak[j] = AutoCSA.NOISE_PEAK;
            }
            if( startedMatching ) {
              sall = 0;             // redefine search limits
              saul = maxSearchLen;  // for efficiency
              prevScan=prev_peak[prev_peak[0][0]] [3]; // scan of last base
              matchSeqLenCrit=12;     // relax required consecutive bases
              if( refSearchEnd-position+1 < matchSeqLenCrit ) {
                matchSeqLenCrit=refSearchEnd-position+1;
                if( matchSeqLenCrit < 6 ) { break whileloop; }
              }
            } else {
              if( adjustmentType == 1 ) {
                prevScan=startOffset+startInc*(position-startBaseNo);
              }
            }
//if(log.isInfoEnabled()) log.info("Search (L2) base no "+position+" of "+refLen+" i: "+i+" prevScan: "+prevScan);
            matchSequence(position, sall, saul, prev_peak, prevScan, matchSeqLenCrit,true);
//  save data if we've found >= matchSeqLenCrit further bases
            if ( prev_peak[0][ 0] - position + 1 >= matchSeqLenCrit ) {
                startedMatching=true;
                if( startBase == 0 ) { startBase=position; }
                for( j = position ; j<= prev_peak[0][ 0] ; ++j) {
                    if ( prev_peak[j][ 2] > 0 ) {
                        index = prev_peak[j][ 2];
                        fork[index] = prev_peak[j][ 1];
                        realPeak[index] = AutoCSA.CALLED_PEAK;
                        bpPos[index] = j;
                    }
                }
                position = prev_peak[0][ 0] + 1;
                continue whileloop;
            }
            ++position;
// must check that position is still < refSearchEnd
            if( position >= refSearchEnd ) {break; }
        }
        break;    // break out of while loop
    }

    if( startedMatching ) { break; }
    }  // end of mainloop

    if( startedMatching ) { matchedSomeRefSeq=true; }
    if( !startedMatching ) { return; }
//  now perform 'clean up' to try to fill in any holes with a
//  targeted matching of exactly the required lengths
//  n.b. only if adjacent scans indicate we haven't a HOM DEL
    int refSearchEnd_orig=refSearchEnd;
    int iPrevBase=1;
    sall=0;
    for( i = 2 ; i< nPeaks  ; ++i) {
      if( bpPos[i] == -1 ) { continue; }
      if( iPrevBase == 1 ) { iPrevBase=i; }
      int lenHole=bpPos[i]-bpPos[iPrevBase];
      int scanDiff=curScan[i]-curScan[iPrevBase];
// check for a genuine (i.e. not a HOM DEL) trace hole > 3
      if( lenHole > 4 && scanDiff > Constants.AVERAGE_BASE_SPACING+3 ) {
            position=bpPos[iPrevBase]+2;
            prevScan=curScan[iPrevBase]+Constants.AVERAGE_BASE_SPACING;
            iPrevBase=i;  // reassign here to save every time we continue
            saul=curScan[i]-2*(Constants.AVERAGE_BASE_SPACING-2)-prevScan;
            if ( saul < Constants.AVERAGE_BASE_SPACING ) { continue; }
            matchSeqLenCrit=lenHole-3;     // set to exact len
            if( lenHole > 5 ) {
              matchSeqLenCrit=Math.min(3,matchSeqLenCrit);  // relax a bit
            }
            refSearchEnd=bpPos[i]-2;       // n.b. need to reset
            prev_peak[0][ 0]=position;     // reset to curr base
            matchSequence(position, sall, saul, prev_peak, prevScan,2,true);
            int nFound=prev_peak[0][ 0] - position + 1;
            if ( nFound >= matchSeqLenCrit ) {
              if ( curScan[prev_peak[prev_peak[0][ 0]][ 2]] < curScan[i]) {
                for( j = position ; j<= prev_peak[0][ 0] ; ++j) {
                    if ( prev_peak[j][ 2] > 0 ) {
                        index = prev_peak[j][ 2];
                        fork[index] = prev_peak[j][ 1];
                        realPeak[index] = AutoCSA.CALLED_PEAK;
                        bpPos[index] = j;
                    }
                }
              } else {
                continue;  // local match gone past end of hole - bin it
              }
            } else {
              index=0;
              if( nFound > 1 ) {
                if ( curScan[prev_peak[prev_peak[0][ 0]][ 2]] < curScan[i] ) {
                  for( j = position ; j<= prev_peak[0][ 0] ; ++j) {
                    if ( prev_peak[j][ 2] > 0 ) {
                        index = prev_peak[j][ 2];
                        fork[index] = prev_peak[j][ 1];
                        realPeak[index] = AutoCSA.CALLED_PEAK;
                        bpPos[index] = j;
                    }
                  }
                } else {
                  continue;  // local match gone past end of hole - bin it
                }
              }
              if( lenHole > (5 + nFound) ) {
                position+=(nFound+1);
                if( index > 0 ) {
                  prevScan=curScan[index];   // use last saved index
                } else {
                  prevScan+=(nFound+1)*(Constants.AVERAGE_BASE_SPACING-2);
                }
                saul=curScan[i]-2*(Constants.AVERAGE_BASE_SPACING-2)-prevScan;
                if ( saul < Constants.AVERAGE_BASE_SPACING ) { continue; }
                matchSeqLenCrit-=(nFound+1);  // re-set
                matchSeqLenCrit=Math.max(2,matchSeqLenCrit);
                matchSeqLenCrit=Math.min(3,matchSeqLenCrit);
                prev_peak[0][ 0]=position;     // reset to curr base
                matchSequence(position, sall, saul, prev_peak, prevScan, matchSeqLenCrit,true);
                if ( prev_peak[0][ 0] - position + 1 >= matchSeqLenCrit) {
                  if ( curScan[prev_peak[prev_peak[0][ 0]][ 2]] < curScan[i]) {
                    for( j = position ; j<= prev_peak[0][ 0] ; ++j) {
                      if ( prev_peak[j][ 2] > 0 ) {
                          index = prev_peak[j][ 2];
                          fork[index] = prev_peak[j][ 1];
                          realPeak[index] = AutoCSA.CALLED_PEAK;
                          bpPos[index] = j;
                      }
                    }
                  }
                }
              }
            }
      }
      iPrevBase=i;
    }
    refSearchEnd=refSearchEnd_orig;  // reset original
// finally check if we are > 5 bases from the refSearchStart
// if so then try a reverse match from match start
    if( startBase > refSearchStart+5 ) {
                position=startBase-2;   // skip previous base
                index=getAnalysisIndex(startBase);
                prevScan=curScan[index];
                saul=0;
                sall=5*Constants.AVERAGE_BASE_SPACING;
                matchSeqLenCrit=5;
                prev_peak[0][ 0]=position;     // reset to curr base
                matchSequence(position, sall, saul, prev_peak, prevScan, matchSeqLenCrit,false);
                if ( position - prev_peak[0][ 0] + 1 >= matchSeqLenCrit)
{
                    for( j = position ; j>= prev_peak[0][ 0] ; --j) {
                      if ( prev_peak[j][ 2] > 0 ) {
                          index = prev_peak[j][ 2];
                          fork[index] = prev_peak[j][ 1];
                          realPeak[index] = AutoCSA.CALLED_PEAK;
                          bpPos[index] = j;
                      }
                    }
                }
    }
// finally sweep for single trace holes to check if we can find a peak
// in the gap - probably was lost due to adverse mobility, i.e. too
// close to one of it's neighbours also peaks lost due to Underloading
// N.B. tests to make sure hole is not due to a single base Del
    fillSingleTraceHoles();

  }

  private void matchSequence(int position, int sall, int saul, int[][] prev_peak, int lastBaseScan, int matchSeqLenCrit, boolean forwardSearch) {
    
    int refLen = refSeq.length();
    int peak[][] = new int[refLen+1][6];
    
//  set trace parameterisation coefficients
//  n.b. these are obtained from non-mobility corrected 3730 POP7 data
    float curve,slope,intercept;
// must use trace parameterisation if data is not mobility corrected
    if( Adjustment.getAdjustmentType() == 0 ) {
      curve = 0.00318f;
      slope = 9.56794f;
      intercept = 1171.146f;
    } else {
      curve = 0.0f;
      slope = 0.0f;
      intercept = 0.0f;
    }
   
    int direction,refSearchLimit;
    if ( forwardSearch ) {
      direction=1;
      refSearchLimit=refSearchEnd;
    } else {
      direction=-1;
      refSearchLimit=refSearchStart;
    }
    float pos = (float) position;
    float sac;
    if( lastBaseScan > 0 ) {
      sac= (float) (lastBaseScan+direction*1);   // inc 1 so don't pick up last base
    } else {
      sac = curve*pos*pos + pos*slope + intercept;
      if( sac == 0.0f ) {
        if(log.isWarnEnabled()) log.warn("WARNING: Quadratic Trace parameterisation should NOT be used with zero coefficients");
      }
    }
// define upper and lower search limits
    float sal = Math.max(50.0f, sac - (float) sall);
    float sau = sac + (float) saul;
    
    int i,j,k,l,n;
    int branch,maxPos,scan,bp;
    int index,prevIndex,searchType,searchFilter;
    int resPeak,resScan,resIntensity;
    String base1,base2,base3,searchBase;
    float bin;
    float searchLimit1,searchLimit2,baseOffset;
    boolean selectMax;

//  define index range that fall in search range
//  n.b. curScan[y] is int, sal,sau are float
    int y = 1;
    while( y < nPeaks-1 && (float) curScan[y] < sal ) {
        ++y;
    }
    int lll = y;       // lower line limit  (was y-1)
    
    while( y < nPeaks && (float) curScan[y] < sau ) {
        ++y;
    }
    int ull = y;   // upper line limit
    
    int trialIndex;
    if ( forwardSearch ) { trialIndex = lll; } else { trialIndex = ull; }
    String startBase = refSeq.substring(position-1, position);
    int numTrialBases = countChars(curBase,startBase.charAt(0),lll,ull);
    
    int forkValue[] = new int[nPeaks+1];

// main loop over No. of bases (of type startBase) in range [lll,ull]
  iloop:
    for(i = 1 ; i<=numTrialBases ; ++i) {
        for(j = 0 ; j<=nPeaks ; ++j) {
          forkValue[j] = 0;
        }
        for(j = refSearchStart ; j<=refSearchEnd ; ++j) {
            for(k = 1 ; k<=5 ; ++k) {
                peak[j][ k] = 0;
            }
        }
        branch = 0;
        maxPos = 0;
        index  = 0;
       
// find start base i.e next base of required type in refseq portion
// n.b. all peaks are default noise peaks until matched as called peaks
        if ( forwardSearch ) {
          if( trialIndex > nPeaks ) { return; }
          while( trialIndex < nPeaks && (curBase[trialIndex] != startBase.charAt(0) || filter[trialIndex] != 1 || realPeak[trialIndex] != AutoCSA.NOISE_PEAK) ) {
            ++trialIndex;
          } 
        } else {
          if( trialIndex <= 0 ) { return; } 
          while( trialIndex > 0 && (curBase[trialIndex] != startBase.charAt(0) || filter[trialIndex] != 1 || realPeak[trialIndex] != AutoCSA.NOISE_PEAK) ) {
            --trialIndex;
          }
        }
        
//      set start base info
        scan = curScan[trialIndex];
        bp = position;
        
        peak[bp][ 2] = trialIndex;
        peak[bp][ 3] = scan;
        peak[bp][ 4] = curIntensity[trialIndex];
                
//      search for seq at selected start base
        while ( true) {
  outerwhile:
        while( keepSearching(bp,forwardSearch,refSearchLimit) ) {
            bin = peakSearchBin;
                        
// search (in window) either side of estimated base position (in scans)
            
            if ( forwardSearch ) {
              base1 = refSeq.substring(bp-1, bp);
              base2 = refSeq.substring(bp, bp+1);
              base3 = refSeq.substring(bp+1, bp+2);
            } else {
              base1 = refSeq.substring(bp-1, bp);
              base2 = refSeq.substring(bp-2, bp-1);
              base3 = refSeq.substring(bp-3, bp-2);
            }
            baseOffset=Adjustment.getAdjustmentOffset(base1,base2,bp,bp+1);
            baseOffset*=direction;   // account for fwd or rev search
//          setup criteria for array searching
            
            if ( forkValue[peak[bp][ 2]] > 0 ) {
                peak[bp][ 1] = forkValue[peak[bp][ 2]];
                peak[bp][ 5] = 1;
            }
            
// inner while
            selectMax=false;
            while( bin > 0.0f ) {  // allows smallest bin=0.5
              if ( peak[bp][ 1] == 0 ) {
                searchLimit1 = (float) scan + baseOffset - bin;
                searchLimit2 = (float) scan + baseOffset + bin;
                searchType=1;    // > lim1 && < lim2
              } else if ( peak[bp][ 1] == 1 ) {
                searchLimit1 = (float) scan + baseOffset - bin;
                searchLimit2 = (float) scan + baseOffset;
                searchType=2;    // > lim1 && <= lim2
              } else {
                searchLimit1 = (float) scan + baseOffset;
                searchLimit2 = (float) scan + baseOffset + bin;
                searchType=3;    // >= lim1 && < lim2
              }
              searchFilter = 1;
              searchBase = base2;

              prevIndex=index;
              index=searchData(searchType,searchLimit1,searchLimit2,searchFilter,searchBase,selectMax);
            
//          examine results of peak search
              if ( index < 0 ) {
                if ( index < -1 ) {    // >1 peak found - reduce bin
                //  if ( peak[bp][ 1] == 0 ) {
                //      peak[bp][ 1] = 1;
                //      continue;
                //  } else {
                //      if ( peak[bp][ 1] == 2 ) {
                //          peak[bp][ 1] = 1;
                //      }
                //      bin = bin - 0.5f;
                //      continue;
                //  }
if( index < -2 ) {
//if(log.isInfoEnabled()) log.info("Search3: "+i+" "+sal+" "+sau+" "+searchBase+ " "+searchLimit1+" "+searchLimit2+" "+index);
}
                    bin = bin - 1.0f;  // decrease bin and search again
if( bin == 0.0 ) {
if(log.isInfoEnabled()) log.info("Returning: bp: "+bp+" Retval: "+index);
        if( maxPos == 0 ) {return; }
        if ( direction*(maxPos - position) + 1 >= matchSeqLenCrit ) {
                prev_peak[0][ 0] = maxPos;
                prev_peak[0][ 1] = Math.abs(maxPos - position) + 1;
                for(j = refSearchStart ; j<=refSearchEnd ; ++j) {
                    for(k = 1 ; k<=5 ; ++k) {
                        prev_peak[j][ k] = peak[j][ k];
                    }
                }
        }
        return;   // return here prevents infinite loop if bin=0.0
}
                        continue;
                } else {                // 0 peaks found
                    if( prevIndex == -2 && ! base3.equals(base2) ) {
                      selectMax=true;
                      bin = bin + 1.0f;    // reset previous bin
                      continue;
                    } else {
                      if ( forwardSearch ) { ++bp; } else { --bp; }
                      break outerwhile;   // end current matching process
                    }
                }
              } else {                  // 1 peak found
                resPeak = tracePeakNo[index];
                resScan = curScan[index];
                resIntensity = curIntensity[index];
                if ( branch != 0 ) {
if(log.isInfoEnabled()) log.info("Branch > 0");
                    if ( resPeak != peak[bp + 1][ 2] ) {
                        for(k = bp + 1 ; k<= prev_peak[0][ 0] ; ++k ) {
                            for(l = 1 ; l <=5 ; ++l ) {
                                peak[k][ l] = 0;
                            }
                        }
                    }
                }
                if ( forwardSearch ) { ++bp; } else { --bp; }
                peak[bp][ 2] = resPeak;
                peak[bp][ 3] = resScan;
                peak[bp][ 4] = resIntensity;
                scan = resScan;
                maxPos = bp;
                break;  // break out of inner while i.e. continue matching
              }
            }
        }
            
//      store data for longest sequence found (may need to be >20)
//      if ( maxPos - position + 1 > matchSeqLenCrit && maxPos - position + 1 >= prev_peak[0][ 1] ) {
        if ( maxPos > 0 && direction*(maxPos - position) + 1 >= matchSeqLenCrit ) {
//          if ( maxPos - position + 1 > prev_peak[0][ 1] ) {
            if ( true ) {
                prev_peak[0][ 0] = maxPos;
                prev_peak[0][ 1] = Math.abs(maxPos - position) + 1;
                for(j = refSearchStart ; j<=refSearchEnd ; ++j) {
                    for(k = 1 ; k<=5 ; ++k) {
                        prev_peak[j][ k] = peak[j][ k];
                    }
                }
                return;
            } else {
                if ( maxPos - position + 1 == prev_peak[0][ 1] ) {
                    if ( peak[branch + 1][4] >= prev_peak[branch + 1][4] ) {
                        prev_peak[0][ 0] = maxPos;
                        prev_peak[0][ 1] = maxPos - position + 1;
                        for(j = refSearchStart ; j<=refSearchEnd ; ++j) {
                            for(k = 1 ; k<=5 ; ++k) {
                                prev_peak[j][ k] = peak[j][ k];
                            }
                        }
                        return;
                    }
                }
            }
 //         if ( peak[branch][ 1] >= 2 && peak[branch][ 5] == 0 ) {
 //             peak[branch][ 5] = 1;
 //         }
        }
//      if ( maxPos == refSearchEnd) { return; }  // return if matched to end
        if ( true ) {
              if( forwardSearch ) { ++trialIndex; } else { --trialIndex; }
              continue iloop;
        }
        
        if ( peak[branch][ 1] == 1 && peak[branch][ 5] == 0 ) {
            peak[branch][ 1] = 2;
        } else {
            if ( peak[branch][ 1] == 2 && peak[branch][ 5] == 0 ) {
                peak[branch][ 1] = 1;
                peak[branch][ 5] = 1;
            }
        }
        
        for(j = branch ; j<=maxPos ; ++j ) {
            if ( peak[j][ 1] > 0 && peak[j][ 5] == 1 ) {
                forkValue[peak[j][ 2]] = peak[j][ 1];
            }
        }
        
        while( bp >= 0 ) {
            if( peak[bp][ 1] == 1 && peak[bp][ 5] == 0 ) { break; }
            bp = bp - 1;
            if ( bp == position - 1 ) {
                ++trialIndex;
                continue iloop;
            }
        }
        
        peak[bp][ 1] = 2;
        branch = bp;
        scan = peak[bp][ 3];
        
//  at this point control transfers to start of while(true) loop
//  so we can resume searching for refseq
      }
    }
  }
  private int searchData(int searchType,float searchLimit1,float searchLimit2,int searchFilter,String searchBase, boolean selectMax) {

    int index=0;
    int count=0;
    int curMax=0;
    float scan=0.0f;
    for( int i=1; i<= nPeaks ; ++i) {
      if ( realPeak[i] == AutoCSA.CALLED_PEAK ) {continue; }
      if ( filter[i] != searchFilter ) {continue; }
      if ( curBase[i] != searchBase.charAt(0) ) {continue; }
      scan = (float) curScan[i];
      switch(searchType) {
      case 1: {
//        if ( scan > searchLimit1 && scan < searchLimit2 ) {
          if ( scan >= searchLimit1 && scan <= searchLimit2 ) {
            if( selectMax ) {
              if( curIntensity[i] > curMax ) {
                curMax=curIntensity[i];
                index = i;
              }
              count=1;
            } else {
              index = i;
              ++count;
            }
          }
          break;
        }
      case 2: {
//        if ( scan > searchLimit1 && scan <= searchLimit2 ) {
          if ( scan >= searchLimit1 && scan <= searchLimit2 ) {
            index = i;
            ++count;
          }
          break;
        }
      case 3: {
//        if ( scan >= searchLimit1 && scan < searchLimit2 ) {
          if ( scan >= searchLimit1 && scan <= searchLimit2 ) {
            index = i;
            ++count;
          }
          break;
        }
      }
    }

    switch(count) {
      case 0: {
          return -1;
        }
      case 1: {
          return index;
        }
      default: {
          return -count;
        }
    }
  }
  private int countChars(String source, String match) {

    int count=0;
    for(int i=0 ; i< source.length() ; ++i) {
      if ( (source.substring(i,i+1)).equals(match) ) {
        ++count;
      }
    }
    return count;
  }
  private int countChars(char[] source, char match, int lim1, int lim2) {

    int count=0;
    if ( lim2 >= source.length ) { return -1; }
    for(int i=lim1 ; i<= lim2 ; ++i) {
      if ( source[i] == match && filter[i] == 1 && realPeak[i] == AutoCSA.NOISE_PEAK) {
        ++count;
      }
    }
    return count;
  }
/**
 * Creates a text file of analysis parameters in tabular form
 * @param outFile the name of the output text file
 */
  public void createAnalysisTable(String outFile) {

    if(log.isInfoEnabled()) log.info("CSAFile : "+outFile);
    String heading="|\tBP\tBase\tScan\tHt\tPeak\tReal\tFilter\tQuality";
    File oFile=new File(outFile);
    FileWriter out=null;
    String brow=new String("|\t");
    try {
      oFile.createNewFile();
      out = new FileWriter(oFile);
      out.write(heading);
      out.write('\n');
      for(int i=1 ; i<nPeaks+1 ; ++i ) {
        String row=brow+bpPos[i]+"\t"+curBase[i]+"\t"+curScan[i]+"\t"+curIntensity[i]+"\t"+tracePeakNo[i]+"\t"+realPeak[i]+"\t"+filter[i]+"\t"+quality[i];
        out.write(row);
        out.write('\n');
      }
      out.close();
    } catch ( IOException e ) {
      if(log.isWarnEnabled()) log.warn("Caught IOException: "+e.getMessage());
    }

  }
/**
 * Reports to stdout a sumary of the whole and ROI (if set) coverage
 * @param roiStartCoord start base of ROI
 * @param roiEndCoord end base of ROI
 * @param refSeqLength length of amplimer sequence
 */
  public void reportCoverage(int roiStartCoord, int roiEndCoord, int refSeqLength) {
    // find extent of ref seq analysed
    int i;
    float roiCoverage= 0.0f;
    float stsCoverage= 0.0f;
    for( i=0; i < nPeaks+1; ++i ) {
      if( bpPos[i] > -1 ) { break; }
    }
    if( i <= nPeaks ) {
      int bp1=bpPos[i];
      int bp2=bpPos[nPeaks];
      int lim1=Math.min(bp1,roiStartCoord);
      int lim2=Math.min(bp2,roiEndCoord);
      int ourRange=lim2-lim1+1;
      int roiRange=roiEndCoord-roiStartCoord+1;
      roiCoverage= 100.0f* (float) ourRange/ (float) roiRange;
      roiCoverage=Math.min(100.0f,roiCoverage);
      stsCoverage= 100.0f* (float) (bp2-bp1+1)/ (float) refSeqLength;
    }
    if(log.isInfoEnabled()) log.info("ROI Coverage: "+roiCoverage+" %");
    if(log.isInfoEnabled()) log.info("STS Coverage: "+stsCoverage+" %");
  }
/**
 * <p>Sets the individual base trace quality array.</p>
 * <p>Note the analysis datasets must have been sorted before calling this method.</p>
 * <p>Method uses a per base Signal to Noise ratio.</p>
 * @param peakRange Deprecated not required
 */
  public void setTraceQuality(int peakRange) {
// n.b datasets must have been sorted wrt bpPos, this puts noise peaks
// at start of array so cannot use data until we've reached the first
// index corresponding to a real peak
// Trace base quality q is defined as signal/noise ratio at the base
    int i;
    int peakIncr=(peakRange-1)/2;
    for( i=0; i < nPeaks+1; ++i ) {
      if( bpPos[i] > -1 ) { break; }
    }
    int realPeakStart=i;
    float maxQ=maxAllowableQuality;
    for( i=0; i < nPeaks+1; ++i ) {
      if( realPeak[i] != AutoCSA.CALLED_PEAK ) {
        quality[i]=-1.0f;
        continue;
      }
      // max called peak is actually intensity of the base itself
      int maxCalledPeak=curIntensity[i];
      // find max uncalled peak in window of +-5 around base
      int peakScan=curScan[i];
      int scanLim1=peakScan-5;
      int scanLim2=peakScan+5;
      int maxUnCalledPeak=getMaxNoisePeak(scanLim1,scanLim2);
      if( maxUnCalledPeak > 0 ) {
        quality[i]= Math.min(maxQ,(float) maxCalledPeak/ (float) maxUnCalledPeak);
//if(log.isInfoEnabled()) log.info("Q: "+quality[i]+" "+maxCalledPeak+" "+maxUnCalledPeak);
      } else {
        int chan=getTraceObj().getChannelIndex(String.valueOf(curBase[i]));
        int maxNoise=-30000;
        for (int j=0; j< 4; ++j ) {
          if( j == chan ) { continue; }
          int noise= getTraceObj().getChannel(j).getMaxIntensity(peakScan-1,peakScan+1);
          maxNoise= Math.max(maxNoise,noise);
        }
        if( maxNoise < 1 ) { maxNoise=1; }  // prevents divide by 0
        quality[i]= Math.min(maxQ,(float) maxCalledPeak/ (float) maxNoise);
      }
    }
  }
/**
 * <p>Returns a trace quality array calculated using peak data only.</p>
 * <p>Note the analysis datasets must have been sorted before calling
this method.</p>
 * <p>Method uses a Peak Signal to Peak Noise ratio at a base location.</p>
 * @param maxQual the max allowable quality value
 * @param bin the half length of the search window for noise peaks
 * @return float array of quality values
 */
  public float[] getTraceQualityUsingPeaksOnly(float maxQual, int bin) {
// n.b datasets must have been sorted wrt bpPos, this puts noise peaks
// at start of array so cannot use data until we've reached the first
// index corresponding to a real peak
// Trace base quality q is defined as signal/noise ratio at the base
// n.b. this method returns an alternative q array calculated with a
// noise peaks bases on peaks only - used for subtle het detection
    int i;
    float[] qual = new float[nPeaks+1];
    for( i=0; i < nPeaks+1; ++i ) {
      if( realPeak[i] != AutoCSA.CALLED_PEAK ) {
        qual[i]=-1.0f;
        continue;
      }
      // max called peak is actually intensity of the base itself
      int maxCalledPeak=curIntensity[i];
      // find max uncalled peak in window of +- bin around base
      int peakScan=curScan[i];
      int scanLim1=peakScan-bin;
      int scanLim2=peakScan+bin;
      int maxUnCalledPeak=getMaxNoisePeak(scanLim1,scanLim2);
      if( maxUnCalledPeak > 0 ) {
        qual[i]= Math.min(maxQual,(float) maxCalledPeak/ (float) maxUnCalledPeak);
      } else {
        qual[i]= maxQual;
      }
    }
    return qual;
  }
/**
 * <p>Returns a trace quality array calculated over a window.</p>
 * <p>Note the analysis datasets must have been sorted before calling this method.</p>
 * <p>Method uses a Signal to Noise ratio over a window of bases.</p>
 * @param peakRange the window size (bases)
 * @param maxQual the max allowable quality value
 * @return float array of quality values
 */
  public float[] getTraceQualityOverWindow(int peakRange, float maxQual) {
// n.b datasets must have been sorted wrt bpPos, this puts noise peaks
// at start of array so cannot use data until we've reached the first
// index corresponding to a real peak
// n.b. this method returns an alternative q array calculated with a
// signal/noise ratio over a window of bases - used for indel Pos detection
    int i;
    int peakIncr=(peakRange-1)/2;
    for( i=0; i < nPeaks+1; ++i ) {
      if( bpPos[i] > -1 ) { break; }
    }
    int realPeakStart=i;
    float[] qual = new float[nPeaks+1];
    for( i=0; i < nPeaks+1; ++i ) {
      if( realPeak[i] != AutoCSA.CALLED_PEAK && realPeak[i] != AutoCSA.DYE_BLOB_PEAK ) {
        qual[i]=-1.0f;
        continue;
      }
      int lim1=Math.max(i-peakIncr,realPeakStart);
      int lim2=Math.min(i+peakIncr,nPeaks);
      // find max called peak in window of peakRange peaks centred on current
      int maxCalledPeak=0;
      int maxIndex=0;
      for( int j=lim1; j <= lim2; ++j ) {
        maxCalledPeak=Math.max(maxCalledPeak,curIntensity[j]);
        if ( maxCalledPeak == curIntensity[j] ) {maxIndex=j; }
      }
      // find max uncalled peak in same window (add a bit extra)
      int scanLim1=curScan[lim1]-5;
      int scanLim2=curScan[lim2]+5;
      int maxUnCalledPeak=getMaxNoisePeak(scanLim1,scanLim2);
      if( maxUnCalledPeak > 0 ) {
        qual[i]= Math.min(maxQual,(float) maxCalledPeak/ (float) maxUnCalledPeak);
      } else {
        int chan=getTraceObj().getChannelIndex(String.valueOf(curBase[maxIndex]));
        int maxNoise=-30000;
        for (int j=0; j< 4; ++j ) {
          if( j == chan ) { continue; }
          maxNoise= Math.max(maxNoise,getTraceObj().getChannel(j).getMaxIntensity(curScan[maxIndex],curScan[maxIndex]));
        }
        if( maxNoise < 1 ) { maxNoise=1; }  // prevents divide by 0
//  if(log.isInfoEnabled()) log.info("Max Noise: "+i+" "+bpPos[i]+" "+curBase[i]+" "+maxNoise);
        qual[i]= Math.min(maxQual,(float) maxCalledPeak/ (float) maxNoise);
      }
    }
    return qual;
  }
/**
 * Returns the trace quality array.
 * @return float array of quality values
 */
  public float[] getTraceQuality() {
    return quality;
  }
/**
 * <p>Returns the average trace quality.</p>
 * <p>Note calculated over the ROI if set.</p>
 * @return the average trace quality
 */
  public float getAverageTraceQuality() {
    int lim1,lim2,count=0;
    if( ROIStartCoord > 0 ) {
      int[] lims=getROIArrayIndices();
      if ( lims[0] == 0 ) { return 0.0f; } // no ROR covered
      lim1=lims[0];
      lim2=lims[1]+1;
    } else {
      lim1=0;
      lim2=nPeaks+1;
    }
    float avQuality=0.0f;
    for( int i=lim1; i < lim2; ++i ) {
      if( quality[i] > 0.0f ) {
        avQuality+=quality[i];
        ++count;
      }
    }
    if ( count == 0 ) {
      return 0.0f;
    } else {
      return avQuality/ (float) count;
    }
  }
/**
 * <p>Returns the % of bases with quality above a critical value.</p>
 * <p>Note calculated over the ROI if set.</p>
 * @param critQuality the critical quality value
 * @return the % of bases with q > critQuality
 */
  public float getQualityProfile(float critQuality) {
    int lim1,lim2,count=0;
    if( ROIStartCoord > 0 ) {
      int[] lims=getROIArrayIndices();
      if ( lims[0] == 0 ) { return 0.0f; } // no ROR covered
      lim1=lims[0];
      lim2=lims[1]+1;
    } else {
      lim1=0;
      lim2=nPeaks+1;
    }
    float avQuality=0.0f;
    for( int i=lim1; i < lim2; ++i ) {
      if( quality[i] > 0.0f ) {
        if( quality[i] > critQuality ) { avQuality+=1.0f; }
        ++count;
      }
    }
// returns the % of bases with q > critQuality
    if ( count == 0 ) {
      return 0.0f;
    } else {
      return 100.0f * avQuality/ (float) count;
    }
  }
  private int getMaxNoisePeak(int limit1, int limit2) {
    int maxNoisePeak=0;
    for( int i=1; i < nPeaks+1; ++i ) {
      if( realPeak[i] != AutoCSA.NOISE_PEAK ) { continue; }  // selects only unassigned
      if( curScan[i] >= limit1 && curScan[i] <= limit2 ) {
        maxNoisePeak=Math.max(maxNoisePeak,curIntensity[i]);
      }
    }
// returns the max noise peak over a scan range
    return maxNoisePeak;
  }
  private int getMaxNoisePeak(int limit1, int limit2, String base) {
    int maxNoisePeak=0;
    for( int i=1; i < nPeaks+1; ++i ) {
      if( realPeak[i] != AutoCSA.NOISE_PEAK ) { continue; }  // selects only unassigned
      if( ! base.equals(String.valueOf(curBase[i])) ) {continue; }
      if( curScan[i] >= limit1 && curScan[i] <= limit2 ) {
        maxNoisePeak=Math.max(maxNoisePeak,curIntensity[i]);
      }
    }
// returns the max noise peak (of a specified base) over a scan range
    return maxNoisePeak;
  }
/**
 * <p>Returns analysis array index of the specified base.</p>
 * <p>Note if the base is not in the analysis then returns -1.</p>
 * @param targetPos the specified base
 * @return the analysis array index
 */
  public int getAnalysisIndex(int targetPos) {
    int i;
    for( i=0; i < nPeaks+1; ++i ) {
      if( bpPos[i] == targetPos ) { break; }
    }
// returns index of a base position else -1 if position not analysed
    return i < nPeaks+1 ? i : -1;
  }
/**
 * <p>Returns an array containing the analysis limits.</p>
 * <p>Note if the base is not in the analysis then returns -1.</p>
 * @return the int[2] array containing the start and end analysis base No
 */
  public int[] getAnalysisLimits() {
    int i;
    int[] limits= new int[2];
    for( i=0; i < nPeaks+1; ++i ) {
      if( bpPos[i] > -1 ) { break; }
    }
    if( i > nPeaks ) {
      limits[0]=0;
      limits[1]=0;
    } else {
      limits[0]=bpPos[i];
      limits[1]=bpPos[nPeaks];
    }
    return limits;
  }
/**
 * <p>Setter of ArrayList of trace holes.</p>
 * <p>Note trace holes are amplimer bases not matched with trace peaks.</p>
 * @param holes the ArrayList of trace holes
 */
  public void setTraceHoles(ArrayList holes) {
    traceHoles=holes;
  }
/**
 * <p>Getter of ArrayList of trace holes.</p>
 * <p>Note trace holes are amplimer bases not matched with trace peaks.</p>
 * @return the ArrayList of trace holes
 */
  public ArrayList getTraceHoles() {
    return traceHoles;
  }
/**
 * <p>Setter of class analysis parameters.</p>
 * @param bin the scan length to use for the window for peak searching
 * @param refstart the start base No for amplimer matching
 * @param refend the end base No for amplimer matching
 * @param maxmiss the max allowable No of bases to miss in the amplimer
 * @param mindist Deprecated has no function
 */
  public void setAnalysisParams(float bin, int refstart, int refend, int refinc, int maxmiss, int mindist) {
    peakSearchBin=bin;
    refSearchStart=refstart;
    refSearchEnd=refend;
    refSearchStartInc=refinc;
    maxBasesMissed=maxmiss;
    minPeakSpacing=mindist;
  }
/**
 * <p>Calculate and get of array indices corresponding to the set ROI.</p>
 * <p>Note returns [0,0] if the ROI is not set.</p>
 * @return the array indices of the ROI coordinates
 */
  public int[] getROIArrayIndices() {
    int i,j,index=0;
    int[] indices = new int[] {0,0};
    int[] roi = new int[] {ROIStartCoord,ROIEndCoord};
    if ( ROIStartCoord == 0 ) { return indices; }
//  check we have covered at least some of ROI
    int[] lims=getAnalysisLimits();
    if ( ROIStartCoord > lims[1] ) { return indices; }
    for( j=0; j < 2; ++j ) {
      for( i=roi[j]; i > 0; --i ) {
        index=getAnalysisIndex(i);
        if( index > -1 ) { break; }
      }
      indices[j]=index;
      if( index == -1 ) {
        for( i=roi[j]+1; i < 1000; ++i ) {
          index=getAnalysisIndex(i);
          if( index > -1 ) { break; }
        }
        indices[j]=index;
      }
      if( indices[j] == -1 ) { indices[j]=0; }
    }
// return indices of ROR or nearest analysed bases
    return indices;
  }
/**
 * <p>Return an integer code for the trace status.</p>
 * <p>Note any trace status is based on the ROI if set.</p>
 * <p>Note the trace status integer codes are set in {@link Constants}.</p>
 * @return integer code trace status
 */
  public int getTraceStatus() {
// n.b this routine uses SNP database codes from TRACE_STATUS_DICT
//            1 Ready for rerun analysis
//            6 Rerun - Underloaded
//            7 Rerun - Overloaded
//            8 Rerun - Failed to locate parent seq
    int code=Constants.ORA_ID_RERUN_ANALYSIS;    // indicates trace is OK
    int over=0;
    int under=0;
    int count=0;
    int lim1, lim2;
    if( ROIStartCoord > 0 ) {
      int[] lims=getROIArrayIndices();
      lim1=lims[0];
      lim2=lims[1]+1;
    } else {
      lim1=0;
      lim2=nPeaks+1;
    }
//  int critULIntensity=getUnderloadIntensity();
    for( int i=lim1; i < lim2; ++i ) {
      if( realPeak[i] == 1 ) {
        if( curIntensity[i] > Constants.INTENSITY_OVERLOADED ) { ++over; }
//      if( curIntensity[i] < critULIntensity ) { ++under; }
        if( isBaseUnderloadedComplex(curIntensity[i],String.valueOf(curBase[i])) ) ++under;
        ++count;
      }
    }
    float critOverloadPC,critUnderloadPC;
    if( isTumourSample ) {
      critOverloadPC=Constants.CRIT_OVERLOAD_PC_TUMOUR;
      critUnderloadPC=Constants.CRIT_UNDERLOAD_PC_TUMOUR;
    } else {
      critOverloadPC=Constants.CRIT_OVERLOAD_PC_NORMAL;
      critUnderloadPC=Constants.CRIT_UNDERLOAD_PC_NORMAL;
    }
// n.b. flagging under/over loaded more important than locating seq
    if( count == 0 ) {
      if( matchedSomeRefSeq ) {
        code=Constants.ORA_ID_RERUN_NO_ROR_COVERED;
      } else {
        code=Constants.ORA_ID_RERUN_MATCH_FAILED;
      }
    } else {
      if( 100.0* (float) under/ (float) (count) > critUnderloadPC ) {
        code=Constants.ORA_ID_RERUN_UNDERLOADED;
      }
      if( 100.0* (float) over/ (float) (count) > critOverloadPC ) {
        code=Constants.ORA_ID_RERUN_OVERLOADED;
      }
    }
    return code;  // return trace status code (n.b. see Constants.java)
  }
  private boolean isBaseUnderloadedComplex(int intensity, String base) {
    int crit=Constants.INTENSITY_UNDERLOADED;
    int chan=getTraceObj().getChannelIndex(base);
    if( minPeakIntensity[chan] < Constants.STRICT_MIN_PEAK ) {
      crit=(int) (Constants.UL_DECREASE_FACTOR * (float) Constants.INTENSITY_UNDERLOADED);
    }
    if( intensity < crit ) return true;
    return false;
  }
  private int getUnderloadIntensity(){
// n.b datasets must have been sorted wrt bpPos
    int i,critIntensity;
    float avIntensity=0.0f;
    int count=0;
    for( i=0; i < nPeaks+1; ++i ) {
      if( bpPos[i] == -1 ) { continue; }
      avIntensity+=(float) curIntensity[i];
      ++count;
    }
    avIntensity/=(float) count;
    if( avIntensity < 5000.0f ) {
      critIntensity=(int) (0.5f*(float)Constants.INTENSITY_UNDERLOADED);
    } else if( avIntensity < 10000.0f ) {
      critIntensity=(int) (0.75f*(float)Constants.INTENSITY_UNDERLOADED);
    } else {
      critIntensity=Constants.INTENSITY_UNDERLOADED;
    }
    critIntensity=Constants.INTENSITY_UNDERLOADED; // 1 value used
    return critIntensity;
  }
/**
 * <p>Setter for the Trace object associated with this class.</p>
 * @param trace the trace object to associate
 */
  public void setTraceObj(Trace trace){
    traceObj=trace;
  }
/**
 * <p>Getter for the Trace object associated with this class.</p>
 * @return the associated trace object
 */
  public Trace getTraceObj(){
    return traceObj;
  }
/**
 * <p>Getter for the int array of scan indices of base positions.</p>
 * <p>Note the class array bpPos must have been sorted.</p>
 * @return the int array of base scan indices
 */
  public int[] getBaseScans(){
// n.b datasets must have been sorted wrt bpPos
    int i;
    for( i=0; i < nPeaks+1; ++i ) {
      if( bpPos[i] > -1 ) { break; }
    }
    int num=nPeaks+1-i;
    int[] scan  = new int[num];
    System.arraycopy(curScan,i,scan,0,num);
    return scan;
  }
/**
 * <p>Calculate and return a value for the average base spacing over the specified number of proceeding bases.</p>
 * <p>Note returns -1.0 if calculation tries to use unanalysed regions.</p>
 * @param pos the base position
 * @param nBases the number of proceeding bases
 * @return the average base spacing
 */
  public float getAverageBaseSpacing(int pos, int nBases) {
// n.b datasets must have been sorted wrt bpPos
// n.b can't be gaps to interfere with region [pos-nBases+1,pos]
// as must have a region of at least 15 consecutive bases before gap
// n.b will also never creep into unanalysed data if nBases < 20
    int iStart=getAnalysisIndex(pos);
    if( bpPos[iStart-nBases+1] == -1 ) {
      return -1.0f;      // crept into unanalysed data
    }
    float spacing= (float) (curScan[iStart]-curScan[iStart-nBases+1])
                   / (float) (nBases-1);
    return spacing;
  }
/**
 * <p>Calculate and return the number of Viable peaks in a scan index range .</p>
 * <p>Note A viable peak is defined as a peak deemed by CSA to be potentially associated with an amplimer base, and a peak in the scan range must be > 20 % of the max peak in the same range.</p>
 * @param scan1 the start scan index
 * @param scan2 the end scan index
 * @param viableCrit Depracated Not currently used
 * @return the number of viable peaks
 */
  public int getNumViablePeaks(int scan1, int scan2, int viableCrit) {
    int count=0;
    int i;
    ArrayList peaks = new ArrayList();
    for( i=1; i < nPeaks+1; ++i ) {
// selects only unassigned and peaks passed filtering
      if( realPeak[i] != AutoCSA.NOISE_PEAK || filter[i] != 1) { continue; }
      if( curScan[i] > scan1 && curScan[i] < scan2 ) {
        ++count;
        peaks.add(new Integer(curIntensity[i]));
      }
    }
// sort peakHt array and test ratios to estimate which are real
    if( count > 1 ) {
      int[] peakHt = new int[count];
      for( i=0; i < count; ++i ) {
        peakHt[i]=((Integer) peaks.get(i)).intValue();
      }
      Arrays.sort(peakHt);
      int total=count;
      for( i=0; i < total-1; ++i ) {
        if( (float) peakHt[total-1]/ (float) peakHt[i] > 5.0f ) {
          --count;
        }
      }
    }
// return No. of peaks in scan range which could be assumed real
    return count;
  }
/**
 * <p>Gets the quality value for a specified base number.</p>
 * <p>Note returns -1 if base has not been analysed.</p>
 * @param pos the base position
 * @return the base quality value (Signal : Noise ratio)
 */
  public float getBaseQuality(int pos) {
    int index=getAnalysisIndex(pos);
    return index > -1 ? quality[index] : -1.0f;
  }
/**
 * <p>Gets the scan index for a specified base number.</p>
 * <p>Note returns -1 if base has not been analysed.</p>
 * @param pos the base position
 * @return the base scan index value
 */
  public int getScanAtBase(int pos) {
    int index=getAnalysisIndex(pos);
    return index > -1 ? curScan[index] : -1;
  }
/**
 * <p>Sets the ROI coordinates for this class.</p>
 * @param start the starting ROI coordinate
 * @param end the ending ROI coordinate
 */
  public void setROICoords(int start, int end) {
      ROIStartCoord=start;
      ROIEndCoord=end;
  }
/**
 * <p>Calculates and gets the true local quality for a specified base
type and region of a trace.</p>
 * <p>Note here true local quality refers to the fact that we ignore any mutant (het) peaks in the range and consider only genuine noise.</p>
 * @param pos the base position
 * @param base the base type (channel) to consider
 * @param windowHalfLen the half length of the window around the base
 * @return the true local quality
 */
  public float getLocalTrueChannelQuality(int pos, String base, int windowHalfLen) {
// This method calculates the local trace quality over a base window
// for a specified base type ignoring any possible mutant peaks
//  - hence gives a genuine trace quality (for a channel)
    float quality;
    int index=getAnalysisIndex(pos);
    int lim1=pos;
    int lim2=pos;
    int curpos=pos;
    for( int i=pos+1; i < pos+windowHalfLen+1; ++i ) {
      if( bpPos[Math.min(nPeaks,++index)] == ++curpos ) { ++lim2; }
    }
    index=getAnalysisIndex(pos);
    curpos=pos;
    for( int i=pos-1; i > pos-windowHalfLen-1; --i ) {
      if( bpPos[Math.max(0,--index)] == --curpos ) { --lim1; }
    }
    lim1=getAnalysisIndex(lim1);
    lim2=getAnalysisIndex(lim2);
    float maxQ=maxAllowableQuality;
      // find max called peak in window of peakRange peaks centred on current
      int maxCalledPeak=0;
      int maxIndex=0;
      for( int j=lim1; j <= lim2; ++j ) {
        if( ! base.equals(String.valueOf(curBase[j])) ) {continue; }
        maxCalledPeak=Math.max(maxCalledPeak,curIntensity[j]);
        if ( maxCalledPeak == curIntensity[j] ) {maxIndex=j; }
      }
      // find max uncalled peak in same window
      int scanLim1=curScan[lim1];
      int scanLim2=curScan[lim2];
      int maxUnCalledPeak=getMaxNoisePeak(scanLim1,scanLim2,base);
//if(log.isInfoEnabled()) log.info("True Quality: Pos: "+pos+" "+base+" "+maxCalledPeak+" "+maxUnCalledPeak);
//    if( maxUnCalledPeak == 0 ) {
//      return maxQ;     // no noise in channel so return high q
//    }
      if( maxCalledPeak == 0 ) {
        return 100.0f;     // no like peak (to base) found in window
      }
      if( maxUnCalledPeak > 0 ) {
        quality= Math.min(maxQ,(float) maxCalledPeak/ (float) maxUnCalledPeak);
      } else {
//      int chan=getTraceObj().getChannelIndex(base);
        int maxNoise=-30000;
//      maxNoise= Math.max(maxNoise,getTraceObj().getChannel(chan).getMaxIntensity(scanLim1,scanLim2));
        maxNoise=getMinPeakIntensity(base);
        if( maxNoise < 1 ) { maxNoise=1; }  // prevents divide by 0
        quality= Math.min(maxQ,(float) maxCalledPeak/ (float) maxNoise);
      }
//if(log.isInfoEnabled()) log.info("True Quality: Pos: "+pos+" "+quality);
    return quality;
  }
  public int[] searchChannelForPeak(String searchBase, int searchLimit1,int searchLimit2, float intensityFactor) {
      int chan=getTraceObj().getChannelIndex(searchBase);
      int[] peakInfo = getTraceObj().getChannel(chan).searchForPeak(searchLimit1,searchLimit2,intensityFactor);
      if ( peakInfo[0] > -1 ) {
        return peakInfo;
      }
      return new int[] {-1,-1};
  }
/**
 * <p>Searches for Dye Blobs over a specified base range, and if found then resets the base q to zero, returns the number of bases affected.</p>
 * @param windowStart the start base to search for dye blobs
 * @param windowWidth the window width (added to windowStart to define a
window) to search for dye blobs (bases)
 * @param critIntensityRatio multiplier for the average channel intensity to define a critical intensity for a dye blob
 * @return the number of bases affected by dye blob
 */
  public int searchForDyeBlobs(int windowStart, int windowWidth, float critIntensityRatio) {
    int count=0;
    int i,index1=-1,index2=-1;
// find analysis indices of dye blob extremities
// n.b. must allow for trace holes at the extremities
    for( i=windowStart; i > 0; --i ) {
      index1=getAnalysisIndex(i);
      if( index1 > -1 ) { break; }
    }
    for( i=windowStart+windowWidth; i < 1000; ++i ) {
      index2=getAnalysisIndex(i);
      if( index2 > -1 ) { break; }
    }
//  return if we can't identify dye blob search region
    if( index1 == -1 || index2 == -1 ) { return -1; }
// derive average peak intensities in region after DYE_BLOB_WINDOW
// for each base - if failure for a base can't test for that base type
// n.b. avearge is over 20 bases - be carefull to avoid next dye blob
    float avIntensityA=getAveragePeakIntensity('A',index2,index2+20);
    float avIntensityC=getAveragePeakIntensity('C',index2,index2+20);
    float avIntensityG=getAveragePeakIntensity('G',index2,index2+20);
    float avIntensityT=getAveragePeakIntensity('T',index2,index2+20);
// identify dye blobs over search area
    for( i=index1; i <= index2; ++i ) {
      switch (curBase[i]) {
      case 'A': {
        quality[i]=getModifiedBaseQuality(i,critIntensityRatio*avIntensityA);
        break;
      }
      case 'C': {
        quality[i]=getModifiedBaseQuality(i,critIntensityRatio*avIntensityC);
        break;
      }
      case 'G': {    //  n.b. allow for small G's i.e. raise crit
        quality[i]=getModifiedBaseQuality(i,1.4f*critIntensityRatio*avIntensityG);
        break;
      }
      case 'T': {
        quality[i]=getModifiedBaseQuality(i,critIntensityRatio*avIntensityT);
        break;
      }
      }
      if ( quality[i] == 0.0f ) {
        ++count;
        realPeak[i]=AutoCSA.DYE_BLOB_PEAK;
      }
    }
    return count;
  }
  private float getModifiedBaseQuality(int index,float critIntensity ) {
    int i=index;
    float originalQ=quality[i];
    if( critIntensity < 0.0f ) { return originalQ; }
// if neighbours are both high quality (q > 3) and are properly
// adjacent in sequence then return original Q (as isolated high peak)
    if( quality[i-1] > 3.0f && quality[i+1] > 3.0f ) {
      if( bpPos[i]-bpPos[i-1] == 1 && bpPos[i+1]-bpPos[i] == 1 ) {
        return originalQ;
      }
    }
    if( curIntensity[i] > critIntensity ) {
      return 0.0f;   // set base q to zero as a suspected dye blob
    }
    return originalQ;
  }
  private float getAveragePeakIntensity(char base, int index1, int index2) {
    int count=0, total=0;
    float nullAverage=-1.0f;
    index2=Math.min(index2,nPeaks);  // safety check on final index
    for( int i=index1; i <= index2; ++i ) {
      if (curBase[i] != base ) { continue; }
      total+=curIntensity[i];
      ++count;
    }
// return average peak intensity over index range (if > 2 peaks)
    if (count < 3 ) {
      return nullAverage;
    } else {
      return (float) total/(float) count;
    }
  }
/**
 * Sets the base quality value to zero at every Underloaded base
 */
  public void setModifiedQwhereUnderloaded() {
    int critIntensity=getUnderloadIntensity();
    for( int i=0; i < nPeaks+1; ++i ) {
      if( bpPos[i] == -1 ) { continue; }
      if( curIntensity[i] < critIntensity ) {
        quality[i]=0.0f;  // modify q if base is Underloaded
      }
    }
  }
/**
 * Sets the base quality value to zero at every Overloaded base
 */
  public void setModifiedQwhereOverloaded() {
    int critIntensity=Constants.INTENSITY_OVERLOADED;
    for( int i=0; i < nPeaks+1; ++i ) {
      if( bpPos[i] == -1 ) { continue; }
      if( curIntensity[i] > critIntensity ) {
        quality[i]=0.0f;  // modify q if base is Overloaded
      }
    }
  }
/**
 * Returns a boolean to indicate whether base is affected by Dye Blob
 * @param baseNum the base number
 * @return boolean indicating if base is Dye Blob affected
 */
  public boolean isBaseDyeBlob( int baseNum) {
    int i=getAnalysisIndex(baseNum);  // base definately exists
    if ( i == -1 ) return false;      // safety net
    if( realPeak[i] == AutoCSA.DYE_BLOB_PEAK ) {
      return true;
    } else {
      return false;
    }
  }
/**
 * Returns a boolean to indicate whether base is Underloaded.
 * Note If base is not analysed then returns false.
 * @param targetPos the base number
 * @return boolean indicating if base is Underloaded
 */
  public boolean isBaseUnderloaded(int targetPos) {
    int index=getAnalysisIndex(targetPos);
    if( index == -1 ) { return false; }
    if( curIntensity[index] < getUnderloadIntensity() ) {
      return true;
    } else {
      return false;
    }
  }
/**
 * Returns a boolean to indicate whether base is Overloaded.
 * Note If base is not analysed then returns false.
 * @param targetPos the base number
 * @param critIntensity if this value > 0 then this is used as the critical intensity
 * @return boolean indicating if base is Overloaded
 */
  public boolean isBaseOverloaded(int targetPos, int critIntensity) {
    int index=getAnalysisIndex(targetPos);
    if( index == -1 ) { return false; }
    int critical=Constants.INTENSITY_OVERLOADED;
    if ( critIntensity > 0 ) critical=critIntensity; // reset if supplied
    if( curIntensity[index] > critical ) {
      return true;
    } else {
      return false;
    }
  }
/**
 * Returns an indication of whether there is an Het Indel in the trace.
 * Method uses a ratio of the average quality over left and right windows.
 * Note If there is not enough data points to work with then returns zero.
 * @param roiOnly boolean specifying whether to consider ROI region only
 * @return position of start of suspected Indel, otherwise zero
 */
  public int checkQualityProfileForInDel(boolean roiOnly) {
    int i,j,jj;
    int pos=0;
    int startLhsWindow=10;
    int windowLim=30;
    float[] qual = getTraceQualityOverWindow(5,15.0f);
    int[] lims=getAnalysisLimits();
    int startBP=lims[0];
    if( startBP == 0 ) { return 0; }
//  if( windowLim > 0 && lims[1]-lims[0]+1 < 2*windowLim ) { return 0; }
    for( i=0; i < qual.length; ++i ) {
      if( qual[i] > 0.0f ) { break; }
    }
    int start=i;
    int len=qual.length-start;
    if( len < 2*startLhsWindow+10 ) { return 0; }
    float[] wQual = new float[len];
    float[] qtemp = new float[5];
    for( i=start; i < start+len; ++i ) {
        wQual[pos++]=qual[i];
    }
    int nData=len;
    float avQualityLhs=0.0f;
    float avQualityRhs=0.0f;
    int numLhs=startLhsWindow;
    int numRhs=nData-startLhsWindow;
    int nQualComparisons=nData-2*startLhsWindow;
//  if( windowLim > 0 ) { nData=startLhsWindow+windowLim; }
    for( i=0; i < startLhsWindow; ++i ) {
        avQualityLhs+=wQual[i];
    }
    for( i=startLhsWindow; i < Math.min(nData,startLhsWindow+windowLim); ++i ) {
        avQualityRhs+=wQual[i];
    }
    numRhs=Math.min(numRhs,windowLim);
    avQualityLhs/= (float) numLhs;
    avQualityRhs/= (float) numRhs;
    float firstQualL=avQualityLhs;
    float firstQualR=avQualityRhs;
//  derive a critical Qav ratio based on Av Q in rhs window
    float critRatio=2.0f;

    int curIndex=startLhsWindow;
    pos=startBP+startLhsWindow;
//if(log.isInfoEnabled()) log.info("Q Pos: "+pos+" L "+avQualityLhs+" R: "+avQualityRhs);
    float maxRatio=0.1f;    // initialise just > 0.0
    int maxPos=0;
    int critLen=5000;
    for( i=0; i < nQualComparisons ; ++i ) {
      if( i >= nQualComparisons ) break;    // fix for glitch in JVM
      if( maxPos == 0 ) {
        critRatio=2.0f;
        if( avQualityLhs > 5.0f ) critRatio=2.5f;
        if( avQualityLhs > 7.5f ) critRatio=3.0f;
        if( avQualityLhs > 10.0f ) critRatio=3.5f;
      }
// must take account of trace holes by using base No array bpPos
      try {
        pos=bpPos[start+curIndex+1];
      } catch (ArrayIndexOutOfBoundsException e ) {
String out="pos assign: bpPos.len "+bpPos.length+
" start: "+start+
" nQualComparisons: "+nQualComparisons+
" curIndex: "+curIndex+
" windowLim: "+windowLim+
" startLhsWindow: "+startLhsWindow+
" qual.length: "+qual.length+
" nPeaks: "+nPeaks+
" nData: "+nData+
" len: "+len+
" i: "+i;
        throw new CheckedRuntimeCSAException(out,e);
      }
      avQualityLhs=(avQualityLhs*(float) numLhs)+wQual[curIndex];
      avQualityRhs=(avQualityRhs*(float) numRhs)-wQual[curIndex];
      if( windowLim > 0 ) {
        ++numLhs;
        numLhs=Math.min(numLhs,windowLim);
        if( curIndex-windowLim >= 0 ) {
          avQualityLhs-=wQual[curIndex-windowLim];
        }
        if( curIndex+windowLim > nData-1 ) {
          --numRhs;
        } else {
          avQualityRhs+=wQual[curIndex+windowLim];
        }
        avQualityLhs/= (float) numLhs;
        avQualityRhs/= (float) numRhs;
      } else {
        avQualityLhs/= (float) ++numLhs;
        avQualityRhs/= (float) --numRhs;
      }
      ++curIndex;
      float ratio=avQualityLhs/avQualityRhs;
//if(log.isInfoEnabled()) log.info("Q ratios: Pos: "+pos+" Ratio: "+ratio+" L "+avQualityLhs+" R: "+avQualityRhs+"Q:"+wQual[curIndex-1]);
//  skip if we're before start of ROI (if set)
      if( roiOnly && pos < ROIStartCoord ) {
        ratio=0.0f;
        continue;
      }
      if ( ratio > critRatio ) {
// n.b. don't assign a max if > critLen bases from prev max
        if ( ratio > maxRatio && pos-maxPos < critLen ) {
          critLen=5;
          maxRatio=ratio;
          maxPos=pos;
        }
      }
    }
    if ( maxPos == pos ) { maxPos=0; }  // zero if max at end

    if( maxPos > 0 ) {
if(log.isInfoEnabled()) log.info("INDEL: "+maxPos+" "+maxRatio+" "+firstQualL+" "+firstQualR);
    }
    if( maxPos > 0 ) {
      maxPos=filterIndelPos(maxPos,maxRatio,critRatio,roiOnly);
    }
    return maxPos;
  }
  private int filterIndelPos(int pos, float maxQRatio, float critRatio, boolean roiOnly) {
// n.b. the suspect indel position calculated from the Q ratios
// can be false - so filter out if any of the following
// a) outside of ROI
// b) match of refSeq has stopped > 25% from end (probably UL)
// c) occurs near (<5) bases from a 'poly base'
// d) severe intensity drop off
    if( roiOnly && ROIStartCoord > 0 ) {
      if( pos < ROIStartCoord || pos > ROIEndCoord ) {
if(log.isInfoEnabled()) log.info("INDEL: ROR "+pos);
        return 0;
      }
    }
    int[] lims=getAnalysisLimits();
// only disable pos if we have < 100 bases after pos
    if( maxQRatio < 1.10f*critRatio ) {
      if( lims[1]-pos < 100 ) {
        if( (float) (lims[1]-lims[0])/(float) (refSearchEnd-lims[0]) < 0.75f) {
if(log.isInfoEnabled()) log.info("INDEL: COV SHORT "+lims[1]);
          return 0;
        }
      }
    }
    String repStr[]={"AAAAAAAA","CCCCCCCC","GGGGGGGG","TTTTTTTT"};
    StringBuffer buf=new StringBuffer(refSeq);
    String revRefSeq=new String(buf.reverse());
    for( int i=0; i < 4 ; ++i ) {
      int match=revRefSeq.indexOf(repStr[i],refSeq.length()-pos);
      if( match > -1 && match-(refSeq.length()-pos) < 5 ) {
if(log.isInfoEnabled()) log.info("INDEL: POLY "+repStr[i]+" "+pos);
//      return 0;
      }
    }
// test for a marked intensity drop off
    if( maxQRatio < 1.10f*critRatio ) {
// test we have decent No. to work with & we're past halfway
    if( lims[1]-pos > 60 &&
      (float) (pos-lims[0])/(float) (lims[1]-lims[0]) > 0.5f) {
    int index=1+getAnalysisIndex(pos);
    float total1=0.0f,total2=0.0f;
    int num=10;
    for( int i=index; i < index+num ; ++i ) {
      total1+=curIntensity[i];
    }
    for( int i=nPeaks+1-num; i < nPeaks+1 ; ++i ) {
      total2+=curIntensity[i];
    }
    if( total1/total2 > 2.0f ) {
if(log.isInfoEnabled()) log.info("INDEL: IDROP "+total1+" "+total2);
        return 0;
      }
    }
    }
    return pos;
  }
  private boolean keepSearching(int pos, boolean forwardSearch, int refSearchLimit) {
    boolean resumeSearch=false;
    if ( forwardSearch ) {
      if( pos < refSearchLimit ) { resumeSearch=true; }
    } else {
      if( pos > refSearchLimit ) { resumeSearch=true; }
    }
    return resumeSearch;
  }
  private void fillSingleTraceHoles() {
    int searchType=1;
    int searchFilter = 1;
    int i,j,index,saveIndex,nUnder;
    int iPrevBase=1;
    int scan1,scan2;
    String baseList[]={"A","C","G","T"};
    boolean selectMax=true;   // search for max peak in region
    boolean keep=true;
    int[] vals=null;
    for( i = 2 ; i< nPeaks  ; ++i) {
      if( bpPos[i] == -1 ) { continue; }
      if( iPrevBase == 1 ) { iPrevBase=i; }
      if ( bpPos[i]-bpPos[iPrevBase] == 2 ) {

        int scanDif=curScan[i]-curScan[iPrevBase];
        if( scanDif < Constants.AVERAGE_BASE_SPACING+2 ) continue; // prob Del
        int missingBaseNo=bpPos[i]-1;
        String searchBase = refSeq.substring(missingBaseNo-1,missingBaseNo);
// can't work with any non searchable bases like N, R etc.
        if( nonSearchableBases.indexOf(searchBase) > -1 ) continue;

// test if we're in an Underloaded region n.b. data hasn't been sorted yet so real/noise peaks are mixed
        nUnder=0;
        int found=0;
        for( j = iPrevBase ; j > iPrevBase-10 ; --j ) {
          if( realPeak[j] != AutoCSA.CALLED_PEAK ) { continue; }
          if( curIntensity[j] < Constants.INTENSITY_UNDERLOADED ) nUnder++;
          if( ++found == 2 ) { break; }
        }
        found=0;
        for( j = i ; j < i+10 ; ++j ) {
          if( realPeak[j] != AutoCSA.CALLED_PEAK ) { continue; }
          if( curIntensity[j] < Constants.INTENSITY_UNDERLOADED ) nUnder++;
          if( ++found == 2 ) { break; }
        }
        if( nUnder >= 2 ) {
          if( scanDif > 30 ) continue;
          float mid=(float) (curScan[iPrevBase]+0.5f*(curScan[i]-curScan[iPrevBase]));
          scan1= (int) Math.floor(mid-peakSearchBin);
          scan2= (int) Math.ceil(mid+peakSearchBin);
          scan1= Math.max(scan1,curScan[iPrevBase]+5);
          scan2= Math.min(scan2,curScan[i]-5);
          if( scan1 >= scan2 ) continue;
          vals=searchChannelForPeak(searchBase,scan1,scan2,0.5f);
          if( vals[0] == -1 ) { iPrevBase=i; continue; }  // no peak found
          saveIndex=++reservePeakIndex;
          if( saveIndex > totalReservedPeaks ) { iPrevBase=i; continue; }  // reserved peaks all used
          curScan[saveIndex]=vals[0];
          curIntensity[saveIndex]=vals[1];
        } else {
          scan1=curScan[iPrevBase]+3;
          scan2=curScan[i]-3;
          saveIndex=searchData(searchType,scan1,scan2,searchFilter,searchBase,selectMax);
        }
        keep=true;   // reset to true every time as test below will set to false
        if( saveIndex > 0 ) {
          for( j = 0 ; j< 4 ; ++j ) {  //  cycle bases
            if( searchBase.equals(baseList[j]) ) { continue; }
            index=searchData(searchType,scan1,scan2,searchFilter,baseList[j],selectMax); // turn off result if we find a larger alternative peak (i.e is probably a Hom Sub)
            if( index > 0 && curIntensity[index] > curIntensity[saveIndex] ) { keep=false; }
          }
// n.b. below here nUnder is used to signify if we have added a new peak
          if( keep ) {
             fork[saveIndex] = 0;
             realPeak[saveIndex] = AutoCSA.CALLED_PEAK;
             bpPos[saveIndex] = missingBaseNo;
             if( nUnder >= 2 ) {
               curBase[saveIndex]=searchBase.charAt(0);
               filter[saveIndex]=1;
             }
          } else {
             if( nUnder >= 2 ) {
               curScan[saveIndex]=reservePeakIndex;
               curIntensity[saveIndex]=1;
               --reservePeakIndex;   // peak not needed
             }
          }
        }
      }
      iPrevBase=i;
    }
  }
  public int getTotalAvailablePeaksInRange(int base1, int base2) {
    int index1=getAnalysisIndex(base1);
    if ( index1 == -1 ) return -1;
    int peakNum1=tracePeakNo[index1];
    int index2=getAnalysisIndex(base2);
    if ( index2 == -1 ) return -1;
    int peakNum2=tracePeakNo[index2];
    int numPeaks=0;
// search for peaks over a "peak number" range
    for( int i = 2 ; i < nPeaks ; ++i ) {
      if( tracePeakNo[i] <= peakNum1 ) continue;
      if( tracePeakNo[i] >= peakNum2 ) break;
      if( filter[i] == 1 ) ++numPeaks;
    }
    return numPeaks;
  }
  public int addPeakToAnalysis(String newBase, int newScan, int newIntensity) {
// firstly check if peak has been added before - if so return
    int saveIndex=++reservePeakIndex;
    if( saveIndex > 1 ) {
      if( curScan[saveIndex-1] == newScan ) {
        --reservePeakIndex;
        return -1;
      }
    }
// check if all reserved peaks already used
    if( saveIndex > totalReservedPeaks ) {
      --reservePeakIndex;
      return -1;
    }
    curScan[saveIndex]=newScan;
    curIntensity[saveIndex]=newIntensity;
    curBase[saveIndex]=newBase.charAt(0);
    fork[saveIndex] = 0;
    realPeak[saveIndex] = AutoCSA.NOISE_PEAK;
    bpPos[saveIndex] = -1;
    filter[saveIndex]=1;
    return 1;
  }
/** 
 * <p>Finds the int array of scan indices of Inserted base positions.</p>
 * <p>Note the bases must have been set to "Inserted" type to be recognised.</p>
 * @param bases the String array of Inserted bases
 * @param StartScan the scan from which to start the search
 * @return the int array of Inserted base scan indices
 */
  public int[] getScansOfInsertedBases(String bases, int StartScan) {
    int len=bases.length();
    int[] scans=new int[len];
    int count=0;
    for( int i=0; i < nPeaks+1; ++i ) {
      if( curScan[i] < StartScan ) continue;
      if( realPeak[i] == AutoCSA.ASSIGNED_INSERTION ) {
        scans[count++]=curScan[i];
        if( count >= len ) return scans;
      }
    }
    return scans;
  }
}
