package uk.ac.sanger.cgp.autocsa.util;

import java.util.*;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/**
 *<p> Class for trace channel operations, for example peak finding.</p>
 *
 *Original author:  emd
 *@author $Author$
 *@version $Revision$
 */
public class Channel {
	
	protected static Log log = LogFactory.getLog(Channel.class.getName());

  public String name;
  private int[] peaksPos;
  private int[] peaks;

// set ABI file codes for raw data channels
  public static final int[] CHANNEL_NAME_CODES={10,11,12,13,14};
  public static final int[] CHANNEL_DATA_CODES={5,6,7,8,9};
// set int flags for adjacent peak indicators
  private static final int PEAK_TO_LEFT=-1;
  private static final int PEAK_TO_RIGHT=1;

  private int[] dataPoints; 

/**
 * null constructor
 */
  public Channel() {
  }
/**
 * Sets the data points for the channel
 * @param points array of channel data
 */
  public void setDataPoints(int[] points){
	dataPoints = points;
  }

/**
 * Gets the data points for the channel
 * @return the array of channel intensities
 */
  public int[] getDataPoints(){
	return dataPoints;
  }

/**
 * Gets the length of the channel
 * @return the length of the channel intensities array
 */
  public int getChannelLength(){
       return dataPoints.length;
  }

/**
 * Gets the minimum intensity of the channel
 * @return the minimum of the channel intensities array
 */
  public int getMinIntensity(){
    int min=30000;
    for (int i=0; i< dataPoints.length; ++i ) {
      min= Math.min(min,dataPoints[i]);
    }
    return min;
  }

/**
 * Gets the maximum intensity of the channel
 * @return the maximum of the channel intensities array
 */
  public int getMaxIntensity(){
    int max=-30000;
    for (int i=0; i< dataPoints.length; ++i ) {
      max= Math.max(max,dataPoints[i]);
    }
    return max;
  }
/**
 * Gets the minimum intensity over a specified range of the channel
 * @param index start array index
 * @param index2 end array index
 * @return the minimum of the channel intensities array over the range
 */
  public int getMinIntensity(int index, int index2){
    int min=30000;
    index=Math.max(index,0);
    index2=Math.min(index2,dataPoints.length-1);
    for (int i=index; i <= index2; ++i ) {
      min= Math.min(min,dataPoints[i]);
    }
    return min;
  }

/**
 * Gets the maximum intensity over a specified range of the channel
 * @param index start array index
 * @param index2 end array index
 * @return the maximum of the channel intensities array over the range
 */
  public int getMaxIntensity(int index, int index2){
    int max=-30000;
    index=Math.max(index,0);
    index2=Math.min(index2,dataPoints.length-1);
    for (int i=index; i <= index2; ++i ) {
      max= Math.max(max,dataPoints[i]);
    }
    return max;
  }
/**
 * Gets the array of derived peak intensities
 * @return the array of derived peak intensities
 */
  public int[] getPeaks () {
    return peaks;
  }
/**
 * Gets the array of derived peak positions
 * @return the array of derived peak positions (scan indices)
 */
  public int[] getPeaksPos () {
    return peaksPos;
  }
/**
 * Calculates the set of peak intensities and positions of the channel
 * @param minAmp Minimum allowed amplitude for peaks
 * @param minShoulder Minimum allowed amplitude for shoulder peaks
 * @param minRelaxedShoulder Relaxed Minimum allowed amplitude for shoulder peaks. <br/> N.B. this is only applied after a certain scan index
 * @param minDist Minimum allowed distance between peaks
 */
  public void findPeaks (int minAmp, int minShoulder, int minRelaxedShoulder, int minDist) {

  int curPos,maxPos,maxValue,peakIndex,dlhs1,drhs1,dlhs2,drhs2;
  int n,num,nShift,count,nStart;

  int nPoints=dataPoints.length;
  int descData[] = new int[nPoints];
  
  List dpList = new ArrayList();

// use DataPoint class to pair up data values with array indices
  for(n=0 ; n < nPoints ; ++n ) {
    dpList.add(new DataPoint(dataPoints[n],n));
  }

// form an ordered descending list of data values
  Collections.sort(dpList);
  Collections.reverse(dpList);
  int scanForRelaxedShoulders=Constants.SCAN_SHOULDER_RELAXATION;

  ArrayList index= new ArrayList();
  for(n=0 ; n < nPoints ; ++n ) {
    descData[n]=((DataPoint) dpList.get(n)).getDataPoint();
// only add elements to index where descData < minAmp
    if( descData[n] >= minAmp ) {
      index.add(new Integer(((DataPoint) dpList.get(n)).getIndex()));
    }
  }

// if first element < minAmp then all elements are
  if ( descData[0] < minAmp ) return;

  // Initialise peaks info arrayLists
  ArrayList peakList= new ArrayList();
  ArrayList peakPosList= new ArrayList();
  
  // Repeat the following procedure until index is empty
  while (index.size() > 0) {  
     // save index of current peak wrt original unsorted array
     peakIndex=((Integer) index.get(0)).intValue();
     maxPos = peakIndex;
     maxValue = dataPoints[peakIndex];
     
     // Store peaks info in a separate array
     // However Discard peaks that are not within a maximum
     // also get rid of isolated spikes
     if ( maxPos > 1 && maxPos < nPoints-2 ) {
       dlhs1 = dataPoints[ maxPos] - dataPoints[ maxPos - 2];
       drhs1 = dataPoints[ maxPos] - dataPoints[ maxPos + 2];
       dlhs2 = dataPoints[ maxPos-1] - dataPoints[ maxPos - 2];
       drhs2 = dataPoints[ maxPos+1] - dataPoints[ maxPos + 2];
     } else {
       dlhs1 = -1;
       drhs1 = -1;
       dlhs2 = -1;
       drhs2 = -1;
     }

// A condition GT produces less FALSE POSITIVE bases than just using GE
// Tested with cross_match.
//   if ( dlhs1 > 0 && drhs1 > 0 && dlhs2 > 0 && drhs2 > 0 ) {
     boolean foundPeak=false;
// if we fail to identify a peak with our conventional model of a
// peak i.e. a maxima, look for a shoulder located within a base
// spacing of a previously found peak (i.e. shoulders will always
// be of lower intensity and we work from largest intensity first)
     if ( dlhs1 >= 0 && drhs1 >= 0 ) {
       foundPeak=true;
     } else {
       int df1=0,df2=0,df3=0;
       int incr=0;
       if( isAdjacentPeak(peakPosList,maxPos,PEAK_TO_LEFT) ) {
         if ( maxPos+3 < nPoints ) {
           df1 = dataPoints[ maxPos    ] - dataPoints[ maxPos + 1];
           df2 = dataPoints[ maxPos + 1] - dataPoints[ maxPos + 2];
           df3 = dataPoints[ maxPos + 2] - dataPoints[ maxPos + 3];
         }
         incr=1;
       } else if( isAdjacentPeak(peakPosList,maxPos,PEAK_TO_RIGHT) ) {
         if ( maxPos-3 >= 0 ) {
           df1 = dataPoints[ maxPos    ] - dataPoints[ maxPos - 1];
           df2 = dataPoints[ maxPos - 1] - dataPoints[ maxPos - 2];
           df3 = dataPoints[ maxPos - 2] - dataPoints[ maxPos - 3];
         }
         incr=-1;
       }
// model of inflexion point where gradients all same sign
       if( df1 >= 0 && df2 > 0 && df3 > 0
                    && (float)df2 >= 2.5f*(float)df1 ) { foundPeak=true; }
// model of inflexion where 1 -ve gradient caused by local upslope
// can relax gradient ratio here as we do have a part peak
       if( df1 < 0 && df2 > 0 && df3 > 0
                   && df2 >= Math.abs(df1) ) { foundPeak=true; }
// n.b. need to corect maxPos as actual peak is at adjacent point
       if ( foundPeak && df1 < 0 ) {
         maxPos += incr;
         maxValue = dataPoints[maxPos];
       }
// reject shoulder if it falls below min value allowed
       int minShoulderForScan=minShoulder;
       if( maxPos > scanForRelaxedShoulders ) minShoulderForScan=minRelaxedShoulder;
       if( foundPeak && maxValue < minShoulderForScan ) foundPeak=false;
     }
// if a valid peak add intensity and position to array lists
     if ( foundPeak ) {
       peakPosList.add(new Integer(maxPos));
       peakList.add(new Integer(maxValue));
     
// remove all values at > distance than minDist to the current value
       nShift=0;
       for(n=0; n < index.size() ; ++n ) {
         curPos=((Integer) index.get(n-nShift)).intValue();
         if( Math.abs(maxPos - curPos) <= minDist ) {
           index.remove(n-nShift++);    // list is left shifted on remove
         }
       }
     } else {
       index.remove(0);    // remove current (failed as a peak)
     }
  }
     
// deal with case where no peaks found
  if( peakPosList.size() == 0 ) {
    peaksPos=null;
    peaks=null;
  } else {
    Integer[] peakPos=(Integer[]) peakPosList.toArray(new Integer[1]);
    Integer[] peak=(Integer[]) peakList.toArray(new Integer[1]);

    TreeMap tmap= new TreeMap();

    for(n=0 ; n < peakPos.length ; ++n ) {
      tmap.put(peakPos[n],peak[n]);
    }
    peaksPos=new int[peakPosList.size()];
    peaks=new int[peakPosList.size()];

    Set set=tmap.keySet();
    Iterator iter = set.iterator();
    count=0;
// save peak & pos into arrays
    while(iter.hasNext() ) {
      Integer key =(Integer) iter.next();
      peaksPos[count]=key.intValue(); 
      peaks[count++]=((Integer) tmap.get(key)).intValue();
    }

    tmap.clear();
    index.clear();
    peakPosList.clear();
    peakList.clear();
  }

  }        // end findPeaks

/**
 * Searches for a single peak over a specified range
 * @param searchLimit1 start array index
 * @param searchLimit2 end array index
 * @return a 2 element array of [peakIndex,peakIntensity], if no peak found then returns [-1,-1].
 */
  public int[] searchForPeak(int searchLimit1, int searchLimit2, float intensityFactor) {
      searchLimit1=Math.max(searchLimit1,0);
      searchLimit2=Math.min(searchLimit2,dataPoints.length-1);
      int maxData=getMaxIntensity(searchLimit1,searchLimit2);
      if( maxData < Constants.STRICT_MIN_PEAK * intensityFactor ) {
        return new int[] {-1,-1};
      }
      int indexAtMax=0;
      for( int i=searchLimit1; i <= searchLimit2; ++i ) {
        if( dataPoints[i] == maxData ) { indexAtMax=i; }
      }
// check if we have a peak structure at this index
      int nChanPoints=dataPoints.length;
      int dlhs1,drhs1;
      if ( indexAtMax > 1 && indexAtMax < nChanPoints-2 ) {
        dlhs1 = dataPoints[ indexAtMax] - dataPoints[ indexAtMax - 2];
        drhs1 = dataPoints[ indexAtMax] - dataPoints[ indexAtMax + 2];
      } else {
        dlhs1 = -1;
        drhs1 = -1;
      }
      if ( dlhs1 >= 0 && drhs1 >= 0 ) {
        return new int[] {indexAtMax,maxData};
      }
      return new int[] {-1,-1};
  }
/**
 * Asks if there is a close peak on either the lhs or rhs contained in
the input list, wrt an input peak position.
 * n.b. spacing < Constants.AVERAGE_BASE_SPACING-2
 * @param peakPosList ArrayList of peak positions
 * @param position scan position of peak to query
 * @param direction direction indicator for peak (PEAK_TO_LEFT or PEAK_TO_RIGHT)
 * @return boolean result
 */
  private boolean isAdjacentPeak(ArrayList peakPosList, int position, int direction) {
// n.b. shoulders will be relatively closely spaced so relax spacing
    int maxSpacing=Constants.AVERAGE_BASE_SPACING-2;
// search through list of peaks found so far
    Iterator iter = peakPosList.iterator();
    while(iter.hasNext() ) {
      int pos=((Integer) iter.next()).intValue();
      if ( Math.abs(pos - position) < maxSpacing) {
        if( direction == PEAK_TO_LEFT && pos < position) {return true; }
        if( direction == PEAK_TO_RIGHT && pos > position) {return true; }
      }
    }
    return false;
  }
/**
 * Outputs to stdout a print out of channel intensities over a range
 * @param startIndex start array index
 * @param numPoints No of points to output
 */
  public void printVals(int startIndex, int numPoints){
    for (int i=startIndex; i< startIndex+numPoints; ++i ) {
      if(log.isInfoEnabled()) log.info("Channel: Trace at Index: "+i+" "+dataPoints[i]);
    }
  }

}
