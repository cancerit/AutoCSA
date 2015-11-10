package uk.ac.sanger.cgp.autocsa.beans;

import uk.ac.sanger.cgp.autocsa.interfaces.CsaTraceInterface;
import uk.ac.sanger.cgp.autocsa.util.UserConfigHelper;

public class CSATrace implements java.io.Serializable, CSAExtractableParameters, CsaTraceInterface {
  
  private static final int refSearchEnd = UserConfigHelper.getRefSearchEnd();
  private static final int refSearchStart = UserConfigHelper.getRefSearchStart();
  private static final int minPeakSpacing = UserConfigHelper.getMinPeakSpacing();
  private static final float critMutRatio = UserConfigHelper.getCritMutRatio();
  
  private static final float peakSearchBin = 6.0f;
  private static final int refSearchStartInc = 40;
  private static final int maxBasesMissed = 50;
  private static final int mobilityCorrection = 1;
  

  private TraceDetails traceDetails;
  private int minPeakIntensity;
  private String refSeq;

  /**
   * Creates a new instance of CSATrace
   */
  public CSATrace() {
  }

  /**
   * Sets the value of refSeq
   * @param refSeq New value for refSeq
   */
  public void setRefSeq(String refSeq) {
    this.refSeq=refSeq;
  }

  /**
   * Gets the current value of refSeq
   * @return Current value of refSeq
   */
  public String getRefSeq() {
    return refSeq;
  }
  /**
   * Gets the current value of traceDetails
   * @return Current value of traceDetails
   */
  public TraceDetails getTraceDetails() {
    return traceDetails;
  }

  /**
   * Sets the value of traceDetails
   * @param traceDetails New value for traceDetails
   */
  public void setTraceDetails(TraceDetails traceDetails) {
    this.traceDetails=traceDetails;
  }

  /**
   * Gets the current value of peakSearchBin
   * @return Current value of peakSearchBin
   */
  public float getPeakSearchBin() {
    return peakSearchBin;
  }

  /**
   * Gets the current value of refSearchStartInc
   * @return Current value of refSearchStartInc
   */
  public int getRefSearchStartInc() {
    return refSearchStartInc;
  }

  /**
   * Gets the current value of maxBasesMissed
   * @return Current value of maxBasesMissed
   */
  public int getMaxBasesMissed() {
    return maxBasesMissed;
  }

  /**
   * Gets the current value of mobilityCorrection
   * @return Current value of mobilityCorrection
   */
  public int getMobilityCorrection() {
    return mobilityCorrection;
  }

  /**
   * Gets the current value of refSearchEnd
   * @return Current value of refSearchEnd
   */
  public int getRefSearchEnd() {
    return refSearchEnd;
  }

  /**
   * Gets the current value of refSearchStart
   * @return Current value of refSearchStart
   */
  public int getRefSearchStart() {
    return refSearchStart;
  }

  /**
   * Gets the current value of minPeakSpacing
   * @return Current value of minPeakSpacing
   */
  public int getMinPeakSpacing() {
    return minPeakSpacing;
  }

  /**
   * Gets the current value of minPeakIntensity
   * @return Current value of minPeakIntensity
   */
  public int getMinPeakIntensity() {
    return minPeakIntensity;
  }

  /**
   * Sets the value of minPeakIntensity
   * @param minPeakIntensity New value for minPeakIntensity
   */
  public void setMinPeakIntensity(int minPeakIntensity) {
    this.minPeakIntensity=minPeakIntensity;
  }
	
	public String toString() {
     String nl = System.getProperty("line.separator");
     StringBuffer sb = new StringBuffer();

     String dec = "========";

     sb.append(dec).append(" [S] uk.ac.sanger.cgp.autocsa.beans.CSATrace ").append(dec).append(nl);
     sb.append("traceDetails=").append(traceDetails).append(nl);
     sb.append("peakSearchBin=").append(peakSearchBin).append(nl);
     sb.append("refSearchStartInc=").append(refSearchStartInc).append(nl);
     sb.append("maxBasesMissed=").append(maxBasesMissed).append(nl);
     sb.append("mobilityCorrection=").append(mobilityCorrection).append(nl);
     sb.append("refSearchEnd=").append(refSearchEnd).append(nl);
     sb.append("refSearchStart=").append(refSearchStart).append(nl);
     sb.append("minPeakSpacing=").append(minPeakSpacing).append(nl);
     sb.append("minPeakIntensity=").append(minPeakIntensity).append(nl);
     sb.append("refSeq=").append(refSeq).append(nl);

     sb.append(dec).append(" [E] uk.ac.sanger.cgp.autocsa.beans.CSATrace ").append(dec);

     return sb.toString();
 }
}
