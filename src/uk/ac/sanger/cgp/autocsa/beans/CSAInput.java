package uk.ac.sanger.cgp.autocsa.beans;
import uk.ac.sanger.cgp.autocsa.beans.*;

import java.io.Serializable;
import uk.ac.sanger.cgp.autocsa.interfaces.CsaInputInterface;
import uk.ac.sanger.cgp.autocsa.util.UserConfigHelper;

public class CSAInput implements java.io.Serializable, CSAExtractableParameters, CsaInputInterface {
  
  private static final float peakSearchBin = 6.0f;
  private static final int refSearchStartInc = 40;
  private static final int maxBasesMissed = 50;
  private static final int mobilityCorrection = 1;
  
  private static final int refSearchEnd = UserConfigHelper.getRefSearchEnd();
  private static final int refSearchStart = UserConfigHelper.getRefSearchStart();
  private static final int minPeakSpacing = UserConfigHelper.getMinPeakSpacing();
  private static final float critMutRatio = UserConfigHelper.getCritMutRatio();
  
  private TraceDetails traceMut;
  private TraceDetails traceWt;
  private String refSeq;
  private boolean isNormalNormalComp;
  
  // not used
  private float minAllowedQuality;
  private int minPeakIntensity;

  /**
   * Creates a new instance of CSAInput
   */
  public CSAInput() {
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
   * Gets the current value of traceMut
   * @return Current value of traceMut
   */
  public TraceDetails getMutTraceDetails() {
    return traceMut;
  }

  /**
   * Sets the value of traceMut
   * @param traceMut New value for traceMut
   */
  public void setMutTraceDetails(TraceDetails traceMut) {
    this.traceMut=traceMut;
  }

  /**
   * Gets the current value of traceWt
   * @return Current value of traceWt
   */
  public TraceDetails getWtTraceDetails() {
    return traceWt;
  }

  /**
   * Sets the value of traceWt
   * @param traceWt New value for traceWt
   */
  public void setWtTraceDetails(TraceDetails traceWt) {
    this.traceWt=traceWt;
  }

  /**
   * Gets the current value of critMutRatio
   * @return Current value of critMutRatio
   */
  public float getCritMutRatio() {
    return critMutRatio;
  }

  /**
   * Gets the current value of minAllowedQuality
   * @return Current value of minAllowedQuality
   */
  public float getMinAllowedQuality() {
    return minAllowedQuality;
  }

  /**
   * Sets the value of minAllowedQuality
   * @param minAllowedQuality New value for minAllowedQuality
   */
  public void setMinAllowedQuality(float minAllowedQuality) {
    this.minAllowedQuality=minAllowedQuality;
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

  /**
   * Gets the current value of isNormalNormalComp
   * @return Current value of isNormalNormalComp
   */
  public boolean getNormalNormalComp() {
    return isNormalNormalComp;
  }
  /**
   * Sets the value of isNormalNormalComp
   * @param isNormalNormalComp New value for isNormalNormalComp
   */
  public void setNormalNormalComp(boolean isNormalNormalComp) {
    this.isNormalNormalComp=isNormalNormalComp;
  }
	
	public String toString() {
     String nl = System.getProperty("line.separator");
     StringBuffer sb = new StringBuffer();

     String dec = "========";

     sb.append(dec).append(" [S] uk.ac.sanger.cgp.autocsa.beans.CSAInput ").append(dec).append(nl);
     sb.append("traceMut=").append(traceMut).append(nl);
     sb.append("traceWt=").append(traceWt).append(nl);
     sb.append("critMutRatio=").append(critMutRatio).append(nl);
     sb.append("minAllowedQuality=").append(minAllowedQuality).append(nl);
     sb.append("peakSearchBin=").append(peakSearchBin).append(nl);
     sb.append("refSearchStartInc=").append(refSearchStartInc).append(nl);
     sb.append("maxBasesMissed=").append(maxBasesMissed).append(nl);
     sb.append("mobilityCorrection=").append(mobilityCorrection).append(nl);
     sb.append("refSearchEnd=").append(refSearchEnd).append(nl);
     sb.append("refSearchStart=").append(refSearchStart).append(nl);
     sb.append("minPeakSpacing=").append(minPeakSpacing).append(nl);
     sb.append("minPeakIntensity=").append(minPeakIntensity).append(nl);
     sb.append("refSeq=").append(refSeq).append(nl);
     sb.append("isNormalNormalComp=").append(isNormalNormalComp).append(nl);

     sb.append(dec).append(" [E] uk.ac.sanger.cgp.autocsa.beans.CSAInput ").append(dec);

     return sb.toString();
 }

}
