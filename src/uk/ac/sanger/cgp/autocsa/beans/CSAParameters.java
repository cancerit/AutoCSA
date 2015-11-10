package uk.ac.sanger.cgp.autocsa.beans;

import java.io.Serializable;

import uk.ac.sanger.cgp.autocsa.beans.*;

public class CSAParameters implements java.io.Serializable {

  private float critMutRatio;
  private float minAllowedQuality;
  private float peakSearchBin;
  private int refSearchStartInc;
  private int maxBasesMissed;
  private int mobilityCorrection;
  private int refSearchEnd;
  private int refSearchStart;
  private int minPeakSpacing;
  private int minPeakIntensity;

  /**
   * Creates a new instance of CSAParameters
   */
  public CSAParameters() {
  }

  /**
   * Gets the current value of critMutRatio
   * @return Current value of critMutRatio
   */
  public float getCritMutRatio() {
    return critMutRatio;
  }

  /**
   * Sets the value of critMutRatio
   * @param critMutRatio New value for critMutRatio
   */
  public void setCritMutRatio(float critMutRatio) {
    this.critMutRatio=critMutRatio;
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
   * Sets the value of peakSearchBin
   * @param peakSearchBin New value for peakSearchBin
   */
  public void setPeakSearchBin(float peakSearchBin) {
    this.peakSearchBin=peakSearchBin;
  }

  /**
   * Gets the current value of refSearchStartInc
   * @return Current value of refSearchStartInc
   */
  public int getRefSearchStartInc() {
    return refSearchStartInc;
  }

  /**
   * Sets the value of refSearchStartInc
   * @param refSearchStartInc New value for refSearchStartInc
   */
  public void setRefSearchStartInc(int refSearchStartInc) {
    this.refSearchStartInc=refSearchStartInc;
  }

  /**
   * Gets the current value of maxBasesMissed
   * @return Current value of maxBasesMissed
   */
  public int getMaxBasesMissed() {
    return maxBasesMissed;
  }

  /**
   * Sets the value of maxBasesMissed
   * @param maxBasesMissed New value for maxBasesMissed
   */
  public void setMaxBasesMissed(int maxBasesMissed) {
    this.maxBasesMissed=maxBasesMissed;
  }

  /**
   * Gets the current value of mobilityCorrection
   * @return Current value of mobilityCorrection
   */
  public int getMobilityCorrection() {
    return mobilityCorrection;
  }

  /**
   * Sets the value of mobilityCorrection
   * @param mobilityCorrection New value for mobilityCorrection
   */
  public void setMobilityCorrection(int mobilityCorrection) {
    this.mobilityCorrection=mobilityCorrection;
  }

  /**
   * Gets the current value of refSearchEnd
   * @return Current value of refSearchEnd
   */
  public int getRefSearchEnd() {
    return refSearchEnd;
  }

  /**
   * Sets the value of refSearchEnd
   * @param refSearchEnd New value for refSearchEnd
   */
  public void setRefSearchEnd(int refSearchEnd) {
    this.refSearchEnd=refSearchEnd;
  }

  /**
   * Gets the current value of refSearchStart
   * @return Current value of refSearchStart
   */
  public int getRefSearchStart() {
    return refSearchStart;
  }

  /**
   * Sets the value of refSearchStart
   * @param refSearchStart New value for refSearchStart
   */
  public void setRefSearchStart(int refSearchStart) {
    this.refSearchStart=refSearchStart;
  }

  /**
   * Gets the current value of minPeakSpacing
   * @return Current value of minPeakSpacing
   */
  public int getMinPeakSpacing() {
    return minPeakSpacing;
  }

  /**
   * Sets the value of minPeakSpacing
   * @param minPeakSpacing New value for minPeakSpacing
   */
  public void setMinPeakSpacing(int minPeakSpacing) {
    this.minPeakSpacing=minPeakSpacing;
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

}
