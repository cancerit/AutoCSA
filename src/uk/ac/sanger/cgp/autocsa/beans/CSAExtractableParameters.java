package uk.ac.sanger.cgp.autocsa.beans;

import java.io.Serializable;

import uk.ac.sanger.cgp.autocsa.beans.*;

public interface CSAExtractableParameters extends java.io.Serializable {

  /**
   * Gets the current value of peakSearchBin
   * @return Current value of peakSearchBin
   */
  public float getPeakSearchBin() ;

  /**
   * Gets the current value of refSearchStartInc
   * @return Current value of refSearchStartInc
   */
  public int getRefSearchStartInc() ;

  /**
   * Gets the current value of maxBasesMissed
   * @return Current value of maxBasesMissed
   */
  public int getMaxBasesMissed() ;

  /**
   * Gets the current value of mobilityCorrection
   * @return Current value of mobilityCorrection
   */
  public int getMobilityCorrection() ;

  /**
   * Gets the current value of refSearchEnd
   * @return Current value of refSearchEnd
   */
  public int getRefSearchEnd() ;

  /**
   * Gets the current value of refSearchStart
   * @return Current value of refSearchStart
   */
  public int getRefSearchStart() ;

  /**
   * Gets the current value of minPeakSpacing
   * @return Current value of minPeakSpacing
   */
  public int getMinPeakSpacing() ;

  /**
   * Gets the current value of minPeakIntensity
   * @return Current value of minPeakIntensity
   */
  public int getMinPeakIntensity() ;

}
