/*
 * CsaTraceInterface.java
 *
 * Created on October 10, 2006, 10:41 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package uk.ac.sanger.cgp.autocsa.interfaces;

import uk.ac.sanger.cgp.autocsa.beans.*;

/**
 *
 * @author kr2
 */
public interface CsaTraceInterface {
	/**
	 * Gets the current value of maxBasesMissed
	 * 
	 * @return Current value of maxBasesMissed
	 */
	int getMaxBasesMissed();

	/**
	 * Gets the current value of minPeakIntensity
	 * 
	 * @return Current value of minPeakIntensity
	 */
	int getMinPeakIntensity();

	/**
	 * Gets the current value of minPeakSpacing
	 * 
	 * @return Current value of minPeakSpacing
	 */
	int getMinPeakSpacing();

	/**
	 * Gets the current value of mobilityCorrection
	 * 
	 * @return Current value of mobilityCorrection
	 */
	int getMobilityCorrection();

	/**
	 * Gets the current value of peakSearchBin
	 * 
	 * @return Current value of peakSearchBin
	 */
	float getPeakSearchBin();

	/**
	 * Gets the current value of refSearchEnd
	 * 
	 * @return Current value of refSearchEnd
	 */
	int getRefSearchEnd();

	/**
	 * Gets the current value of refSearchStart
	 * 
	 * @return Current value of refSearchStart
	 */
	int getRefSearchStart();

	/**
	 * Gets the current value of refSearchStartInc
	 * 
	 * @return Current value of refSearchStartInc
	 */
	int getRefSearchStartInc();

	/**
	 * Gets the current value of refSeq
	 * 
	 * @return Current value of refSeq
	 */
	String getRefSeq();

	/**
	 * Gets the current value of traceDetails
	 * 
	 * @return Current value of traceDetails
	 */
	TraceDetails getTraceDetails();

	/**
	 * Sets the value of minPeakIntensity
	 * 
	 * @param minPeakIntensity New value for minPeakIntensity
	 */
	void setMinPeakIntensity(int minPeakIntensity);

	/**
	 * Sets the value of traceDetails
	 * 
	 * @param traceDetails New value for traceDetails
	 */
	void setTraceDetails(TraceDetails traceDetails);
  
  /**
	 * Sets the value of traceDetails
	 * 
	 * @return refSeq New value for traceDetails
	 */
	void setRefSeq(String refSeq);
	
	String toString();
	
}
