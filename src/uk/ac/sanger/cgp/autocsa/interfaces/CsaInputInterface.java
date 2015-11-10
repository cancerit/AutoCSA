/*
 * NewInterface.java
 *
 * Created on October 10, 2006, 10:58 AM
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
public interface CsaInputInterface {
	/**
	 * Gets the current value of critMutRatio
	 * 
	 * @return Current value of critMutRatio
	 */
	float getCritMutRatio();

	/**
	 * Gets the current value of maxBasesMissed
	 * 
	 * @return Current value of maxBasesMissed
	 */
	int getMaxBasesMissed();

	/**
	 * Gets the current value of minAllowedQuality
	 * 
	 * @return Current value of minAllowedQuality
	 */
	float getMinAllowedQuality();

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
	 * Gets the current value of traceMut
	 * 
	 * @return Current value of traceMut
	 */
	TraceDetails getMutTraceDetails();

	/**
	 * Gets the current value of isNormalNormalComp
	 * 
	 * @return Current value of isNormalNormalComp
	 */
	boolean getNormalNormalComp();

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
	 * Gets the current value of traceWt
	 * 
	 * @return Current value of traceWt
	 */
	TraceDetails getWtTraceDetails();

	/**
	 * Sets the value of minAllowedQuality
	 * 
	 * @param minAllowedQuality New value for minAllowedQuality
	 */
	void setMinAllowedQuality(float minAllowedQuality);

	/**
	 * Sets the value of minPeakIntensity
	 * 
	 * @param minPeakIntensity New value for minPeakIntensity
	 */
	void setMinPeakIntensity(int minPeakIntensity);

	/**
	 * Sets the value of traceMut
	 * 
	 * @param traceMut New value for traceMut
	 */
	void setMutTraceDetails(TraceDetails traceMut);

	/**
	 * Sets the value of isNormalNormalComp
	 * 
	 * @param isNormalNormalComp New value for isNormalNormalComp
	 */
	void setNormalNormalComp(boolean isNormalNormalComp);

	/**
	 * Sets the value of refSeq
	 * 
	 * @param refSeq New value for refSeq
	 */
	void setRefSeq(String refSeq);

	/**
	 * Sets the value of traceWt
	 * 
	 * @param traceWt New value for traceWt
	 */
	void setWtTraceDetails(TraceDetails traceWt);

	String toString();
	
}
