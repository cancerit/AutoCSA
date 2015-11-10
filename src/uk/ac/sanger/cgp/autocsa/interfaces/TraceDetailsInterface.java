package uk.ac.sanger.cgp.autocsa.interfaces;

import java.io.File;
import java.util.ArrayList;
import org.biojava.bio.chromatogram.Chromatogram;

/**
 *@author  $Author$
 * @version $Revision$
 * @author kr2
 */
public interface TraceDetailsInterface {
	Chromatogram getChromatogram();

	/**
	 * Gets the current value of csa_version
	 * 
	 * @return Current value of csa_version
	 */
	String getCsaVersion();

	/**
	 * Gets the current value of dnaType
	 * 
	 * @return Current value of dnaType
	 */
	int getDnaType();

	/**
	 * Gets the current value of forward_or_reverse_strand
	 * 
	 * @return Current value of forward_or_reverse_strand
	 */
	char getForwardOrReverseStrand();

	/**
	 * Gets the current value of maxScanIndex
	 * 
	 * @return Current value of maxScanIndex
	 */
	int getMaxScanIndex();

	/**
	 * Gets the current value of overall_trace_qual
	 * 
	 * @return Current value of overall_trace_qual
	 */
	Float getOverallTraceQuality();

	/**
	 * Gets the current value of roiCoords
	 * 
	 * @return Current value of roiCoords
	 */
	int[] getROICoords();

	/**
	 * Gets the current value of stsSequence
	 * 
	 * @return Current value of stsSequence
	 */
	String getRefSeq();
	
	/**
   * Gets the current value of file_position_from_root
   * @return Current value of file_position_from_root
   */
  File getFilePositionFromRoot();

	/**
	 * Gets the current value of sts_coverage_start
	 * 
	 * @return Current value of sts_coverage_start
	 */
	Integer getStsCoverageStart();

	/**
	 * Gets the current value of sts_coverage_stop
	 * 
	 * @return Current value of sts_coverage_stop
	 */
	Integer getStsCoverageStop();

	/**
	 * Gets the current value of traceBaseParam
	 * 
	 * @return Current value of traceBaseParam
	 */
	ArrayList getTraceBaseParams();

	/**
	 * Gets the current value of traceHoles
	 * 
	 * @return Current value of traceHoles
	 */
	ArrayList getTraceHoles();

	/**
	 * Gets the current value of traceQualityPercentage
	 * 
	 * @return Current value of traceQualityPercentage
	 */
	Float getTraceQualityPercentage();

	/**
	 * Gets the current value of trace_status
	 * 
	 * @return Current value of trace_status
	 */
	int getTraceStatus();

	/**
	 * Gets the current value of normal
	 * 
	 * @return Current value of normal
	 */
	boolean isNormal();

	void setChromatogram(Chromatogram chromatogram);

	/**
	 * Sets the value of csa_version
	 * 
	 * @param csa_version New value for csa_version
	 */
	void setCsaVersion(String csa_version);

	/**
	 * Sets the value of dnaType
	 * 
	 * @param dnaType New value for dnaType
	 */
	void setDnaType(int dnaType);

	/**
	 * Sets the value of forward_or_reverse_strand to 'f'
	 */
	void setForwardStrand();
	
	/**
   * Sets the value of file_position_from_root
   * @param file_position_from_root New value for file_position_from_root
   */
  void setFilePositionFromRoot(File file_position_from_root);

	/**
	 * Sets the value of maxScanIndex
	 * 
	 * @param maxScanIndex New value for maxScanIndex
	 */
	void setMaxScanIndex(int maxScanIndex);

	/**
	 * Sets the value of normal
	 * 
	 * @param normal New value for normal
	 */
	void setNormal(boolean normal);

	/**
	 * Sets the value of overall_trace_qual
	 * 
	 * @param overall_trace_qual New value for overall_trace_qual
	 */
	void setOverallTraceQuality(Float overall_trace_qual);

	/**
	 * Sets the value of roiCoords
	 * 
	 * @param roiCoords New value for roiCoords
	 */
	void setROICoords(int[] roiCoords);

	/**
	 * Sets the value of stsSequence
	 * 
	 * @param stsSequence New value for stsSequence
	 */
	void setRefSeq(String stsSequence);

	/**
	 * Sets the value of forward_or_reverse_strand to 'r'
	 */
	void setReverseStrand();

	/**
	 * Sets the value of seqread_name
	 * 
	 * @param seqread_name New value for seqread_name
	 */
	void setSeqReadName(String seqread_name);

	void setStrand(char strand);

	/**
	 * Sets the value of sts_coverage_start
	 * 
	 * @param sts_coverage_start New value for sts_coverage_start
	 */
	void setStsCoverageStart(Integer sts_coverage_start);

	/**
	 * Sets the value of sts_coverage_stop
	 * 
	 * @param sts_coverage_stop New value for sts_coverage_stop
	 */
	void setStsCoverageStop(Integer sts_coverage_stop);

	/**
	 * Sets the value of traceBaseParam
	 * 
	 * @param traceBaseParam New value for traceBaseParam
	 */
	void setTraceBaseParams(ArrayList traceBaseParam);

	/**
	 * Sets the value of traceHoles
	 * 
	 * @param traceHoles New value for traceHoles
	 */
	void setTraceHoles(ArrayList traceHoles);

	/**
	 * Sets the value of traceQualityPercentage
	 * 
	 * @param traceQualityPercentage New value for traceQualityPercentage
	 */
	void setTraceQualityPercentage(Float traceQualityPercentage);

	/**
	 * Sets the value of trace_status
	 * 
	 * @param trace_status New value for trace_status
	 */
	void setTraceStatus(int trace_status);

	/**
	 * Overides the toString() object method
	 */
	String toString();

	/**
	 * Checks the internal objects chromatogram and file_position_from_root
	 * for null values and returns a boolean indicating if it is possible to
	 * use the chromatogram object or not
	 * 
	 * @return Indicates if you should use the internal Chromatogram object
	 * or you should retrieve the chromatogram from {@link #getFilePositionFromRoot()}
	 */
	boolean useChromatogram();
	
}
