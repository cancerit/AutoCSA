package uk.ac.sanger.cgp.autocsa.beans;

/*
 * Integer objects replace primitive type int.  This allows this bean to hold null values for these fields.  The primitive int value will need to be extracted before
 * it can be used in autoCSA - Ed will need to add that to his code as well as checking if the Integer contains null.
 *
 * This is the bean wrapper for Trace.  Originally started by Ed but added to by Ken - 8th Sept 2003.
 *
 */

import java.io.File;
import java.util.ArrayList;
import org.biojava.bio.chromatogram.Chromatogram;
import uk.ac.sanger.cgp.autocsa.interfaces.TraceDetailsInterface;

/**
 *<p> Class for holding trace parameters.</p>
 *
 *Original author:  emd
 *@author $Author$
 *@version $Revision$
 */
public class TraceDetails implements java.io.Serializable, TraceDetailsInterface {

  private static final char FORWARD = 'f';
  private static final char REVERSE = 'r';

    // TRACE table fields that will already be populated on the trace table before autoCSA runs and that will not be changed
  private Long id_trace=null; // this is Long as that's the maximum value being used.
  private Long id_pcrdes=null;
  private char forward_or_reverse_strand='a'; // this will be either f or r and need to convert F or R to lc and throw exception if not.  Doesn't need to be a java object as can't be null.
  private File file_position_from_root=null;
  private String seqread_name=null;

    //  TRACE table fields populated / updated by autoCSA
  private String csa_version = null;
  private Long id_ind = null;
  private Long id_sts = null;
  private Float overall_trace_qual = null;
  private Integer sts_coverage_start = null;
  private Integer sts_coverage_stop = null;
  private ArrayList traceBaseParam = null;
  private ArrayList traceHoles = null;
  private Float traceQualityPercentage = null;
  private int trace_status = 0;

  // these are not on the TRACE table currently (1/04/05)
  private Long id_desexpt = null;
  private String stsSequence = null;
  private int maxScanIndex = 0;
  private boolean normal = false;
  private int[] roiCoords = null;
  private int dnaType = -1;
	private String stsName;

  //New variable for storing a parsed chromatogram
  private Chromatogram chromatogram = null;
	
  /**
   * Creates a new instance of TraceDetails
   */
  public TraceDetails() {
  }

  /**
   * Shortcut constructor, creates a new instance of TraceDetails
   * @param start start base of trace coverage
   * @param stop stop base of trace coverage
   * @param quality overall average trace quality
   */
  public TraceDetails(Integer start, Integer stop, Float quality)
  {
      sts_coverage_start=start;
      sts_coverage_stop=stop;
      overall_trace_qual=quality;
  }

  /**
   * Gets the current value of csa_version
   * @return Current value of csa_version
   */
  public String getCsaVersion() {
    return csa_version;
  }

  /**
   * Sets the value of csa_version
   * @param csa_version New value for csa_version
   */
  public void setCsaVersion(String csa_version) {
    this.csa_version=csa_version;
  }

  /**
   * Gets the current value of file_position_from_root
   * @return Current value of file_position_from_root
   */
  public File getFilePositionFromRoot() {
    return file_position_from_root;
  }

  /**
   * Sets the value of file_position_from_root
   * @param file_position_from_root New value for file_position_from_root
   */
  public void setFilePositionFromRoot(File file_position_from_root) {
    this.file_position_from_root=file_position_from_root;
  }

  /**
   * Gets the current value of forward_or_reverse_strand
   * @return Current value of forward_or_reverse_strand
   */
  public char getForwardOrReverseStrand() {
    return forward_or_reverse_strand;
  }
  
  public void setStrand(char strand) {
		if(strand == FORWARD) {
			setForwardStrand();
		}
		else if(strand == REVERSE) {
			setReverseStrand();
		}
	}

  /**
   * Sets the value of forward_or_reverse_strand to 'f'
   */
  public void setForwardStrand()
  {
      forward_or_reverse_strand = FORWARD;
  }

  /**
   * Sets the value of forward_or_reverse_strand to 'r'
   */
  public void setReverseStrand()
  {
      forward_or_reverse_strand = REVERSE;
  }

  /**
   * Gets the current value of id_desexpt
   * @return Current value of id_desexpt
   */
  public Long getDesExpt() {
    return id_desexpt;
  }

  /**
   * Sets the value of id_desexpt
   * @param id_desexpt New value for id_desexpt
   */
  public void setDesExpt(Long id_desexpt) {
    this.id_desexpt=id_desexpt;
  }

  /**
   * Gets the current value of id_ind
   * @return Current value of id_ind
   */
  public Long getIdInd() {
    return id_ind;
  }

  /**
   * Sets the value of id_ind
   * @param id_ind New value for id_ind
   */
  public void setIdInd(Long id_ind) {
    this.id_ind=id_ind;
  }

  /**
   * Gets the current value of id_pcrdes
   * @return Current value of id_pcrdes
   */
  public Long getIdPcrDes() {
    return id_pcrdes;
  }

  /**
   * Sets the value of id_pcrdes
   * @param id_pcrdes New value for id_pcrdes
   */
  public void setIdPcrDes(Long id_pcrdes) {
    this.id_pcrdes=id_pcrdes;
  }

  /**
   * Gets the current value of id_sts
   * @return Current value of id_sts
   */
  public Long getIdSts() {
    return id_sts;
  }

  /**
   * Sets the value of id_sts
   * @param id_sts New value for id_sts
   */
  public void setIdSts(Long id_sts) {
    this.id_sts=id_sts;
		setStsName(id_sts.toString());
  }

  /**
   * Gets the current value of id_trace
   * @return Current value of id_trace
   */
  public Long getIdTrace() {
    return id_trace;
  }

  /**
   * Sets the value of id_trace
   * @param id_trace New value for id_trace
   */
  public void setIdTrace(Long id_trace) {
    this.id_trace=id_trace;
  }

  /**
   * Gets the current value of maxScanIndex
   * @return Current value of maxScanIndex
   */
  public int getMaxScanIndex() {
    return maxScanIndex;
  }

  /**
   * Sets the value of maxScanIndex
   * @param maxScanIndex New value for maxScanIndex
   */
  public void setMaxScanIndex(int maxScanIndex) {
    this.maxScanIndex=maxScanIndex;
  }

  /**
   * Gets the current value of normal
   * @return Current value of normal
   */
  public boolean isNormal() {
    return normal;
  }

  /**
   * Sets the value of normal
   * @param normal New value for normal
   */
  public void setNormal(boolean normal) {
    this.normal=normal;
  }

  /**
   * Gets the current value of overall_trace_qual
   * @return Current value of overall_trace_qual
   */
  public Float getOverallTraceQuality() {
    return overall_trace_qual;
  }

  /**
   * Sets the value of overall_trace_qual
   * @param overall_trace_qual New value for overall_trace_qual
   */
  public void setOverallTraceQuality(Float overall_trace_qual) {
    this.overall_trace_qual=overall_trace_qual;
  }

  /**
   * Gets the current value of roiCoords
   * @return Current value of roiCoords
   */
  public int[] getROICoords() {
    return roiCoords;
  }

  /**
   * Sets the value of roiCoords
   * @param roiCoords New value for roiCoords
   */
  public void setROICoords(int[] roiCoords) {
    this.roiCoords=roiCoords;
  }

  /**
   * Gets the current value of seqread_name
   * @return Current value of seqread_name
   */
  public String getSeqReadName() {
    return seqread_name;
  }

  /**
   * Sets the value of seqread_name
   * @param seqread_name New value for seqread_name
   */
  public void setSeqReadName(String seqread_name) {
    this.seqread_name=seqread_name;
  }

  /**
   * Gets the current value of stsSequence
   * @return Current value of stsSequence
   */
  public String getRefSeq() {
    return stsSequence;
  }

  /**
   * Sets the value of stsSequence
   * @param stsSequence New value for stsSequence
   */
  public void setRefSeq(String stsSequence) {
    this.stsSequence=stsSequence;
  }

  /**
   * Gets the current value of sts_coverage_start
   * @return Current value of sts_coverage_start
   */
  public Integer getStsCoverageStart() {
    return sts_coverage_start;
  }

  /**
   * Sets the value of sts_coverage_start
   * @param sts_coverage_start New value for sts_coverage_start
   */
  public void setStsCoverageStart(Integer sts_coverage_start) {
    this.sts_coverage_start=sts_coverage_start;
  }

  /**
   * Gets the current value of sts_coverage_stop
   * @return Current value of sts_coverage_stop
   */
  public Integer getStsCoverageStop() {
    return sts_coverage_stop;
  }

  /**
   * Sets the value of sts_coverage_stop
   * @param sts_coverage_stop New value for sts_coverage_stop
   */
  public void setStsCoverageStop(Integer sts_coverage_stop) {
    this.sts_coverage_stop=sts_coverage_stop;
  }

  /**
   * Gets the current value of traceBaseParam
   * @return Current value of traceBaseParam
   */
  public ArrayList getTraceBaseParams() {
    return traceBaseParam;
  }

  /**
   * Sets the value of traceBaseParam
   * @param traceBaseParam New value for traceBaseParam
   */
  public void setTraceBaseParams(ArrayList traceBaseParam) {
    this.traceBaseParam=traceBaseParam;
  }

  /**
   * Gets the current value of traceHoles
   * @return Current value of traceHoles
   */
  public ArrayList getTraceHoles() {
    return traceHoles;
  }

  /**
   * Sets the value of traceHoles
   * @param traceHoles New value for traceHoles
   */
  public void setTraceHoles(ArrayList traceHoles) {
    this.traceHoles=traceHoles;
  }

  /**
   * Gets the current value of traceQualityPercentage
   * @return Current value of traceQualityPercentage
   */
  public Float getTraceQualityPercentage() {
    return traceQualityPercentage;
  }

  /**
   * Sets the value of traceQualityPercentage
   * @param traceQualityPercentage New value for traceQualityPercentage
   */
  public void setTraceQualityPercentage(Float traceQualityPercentage) {
    this.traceQualityPercentage=traceQualityPercentage;
  }

  /**
   * Gets the current value of trace_status
   * @return Current value of trace_status
   */
  public int getTraceStatus() {
    return trace_status;
  }

  /**
   * Sets the value of trace_status
   * @param trace_status New value for trace_status
   */
  public void setTraceStatus(int trace_status) {
    this.trace_status=trace_status;
  }

  /**
   * Gets the current value of dnaType
   * @return Current value of dnaType
   */
  public int getDnaType() {
    return dnaType;
  }

  /**
   * Sets the value of dnaType
   * @param dnaType New value for dnaType
   */
  public void setDnaType(int dnaType) {
    this.dnaType=dnaType;
  }
	
  public Chromatogram getChromatogram() {
		return chromatogram;
	}

	public void setChromatogram(Chromatogram chromatogram) {
		this.chromatogram = chromatogram;
	}

	/**
	 * Checks the internal objects chromatogram and file_position_from_root
	 * for null values and returns a boolean indicating if it is possible to
	 * use the chromatogram object or not
	 *
	 * @return Indicates if you should use the internal Chromatogram object
	 * or you should retrieve the chromatogram from {@link #getFilePositionFromRoot()}
	 */
	public boolean useChromatogram() {
		boolean useChromatogram = false;
		if(getChromatogram() != null) {
			useChromatogram = true;
		}
		return useChromatogram;
	}
	
  /**
   * Overides the toString() object method
   */
  public String toString() {
    
    StringBuffer sb = new StringBuffer();
    
    sb.append(this.getClass().getName()+":\n");
    
    if(csa_version==null) {
      sb.append("csa_version: null\n");
    }
    else {
      sb.append("csa_version: "+csa_version+"\n");
    }
    
    if(file_position_from_root==null) {
      sb.append("file_position_from_root: null\n");
    }
    else {
      sb.append("file_position_from_root: "+file_position_from_root+"\n");
    }
    
    sb.append("forward_or_reverse_strand: "+forward_or_reverse_strand+"\n");
    
    if(id_desexpt==null) {
      sb.append("id_desexpt: null\n");
    }
    else {
      sb.append("id_desexpt: "+id_desexpt+"\n");
    }
    
    if(id_ind==null) {
      sb.append("id_ind: null\n");
    }
    else {
      sb.append("id_ind: "+id_ind+"\n");
    }
    
    if(id_pcrdes==null) {
      sb.append("id_pcrdes: null\n");
    }
    else {
      sb.append("id_pcrdes: "+id_pcrdes+"\n");
    }
    
    if(id_sts==null) {
      sb.append("id_sts: null\n");
    }
    else {
      sb.append("id_sts: "+id_sts+"\n");
    }
    
    if(id_trace==null) {
      sb.append("id_trace: null\n");
    }
    else {
      sb.append("id_trace: "+id_trace+"\n");
    }
    
    sb.append("maxScanIndex: "+maxScanIndex+"\n");
    
    sb.append("normal: "+normal+"\n");
    
    if(overall_trace_qual==null) {
      sb.append("overall_trace_qual: null\n");
    }
    else {
      sb.append("overall_trace_qual: "+overall_trace_qual+"\n");
    }
    
    if(roiCoords==null) {
      sb.append("roiCoords: null\n");
    }
    else {
      sb.append("roiCoords: "+roiCoords[0]+","+roiCoords[1]+"\n");
    }
    
    if(seqread_name==null) {
      sb.append("seqread_name: null\n");
    }
    else {
      sb.append("seqread_name: "+seqread_name+"\n");
    }
    
    if(stsSequence==null) {
      sb.append("stsSequence: null\n");
    }
    else {
      sb.append("stsSequence: "+stsSequence+"\n");
    }
    
    if(sts_coverage_start==null) {
      sb.append("sts_coverage_start: null\n");
    }
    else {
      sb.append("sts_coverage_start: "+sts_coverage_start+"\n");
    }
    
    if(sts_coverage_stop==null) {
      sb.append("sts_coverage_stop: null\n");
    }
    else {
      sb.append("sts_coverage_stop: "+sts_coverage_stop+"\n");
    }
    
    if(traceBaseParam==null) {
      sb.append("traceBaseParam: null\n");
    }
    else {
      sb.append("traceBaseParam: "+traceBaseParam+"\n");
    }
    
    if(traceHoles==null) {
      sb.append("traceHoles: null\n");
    }
    else {
      sb.append("traceHoles: "+traceHoles+"\n");
    }
    
    if(traceQualityPercentage==null) {
      sb.append("traceQualityPercentage: null\n");
    }
    else {
      sb.append("traceQualityPercentage: "+traceQualityPercentage+"\n");
    }
    
    sb.append("trace_status: "+trace_status+"\n");
    
    return sb.toString();
  }

	public String getStsName() {
		return stsName;
	}

	public void setStsName(String stsName) {
		this.stsName = stsName;
	}

}
