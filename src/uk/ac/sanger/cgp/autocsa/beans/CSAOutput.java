package uk.ac.sanger.cgp.autocsa.beans;

import java.util.ArrayList;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import uk.ac.sanger.cgp.autocsa.analysis.*;
import uk.ac.sanger.cgp.autocsa.util.*;

/**
 *<p> Class for holding objects associated with CSA output.</p>
 
 *Original author:  emd
 *@author $Author$
 *@version $Revision$
 */
public class CSAOutput implements java.io.Serializable {
	
	protected static Log log = LogFactory.getLog(CSAOutput.class.getName());

  private SeqTraceComparison comparison;
  private TraceDetails traceMut;
  private TraceDetails traceWt;
  private ComparisonDetails compDetails;
  private ArrayList mutationList=null;
  private ArrayList mutationDetailsList=null;
  private boolean isNormalNormalComp;

/**
 * null constructor
 */
  public CSAOutput () {
  }
/**
 * Gets the current ComparisonDetails object
 * @return the current value of compDetails
 */
  public ComparisonDetails getComparisonDetails() {
    return compDetails;
  }
/**
 * Gets the current mutation list ArrayList
 * @return the current value of mutationList
 */
  public ArrayList getMutationList() {
		if(mutationList == null) {
			mutationList = new ArrayList();
		}
    return mutationList;
  }
/**
 * Gets the current mutation details list ArrayList
 * @return the current value of mutationDetailsList
 */
  public ArrayList getMutationDetailsList() {
    return mutationDetailsList;
  }
/**
 * Gets the current Wild Type TraceDetails bean
 * @return the current value of traceWt
 */
  public TraceDetails getWtTraceDetails() {
    return traceWt;
  }
/**
 * Gets the current Tumour Type TraceDetails bean
 * @return the current value of traceMut
 */
  public TraceDetails getMutTraceDetails() {
    return traceMut;
  }
/**
 * Updates the Tumour Type TraceDetails bean, traceMut.
 * <p> Note, traceMut must have been set before calling this method.
 * @param traceAnal the SeqTraceAnalysis object holding the trace analysis
 */
  public void setCSAMutTraceAnalInfo(SeqTraceAnalysis traceAnal) {
    int[] limits={0, 0};
    float quality=0.0f;
    int statusCode;
    if ( traceAnal != null ) {
      traceAnal.setTraceQuality(Constants.QUALITY_WINDOW_SIZE);
      quality=traceAnal.getAverageTraceQuality();
      limits=traceAnal.getAnalysisLimits();
      statusCode=traceAnal.getTraceStatus();
    } else {
// the Analysis has failed only if traceAnal == null
      statusCode=Constants.ORA_ID_RERUN_MATCH_FAILED;     // code for failure to locate parent seq
    }
    traceMut.setOverallTraceQuality(new Float(quality));
    traceMut.setStsCoverageStart(new Integer(limits[0]));
    traceMut.setStsCoverageStop(new Integer(limits[1]));
    traceMut.setTraceStatus(statusCode);
  }
/**
 * Updates the Wild Type TraceDetails bean, traceWt.
 * <p> Note, traceWt must have been set before calling this method.
 * @param traceAnal the SeqTraceAnalysis object holding the trace analysis
 */
  public void setCSAWtTraceAnalInfo(SeqTraceAnalysis traceAnal) {
    int[] limits={0, 0};
    float quality=0.0f;
    int statusCode;
    if ( traceAnal != null ) {
      traceAnal.setTraceQuality(Constants.QUALITY_WINDOW_SIZE);
      quality=traceAnal.getAverageTraceQuality();
      limits=traceAnal.getAnalysisLimits();
      statusCode=traceAnal.getTraceStatus();
    } else {
// the Analysis has failed only if traceAnal == null
      statusCode=Constants.ORA_ID_RERUN_MATCH_FAILED;     // code for failure to locate parent seq
    }
    traceWt.setOverallTraceQuality(new Float(quality));
    traceWt.setStsCoverageStart(new Integer(limits[0]));
    traceWt.setStsCoverageStop(new Integer(limits[1]));
    traceWt.setTraceStatus(statusCode);
  }
/**
 * Sets the current ComparisonDetails bean, derived from inputs
 * @param comp the input SeqTraceComparison object
 * @param CSAVersion the AutoCSA version
 */
  public void setComparisonInfo(SeqTraceComparison comp, String CSAVersion) {
    comparison=comp;
    int[] limits={0, 0};
    int quality=0;
    int statusCode;
    if ( comp != null ) {
      quality=comp.getComparisonQuality();
      limits=comp.getComparisonLimits();
      statusCode=comp.getStatusCode();
    } else {
// the comparison has failed only if comp == null
      statusCode=Constants.ORA_ID_FAILED_IN_CSA;    // code for csa failure
    }
    compDetails= new ComparisonDetails();
    compDetails.setIdCompStatus(statusCode);
    compDetails.setComparisonScore(quality);
    compDetails.setCompCoverageStart(limits[0]);
    compDetails.setCompCoverageStop(limits[1]);
    compDetails.setCsaVersion(CSAVersion);
  }
/**
 * Sets the current Tumour Type trace hole ArrayList
 * @param holes the Tumour Type trace holes
 */
  public void setMutTraceHoles(ArrayList holes) {
    traceMut.setTraceHoles(holes);
  }
/**
 * Sets the current Wild Type trace hole ArrayList
 * @param holes the Wild Type trace holes
 */
  public void setWtTraceHoles(ArrayList holes) {
    traceWt.setTraceHoles(holes);
  }
/**
 * Sets the current Tumour Type TraceDetails bean
 * @param traceDet the Tumour Type TraceDetails bean
 */
  public void setMutTraceDetails(TraceDetails traceDet) {
    traceMut=traceDet;
  }
/**
 * Sets the current Wild Type TraceDetails bean
 * @param traceDet the Wild Type TraceDetails bean
 */
  public void setWtTraceDetails(TraceDetails traceDet) {
    traceWt=traceDet;
  }
/**
 * Adds a Mutation and MutationDetails bean pair to the mutationList and
 * mutationDetailsList ArrayLists.
 * <p> Note, the ArrayLists are created if they don't already exist.
 * @param mut the Mutation ArrayList
 * @param mutDetails the MutationDetails ArrayList
 */
  public void addMutation(Mutation mut, MutationDetails mutDetails) {
    if( mutationList == null || mutationDetailsList == null ) {
      mutationList=new ArrayList();
      mutationDetailsList=new ArrayList();
    }
    mutationList.add(mut);
    mutationDetailsList.add(mutDetails);
  }
/**
 * Modifies the mutation type code for a specified mutation in the ArrayList of mutations
 * The mutation to modify is speciified by position and type code
 * <p> For example inputs of (100,1,10) changes the type code from 1 to
10 for a mutation at position 100
 * @param mutPos the mutation position on the sts
 * @param mutType the code indicating the mutation type to modify
 * @param newMutType the new code to use for the mutation
 */
  public void modifyMutationType(int mutPos, int mutType, int newMutType) {
    for(int i=0 ; i < mutationList.size() ; ++i ) {
      Mutation mutObj= (Mutation) mutationList.get(i);
      if ( mutObj.getStsMutationStart() == mutPos &&
           mutObj.getMutationType() == mutType ) {
        mutObj.setMutationType(newMutType);
        return;
      }
    }
  }
/**
 * Removes a specified number of Mutation and MutationDetails bean pairs from the mutationList
 * and mutationDetailsList ArrayLists.
 * <p> For example inputs of (3,2) remove mutations n-2,n-3,n-4.
 * @param numMutations the number of mutations to remove
 * @param offset the list offset at which to start removing
 */
  public int removeMutations(int numMutations, int offset) {
    int total=mutationList.size();
    if( total < numMutations+offset ) {
      return -1;
    }
    int count=0;
    for(int i=total-1-offset ; i>total-1-offset-numMutations ; --i ) {
      mutationList.remove(i);
      mutationDetailsList.remove(i);
      ++count;
    }
    return count;
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
  /**
   * Writes to STDOUT basic information about mutations stored in the
   * arraylist mutationList held in the class. The objects in the 
   * arraylist are type {@link Mutation}. (Mainly used for debugging)
   * @param mutType - DB code for mutation type to write information for
   */
  public void writeMutations(int mutType) {
    if(log.isInfoEnabled()) log.info("Type: "+mutType+" Mutations;");
    for(int i=0 ; i < mutationList.size() ; ++i ) {
      Mutation mutObj= (Mutation) mutationList.get(i);
      int mType=mutObj.getMutationType();
      if( mType == mutType ) {
        int mPos=mutObj.getStsMutationStart();
        int mZyg=mutObj.getZygosity();
        String mTT=mutObj.getMutationSeq();
        String mWT=mutObj.getWildTypeSeq();
        if(log.isInfoEnabled()) log.info("Pos: "+mPos+" Zyg: "+mZyg+" WT: "+mWT+" TT: "+mTT);
      }
    }
  }
	
	public String toString() {
     String nl = System.getProperty("line.separator");
     StringBuffer sb = new StringBuffer();

     String dec = "========";

     sb.append(dec).append(" [S] uk.ac.sanger.cgp.autocsa.beans.CSAOutput ").append(dec).append(nl);
     sb.append("comparison=").append(comparison).append(nl);
     sb.append("traceMut=").append(traceMut).append(nl);
     sb.append("traceWt=").append(traceWt).append(nl);
     sb.append("compDetails=").append(compDetails).append(nl);
     sb.append("mutationList=").append(mutationList).append(nl);
     sb.append("mutationDetailsList=").append(mutationDetailsList).append(nl);
     sb.append("isNormalNormalComp=").append(isNormalNormalComp).append(nl);

     sb.append(dec).append(" [E] uk.ac.sanger.cgp.autocsa.beans.CSAOutput ").append(dec);

     return sb.toString();
 }

}
