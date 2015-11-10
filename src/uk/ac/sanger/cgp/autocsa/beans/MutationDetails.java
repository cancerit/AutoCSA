package uk.ac.sanger.cgp.autocsa.beans;

import java.io.Serializable;

public class MutationDetails implements java.io.Serializable, Comparable {

  private int mut_scan_index_start;
  private int mut_scan_index_stop;
  private int wt_scan_index_start;
  private int wt_scan_index_stop;
  private int mut_confidence;
  private float wt_peak_decrease;
  private boolean unique_peak;
  private boolean has_contextual_change;
  private int overall_mut_confidence;

  public MutationDetails(int mut_start, int mut_stop,int wt_start, int wt_stop, int confidence, float peak_dec, int overall_confidence, boolean unique_pk, boolean has_context) {
    mut_scan_index_start=mut_start;
    mut_scan_index_stop=mut_stop;
    wt_scan_index_start=wt_start;
    wt_scan_index_stop=wt_stop;
    mut_confidence=confidence;
    overall_mut_confidence=overall_confidence;
    wt_peak_decrease=peak_dec;
    unique_peak=unique_pk;
    has_contextual_change=has_context;
  }
  public MutationDetails() {
  }
  public void setMutScanStart(int start) {
    mut_scan_index_start=start;
  }
  public int getMutScanStart() {
    return mut_scan_index_start;
  }
  public void setMutScanStop(int stop) {
    mut_scan_index_stop=stop;
  }
  public int getMutScanStop() {
    return mut_scan_index_stop;
  }
  public void setWtScanStart(int start) {
    wt_scan_index_start=start;
  }
  public int getWtScanStart() {
    return wt_scan_index_start;
  }
  public void setWtScanStop(int stop) {
    wt_scan_index_stop=stop;
  }
  public int getWtScanStop() {
    return wt_scan_index_stop;
  }
  public void setMutConfidence(int confidence) {
    mut_confidence=confidence;
  }
  public int getMutConfidence() {
    return mut_confidence;
  }
  public void setOverallMutConfidence(int confidence) {
    overall_mut_confidence=confidence;
  }
  public int getOverallMutConfidence() {
    return overall_mut_confidence;
  }
  public void setWtPeakDecrease(float peak_dec) {
    wt_peak_decrease=peak_dec;
  }
  public float getWtPeakDecrease() {
    return wt_peak_decrease;
  }
  public void setUniquePeak(boolean unique_pk) {
    unique_peak=unique_pk;
  }
  public boolean getUniquePeak() {
    return unique_peak;
  }
  public void setContextChange(boolean has_context) {
    has_contextual_change=has_context;
  }
  public boolean getContextChange() {
    return has_contextual_change;
  }
	
	public String toString() {
     String nl = System.getProperty("line.separator");
     StringBuffer sb = new StringBuffer();

     String dec = "========";

     sb.append(dec).append(" [S] uk.ac.sanger.cgp.autocsa.beans.MutationDetails ").append(dec).append(nl);
     sb.append("mut_scan_index_start=").append(mut_scan_index_start).append(nl);
     sb.append("mut_scan_index_stop=").append(mut_scan_index_stop).append(nl);
     sb.append("wt_scan_index_start=").append(wt_scan_index_start).append(nl);
     sb.append("wt_scan_index_stop=").append(wt_scan_index_stop).append(nl);
     sb.append("mut_confidence=").append(mut_confidence).append(nl);
     sb.append("wt_peak_decrease=").append(wt_peak_decrease).append(nl);
     sb.append("unique_peak=").append(unique_peak).append(nl);
     sb.append("has_contextual_change=").append(has_contextual_change).append(nl);
     sb.append("overall_mut_confidence=").append(overall_mut_confidence).append(nl);

     sb.append(dec).append(" [E] uk.ac.sanger.cgp.autocsa.beans.MutationDetails ").append(dec);

     return sb.toString();
 }

  public int compareTo(Object o) {
		MutationDetails mut = (MutationDetails) o;
    int result = 0;
		if(this.getMutScanStart() > mut.getMutScanStart()) {
      result = 1;
    }
    else {
      result = -1;
    }
    return result;
	}
}
