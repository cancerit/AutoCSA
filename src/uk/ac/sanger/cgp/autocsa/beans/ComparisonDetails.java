package uk.ac.sanger.cgp.autocsa.beans;

import java.io.Serializable;

public class ComparisonDetails implements java.io.Serializable {

  private int compCoverageStart = 0;
  private int compCoverageStop = 0;
  private int comparisonScore = 0;
  private String csaVersion = null;
  private Long idComparison = null;
  private int idCompStatus = 0;

  /**
   * Creates a new instance of ComparisonDetails
   */
  public ComparisonDetails() {
  }

  /**
   * Gets the current value of compCoverageStart
   * @return Current value of compCoverageStart
   */
  public int getCompCoverageStart() {
    return compCoverageStart;
  }

  /**
   * Sets the value of compCoverageStart
   * @param compCoverageStart New value for compCoverageStart
   */
  public void setCompCoverageStart(int compCoverageStart) {
    this.compCoverageStart=compCoverageStart;
  }

  /**
   * Gets the current value of compCoverageStop
   * @return Current value of compCoverageStop
   */
  public int getCompCoverageStop() {
    return compCoverageStop;
  }

  /**
   * Sets the value of compCoverageStop
   * @param compCoverageStop New value for compCoverageStop
   */
  public void setCompCoverageStop(int compCoverageStop) {
    this.compCoverageStop=compCoverageStop;
  }

  /**
   * Gets the current value of comparisonScore
   * @return Current value of comparisonScore
   */
  public int getComparisonScore() {
    return comparisonScore;
  }

  /**
   * Sets the value of comparisonScore
   * @param comparisonScore New value for comparisonScore
   */
  public void setComparisonScore(int comparisonScore) {
    this.comparisonScore=comparisonScore;
  }

  /**
   * Gets the current value of csaVersion
   * @return Current value of csaVersion
   */
  public String getCsaVersion() {
    return csaVersion;
  }

  /**
   * Sets the value of csaVersion
   * @param csaVersion New value for csaVersion
   */
  public void setCsaVersion(String csaVersion) {
    this.csaVersion=csaVersion;
  }

  /**
   * Gets the current value of idComparison
   * @return Current value of idComparison
   */
  public Long getIdComparison() {
    return idComparison;
  }

  /**
   * Sets the value of idComparison
   * @param idComparison New value for idComparison
   */
  public void setIdComparison(Long idComparison) {
    this.idComparison=idComparison;
  }

  /**
   * Gets the current value of idCompStatus
   * @return Current value of idCompStatus
   */
  public int getIdCompStatus() {
    return idCompStatus;
  }

  /**
   * Sets the value of idCompStatus
   * @param idCompStatus New value for idCompStatus
   */
  public void setIdCompStatus(int idCompStatus) {
    this.idCompStatus=idCompStatus;
  }

  /**
   * Overides the toString() object method
   */
  public String toString() {
    
    StringBuffer sb = new StringBuffer();
    
    sb.append(this.getClass().getName()+":\n");
    
    sb.append("compCoverageStart: "+compCoverageStart+"\n");
    
    sb.append("compCoverageStop: "+compCoverageStop+"\n");
    
    sb.append("comparisonScore: "+comparisonScore+"\n");
    
    if(csaVersion==null) {
      sb.append("csaVersion: null\n");
    }
    else {
      sb.append("csaVersion: "+csaVersion+"\n");
    }
    
    if(idComparison==null) {
      sb.append("idComparison: null\n");
    }
    else {
      sb.append("idComparison: "+idComparison+"\n");
    }
    
    sb.append("idCompStatus: "+idCompStatus+"\n");
    
    return sb.toString();
  }

}
