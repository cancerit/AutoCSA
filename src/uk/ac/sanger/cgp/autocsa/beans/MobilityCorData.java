package uk.ac.sanger.cgp.autocsa.beans;

public class MobilityCorData implements java.io.Serializable {

  private int doCorrection = 0;
  private int numShift = 0;
  private int[] scan = null;
  private int[] shift = null;

  /**
   * Creates a new instance of MobilityCorData
   */
  public MobilityCorData() {
  }

  /**
   * Gets the current value of doCorrection
   * @return Current value of doCorrection
   */
  public int getDoCorrection() {
    return doCorrection;
  }

  /**
   * Sets the value of doCorrection
   * @param doCorrection New value for doCorrection
   */
  public void setDoCorrection(int doCorrection) {
    this.doCorrection=doCorrection;
  }

  /**
   * Gets the current value of numShift
   * @return Current value of numShift
   */
  public int getNumShift() {
    return numShift;
  }

  /**
   * Sets the value of numShift
   * @param numShift New value for numShift
   */
  public void setNumShift(int numShift) {
    this.numShift=numShift;
  }

  /**
   * Gets the current value of scan
   * @return Current value of scan
   */
  public int[] getScan() {
    return scan;
  }

  /**
   * Sets the value of scan
   * @param scan New value for scan
   */
  public void setScan(int[] scan) {
    this.scan=scan;
  }

  /**
   * Gets the current value of shift
   * @return Current value of shift
   */
  public int[] getShift() {
    return shift;
  }

  /**
   * Sets the value of shift
   * @param shift New value for shift
   */
  public void setShift(int[] shift) {
    this.shift=shift;
  }

  /**
   * Overides the toString() object method
   */
  public String toString() {
    
    StringBuffer sb = new StringBuffer();
    
    sb.append(this.getClass().getName()+":\n");
    
    sb.append("doCorrection: "+doCorrection+"\n");
    
    sb.append("numShift: "+numShift+"\n");
    
    if(scan==null) {
      sb.append("scan: null\n");
    }
    else {
      sb.append("scan: "+scan+"\n");
    }
    
    if(shift==null) {
      sb.append("shift: null\n");
    }
    else {
      sb.append("shift: "+shift+"\n");
    }
    
    return sb.toString();
  }

}
