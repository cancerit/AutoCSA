package uk.ac.sanger.cgp.autocsa.util;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/*
 * DataPoint.java
 *
 * Created on 15 April 2003, 18:43
 */

/** Datapoint object for storing information from traces
 * @author ady
 */
public class DataPoint implements Comparable {
	
	protected static Log log = LogFactory.getLog(DataPoint.class.getName());
  
  /** Instance of property dataPoint */
  private int dataPoint;
  
  /** Instance of property index */
  private int index;
  
  /** Creates a new instance of DataPoint */
  public DataPoint() {
  }
  
  /** Working constructor for DataPoint
   * @param dataPoint Value of property dataPoint
   * @param index Value of property index
   */
  public DataPoint(int dataPoint, int index) {
    this.dataPoint=dataPoint;
    this.index=index;
  }
  
  /** Over-ridden method of interface Comparable
   * @return Returns an int based upon the results of this method
   * @param obj The object to be compared to
   */  
  public int compareTo(Object obj) {
    //throws ClassCastException
    try {		
      return(compareTo((DataPoint)obj));		
    }		
    
    catch(ClassCastException e){
      if(log.isInfoEnabled()) log.info("Cannot cast Object x to myInteger");

      //throw e;		
    }		
    return(1000);	
  }
  
  /** Over-ridden method of interface Comparable
   * @return Returns an int based upon the results of this method
   * @param dp The <CODE>DataPoint</CODE> to be compared to
   */  
  public int compareTo(DataPoint dp) {
    if (this.dataPoint < dp.dataPoint) return(-1);		
    if (this.dataPoint == dp.dataPoint) return(0);		 
    return(1);
  }
  
  /** Over-ridden method from <CODE>java.lang.Object</CODE>
   * @return String representation of the object
   */  
  public String toString() {
    StringBuffer sb = new StringBuffer();
    sb.append("[dataPoint: ");
    sb.append(dataPoint);
    sb.append(", index: ");
    sb.append(index);
    sb.append("]");
    return sb.toString();
  }
  
  /** Getter for property dataPoint.
   * @return Value of property dataPoint.
   */
  public int getDataPoint() {
    return this.dataPoint;
  }
  
  /** Setter for property dataPoint.
   * @param dataPoint New value of property dataPoint.
   */
  public void setDataPoint(int dataPoint) {
    this.dataPoint = dataPoint;
  }
  
  /** Getter for property index.
   * @return Value of property index.
   */
  public int getIndex() {
    return this.index;
  }
  
  /** Setter for property index.
   * @param index New value of property index.
   */
  public void setIndex(int index) {
    this.index = index;
  }
}
