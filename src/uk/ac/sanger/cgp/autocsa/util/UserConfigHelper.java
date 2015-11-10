package uk.ac.sanger.cgp.autocsa.util;

import java.util.Properties;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import uk.ac.sanger.cgp.autocsa.exceptions.CheckedRuntimeCSAException;

/**
 * A class which encapsulates common acitivties against the user configurable parameters
 * in the properties file './resources/csa_analysis.properties'
 *
 * @author Original: kr2
 * @author  $Author$
 * @version $Revision$
 */
public class UserConfigHelper {
  
  protected static final Log log = LogFactory.getLog(UserConfigHelper.class.getName());
  
  //protected static final Log log = LogFactory.getLog(UserConfigHelper.class.getName());
  
  private static final String ANALYSIS_PROP_LOC = "./resources/csa_analysis.properties";
  private static final Properties ANALYSIS_PROPERTIES = PropertyUtils.getProperties(ANALYSIS_PROP_LOC);
  private static float critMutRatio = -1f;
  private static int refSearchEnd = Integer.MIN_VALUE;
  private static int refSearchStart = Integer.MIN_VALUE;
  private static int minPeakSpacing = Integer.MIN_VALUE;
  
  /** Creates a new instance of UserConfigHelper */
  private UserConfigHelper() {
  }
  
  public static float getCritMutRatio() {
    if(critMutRatio < 0.0f) {
      if(log.isInfoEnabled()) log.info("getting critMutRatio from "+ ANALYSIS_PROP_LOC +" file");
      String prop = ANALYSIS_PROPERTIES.getProperty("critMutRatio");
      if(prop.endsWith("%")) {
        prop = prop.substring(0, prop.length()-1);
      }
      int propValue = -1;
      try {
        propValue = Integer.parseInt(prop);
      }
      catch(NumberFormatException e) {
        throw new CheckedRuntimeCSAException("The value for configurable property 'critMutRatio' should be of the form '20%'", e);
      }
      critMutRatio = (float)(100 - propValue) / 100;
    }
    return critMutRatio;
  }
  
  private static int getIntProp(String propName) {
    String prop = ANALYSIS_PROPERTIES.getProperty(propName);
    int propValue = -1;
    try {
      propValue = Integer.parseInt(prop);
    }
    catch(NumberFormatException e) {
      throw new CheckedRuntimeCSAException("The value for configurable property '"+ propName +"' should be of the form '[-]20'", e);
    }
    return propValue;
  }
  
  public static int getRefSearchStart() {
    if(refSearchStart == Integer.MIN_VALUE) {
      if(log.isInfoEnabled()) log.info("getting refSearchStart from "+ ANALYSIS_PROP_LOC +" file");
      refSearchStart = getIntProp("refSearchStart");
    }
    return refSearchStart;
  }
  
  public static int getRefSearchEnd() {
    if(refSearchEnd == Integer.MIN_VALUE) {
      if(log.isInfoEnabled()) log.info("getting refSearchEnd from "+ ANALYSIS_PROP_LOC +" file");
      refSearchEnd = getIntProp("refSearchEnd");
    }
    return refSearchEnd;
  }
  
  public static int getMinPeakSpacing() {
    if(minPeakSpacing == Integer.MIN_VALUE) {
      if(log.isInfoEnabled()) log.info("getting minPeakSpacing from "+ ANALYSIS_PROP_LOC +" file");
      minPeakSpacing = getIntProp("minPeakSpacing");
    }
    return minPeakSpacing;
  }
  
  public static void main(String [] args) {
    System.out.println("getCritMutRatio: "+ getCritMutRatio());
    System.out.println("getRefSearchStart: "+ getRefSearchStart());
    System.out.println("getRefSearchEnd: "+ getRefSearchEnd());
    System.out.println("getMinPeakSpacing: "+ getMinPeakSpacing());
    
    
    System.out.println("\n\n\n\n\n");
    
    
    System.out.println("getCritMutRatio: "+ getCritMutRatio());
    System.out.println("getRefSearchStart: "+ getRefSearchStart());
    System.out.println("getRefSearchEnd: "+ getRefSearchEnd());
    System.out.println("getMinPeakSpacing: "+ getMinPeakSpacing());
  }
  
  
}
