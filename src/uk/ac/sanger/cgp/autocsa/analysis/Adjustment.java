package uk.ac.sanger.cgp.autocsa.analysis ;

import java.util.HashMap;

import uk.ac.sanger.cgp.autocsa.beans.*;
import uk.ac.sanger.cgp.autocsa.exceptions.*;
import uk.ac.sanger.cgp.autocsa.util.*;

/**
 *<p> Mobility property/base spacing class for AutoCSA package.</p>
 *<p> Only valid for ABI 3730 POP7 trace files.</p>
 *<p> Default adjustment for pre mobility corrected trace files.</p>
 *
 *Original author:  emd
 *@author $Author$
 *@version $Revision$
 */
public class Adjustment {

// default adjustment type is NO Mobility correction
  private static int adjustmentType=0;
  public static HashMap mobilities=null;

/**
 * dummy constructor as all methods are static
 */
  public Adjustment () {
  }
/**
 * Sets the mobility adjustment type
 * @param type type code of adjustment
 */
  public static void setAdjustmentType(int type) {
    adjustmentType=type;
  }
/**
 * Gets the type code of the adjustment
 * @return type code of the adjustment
 */
  public static int getAdjustmentType() {
    return adjustmentType;
  }
  private static float getAdjustment (String key ) {

// define parameterised mobility correction coefficients (3730 POP7)
    if( mobilities == null ) {
      mobilities = new HashMap(10);
      mobilities.put("A1",new Float(0.0f));
      mobilities.put("A2",new Float(0.0f));
      mobilities.put("C1",new Float(0.5541f));
      mobilities.put("C2",new Float(-3.406f));
      mobilities.put("G1",new Float(1.546f));
      mobilities.put("G2",new Float(-14.252f));
      mobilities.put("T1",new Float(-1.7502f));
      mobilities.put("T2",new Float(13.354f));
      mobilities.put("S1",new Float(1.5458f));
      mobilities.put("S2",new Float(2.6471f));
    }

    return ((Float) mobilities.get(key)).floatValue();
  }
/**
 * Gets the adjustment scan for mutant base searching
 * @param base1 the current base type
 * @param base2 the search base type
 * @param base1_no the current base No
 * @param base1_scan the current base scan index
 * @return adjusted scan for centre of the search base window
 */
  public static float getAdjustedScan(String base1, String base2, 
                                      int base1_no, int base1_scan) {
    float base1_a,base2_a,base1_b,base2_b,base1_adj,base2_adj;

// if mobility corrected data return uncorrected scan
    if ( adjustmentType != 0 ) {
      return (float) base1_scan;
    }
//  lookup slope adj for 1st base
    base1_a = getAdjustment(base1+1);
//  lookup slope adj for 2nd base
    base2_a = getAdjustment(base2+1);
//  lookup intercept adj for 1st base
    base1_b = getAdjustment(base1+2);
//  lookup intercept adj for 2nd base
    base2_b = getAdjustment(base2+2);
    base1_adj = base1_a * (float) Math.log((double) base1_no) + base1_b;
    base2_adj = base2_a * (float) Math.log((double) base1_no) + base2_b;
// return adjusted scan for additional base at same position
    return (float) base1_scan - base1_adj + base2_adj;
  }
/**
 * Gets the adjustment offset for amplimer base searching
 * @param base1 the current base type
 * @param base2 the next base type
 * @param base1_no the current base No
 * @param base2_no the next base No
 * @return spacing (in scans) for next base search
 */
  public static float getAdjustmentOffset(String base1, String base2, 
                                          int base1_no, int base2_no) {
    float base1_a,base2_a,base1_b,base2_b,base1_adj,base2_adj;
    float s_a,s_b,s_offset;

// if mobility corrected data return scan increment for next base search
    if ( adjustmentType != 0 ) {
      if( base1_no > 100 ) {
        return 10.0f;
      } else {
        return 9.0f;
      }
    }
// lookup slope adj for 1st base
    base1_a = getAdjustment(base1+1);
// lookup slope adj for 2nd base
    base2_a = getAdjustment(base2+1);
// lookup intercept adj for 1st base
    base1_b = getAdjustment(base1+2);
// lookup intercept adj for 2nd base
    base2_b = getAdjustment(base2+2);
    base1_adj = base1_a * (float) Math.log((double) base1_no) + base1_b;
    base2_adj = base2_a * (float) Math.log((double) base2_no) + base2_b;
    s_a = getAdjustment("S"+1);
    s_b = getAdjustment("S"+2);
    s_offset = s_a * (float) Math.log((double) base1_no) + s_b;
// return adjusted scan increment for base at neighbouring position
    return s_offset - base1_adj + base2_adj;
  }

}
