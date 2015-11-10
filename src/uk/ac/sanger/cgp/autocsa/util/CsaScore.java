package uk.ac.sanger.cgp.autocsa.util;

import java.math.BigDecimal;
import java.util.ArrayList;
import org.apache.commons.math.util.MathUtils;
import uk.ac.sanger.cgp.autocsa.beans.TraceDetails;

/**
 * @author Original: kr2
 * @author  $Author$
 * @version $Revision$
 */
public class CsaScore {
	
	/** Creates a new instance of CsaScore */
	public CsaScore() {
	}
  
  public static float getCoverageScore(TraceDetails traceBean, boolean reverseCoords) {
		return getCoverageScorewithHoles(traceBean, Integer.MIN_VALUE, reverseCoords);
	}
	
	public static float getCoverageScore(TraceDetails traceBean) {
		return getCoverageScorewithHoles(traceBean, traceBean.getTraceHoles().size(), false);
	}
	
	 /**
	 * AutoCSA normal selection method
	 * @param holes number of trace holes in ROR
	 * @param traceBean contains the details of the normal trace to be given a coverage score.
	 * @return Float object of the coverage assessment value.
	 */
  public static float getCoverageScorewithHoles(TraceDetails traceBean, int holes, boolean reverseCoords)
  {
    //  n.b. TraceStatus is NOT set in traceBean arg but SQL query
    //  (findNormTrace) only selects normals with status=1
    float coverageScore = 0.0f;
    float qFac = 0.0f;
    float hFac = 0.0f;
    int CRIT_N_HOLES = 20;
    int min = CRIT_N_HOLES;
    float qual = MathUtils.round(traceBean.getOverallTraceQuality().floatValue(), 2, BigDecimal.ROUND_HALF_UP);
    //int[] ror=getRORCoordsForAnalysis(traceBean);
		int[] ror=traceBean.getROICoords();
		if(ror[0] == 0 && ror[1] == 0) {
			ror[0] = 1;
			ror[1] = traceBean.getRefSeq().length();
		}
    int covStart = traceBean.getStsCoverageStart().intValue();
    int covStop  = traceBean.getStsCoverageStop().intValue();
    int rorStart = ror[0];
    int rorStop = ror[1];
    if(reverseCoords && (traceBean.getForwardOrReverseStrand() == 'r')) {
      int tmpStart = covStart;
      covStart = (traceBean.getRefSeq().length() - covStop) + 1;
      covStop = (traceBean.getRefSeq().length() - tmpStart) + 1;
      //tmpStart = rorStart;
      //rorStart = (traceBean.getRefSeq().length() - rorStop) + 1;
      //rorStop = (traceBean.getRefSeq().length() - tmpStart) + 1;
    }
    int rorLen = rorStop - rorStart + 1;
    int rorCovLen = Math.min(rorStop,covStop)-Math.max(rorStart,covStart)+1;
    float rorCovPc = 100.0f * (float) rorCovLen/ (float) rorLen;

    // get the % of bases with individual q's > qcrit (3)
    qFac = MathUtils.round(traceBean.getTraceQualityPercentage().floatValue(), 2, BigDecimal.ROUND_HALF_UP);
    
    
    if(holes == Integer.MIN_VALUE) {
      holes = 0;
      ArrayList traceHoles = traceBean.getTraceHoles();
      for(int i=0; i<traceHoles.size(); i++) {
        int currHole = ((Integer)traceHoles.get(i)).intValue();
        if(currHole >= rorStart && currHole <=rorStop) {
          holes++;
        }
      }
    }

    if (min > holes)
    {
      min = holes;
    }
    hFac = (float)(CRIT_N_HOLES - min)/(float)CRIT_N_HOLES;

    coverageScore = rorCovPc * qual * qFac * hFac;

    return coverageScore;
  }
	
}
