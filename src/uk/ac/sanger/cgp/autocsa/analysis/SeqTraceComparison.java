package uk.ac.sanger.cgp.autocsa.analysis ;

import java.util.HashMap;
import java.util.Arrays;
import java.util.ArrayList;
import java.io.IOException;
import java.io.FileWriter;
import java.io.File;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import uk.ac.sanger.cgp.autocsa.beans.*;
import uk.ac.sanger.cgp.autocsa.util.*;

/**
 *<p> Workhorse Class for Sequence Trace Comparison operations.</p>
 *
 *Original author:  Ed Dicks
 *@author $Author$
 *@version $Revision$
 */
public class SeqTraceComparison {
	
	protected static Log log = LogFactory.getLog(SeqTraceComparison.class.getName());

  public int nPoints;
  public String stsName;
  public String well;
  public String analysedRefSeq;
  public String refSeq;
  private char[] baseSeq;
  private int[] baseNo;
  private int[] controlScan;
  private int[] controlHt;
  private int[] sampleScan;
  private int[] sampleHt;
  private float[] factor;
  private float[] scale;
  private float critMutRatio = 0.80f;
  private float mutRatioIncrement = 0.06f;
  private float peakSearchBin=6.0f;
  private ArrayList mutations=null;
  private int ROIStartCoord=0;
  private int ROIEndCoord=0;
  private HashMap normAdj=null;
  private CSAOutput output=null;
  private int[][] missingBases=null;
  private boolean[] isOverloaded;
  private boolean mutationsPresent=false;
  private int indelBasePosFromQ=0;

/**
 * Allocates a SeqTraceComparison object.
 * @param np No of data points in comparison, i.e. No of bases
 */
  public SeqTraceComparison(int np) {
       nPoints=np;
  }
/**
 * Associates a CSAOutput mutation output object with this class
 * @param out the CSAOutput object to associate
 */
  public void setCSAOutput(CSAOutput out) {
       output=out;
  }
/**
 * Returns the associated CSAOutput mutation output object for this class
 * @return the associated CSAOutput object
 */
  public CSAOutput getCSAOutput() {
       return output;
  }
/**
 * Returns the comparison status code
 * @return the comparison status code, i.e. mutant or non-mutant. See {@link Constants}.
 */
  public int getStatusCode() {
// n.b. array missingBases must have been set for this method
// to return a meaningfull status code
       if( getNumHetMutations(0) > 0 || missingBases != null || mutationsPresent) {
         return Constants.ORA_ID_DONE_MUTANT;    // code for mutant
       } else {
         return Constants.ORA_ID_DONE_NON_MUTANT;  // code for non-mutant
       }
  }
/**
 * Sets the array of control information for trace holes
 * @param missing int[][] array of trace hole information
 */
  public void setMissingBases(int[][] missing) {
       missingBases=missing;
  }
/**
 * Sets the scan and intensity class arrays for the Control (Normal) trace.
 * @param scan the control scan array to set
 * @param intensity the control intensity array to set
 */
  public void setControlData(int[] scan, int[] intensity){
	controlScan = scan;
	controlHt = intensity;
  }
/**
 * Sets the scan and intensity class arrays for the Sample (Tumour) trace.
 * @param scan the sample scan array to set
 * @param intensity the sample intensity array to set
 */
  public void setSampleData(int[] scan, int[] intensity){
	sampleScan = scan;
	sampleHt = intensity;
  }
/**
 * Initialises all the other data arrays required for the class.
 * @param startNo the start base of the analyses
 * @param missing the ArrayList of missing base indices (held as Integer Objects).
 */
  public void initialiseData(int startNo, ArrayList missing) {
       baseNo= new int[nPoints];
       factor= new float[nPoints];
       scale= new float[nPoints];
       int bnum=startNo;
       int nmiss=missing.size();
       int index=0;
       for( int i = 0 ; i< nPoints ; ++i ) {
         while ( index < nmiss && bnum == ((Integer) missing.get(index)).intValue() ) {
           ++index;
           ++bnum;
         }
         baseNo[i]=bnum++;
         factor[i]=0.0f;
         scale[i]=0.0f;
       }
  }
/**
 * Gets the number of point Heterozygous substitutions occuring after the specified base number.
 * @param fromBaseNo the base number at which to start counting the Hets
 * @return the number of point Hets occuring after a speciifed base
 */
  public int getNumHetMutations(int fromBaseNo) {
      int count=0;
      for(int i=0 ; i<nPoints ; ++i ) {
        if( baseNo[i] < fromBaseNo) {continue; }
        String mut=(String) mutations.get(i);
        if ( mut != null ) { ++count; }
      }
      return count;
  }
/**
 * Sets the char array of the analysed amplimer sequence.
 * @param startNo the start base of the analyses
 * @param missing the ArrayList of missing base indices (held as Integer Objects).
 */
  public void setAnalysedRefSeqData(int startNo, ArrayList missing) {
      baseSeq = new char[nPoints];
      int bnum=startNo;
      int nmiss=missing.size();
      int index=0;
      for( int i = 0 ; i< nPoints ; ++i ) {
        while ( index < nmiss && bnum == ((Integer) missing.get(index)).intValue() ) {
          ++index;
          ++bnum;
        }
        baseSeq[i]=refSeq.charAt(bnum-1);
        ++bnum;
      }
  }
/**
 * Sets the factor array of intensity ratios. (Sample/Control)
 */
  public void setPeakFactors(){
// n.b. data must be ordered wrt bpPos
      for( int i = 0 ; i< nPoints ; ++i ) {
        if ( controlHt[i] > 0 ) {
            factor[i] = (float) sampleHt[i] / (float) controlHt[i];
        } else {
            factor[i] = 0.0f;
            sampleScan[i] = 0;
            sampleHt[i] = 0;
        }
      }
  }

  private void setBaseOverloads() {

    isOverloaded = new boolean[nPoints];
    for( int i = 0 ; i< nPoints ; ++i ) {
      isOverloaded[i]=false;
      if( sampleHt[i] > Constants.INTENSITY_OVERLOADED ||
         controlHt[i] > Constants.INTENSITY_OVERLOADED ) {
           isOverloaded[i]=true; }
    }
  }
/**
 * Performs the peak intensity ratio Normalisation step of CSA.
 * @param cTrace the Control trace SeqTraceAnalysis object
 * @param sTrace the Sample trace SeqTraceAnalysis object
 */
  public void normalisePeaks(SeqTraceAnalysis cTrace, SeqTraceAnalysis sTrace) {
//  Uses non channel specific peak ratios with base correction factors
//  & uses a median of 5 neighbours to calculate scaling for peak drops
      
    float[] temp = new float[5];
    float scaling,adj;
    char base;
    int i,j,x,lim1,lim2;
            
    float a_fac = 0.0f;
    float c_fac = 0.0f;
    float g_fac = 0.0f;
    float t_fac = 0.0f;
    float a_tot = 0.0f;
    float c_tot = 0.0f;
    float g_tot = 0.0f;
    float t_tot = 0.0f;
    int a_no = 0;
    int c_no = 0;
    int g_no = 0;
    int t_no = 0;
    float a_temp = 0.0f;
    float c_temp = 0.0f;
    float g_temp = 0.0f;
    float t_temp = 0.0f;
    
    setBaseOverloads();       // set boolean overloading switches

// derive region over which to calculate base correction factors
    if( ROIStartCoord > 0 ) {
      int[] lims=getROIArrayIndices();
      lim1=lims[0];
      lim2=lims[1]+1;
// n.b. region for calc of base correction factors must be at least 50
      int len=lim2-lim1;
      if ( len < 50 ) {
        int extra=(50-len+1)/2;
        lim1-=extra;
        if ( lim1 < 0 ) {
          extra+=Math.abs(lim1);
          lim1=0;
        }
        lim2+=extra;
        if ( lim2 > nPoints ) {
          lim1=Math.max(0,lim1-(lim2-nPoints));
          lim2=nPoints;
        }
      }
    } else {
      lim1=0;
      lim2=nPoints;
    }
// if we have suspected Het Indel calculate BCF's on whole coverage
    if ( indelBasePosFromQ > 0 ) {
      lim1=0;
      lim2=nPoints;
    }
    for( i = lim1 ; i < lim2 ; ++i ) {
// don't include Overloaded data in calculation of base correction
// factors as intensity data is not quantitative
        if ( isOverloaded[i] ) { continue; }
// don't include dye blob data as not quantitative
        if ( cTrace.isBaseDyeBlob(baseNo[i]) ) { continue; }
        if ( sTrace.isBaseDyeBlob(baseNo[i]) ) { continue; }
        if( baseSeq[i] == 'A' ) {
            a_temp = factor[i];
            a_tot = a_tot + a_temp;
            a_no = a_no + 1;
        } else {
            if( a_temp > 0.0f ) {
                a_tot = a_tot + a_temp;
                a_no = a_no + 1;
            }
        }
        
        if( baseSeq[i] == 'C' ) {
            c_temp = factor[i];
            c_tot = c_tot + c_temp;
            c_no = c_no + 1;
        } else {
            if( c_temp > 0.0f ) {
                c_tot = c_tot + c_temp;
                c_no = c_no + 1;
            }
        }
        
        if( baseSeq[i] == 'G' ) {
            g_temp = factor[i];
            g_tot = g_tot + g_temp;
            g_no = g_no + 1;
        } else {
            if( g_temp > 0.0f ) {
                g_tot = g_tot + g_temp;
                g_no = g_no + 1;
            }
        }
        
        if( baseSeq[i] == 'T' ) {
            t_temp = factor[i];
            t_tot = t_tot + t_temp;
            t_no = t_no + 1;
        } else {
            if( t_temp > 0.0f ) {
                t_tot = t_tot + t_temp;
                t_no = t_no + 1;
            }
        }
    
    }
    
// allow for zero count of bases of a particular type
    if( a_no > 0 ) { a_fac = a_tot / (float) a_no; }
    if( c_no > 0 ) { c_fac = c_tot / (float) c_no; }
    if( g_no > 0 ) { g_fac = g_tot / (float) g_no; }
    if( t_no > 0 ) { t_fac = t_tot / (float) t_no; }
   
    setBaseCorrectionFactors(a_fac,c_fac,g_fac,t_fac);
    
    for( i = 0 ; i< nPoints ; ++i ) {
        base = baseSeq[i];
        adj = ((Float) normAdj.get(String.valueOf(base))).floatValue();
        factor[i] *= adj;
    }

    System.arraycopy(factor,0,temp,0,5);
    scaling = median(temp); //  median of 5
           
    for( i = 0 ; i< nPoints ; ++i ) {
        //  n.b. always have > 10 analysed bases
            if( i < 2 ) {
                lim1 = 0;
                lim2 = 4;
            } else if( i > nPoints - 3 ) {
                lim1 = nPoints - 5;
                lim2 = nPoints - 1;
            } else {
                lim1 = i - 2;
                lim2 = i + 2;
            }
            
        x = 0;
// select only those points (from the 5) for median calculation
// where the effect of the scaling lies within set bounds
// N.B. This filter is now commented as was originally present when
// normalising wrt the previous base - hence is no loner valid
// Although had effect of calculating proper drops in a region of an
// inDel and for traces with varying intensity ratios as it rejects
// local scaling in favour of upstream scaling.
        for( j = lim1 ; j<= lim2 ; ++j ) {
            if ( isOverloaded[j] ) { continue; }
            if ( indelBasePosFromQ > 0 ) {  // required on whole coverage
              if( Math.abs(1.0f - (scaling / factor[j])) > 1.0f-critMutRatio) {
// only use prev scaling after "settled down" i.e. prevents bad scaling
// at start being used throughout the trace - fixed an example Indel call
// where 0 hets called as bad scaling used right from start of trace
                if( i > 9 ) continue;
              }
            }
            temp[x] = factor[j];  // lim1,lim2 ensure max of temp[4]
            ++x;
        }
           
// if a scaling can't be calculated here then previous value used
        if( x > 0 ) {
            float sc[] = new float[x];
            for( j = 0 ; j < x ; ++j ) {
                sc[j] = temp[j];
            }
            scaling = median(sc);   // median average
        }
        
        base = baseSeq[i];
        adj = ((Float) normAdj.get(String.valueOf(base))).floatValue();
        scale[i] = adj * (float) sampleHt[i] / ((float) controlHt[i] * scaling);
        
    }

  }
  private void setBaseCorrectionFactors(float a_fac,float c_fac,float g_fac,float t_fac) {
    float aRatio=1.0f,cRatio=1.0f,gRatio=1.0f,tRatio=1.0f;
// choose most suitable base for the reference scaling, then if any
// of the other bases values are zero these will have a ?Ratio of 1.0
    if( a_fac > 0.0f ) {
      if( c_fac > 0.0f ) { cRatio=a_fac / c_fac; }
      if( g_fac > 0.0f ) { gRatio=a_fac / g_fac; }
      if( t_fac > 0.0f ) { tRatio=a_fac / t_fac; }
    } else if( c_fac > 0.0f ) {
      if( a_fac > 0.0f ) { aRatio=c_fac / a_fac; }
      if( g_fac > 0.0f ) { gRatio=c_fac / g_fac; }
      if( t_fac > 0.0f ) { tRatio=c_fac / t_fac; }
    } else if( g_fac > 0.0f ) {
      if( a_fac > 0.0f ) { aRatio=g_fac / a_fac; }
      if( c_fac > 0.0f ) { cRatio=g_fac / c_fac; }
      if( t_fac > 0.0f ) { tRatio=g_fac / t_fac; }
    } else if( t_fac > 0.0f ) {
      if( a_fac > 0.0f ) { aRatio=t_fac / a_fac; }
      if( c_fac > 0.0f ) { cRatio=t_fac / c_fac; }
      if( g_fac > 0.0f ) { gRatio=t_fac / g_fac; }
    }
    normAdj = new HashMap(4);
    normAdj.put("A",new Float(aRatio));
    normAdj.put("C",new Float(cRatio));
    normAdj.put("G",new Float(gRatio));
    normAdj.put("T",new Float(tRatio));
    float aa=((Float) normAdj.get("A")).floatValue();
    float ac=((Float) normAdj.get("C")).floatValue();
    float ag=((Float) normAdj.get("G")).floatValue();
    float at=((Float) normAdj.get("G")).floatValue();
    if(log.isInfoEnabled()) log.info("Adj: "+aa+" "+ac+" "+ag+" "+at);
  }
  private float median(float values[] ) {
  
    int len = values.length;
    Arrays.sort(values);
    if ( (int) Math.IEEEremainder( (double) len, 2.0) == 0 ) {
      return 0.5f*(values[(len/2)-1]+values[len/2]);
    } else {
      return values[(len-1)/2];
    }
  }

/**
 * Performs a search of the dataset for point Heterozygous substitutions.
 * @param controlTrace the Control trace SeqTraceAnalysis object
 * @param sampleTrace the Sample trace SeqTraceAnalysis object
 */
  public void hetMutationScan(SeqTraceAnalysis controlTrace, SeqTraceAnalysis sampleTrace) {
    int[] mutation = new int[5];
    int[] mutationPrev = new int[5];
    String baseList[]={"A","C","G","T"};
    int i,j,k;
    float bin,b1_a,b2_a,b1_b,b2_b,b1_adj,b2_adj,m_scan;
    int resPeak,resScan,resIntensity;
    int b1_no,b1_scan;
    float searchLimit1,searchLimit2,baseoffset;
    int searchType,searchFilter;
    String mut_base=null,base1,searchBase;
    char mutationBase='N';
    char mutationBasePrev='N';
    int index=0,lastHetPos=0;
    String nullString=null;
    int mutationScore,curMutationScore,curMutationScorePrev;
    boolean mutPeakUnique=false;
    boolean contextChange=false;
    char uniqueFlag=' ';
    char contextFlag=' ';
    int[] ok = new int[2];
    float binLhs,binRhs;
    
    if ( mutations == null ) {
      mutations=new ArrayList(nPoints);
    }
    //  work over extent of analysed sequence
    
    for( i = 0 ; i< nPoints ; ++i ) {
        
        int numPossHets=0;
        mutations.add(i,nullString);
        if( scale[i] > 0 && scale[i] < critMutRatio ) {
            //  reset mut data
            for( k = 0 ; k < mutation.length ; ++k ) {
                mutation[k] = 0;
            }
            mutationBase= 'N';
            curMutationScore=0;
            curMutationScorePrev=0;
            base1 =String.valueOf(baseSeq[i]);
            
            bin =peakSearchBin;
            float[] bins=setHetSearchBins(bin,i);
            for( j = 0 ; j< 4 ; ++j ) {  //  cycle trial bases for het
                mut_base=baseList[j];

                if( base1.equals(mut_base) ) { continue; }
                binLhs=bins[0];
                binRhs=bins[1];
                b1_no = baseNo[i];
                b1_scan = sampleScan[i];
// get centre of window for search for het peak
                m_scan = Adjustment.getAdjustedScan(base1,mut_base,b1_no,b1_scan);
                
                searchType=1;    // > lim1 && < lim2
                searchFilter = 999;
                searchBase = mut_base;
// allows smallest bin=1.0
                int ntries=0;
                while( Math.max(binLhs,binRhs) > 0.0f ) {
                  ++ntries;
                  searchLimit1 = m_scan - binLhs;
                  searchLimit2 = m_scan + binRhs;

                  index=searchData(sampleTrace,searchType,searchLimit1,searchLimit2,searchFilter,searchBase);
                
                  if( index < 0 ) {
                    if( index == -2 ) { //  >1 peak found - reduce bin
                        binLhs = Math.max(0.0f,binLhs - 1.0f);
                        binRhs = Math.max(0.0f,binRhs - 1.0f);
                        continue;
                    } else { //  no peaks found  - move to next
                        //  however if we have decent drop widen window
                        if( ntries == 1 && scale[i] < 0.9f*critMutRatio ) {
                          binLhs += 1.0f;
                          binRhs += 1.0f;
                          continue;
                        }
                        if( sampleHt[i] < Constants.CRIT_UNDERLOAD_PEAK_FOR_HET ) {
                          int[] peakInfo=sampleTrace.searchChannelForPeak(searchBase,Math.round(searchLimit1),Math.round(searchLimit2),0.5f);
                          if( peakInfo[0] > 0 ) {
// n.b. added peak is put into empty data slots in STA object
// none of the het validation methods used here require the noise
// peaks to be in the current order in the STA object (i.e. ordered by scan)
                            int ok1=sampleTrace.addPeakToAnalysis(searchBase,peakInfo[0],peakInfo[1]);
                            if( ok1 == -1 ) break;
                            continue;  // added peak should now be found
                          }
                        }
                        break;
                    }
                  } else { //  1 peak found - log infomation
                    int[] temp = sampleTrace.getIntensity();
                    resIntensity = temp[index];
                    ++mutation[4];  // No. of possible mutants
                 // if( resIntensity > mutation[2] ) {
// check viability of mut peak i.e. is it comparable to locality
                    ok=checkViableMutant(sampleTrace,baseNo[i],index,mut_base,false);
// if mut peak is viable check its validity
                    mutationScore=0;
                    if( ok[0] == 1 ) {
                      mutationScore=validateMutantPeak(controlTrace,sampleTrace,i,index,mut_base);
                    }
// n.b. if poss mutant peak is not in control then mutationScore < 0
                    if( Math.abs(mutationScore) > curMutationScore ) {
                        ++numPossHets;
                        for( k = 0 ; k < mutation.length ; ++k ) {
                            mutationPrev[k] = mutation[k];
                        }
                        curMutationScorePrev=curMutationScore;
                        mutationBasePrev=mutationBase;
                        curMutationScore=Math.abs(mutationScore);
                        temp = sampleTrace.getScan();
                        resScan = temp[index];
//                      temp = sampleTrace.getTracePeakNo();
//                      resPeak = temp[index];
                        mutation[0] = index;
                        mutation[1] = resScan;
                        mutation[2] = resIntensity;
                        mutation[3] = ok[1];
                        mutationBase= mut_base.charAt(0);
                        if( mutationScore < 0 ) {
                          mutPeakUnique=true;
                          uniqueFlag='U';
                        } else {
                          mutPeakUnique=false;
                          uniqueFlag='-';
                        }
//  try to identify a contextual change in neighbourhood of mutation
                        float contextInc=(1.0f-scale[i])/3.0f;
                        if( i < nPoints-1 && Math.abs(1.0f-scale[i+1]) > contextInc ) {
                          contextChange=true;
                          contextFlag='C';
                        } else {
                          contextChange=false;
                          contextFlag='-';
                        }
                    }
                    break;  // break inner while - continue with next base
                  }
                }
            }
            
            if ( mutation[0] > 0 ) {
              int mutDif=Math.abs(mutation[1]-sampleScan[i]);
              sampleTrace.setRealPeak(mutation[0],AutoCSA.HET_SNP_PEAK);
// check local channel wise true quality if borderline het
              if ( (indelBasePosFromQ > 0 && baseNo[i] < indelBasePosFromQ)
                    || indelBasePosFromQ == 0 ) {
              if( scale[i] > 0.9f*critMutRatio && lastHetPos < i-3 ) {
                float trueQual=sampleTrace.getLocalTrueChannelQuality(baseNo[i],String.valueOf(mutationBase),3);
                if( mutation[2] < Constants.STRICT_MIN_PEAK ) trueQual=5.0f; // turn off test
                if( trueQual < 3.0 ) {  // cancel mutation
                  sampleTrace.setRealPeak(mutation[0],AutoCSA.NOISE_PEAK);
                  if( numPossHets == 1 ) continue;
// as > 1 poss hets return to previous poss het
                  mutationBase=mutationBasePrev;
                  curMutationScore=curMutationScorePrev;
                  for( k = 0 ; k < mutation.length ; ++k ) {
                      mutation[k] = mutationPrev[k];
                  }
                  mutDif=Math.abs(mutation[1]-sampleScan[i]);
                  sampleTrace.setRealPeak(mutation[0],AutoCSA.HET_SNP_PEAK);
                  trueQual=sampleTrace.getLocalTrueChannelQuality(baseNo[i],String.valueOf(mutationBase),3);
                  if( trueQual < 3.0 ) {  // cancel mutation
                    sampleTrace.setRealPeak(mutation[0],AutoCSA.NOISE_PEAK);
                    continue;
                  }
                }
              }
              }
              int conf=getHetOverallConfidence(i,curMutationScore,mutation[3],controlTrace,sampleTrace);
              mutations.set(i,mutationBase+"\t"+mutation[1]+"\t"+mutation[2]+"\t"+curMutationScore+"\t"+uniqueFlag+"\t"+contextFlag);
              lastHetPos=i;
//if(log.isInfoEnabled()) log.info("HET MUT at Bp: "+baseNo[i]+" ("+base1+") Change:" +
//                   mutationBase+" Scan: "+mutation[1]+" Ht: "+mutation[2]+
//                   " Total: "+mutation[4]+" Score: "+curMutationScore+
//                   " "+uniqueFlag+" "+contextFlag+" Drop: "+scale[i]+
//                   " MutPeak: "+mutation[3]+" MutDiff: "+mutDif+
//                   " Confidence: "+conf);
            // create a Mutation & MutationDetails object
              if( output != null ) {
                Mutation mutObj= new Mutation(baseNo[i],baseNo[i],Constants.ORA_ID_SUBSTITUTION,Constants.ORA_ID_MANUAL_REVIEW,getZygosityCode(scale[i]),String.valueOf(mutationBase),base1);
                MutationDetails mutDetailsObj=new MutationDetails(mutation[1],0,controlScan[i],0,curMutationScore,scale[i],conf,mutPeakUnique,contextChange);
                output.addMutation(mutObj,mutDetailsObj);
              }
            }
        }
    }
                    
}
  private int validateMutantPeak(SeqTraceAnalysis cTrace, SeqTraceAnalysis sTrace, int compIndex, int sampIndex, String mut_base ) {

    String base1,searchBase;
    float m_scan;
    int b1_no,b1_scan,index,dScan;

    int[] sPeakHt = sTrace.getIntensity();
    int[] cPeakHt = cTrace.getIntensity();
    int[] sScan = sTrace.getScan();
    int[] cScan = cTrace.getScan();
    int mutsampPeak;
    int mutctrlPeak;
    int mutsampScan;
    int mutctrlScan;
    float mutScale;
//  try to locate trial mutant peak in control trace
    int score=0;
    float binFactor=1.5f;
    base1 =String.valueOf(baseSeq[compIndex]);
    b1_no = baseNo[compIndex];
    b1_scan = controlScan[compIndex];
    m_scan = Adjustment.getAdjustedScan(base1,mut_base,b1_no,b1_scan);
//  use increased bin (using binFactor) as important to locate potential
//  peak matching the trial mutant peak, in control sample
    index=locatePeak(cTrace,m_scan,mut_base,binFactor);
    float scan;
    if( index > 0 ) {
      mutctrlPeak= cPeakHt[index];
      mutctrlScan= cScan[index];
      scan=(float) mutctrlScan;
    } else {
// find noise level in equivalent position in control trace
// n.b can't use a window as would encroach into adjacent peaks
      int s1=(int) m_scan;
      int s2=(int) m_scan;
      int chan=cTrace.getTraceObj().getChannelIndex(mut_base);
      mutctrlPeak= cTrace.getTraceObj().getChannel(chan).getMaxIntensity(s1,s2);
// N.B. only accept this value if it is < min peak value used for
// peak finding for the relevant channel as otherwise we've probably
// caught the shoulder of a real called peak
      int minVal=sTrace.getMinPeakIntensity(mut_base);
      if( mutctrlPeak > minVal ) {
        mutctrlPeak=minVal;
      }
      if( mutctrlPeak < 1 ) { mutctrlPeak=1; }  // prevents divide by 0
      mutctrlScan= 0;
      scan=m_scan;
    }
    mutsampPeak= sPeakHt[sampIndex];  // already know index of
    mutsampScan= sScan[sampIndex];    // poss mutant peak (sampIndex)
    float mRatio= (float) mutsampPeak / ((float) sampleHt[compIndex]);
    mRatio/= (float) mutctrlPeak / ((float) controlHt[compIndex]);
//  if(log.isInfoEnabled()) log.info("Mut in CTRL: "+mut_base+" "+mutctrlPeak+" "+mutctrlScan+" Newratio: "+mRatio);
// derive limits around current location
    int lim1,lim2;
            if( compIndex < 2 ) {
                lim1 = 0;
                lim2 = 4;
            } else if( compIndex > nPoints - 3 ) {
                lim1 = nPoints - 5;
                lim2 = nPoints - 1;
            } else {
                lim1 = compIndex - 2;
                lim2 = compIndex + 2;
            }
    float[] temp=new float[5];
    System.arraycopy(factor,lim1,temp,0,5);
    float scaling = median(temp); //  median of 5
    float adj = ((Float) normAdj.get(String.valueOf(mut_base))).floatValue();
    mutScale=adj* (float) mutsampPeak / (scaling* (float) mutctrlPeak);
    float scoringRatio=mutScale;
//if(log.isInfoEnabled()) log.info("Pos: "+b1_no+" Ratio: "+mut_base+" "+scoringRatio+" Mutpos: "+mutctrlScan+" Scale: "+scaling+" Adj: "+adj);
// must have a het score of > 2 i.e. peak increase factor
    score=Math.round(scoringRatio);
    if ( scoringRatio > 2.0f ) {
      if ( mutctrlScan == 0 ) {
        score+=2;    // increment score if no corresponding peak in control
      }
    } else {
      score=0;
//    if ( mutctrlScan == 0 && scoringRatio > 1.5f ) {
//      // relax ratio test if no corresponding peak in control
//      score+=2;   // increment for passing (reduced ratio) + no ctrl
//    }
    }
    score=Math.min(100,score);
    if ( mutctrlScan == 0 ) {
      return -score;
    } else {
      return score;
    }
  }
  private int locatePeak(SeqTraceAnalysis trace,float scan, String base, float binFactor) {
    float searchLimit1,searchLimit2;
    int searchType,searchFilter,index;
    String searchBase;
    float bin=peakSearchBin * binFactor;
    searchFilter = 999;
    searchBase = base;
    searchType=1;    // > lim1 && < lim2
    while( bin > 0.0f ) {  // allows smallest bin=1.0
      searchLimit1 = scan - bin;
      searchLimit2 = scan + bin;

      index=searchData(trace,searchType,searchLimit1,searchLimit2,searchFilter,searchBase);

      if( index < 0 ) {
        if( index == -2 ) { //  >1 peak found - reduce bin
          bin = bin - 1.0f;
          continue;
        } else { //  no peaks found  - mutant peak is valid
          return 0;
//        break;
        }
      } else { //  1 peak found - log infomation
          return index;
//      break;  // break out of inner while
      }
    }
    return -1;
  }
  private int locateAdjacentPeak(SeqTraceAnalysis trace,float targetScan, String searchBase, int inc, int requiredType) {

    int[] scans = trace.getScan();
    int[] realPeak = trace.getRealPeak();
    char[] base = trace.getBase();
    float scan;
    int index=0;
    int lastIndex=0;
    int i;
    for( i=1; i<= trace.nPeaks ; ++i) {
      if ( realPeak[i] != requiredType ) {continue; }
      if ( base[i] != searchBase.charAt(0) ) {continue; }
      scan = (float) scans[i];
      if ( inc == -1 ) {
        if ( scan >= targetScan ) { 
          index=lastIndex;
          break;
        }
      } else {
        if ( scan > targetScan ) { 
          index=i;
          break;
        }
      }
      lastIndex=i;
    }
    if ( i > trace.nPeaks ) { return 0; }
    return index;
  }
  private int searchData(SeqTraceAnalysis trace,int searchType,float searchLimit1,float searchLimit2,int searchFilter,String searchBase) {

    int index=0;
    int count=0;
    float scan=0.0f;
    int[] scans = trace.getScan();
    int[] realPeak = trace.getRealPeak();
    char[] base = trace.getBase();
    for( int i=1; i<= trace.nPeaks ; ++i) {
// search over rest of peaks i.e. not assigned to bases in parent seq
      if ( realPeak[i] != AutoCSA.NOISE_PEAK ) {continue; }
      if ( base[i] != searchBase.charAt(0) ) {continue; }
      scan = (float) scans[i];
      switch(searchType) {
      case 1: {
          if ( scan > searchLimit1 && scan < searchLimit2 ) {
            index = i;
            ++count;
          }
          break;
        }
      case 2: {
          if ( scan > searchLimit1 && scan <= searchLimit2 ) {
            index = i;
            ++count;
          }
          break;
        }
      case 3: {
          if ( scan >= searchLimit1 && scan < searchLimit2 ) {
            index = i;
            ++count;
          }
          break;
        }
      }
    }

    switch(count) {
      case 0: {
          return -1;
        }
      case 1: {
          return index;
        }
      default: {
          return -2;
        }
    }
  }
/**
 * Creates and saves to an output file the comparison data in tabular form.
 * @param outFile the output filename
 */
  public void createTable (String outFile) {

    if(log.isInfoEnabled()) log.info("CSAFile : "+outFile);
    String heading="No\tBase\tScan\tHt\tScan\tHt\tFactor\tScale\tMutant\tScan\tHt";
    File oFile=new File(outFile);
    FileWriter out=null;
    String brow=" ";
    String mut=null;
    try {
      oFile.createNewFile();
      out = new FileWriter(oFile);
      out.write(heading);
      out.write('\n');
      for(int i=0 ; i<nPoints ; ++i ) {

        String row=brow+baseNo[i]+"\t"+baseSeq[i]+"\t"+controlScan[i]+"\t"+controlHt[i]+"\t"+sampleScan[i]+"\t"+sampleHt[i]+"\t"+factor[i]+"\t"+scale[i];
        if ( mutations != null && mutations.size() > 0 ) {
          mut=(String) mutations.get(i);
        }
        if (mut != null ) {
          row=row+"\t"+mut;
        }
        out.write(row);
        out.write('\n');
      }
      out.close();
    } catch ( IOException e ) {
      if(log.isWarnEnabled()) log.warn("Caught IOException: "+e.getMessage());
    }

  }
/**
 * Sets the Region of Interest (ROI) coordinates for this class.
 * @param start the ROI starting coordinate
 * @param end the ROI ending coordinate
 */
  public void setROICoords(int start, int end) {
      ROIStartCoord=start;
      ROIEndCoord=end;
  }
  private int[] getROIArrayIndices() {
    int i,j,index=0;
    int[] indices = new int[] {0,nPoints-1};
    int[] roi = new int[] {ROIStartCoord,ROIEndCoord};
    if ( ROIStartCoord == 0 ) { return indices; }
//  check we have covered at least some of ROI
    int[] lims=getComparisonLimits();
    if ( ROIStartCoord > lims[1] ) { return indices; }
    for( j=0; j < 2; ++j ) {
      for( i=roi[j]; i > 0; --i ) {
        index=getComparisonIndex(i);
        if( index > -1 ) { break; }
      }
      if( index > -1 ) { indices[j]=index; }
      if( index == -1 ) {
        for( i=roi[j]+1; i < 1000; ++i ) {
          index=getComparisonIndex(i);
          if( index > -1 ) { break; }
        }
        if( index > -1 ) { indices[j]=index; }
      }
    }
// return indices of ROI or nearest analysed coordinates
    return indices;
  }
/**
 * Sets the Reference Amplimer Sequence for the current sts.
 * @param refseq the amplimer sequence
 */
  public void setRefSeq(String refseq) {
      refSeq=refseq;
  }
/**
 * Reports coverage information on the current STS wrt STS and ROI
 * @param ROIStartCoord the starting ROI coordinate
 * @param ROIEndCoord the ending ROI coordinate
 * @param refSeqLength the length of the amplimer sequence
 */
  public void reportCoverage(int ROIStartCoord, int ROIEndCoord, int refSeqLength) {
    int bp1=baseNo[0];
    int bp2=baseNo[nPoints-1];
    int lim1=Math.min(bp1,ROIStartCoord);
    int lim2=Math.min(bp2,ROIEndCoord);
    int ourRange=lim2-lim1+1;
    int ROIRange=ROIEndCoord-ROIStartCoord+1;
    float ROICoverage= 100.0f* (float) ourRange/ (float) ROIRange;
    ROICoverage=Math.min(100.0f,ROICoverage);
    if(log.isInfoEnabled()) log.info("ROI Coverage: "+ROICoverage+" %");
    float stsCoverage= 100.0f* (float) (bp2-bp1+1)/ (float) refSeqLength;
    if(log.isInfoEnabled()) log.info(" STS Coverage: "+stsCoverage+" %");
  }
/**
 * Gets a measure of the Comparison Quality.
 * This is calculated as the average over the traces of the (base corrected) intensity ratios of the two traces, expressed in the range [1,100].
 * @return the Comparison Quality in the range [1,100].
 */
  public int getComparisonQuality() {
// returns an integer in range [0,100] indicating extent to which
// average comparison normalisation factor is close to 1.0
    float total=0.0f;
    for(int i=0 ; i<nPoints ; ++i ) {
      total+=factor[i];
    }
    total= total/(float) nPoints;
    if( total > 1.0 ) {
      return Math.round(100.0f/total);
    } else {
      return Math.round(100.0f*total);
    }
  }
/**
 * Gets an array of Comparison limits of the amplimer.
 * @return the int[2] array of comparison limits (base numbers)
 */
  public int[] getComparisonLimits() {
    int[] limits= new int[2];
    limits[0]=baseNo[0];
    limits[1]=baseNo[nPoints-1];
    return limits;
  }
/**
 * Performs a test for any Heterozygous Insertions or Deletions.
 * @return the status code of the test (1: DEL, 2: INS else: NONE)
 */
  public int testHetIndel() {
// if we have > 10% mutations flagged from a "Mutant Sequence" and
// compare this against the ref seq to asses if we have an indel
    int nmut=getNumHetMutations(indelBasePosFromQ);
    int nmutTot=getNumHetMutations(0);
    if( nmut == 0 ) { return -1; }
    int[] lims=getComparisonLimits();
    float mutFreq;
    if( indelBasePosFromQ > 0 ) {
      mutFreq=(float) nmut/ (float) (lims[1]-indelBasePosFromQ);
    } else {
      mutFreq=(float) nmut/ (float) (lims[1]-lims[0]);
    }
    boolean tooFewHets=false;
if(log.isInfoEnabled()) log.info("HET: Mutation Frequency: "+mutFreq+" Nmuts: "+nmut);
//  if( mutFreq < 0.25 ) { tooFewHets=true; }
    if( mutFreq < 0.10 ) { return -1; }
    String mut=null;
    int[] mutPos= new int[Math.max(nmutTot,1)];
    int indexStartIndel;
    int sScan1,cScan1;
    int count=0;
    int i;
    for(i=0 ; i<nPoints ; ++i ) {
        if ( mutations != null && mutations.size() > 0 ) {
          mut=(String) mutations.get(i);
        }
        if (mut != null ) {
          mutPos[count++]=i;    // assumes no gaps in analysis ??
        }
    }
    int numPos=count-1; // save num elements used
    count=0;
// derive bp increments between mutations
    for(i=0 ; i < numPos-1 ; ++i ) {
      if ( baseNo[mutPos[i]] < indelBasePosFromQ ) { continue; }
      if ( (mutPos[i+1]-mutPos[i]) < 3 ) { ++count; }
    }
if(log.isInfoEnabled()) log.info("HET: No of Close Mutations: "+count);
    if( (float) count < (float) nmut/3.0f ) { tooFewHets=true; }
    if( tooFewHets ) {
      if( indelBasePosFromQ > 0 ) {
        indexStartIndel=5+getComparisonIndex(indelBasePosFromQ);
        indexStartIndel=Math.min(indexStartIndel,nPoints-1); // safety
        sScan1=sampleScan[indexStartIndel];
        cScan1=controlScan[indexStartIndel];
        reportSpeculativeIndel(indelBasePosFromQ+5, cScan1, sScan1,1,2);
      }
      return -1;
    }
    for(i=0 ; i<nmut-2 ; ++i ) {
      if ( (mutPos[i+1]-mutPos[i]) < 5 &&
           (mutPos[i+2]-mutPos[i+1]) < 5 ) {
          if( baseNo[mutPos[i]] > indelBasePosFromQ ) {break; }
      }
    }
//  if the base chosen is > 25 bases from indelBasePosFromQ
//  reset the start het position to the first het after indelBasePosFromQ
    if( baseNo[mutPos[i]]-indelBasePosFromQ > 25 ) {
      for(i=0 ; i<nmut-2 ; ++i ) {
          if( baseNo[mutPos[i]] > indelBasePosFromQ ) {break; }
      }
    }
// start Mutant Sequence at position of grouped mutations
    int totalCloseMuts=nmut-i;
    indexStartIndel=mutPos[i];    // array index of start
    sScan1=sampleScan[indexStartIndel];
    cScan1=controlScan[indexStartIndel];
    int startBP=baseNo[indexStartIndel];
if(log.isInfoEnabled()) log.info("HET INDEL: Starting at Base: "+startBP);
    int curPos=startBP;
    int endBP=baseNo[nPoints-1];
    StringBuffer buf= new StringBuffer();
//  n.b. better to restrict mut seq len to something like 50 rather
//  than 100 as if seq has gone out of phase due to indel then ratio
//  of corrections gets skewed so if we have less seq to correct we
//  have more chance of picking up correct corrections ratio
    int maxMutSeqLen=50;   // only approx if have trace holes
    int nLimit=Math.min(nPoints-indexStartIndel,maxMutSeqLen);
    nLimit+=indexStartIndel;
// further limit nLimit in case when apply the mutSeq shift we
// would go past the end of refSeq
    int maxHetDelLength=Constants.MAX_POSSIBLE_HET_DEL;
    int maxHetInsLength=Constants.MAX_POSSIBLE_HET_INS;
    if( baseNo[nLimit-1]+maxHetDelLength > refSeq.length() ) {
      nLimit-=(baseNo[nLimit-1]+maxHetDelLength-refSeq.length());
    }
// test if we have < 25 bases for the mutated seq - if so flag speculative
// however firstly half maxHetDelLength and test again (with 20)
// n.b. maxHetInsLength unchanged as shifts in opposite direction
    if( nLimit-indexStartIndel < 25 ) {
      maxHetDelLength=(int) (0.5f*maxHetDelLength);
      nLimit+=maxHetDelLength;
      nLimit=Math.min(nLimit,nPoints-1);  // saftey check after addition
      if( nLimit-indexStartIndel < 20 ) {
        if( indelBasePosFromQ > 0 ) {
          int nb=nLimit-indexStartIndel;
if(log.isInfoEnabled()) log.info("HET INDEL: Too few bases for MutSeq: "+nb);
          reportSpeculativeIndel(startBP, cScan1, sScan1,1,2);
        }
        return -1;
      }
    }
    int nMiss=0,nTotalMiss=0;
    for(i=indexStartIndel ; i<nLimit ; ++i ) {
        int pos=baseNo[i];
        if( pos != curPos++ ) {
          nMiss=pos-baseNo[i-1]-1;
          nTotalMiss+=nMiss;
          for(int j=0; j< nMiss; ++j) {
            buf.append("N");    // fill in missing bases with N
          }
          curPos+=(nMiss+1);
        }
        if ( mutations != null && mutations.size() > 0 ) {
          mut=(String) mutations.get(i);
        }
        if (mut != null ) {
          buf.append(mut.substring(0,1));  // mut base is 1st char
        } else {
          buf.append(refSeq.substring(pos-1,pos));  // use base from refSeq
        }
    }
    String mutationSeq=buf.toString();
//if(log.isInfoEnabled()) log.info("MUTSEQ: "+mutationSeq);
//  String subRefSeq=refSeq.substring(startBP-1,startBP);
    StringBuffer origBuf= buf;
    float critRatio=0.50f;
    float ratioOfCorrections=compareMutationSeq(buf,startBP,nTotalMiss,critRatio,maxHetDelLength,maxHetInsLength);
    mutationSeq=buf.toString();
//if(log.isInfoEnabled()) log.info("MUTSEQ: "+mutationSeq);
// check if ratioOfCorrections is acceptable
    int score=25;
    int conf=2;
    if( ratioOfCorrections > critRatio ) {
      reportSpeculativeIndel(startBP, cScan1, sScan1,1,2);
      return -1;
    } else if( ratioOfCorrections < 0.1f ) {
      score=50;
      conf=1;
    }
    int matchPos=refSeq.indexOf(mutationSeq,startBP-1);
//  if we find a match then we have a Deletion
    int nbasesIndel;
    String indelSeq;
    int found=0;
    nbasesIndel=matchPos-startBP+1;
    if ( matchPos == -1 ) { nbasesIndel=1; }  //reset so can use if block
    if ( nbasesIndel > 0 ) {
    if ( matchPos > -1 ) {
      indelSeq=refSeq.substring(startBP-1,matchPos);
//    if(log.isInfoEnabled()) log.info("HET DEL: Pos: "+startBP+" Length: "+nbasesIndel+" DelSeq: "+indelSeq+" Score: "+score+" Confidence: "+conf);
      // create a Mutation & MutationDetails object
      if( output != null ) {
        int sScan2=sampleScan[nPoints-1];
        int cIndex=getComparisonIndex(startBP+nbasesIndel-1);
        if( cIndex == -1 ) { cIndex=nPoints-1; }
        int cScan2=controlScan[cIndex];
        Mutation mutObj= new Mutation(startBP,startBP+nbasesIndel-1,Constants.ORA_ID_DELETION,Constants.ORA_ID_MANUAL_REVIEW,Constants.ORA_ID_HETEROZYGOUS,"-",indelSeq);
        MutationDetails mutDetailsObj=new MutationDetails(sScan1,sScan2,cScan1,cScan2,score,0.0f,0,false,false);
        output.addMutation(mutObj,mutDetailsObj);
      }
      mutationsPresent=true;
      found= 1;         // for DEL
    } else {
//  no match found so try to see if we have an Insertion,
//  successively remove bases from start of mutationSeq and try
//  to match against refSeq, if we find a match then the bases just
//  removed (from orig mutant seq) are the Insertion bases
      buf= new StringBuffer();
      int len=mutationSeq.length();
      for(i=0 ; i<len ; ++i ) {
        if( i == len-1 ) { break; }
// n.b. form ins seq with origBuf as this is uncorrected
        buf.append(origBuf.substring(i,i+1));
        mutationSeq=mutationSeq.substring(1);
// if left with mutationSeq.length < 10 matching is getting easier hence
// an incorrect match is getting more likely - so return
        if( mutationSeq.length() < 10 ) { break; }
        matchPos=refSeq.indexOf(mutationSeq,startBP-1);
//      matchPos=refSeq.indexOf(mutationSeq,0);
//if(log.isInfoEnabled()) log.info("HET INS: "+matchPos+" "+mutationSeq);
        if ( matchPos > -1 ) {
          nbasesIndel=i+1;
          indelSeq=buf.toString();
//        if(log.isInfoEnabled()) log.info("HET INS: Pos: "+startBP+" Length: "+nbasesIndel+" InsSeq: "+indelSeq+" Score: "+score+" Confidence: "+conf);
          // create a Mutation & MutationDetails object
          if( output != null ) {
            int sScan2=sampleScan[nPoints-1];
            int cScan2=controlScan[indexStartIndel+1];
            Mutation mutObj= new Mutation(startBP,startBP+1,Constants.ORA_ID_INSERTION,Constants.ORA_ID_MANUAL_REVIEW,Constants.ORA_ID_HETEROZYGOUS,indelSeq,"-");
            MutationDetails mutDetailsObj=new MutationDetails(sScan1,sScan2,cScan1,cScan2,score,0.0f,0,false,false);
            output.addMutation(mutObj,mutDetailsObj);
          }
          mutationsPresent=true;
          found= 2;         // for INS
          break;
        }
      }
    }
    }
    if ( found == 0 ) {
        reportSpeculativeIndel(startBP, cScan1, sScan1,1,2);
    } else {
//  cancel relevant point mutations if INS/DEL found
//  N.B. currently disabled
//    if( output != null ) {
//      int offset=1;   // set so as not to remove last mutation added
//      output.removeMutations(totalCloseMuts,offset);
//    }
    }
    return found;
  }
  private void reportSpeculativeIndel(int pos, int cScan, int sScan, int length, int type) {
    int zygosity,pos1,pos2;
    if( type == 1 ) {
//    if(log.isInfoEnabled()) log.info("HOM INDEL: Speculative call at Pos: "+pos+" Length: "+length);
      zygosity=Constants.ORA_ID_HOMOZYGOUS;
      pos1=pos;
      pos2=pos+length-1;
    } else if( type == 2 ) {
//    if(log.isInfoEnabled()) log.info("HET INDEL: Speculative call at Pos: "+pos+" Length: "+length);
      zygosity=Constants.ORA_ID_HETEROZYGOUS;
      pos1=pos;
      pos2=pos;
    } else {
      if(log.isInfoEnabled()) log.info("Error reportSpeculativeIndel: Unknown type "+type);
      zygosity=Constants.ORA_ID_HOMOZYGOUS;
      pos1=pos;
      pos2=pos+length-1;
    }
    int score=10;   // lower score to 10%
    if( output != null ) {
      // create a Mutation & MutationDetails object
      Mutation mutObj= new Mutation(pos1,pos2,Constants.ORA_ID_SPECULATIVE_INDEL,Constants.ORA_ID_MANUAL_REVIEW,zygosity,"N","N");
      MutationDetails mutDetailsObj=new MutationDetails(sScan,sScan,cScan,cScan,score,0.0f,0,false,false);
      output.addMutation(mutObj,mutDetailsObj);
    }
    mutationsPresent=true;
  }
/**
 * Performs a test for point Homozygous Substitutions.
 * @param cTrace the Control trace SeqTraceAnalysis object
 * @param sTrace the Sample trace SeqTraceAnalysis object
 * @param missingBases the 2-D information array for trace holes
 */
  public void homMutationScan(SeqTraceAnalysis cTrace, SeqTraceAnalysis sTrace, int[][] missingBases) {
    int[] mutation1 = new int[4];
    int[] mutation2 = new int[4];
    int[] peaks = new int[4];
    String baseList[]={"A","C","G","T"};
    int i,j,pos,index,conf;
    int b1_no,b1_scan;
    int compIndex;
    float m_scan;
    char mutation1Base,mutation2Base;
    String hom_base,base1,missBase;
    int homPeak=0;
    int homScan=0;
    int curMutationScore1=0;
    int curMutationScore2=0;

    int[] sPeakHt = sTrace.getIntensity();
    int[] cPeakHt = cTrace.getIntensity();
    int[] sScan = sTrace.getScan();
    int[] cScan = cTrace.getScan();
    float binFactor=1.0f;
    int[] newSampScan = new int[missingBases[0].length];
    int[] newCtrlScan = new int[missingBases[0].length];
    char[] newSeq = new char[missingBases[0].length];
    for( j = 0 ; j < newSampScan.length ; ++j ) {
      newSampScan[j]=0;
      newCtrlScan[j]=0;
      newSeq[j]='N';
    }
    b1_scan=0;     // scan to LHS of pos (ctrl or samp as appopriate)
    for( i=0; i < missingBases[0].length; ++i ) {
      if ( missingBases[1][i] > AutoCSA.MISSING_IN_BOTH ) {    // skip as not snp's
        continue;
      }
      pos=missingBases[0][i];
      missBase=refSeq.substring(pos-1,pos);
      compIndex=getComparisonIndex(pos-1);
      b1_no = pos-1;
// no comparison point to left (n.b will not happen first time thru)
// n.b. always have point to left if we have isolated single misses
      if( compIndex == -1 ) {
        base1=new String(newSeq,i-1,1);  // n.b. i>0 here
        if ( missingBases[1][i] > AutoCSA.MISSING_IN_NORMAL && newSampScan[i-1] == 0 ) {
          base1=String.valueOf('N');   // undef base1 as can't proceed
        }
        if ( missingBases[1][i] != AutoCSA.MISSING_IN_SAMPLE && newCtrlScan[i-1] == 0 ) {
          base1=String.valueOf('N');   // undef base1 as can't proceed
        }
// test for consecutive missed bases (don't get here if isolated points)
        if( missingBases[0][i-1] != b1_no || base1.equals("N") ) {
          reportMissingHomMutation(pos,0,missBase,0,0,cTrace,sTrace);
          continue;
        }
      } else {
        base1 =String.valueOf(baseSeq[compIndex]);   // base to LHS of pos
      }
      //  reset mut data
      for( j = 0 ; j < mutation1.length ; ++j ) {
        mutation1[j] = 0;
        mutation2[j] = 0;
        peaks[j]= 0;
      }
      mutation1Base= 'N';
      mutation2Base= 'N';
      homPeak=0;
      if ( missingBases[1][i] == AutoCSA.MISSING_IN_SAMPLE ) {    // missing in sample
        if( compIndex > -1 ) {
          b1_scan = sampleScan[compIndex];
        } else {
          b1_scan = newSampScan[i-1];
        }
        for( j = 0 ; j< 4 ; ++j ) {  //  cycle bases
          hom_base=baseList[j];

          if( missBase.equals(hom_base) ) { continue; }
          m_scan=getHomSearchCentre(base1,hom_base,pos,b1_scan,sampleScan);
//  use increased bin (using binFactor) as important to locate potential
//  hom peak n.b. NO! uses too big a window -> wrong homs
          index=locatePeak(sTrace,m_scan,hom_base,binFactor);
          if( index > 0 ) {
// check viability of mut peak i.e. is it comparable to locality
            int[] ok=checkHomMutant(sTrace,pos,index,hom_base);
            if( ok[0] == 1 ) {
              homPeak= sPeakHt[index];
              homScan= sScan[index];
              ++mutation1[3]; // No. of possible mutants
              peaks[j]=homPeak;
            } else {
              homPeak=0;      // potential hom peak too small or dye blob
              homScan= 0;
            }
          }
          if( homPeak > mutation1[2] ) {
            mutation1[0] = index;
            mutation1[1] = homScan;
            mutation1[2] = homPeak;
            mutation1Base= hom_base.charAt(0);
          }
        }
        index = cTrace.getAnalysisIndex(pos);
        if ( index > -1 ) {
          newCtrlScan[i]=cScan[index];
        }
        if( mutation1[2] > 0 ) {
          int minAmp=sTrace.getMinPeakIntensity(String.valueOf(mutation1Base));
          curMutationScore1=getHomPeakScore(peaks,minAmp);
          conf=getHomOverallConfidence(pos,curMutationScore1,mutation1[3],sTrace);
//        if(log.isInfoEnabled()) log.info("HOM MUT at Bp: "+pos+" ("+missBase+") Change:"+
//                mutation1Base+" Scan: "+mutation1[1]+" Ht: "+mutation1[2]+
//                " Total: "+mutation1[3]+" Score: "+curMutationScore1+
//                 " Confidence: "+conf);
          newSampScan[i]=mutation1[1];
          newSeq[i]=mutation1Base;
          sTrace.setRealPeak(mutation1[0],AutoCSA.HOM_SNP_PEAK);
          // create a Mutation & MutationDetails object
          if( output != null ) {
            Mutation mutObj= new Mutation(pos,pos,Constants.ORA_ID_SUBSTITUTION,Constants.ORA_ID_MANUAL_REVIEW,Constants.ORA_ID_HOMOZYGOUS,String.valueOf(mutation1Base),missBase);
            MutationDetails mutDetailsObj=new MutationDetails(mutation1[1],0,newCtrlScan[i],0,curMutationScore1,0.0f,conf,getHomUniqueness(mutation1[3]),getHomContext(pos));
            output.addMutation(mutObj,mutDetailsObj);
          }
        } else {
          reportMissingHomMutation(pos,0,missBase,newCtrlScan[i],0,cTrace,sTrace);
        }
      } else if ( missingBases[1][i] == AutoCSA.MISSING_IN_NORMAL ) {    // missing in control
        if( compIndex > -1 ) {
          b1_scan = controlScan[compIndex];
        } else {
          b1_scan = newCtrlScan[i-1];
        }
        for( j = 0 ; j< 4 ; ++j ) {  //  cycle bases
          hom_base=baseList[j];

          if( missBase.equals(hom_base) ) { continue; }
          m_scan=getHomSearchCentre(base1,hom_base,pos,b1_scan,controlScan);
//  use increased bin (using binFactor) as important to locate potential
//  hom peak
          index=locatePeak(cTrace,m_scan,hom_base,binFactor);
          if( index > 0 ) {
// check viability of mut peak i.e. is it comparable to locality
            int[] ok=checkHomMutant(cTrace,pos,index,hom_base);
            if( ok[0] == 1 ) {
              homPeak= cPeakHt[index];
              homScan= cScan[index];
              ++mutation1[3]; // No. of possible mutants
              peaks[j]=homPeak;
            } else {
              homPeak=0;      // potential hom peak too small or dye blob
              homScan= 0;
            }
          }
          if( homPeak > mutation1[2] ) {
            mutation1[0] = index;
            mutation1[1] = homScan;
            mutation1[2] = homPeak;
            mutation1Base= hom_base.charAt(0);
          }
        }
        index = sTrace.getAnalysisIndex(pos);
        if ( index > -1 ) {
          newSampScan[i]=sScan[index];
        }
        if( mutation1[2] > 0 ) {
          int minAmp=cTrace.getMinPeakIntensity(String.valueOf(mutation1Base));
          curMutationScore1=getHomPeakScore(peaks,minAmp);
          conf=getHomOverallConfidence(pos,curMutationScore1,mutation1[3],cTrace);
//        if(log.isInfoEnabled()) log.info("HOM SNP_WT at Bp: "+pos+" ("+missBase+") Change:"+
//                mutation1Base+" Scan: "+mutation1[1]+" Ht: "+mutation1[2]+
//                 " Total: "+mutation1[3]+" Score: "+curMutationScore1+
//                 " Confidence: "+conf);
          newCtrlScan[i]=mutation1[1];
          newSeq[i]=mutation1Base;
          cTrace.setRealPeak(mutation1[0],AutoCSA.HOM_SNP_PEAK);
          // create a Mutation & MutationDetails object
          if( output != null ) {
            Mutation mutObj= new Mutation(pos,pos,Constants.ORA_ID_REFERENCE_SUB,Constants.ORA_ID_MANUAL_REVIEW,Constants.ORA_ID_HOMOZYGOUS,String.valueOf(mutation1Base),missBase);
            MutationDetails mutDetailsObj=new MutationDetails(newSampScan[i],0,mutation1[1],0,curMutationScore1,0.0f,conf,getHomUniqueness(mutation1[3]),getHomContext(pos));
            output.addMutation(mutObj,mutDetailsObj);
          }
        } else {
          reportMissingHomMutation(pos,0,missBase,0,newSampScan[i],cTrace,sTrace);
        }
      } else if ( missingBases[1][i] == AutoCSA.MISSING_IN_BOTH ) {    // missing in both
        if( compIndex > -1 ) {
          b1_scan = controlScan[compIndex];
        } else {
          b1_scan = newCtrlScan[i-1];
        }
        for( j = 0 ; j< 4 ; ++j ) {  //  cycle bases
          hom_base=baseList[j];

          if( missBase.equals(hom_base) ) { continue; }

          m_scan=getHomSearchCentre(base1,hom_base,pos,b1_scan,controlScan);
//  use increased bin (using binFactor) as important to locate potential
//  hom peak
//        if(log.isInfoEnabled()) log.info("Base: "+base1+" "+pos+" Try: "+hom_base+" "+b1_scan+" "+m_scan);
          index=locatePeak(cTrace,m_scan,hom_base,binFactor);
          if( index > 0 ) {
// check viability of mut peak i.e. is it comparable to locality
            int[] ok=checkHomMutant(cTrace,pos,index,hom_base);
            if( ok[0] == 1 ) {
              homPeak= cPeakHt[index];
              homScan= cScan[index];
              ++mutation1[3]; // No. of possible mutants
              peaks[j]=homPeak;
            } else {
              homPeak=0;      // potential hom peak too small or dye blob
              homScan= 0;
            }
          }
          if( homPeak > mutation1[2] ) {
            mutation1[0] = index;
            mutation1[1] = homScan;
            mutation1[2] = homPeak;
            mutation1Base= hom_base.charAt(0);
          }
        }
        if( mutation1[2] > 0 ) {
          int minAmp=cTrace.getMinPeakIntensity(String.valueOf(mutation1Base));
          curMutationScore1=getHomPeakScore(peaks,minAmp);
          conf=getHomOverallConfidence(pos,curMutationScore1,mutation1[3],cTrace);
//        if(log.isInfoEnabled()) log.info("HOM SNP_WT at Bp: "+pos+" ("+missBase+") Change:"+
//                mutation1Base+" Scan: "+mutation1[1]+" Ht: "+mutation1[2]+
//                 " Total: "+mutation1[3]+" Score: "+curMutationScore1+
//                 " Confidence: "+conf);
          if( mutation1[0] > 0 ) {
            cTrace.setRealPeak(mutation1[0],AutoCSA.HOM_SNP_PEAK);
          }
        }
        if( compIndex > -1 ) {
          b1_scan = sampleScan[compIndex];
        } else {
          b1_scan = newSampScan[i-1];
        }
        homPeak=0;    // reset for next search
        for( j = 0 ; j< 4 ; ++j ) {  //  cycle bases
          hom_base=baseList[j];

          if( missBase.equals(hom_base) ) { continue; }
          m_scan=getHomSearchCentre(base1,hom_base,pos,b1_scan,sampleScan);
//  use increased bin (using binFactor) as important to locate potential
//  hom peak
//        if(log.isInfoEnabled()) log.info("Base: "+base1+" "+pos+" Try: "+hom_base+" "+b1_scan+" "+m_scan);
          index=locatePeak(sTrace,m_scan,hom_base,binFactor);
          if( index > 0 ) {
// check viability of mut peak i.e. is it comparable to locality
            int[] ok=checkHomMutant(sTrace,pos,index,hom_base);
            if( ok[0] == 1 ) {
              homPeak= sPeakHt[index];
              homScan= sScan[index];
              ++mutation2[3]; // No. of possible mutants
              peaks[j]=homPeak;
            } else {
              homPeak=0;      // potential hom peak too small or dye blob
              homScan= 0;
            }
          }
          if( homPeak > mutation2[2] ) {
            mutation2[0] = index;
            mutation2[1] = homScan;
            mutation2[2] = homPeak;
            mutation2Base= hom_base.charAt(0);
          }
        }
        if( mutation2[2] > 0) {
          int minAmp=sTrace.getMinPeakIntensity(String.valueOf(mutation2Base));
          curMutationScore2=getHomPeakScore(peaks,minAmp);
          conf=getHomOverallConfidence(pos,curMutationScore2,mutation2[3],sTrace);
//        if(log.isInfoEnabled()) log.info("HOM SNP_TT at Bp: "+pos+" ("+missBase+") Change:"+
//                mutation2Base+" Scan: "+mutation2[1]+" Ht: "+mutation2[2]+
//                 " Total: "+mutation2[3]+" Score: "+curMutationScore2+
//                 " Confidence: "+conf);
          if( mutation2[0] > 0 ) {
            sTrace.setRealPeak(mutation2[0],AutoCSA.HOM_SNP_PEAK);
          }
        }
        newSampScan[i]=mutation2[1];   // store sample scan
        newCtrlScan[i]=mutation1[1];   // store control scan
// if both mut bases are same then we have a HOM SNP in both samples
// if both mut bases are diff then we have a HOM SNP & a HOM MUT
        if( mutation1[2] > 0 && mutation2[2] > 0) {
          // create a Mutation & MutationDetails object (wrt Ctrl here)
          if( output != null ) {
            int wtMutType;
            if( mutation1Base == mutation2Base ) {
              wtMutType=Constants.ORA_ID_GERMLINE_SUB;
            } else {
              wtMutType=Constants.ORA_ID_REFERENCE_SUB;
            }
            conf=getHomOverallConfidence(pos,curMutationScore1,mutation1[3],cTrace);
            Mutation mutObj= new Mutation(pos,pos,wtMutType,Constants.ORA_ID_MANUAL_REVIEW,Constants.ORA_ID_HOMOZYGOUS,String.valueOf(mutation1Base),missBase);
            MutationDetails mutDetailsObj=new MutationDetails(mutation2[1],0,mutation1[1],0,curMutationScore1,0.0f,conf,getHomUniqueness(mutation1[3]),false); // no context change
            output.addMutation(mutObj,mutDetailsObj);
          }
          if( mutation1Base == mutation2Base ) {
            newSeq[i]=mutation1Base;
          } else {
// also record HOM MUT for sample wrt Parent seq
            if( output != null ) {
              conf=getHomOverallConfidence(pos,curMutationScore2,mutation2[3],sTrace);
              Mutation mutObj= new Mutation(pos,pos,Constants.ORA_ID_SUBSTITUTION,Constants.ORA_ID_MANUAL_REVIEW,Constants.ORA_ID_HOMOZYGOUS,String.valueOf(mutation2Base),missBase);
              MutationDetails mutDetailsObj=new MutationDetails(mutation2[1],0,mutation1[1],0,curMutationScore2,0.0f,conf,getHomUniqueness(mutation2[3]),getHomContext(pos));
              output.addMutation(mutObj,mutDetailsObj);
            }
          }
        } else if( mutation2[2] > 0) {
// i.e. case where alternate base found in Tumour but not Normal
// record HOM MUT for sample wrt Parent seq
          if( output != null ) {
            conf=getHomOverallConfidence(pos,curMutationScore2,mutation2[3],sTrace);
            Mutation mutObj= new Mutation(pos,pos,Constants.ORA_ID_SUBSTITUTION,Constants.ORA_ID_MANUAL_REVIEW,Constants.ORA_ID_HOMOZYGOUS,String.valueOf(mutation2Base),missBase);
            MutationDetails mutDetailsObj=new MutationDetails(mutation2[1],0,mutation1[1],0,curMutationScore2,0.0f,conf,getHomUniqueness(mutation2[3]),getHomContext(pos));
            output.addMutation(mutObj,mutDetailsObj);
          }
// i.e. case where alternate base found in Tumour but not Normal
// record HOM MUT for sample wrt Parent seq
        } else {
// n.b. this case encompasses situations where maybe find an alternate
// base in the Normal but not the Tumour - report as a trace hole wrt Tumour
// and a reference sub (if we find something in the normal)
          reportMissingHomMutation(pos,0,missBase,newCtrlScan[i],newSampScan[i],cTrace,sTrace);
          if( mutation1[2] > 0 ) {
            if( output != null ) {
              conf=getHomOverallConfidence(pos,curMutationScore1,mutation1[3],cTrace);
              Mutation mutObj= new Mutation(pos,pos,Constants.ORA_ID_REFERENCE_SUB,Constants.ORA_ID_MANUAL_REVIEW,Constants.ORA_ID_HOMOZYGOUS,String.valueOf(mutation1Base),missBase);
              MutationDetails mutDetailsObj=new MutationDetails(mutation2[1],0,mutation1[1],0,curMutationScore1,0.0f,conf,getHomUniqueness(mutation1[3]),false);
              output.addMutation(mutObj,mutDetailsObj);
            }
          }
        }
      } else if ( missingBases[1][i] == AutoCSA.ASSIGNED_DELETION ) {
// position has been identified as deletion
      } else {
        if(log.isInfoEnabled()) log.info("homMutationScan: Unknown Case");
      }
    }
  }
  private int getComparisonIndex(int bpPos) {
    int i;
    for( i=0; i < nPoints; ++i ) {
      if( baseNo[i] == bpPos ) { break; }
    }
    return i < nPoints ? i : -1;
  }
  private int getHomPeakScore(int[] peaks, int minAmp) {
    int i;
    int maxval=0;
    int nextMaxval=0;
    for( i=0; i < peaks.length ; ++i ) {
      maxval=Math.max(maxval,peaks[i]);
    }
    for( i=0; i < peaks.length ; ++i ) {
      if( peaks[i] == maxval ) { continue; }
      nextMaxval=Math.max(nextMaxval,peaks[i]);
    }
// define hom score as ratio of hom peak to next highest peak in window
// i.e this is equivalent to a base quality score
    if( nextMaxval == 0 ) {nextMaxval=minAmp; }
    return Math.round( (float) maxval/ (float) nextMaxval);
  }
  private int getZygosityCode(float scale ) {
    if( scale > 0.45 && scale < 0.55 ) { return Constants.ORA_ID_HETEROZYGOUS; }
    if( scale < 0.05 ) { return Constants.ORA_ID_HOMOZYGOUS; }
    return Constants.ORA_ID_AMPLIFICATION;
  }
/**
 * Sets new values for the comparison parameters
 * @param mutratio the critical peak drop factor to use
 * @param bin the window half width to search for Het peaks
 */
  public void setComparisonParams(float mutratio, float bin) {
    critMutRatio = mutratio;
    peakSearchBin=bin;
  }
  private int locateAnalysisUpstream(int pos) {
    int i;
    for( i=pos-1; i > -1; --i ) {
      if( getComparisonIndex(i) > -1 ) { break; }
    }
    return i > -1 ? i : -1;
  }
  private void reportMissingHomMutation(int pos, int mutType, String missBase, int ctrlScan, int sampScan, SeqTraceAnalysis cTrace, SeqTraceAnalysis sTrace) {

      int index;
//    if(log.isInfoEnabled()) log.info("Unable to Determine HOM type at Bp: "+pos);
      // create a Mutation & MutationDetails object
      if( output != null ) {
        Mutation mutObj= new Mutation(pos,pos,Constants.ORA_ID_TRACE_HOLE,Constants.ORA_ID_MANUAL_REVIEW,Constants.ORA_ID_HOMOZYGOUS,String.valueOf('N'),missBase);  // declare as a trace hole
        int posToLeft=locateAnalysisUpstream(pos);
        int estCScan,estSScan;
        estCScan=ctrlScan;
        estSScan=sampScan;
        if ( posToLeft > -1 ) {
          if ( estCScan == 0 ) {
            index = cTrace.getAnalysisIndex(posToLeft);
            if ( index > -1 ) {
              int[] cScan=cTrace.getScan();
              estCScan=cScan[index]+(pos-posToLeft)*Constants.AVERAGE_BASE_SPACING;
            }
          }
          if ( estSScan == 0 ) {
            index = sTrace.getAnalysisIndex(posToLeft);
            if ( index > -1 ) {
              int[] sScan=sTrace.getScan();
              estSScan=sScan[index]+(pos-posToLeft)*Constants.AVERAGE_BASE_SPACING;
            }
          }
        }
        MutationDetails mutDetailsObj=new MutationDetails(estSScan,0,estCScan,0,0,0.0f,0,false,false);
        output.addMutation(mutObj,mutDetailsObj);
      }
  }
  private boolean getHomContext(int pos ) {
    int compIndex=getComparisonIndex(pos+1);
    if( compIndex == -1 ) { return false; }   // no data available
    if( scale[compIndex] > 0.4 && scale[compIndex] < 0.6 ) { return true; }
    if( scale[compIndex] > 1.5 ) { return true; }
    return false;
  }
  private boolean getHomUniqueness(int numHomPeaks ) {
    if ( numHomPeaks == 1 ) { return true; }
    return false;
  }
  private int[] checkViableMutant(SeqTraceAnalysis sTrace, int pos, int mutIndex, String mut_base, boolean hetsInControl ) {
    int[] realPeaks= {0,0,0,0};
    int[] sPeakHt = sTrace.getIntensity();
    int[] sScan = sTrace.getScan();
    int scan,index;
    int toLeft=-1;
    int toRight=1;
    int realPeakCode=AutoCSA.CALLED_PEAK;
    float GPeakCrit,ACTPeakCrit;
    if( hetsInControl ) {
      GPeakCrit=0.30f;
      ACTPeakCrit=0.45f;
    } else {
      GPeakCrit=0.15f;
      ACTPeakCrit=0.23f;
    }
    scan=sScan[mutIndex];
    index=locateAdjacentPeak(sTrace,scan,mut_base,toLeft,realPeakCode);
// can only continue looking if we've found first
    if( index > 0 ) {
      realPeaks[1]=sPeakHt[index];
      scan=sScan[index];
      index=locateAdjacentPeak(sTrace,scan,mut_base,toLeft,realPeakCode);
    }
// assign next peak if search successfull
    if( index > 0 ) {
      realPeaks[0]=sPeakHt[index];
    }
    scan=sScan[mutIndex];
    index=locateAdjacentPeak(sTrace,scan,mut_base,toRight,realPeakCode);
// can only continue looking if we've found first
    if( index > 0 ) {
      realPeaks[2]=sPeakHt[index];
      scan=sScan[index];
      index=locateAdjacentPeak(sTrace,scan,mut_base,toRight,realPeakCode);
    }
// assign next peak if search successfull
    if( index > 0 ) {
      realPeaks[3]=sPeakHt[index];
    }
// check we've found some peaks
    int count=0;
    int i;
    for( i=0; i< realPeaks.length; ++i ) {
      if( realPeaks[i] > 0 ) {++count; }
    }
// can't dismiss potential mutant as no data available
// otherwise take median of peaks that are found (up to 4)
    if( count == 0 ) {
//if(log.isInfoEnabled()) log.info("Real Peaks: Pos: "+pos+" Not Enough Data");
      return new int[] {1,-1};
     }  // can't dismiss
    float[] temp= new float[count];
    count=0;
    for( i=0; i< realPeaks.length; ++i ) {
      if( realPeaks[i] > 0 ) {
        temp[count++]=(float) realPeaks[i];
      }
    }
    float medianPeak=median(temp);
    int pc=(int) (100.0*(float) sPeakHt[mutIndex] /medianPeak);
    float critPeakRatio;
    if ( mut_base.equals("G") ) {
      critPeakRatio=GPeakCrit;   // very small G peaks observed (+muts)
    } else {
      critPeakRatio=ACTPeakCrit;
    }
    if( (float) sPeakHt[mutIndex] /medianPeak > critPeakRatio ) {
      return new int[] {1,pc};     // viable
    } else {
      return new int[] {0,pc};     // dismiss as not-viable
    }
  }
  private int[] checkHomMutant(SeqTraceAnalysis sTrace, int pos, int mutIndex, String hom_base ) {

// check viability of mut peak i.e. is it comparable to locality
    int[] ok=checkViableMutant(sTrace,pos,mutIndex,hom_base,false);
//if(log.isInfoEnabled()) log.info("HOM Peaks: Pos: "+pos+" "+hom_base+" PC: "+ok[1]);
    if( ok[1] < 35 ) {
      return new int[] {0,ok[1]};     // dismiss as not-viable
    }
// examine % multiple of median locality to check for dye blobs
    if( ok[1] > 1000 ) {
      return new int[] {0,ok[1]};     // dismiss as prob dye blob
    }
// test if we're in a dye blob, i.e. if either neighbour is
// dye blob affected then reject
    int[] bpPos=sTrace.getbpPos();    // get analysed base No
    if( sTrace.isBaseDyeBlob(pos-1) || sTrace.isBaseDyeBlob(pos+1) ) {
      return new int[] {0,ok[1]};     // dismiss as prob dye blob
    }
    return new int[] {1,ok[1]};       // viable
  }
  private float[] setHetSearchBins(float maxBin,int index) {
    float halfBaseSpace;
    int pos=baseNo[index];
// set default bins
    float[] bins = new float[] {maxBin,maxBin};
// must use uncentred bins if data is not mobility corrected
    if( Adjustment.getAdjustmentType() == 0 ) {
      return bins;
    }
// check scan to LHS if adjacent base is present
    if( index > 0 && baseNo[index-1] == pos-1 ) {
      halfBaseSpace=(float) Math.ceil(0.5f*(float)(sampleScan[index]-sampleScan[index-1]));
      bins[0]=Math.min(maxBin,halfBaseSpace);
    }
// check scan to RHS if adjacent base is present
    if( index < nPoints-1 && baseNo[index+1] == pos+1 ) {
      halfBaseSpace=(float) Math.ceil(0.5f*(float)(sampleScan[index+1]-sampleScan[index]));
      bins[1]=Math.min(maxBin,halfBaseSpace);
    }
//if(log.isInfoEnabled()) log.info("Bins: "+pos+" "+bins[0]+" "+bins[1]);
// return bins which are limited by half distance to neighbouring bases
    return bins;
  }
  private float getHomSearchCentre(String base1,String hom_base,int pos,float b1_scan, int[] scan) {
    float m_scan;
// must use uncentred bins if data is not mobility corrected
    m_scan=b1_scan+Adjustment.getAdjustmentOffset(base1,hom_base,pos-1,pos);
    if( Adjustment.getAdjustmentType() == 0 ) {
      return m_scan;
    }
// if we have a snp, i.e. data either side then use central scan
    int compIndexLhs=getComparisonIndex(pos-1);
    int compIndexRhs=getComparisonIndex(pos+1);
    if( compIndexLhs > -1 && compIndexRhs > -1 ) {
      float scanLhs = (float) scan[compIndexLhs];
      float scanRhs = (float) scan[compIndexRhs];
      m_scan=(float) Math.ceil(0.5f*(scanLhs+scanRhs));
    }
    return m_scan;
  }
/**
 * Performs a test for Homozygous Deletions.
 * @param cTrace the Control trace SeqTraceAnalysis object
 * @param sTrace the Sample trace SeqTraceAnalysis object
 * @param missingBases the 2-D information array for trace holes
 */
  public void homDeletionScan(SeqTraceAnalysis cTrace, SeqTraceAnalysis sTrace, int[][] missingBases) {
    int i,j,pos,type;

    i=0;
    int startI=0;
// search missing base list for contigs of holes in same trace
    while( i < missingBases[0].length ) {
      startI=i;
      pos=missingBases[0][i];
      type=missingBases[1][i];
      if( type == AutoCSA.MISSING_IN_NORMAL ) {
        ++i;
        continue;
      }
      ++i;
      ++pos;
      while( i < missingBases[0].length && pos == missingBases[0][i] && missingBases[1][i] != AutoCSA.MISSING_IN_NORMAL ) {
        ++i;
        ++pos;
      }
// i index here is 1 after end of hom region (maybe single base)
      int len=i-1-startI+1;
      int compIndex1=getComparisonIndex(missingBases[0][startI]-1);
      int compIndex2=getComparisonIndex(missingBases[0][i-1]+1);
      int sScan1,sScan2;
// n.b must test for missing adjacent base in analysis - will be due
// to adjacent base missing in Normal (otherwise gap in sample bigger)
      if( compIndex1 == -1 ) {
        sScan1=sTrace.getScanAtBase(missingBases[0][startI]-1);
      } else {
        sScan1=sampleScan[compIndex1];
      }
      if( compIndex2 == -1 ) {
        sScan2=sTrace.getScanAtBase(missingBases[0][i-1]+1);
      } else {
        sScan2=sampleScan[compIndex2];
      }
      float totSpacing=(float) (sScan2-sScan1);
      float avSpacing=sTrace.getAverageBaseSpacing(missingBases[0][startI]-1,5);
//  compare our total spacing with average spacing to see if we 
//  have a likely hom deletion - use 25% limits on ratio
//  n.b. if we can't find an avSpacing (i.e =-1) ratio is large -Ve
//  which will pass test for a del - rare cases at start of coverage
      float critRatio=0.25f;
      float ratio=(totSpacing-avSpacing)/avSpacing;
// if > 1 base missing and spacing < 15 then definately flag
      if ( len > 1 && totSpacing < Constants.AVERAGE_BASE_SPACING+3 ) {
        ratio=0.1f;
      }
// test if we have missing bases but scan spacing > 1 base
      if ( ratio > critRatio && len > 1 && totSpacing > Constants.AVERAGE_BASE_SPACING+3 ) {
        int cScan1;
        if( compIndex1 == -1 ) {
          cScan1=cTrace.getScanAtBase(missingBases[0][startI]-1);
        } else {
          cScan1=controlScan[compIndex1];
        }
        int delPos=pos-len;
        String delBases=refSeq.substring(pos-1-len,pos-1);
// check if we have a complex i.e. DEL & INS
        boolean complex=true;
        String insBases="N";
        int insPos=0;
        if( true ) {
          int iLen=Math.round(totSpacing/avSpacing)-1;
          insPos=delPos-1;
          if ( iLen > 0 && iLen <= Constants.MAX_HOM_COMPLEX_INS ) {
            insBases=findInsertBases(sTrace,insPos,iLen,sScan1,sScan2,0);
// reject if any base is N (n.b. findInsertBases may shorten if last is N)
            if ( insBases.indexOf("N") > -1 ) complex=false;
          } else {
            complex=false;
          }
        }
        int check=1;
        if( complex ) check=validateComplexIndelCall(delPos,delBases,insBases,sTrace);
// validation test tries to work out if indel is most likely due to poor
// local trace, if so reject indel call but still turn off trace holes
        if( check == 0 ) {
          if( output != null ) {
            int cScan2;
            if( compIndex2 == -1 ) {
              cScan2=cTrace.getScanAtBase(missingBases[0][i-1]+1);
            } else {
              cScan2=controlScan[compIndex2];
            }
            checkComplexForSimplification(cTrace,sTrace,delPos,len,insBases,delBases,cScan1,cScan2,sScan1,sScan2);
          }
          mutationsPresent=true;
        } else {
          if( check > 0 ) reportSpeculativeIndel(delPos, cScan1, sScan1,len,1);
// turn off single missing points i.e. as found to be CSA matching error
          if( check == -1 ) {
            for( i=startI; i < startI+len; ++i ) {
              missingBases[1][i] = AutoCSA.MATCHING_ERROR;
            }
            continue;
          }
        }
// turn off single missing points i.e. as found to be a spec HOM INDEL
// n.b. don't turn off first in list as could be single snp
        for( i=startI; i < startI+len; ++i ) {
          missingBases[1][i] = AutoCSA.SPECULATIVE_INDEL;   // code spec
        }
        continue;
      }
      int score=25;
      int conf=2;
      int delPos=pos-len;
      String delBases=refSeq.substring(pos-1-len,pos-1);

      if( ratio < critRatio ) {
// check if it's really a case of the base having been missed in the
// matching process, i.e. serach for delBases (single) in [pos-2,pos-1]
        if( len == 1 ) {
          int searchFilter = 999;
          int searchType=1;    // > lim1 && < lim2
          int searchLimit1=sTrace.getScanAtBase(delPos-2) + 2;
          int searchLimit2=sTrace.getScanAtBase(delPos-1) - 2;
          int index=searchData(sTrace,searchType,searchLimit1,searchLimit2,searchFilter,delBases);
          if( index > 0 ) {
// check viability of peak found i.e. is it comparable to locality
              int[] ok=checkHomMutant(sTrace,delPos,index,delBases);
              if( ok[0] == 1 ) {
// turn off single missing points i.e. as found to be CSA matching error
                missingBases[1][startI] = AutoCSA.MATCHING_ERROR;
                continue;
              }
          }
        }
        if( ratio > 0.0f && ratio < 0.125f ) {
          score=50;
          conf=1;
        }
//      if(log.isInfoEnabled()) log.info("HOM DEL: Pos: "+delPos+" Length: "+len+" DelSeq: "+delBases+" Score: "+score+" Confidence: "+conf);
//      if(log.isInfoEnabled()) log.info("HOM DEL: "+len+" "+missingBases[0][startI]+" "+avSpacing+" "+totSpacing);
// turn off single missing points i.e. as found to be a HOM DEL
        for( i=startI; i < startI+len; ++i ) {
          missingBases[1][i] = AutoCSA.ASSIGNED_DELETION;   // code for deletion
        }
        // create a Mutation & MutationDetails object
        if( output != null ) {
          int cScan1,cScan2;
          if( compIndex1 == -1 ) {
            cScan1=cTrace.getScanAtBase(missingBases[0][startI]-1);
          } else {
            cScan1=controlScan[compIndex1];
          }
          if( compIndex2 == -1 ) {
            cScan2=cTrace.getScanAtBase(missingBases[0][i-1]+1);
          } else {
            cScan2=controlScan[compIndex2];
          }
          Mutation mutObj= new Mutation(delPos,delPos+len-1,Constants.ORA_ID_DELETION,Constants.ORA_ID_MANUAL_REVIEW,Constants.ORA_ID_HOMOZYGOUS,"-",delBases);
          MutationDetails mutDetailsObj=new MutationDetails(sScan1,sScan2,cScan1,cScan2,score,0.0f,0,false,false);
          output.addMutation(mutObj,mutDetailsObj);
        }
        mutationsPresent=true;
      }
    }
  }
/**
 * Performs a test for Homozygous Insertions.
 * @param cTrace the Control trace SeqTraceAnalysis object
 * @param sTrace the Sample trace SeqTraceAnalysis object
 * @param missingBases the 2-D information array for trace holes
 */
  public void homInsertionScan(SeqTraceAnalysis cTrace, SeqTraceAnalysis sTrace, int[][] missingBases) {
    int i,j,pos,type;

    int windowLen=5;    // window size for av base spacing
    int[] samplePos=sTrace.getbpPos(); // get analysed base No
    boolean deletionExists;
    i=0;
// loop over all base spacings to check for spacings a certain fraction
// larger than the local upstream average base spacing
// n.b. can't start main loop at 0 as need upstream av base spacing
// reasonable to assume no insertions in first "windowLen" bases
// as have to match 15 consecutive bases at start
    for( i=windowLen; i< nPoints; ++i ) {
// firstly check if a variant (snp) already exists at current location
      if( hasBaseVariant(baseNo[i-1]) ) {
        continue;
      }
// next check if a deletion exists upstream within windowLen
      deletionExists=false;
      if( missingBases != null ) {
        for(j=0; j < missingBases[0].length ; ++j ) {
          if( missingBases[1][j] == AutoCSA.ASSIGNED_DELETION ) {  //test for HOMDEL
            if( missingBases[0][j] >= baseNo[i-1-windowLen+1] && missingBases[0][j] <= baseNo[i-1] ) {
              deletionExists=true;
            }
          }
        }
      }
      int baseDif1=baseNo[i]-baseNo[i-1];
      int baseDif2=baseNo[i-1]-baseNo[i-2];
      if( baseDif1 > 1 ||  baseDif2 > 1 ) {
        continue;    // skip as a trace hole (unreliable matching)
      }
      int scanDif=sampleScan[i]-sampleScan[i-1];
      float avSpacing=sTrace.getAverageBaseSpacing(baseNo[i-1],windowLen);
      if( avSpacing < 0.0f ) {   // check we have avSpacing defined
        continue;
      }
//  check avSpacing by looking for trace holes in the window
      int totBaseDif=0;
      int index=sTrace.getAnalysisIndex(baseNo[i-1]);
      if( index > -1 ) {
        totBaseDif=samplePos[index]-samplePos[index-windowLen+1];
      }
// correct avSpacing if affected by trace holes
// however don't correct if deletion exists upstream within windowLen
      if( ! deletionExists && totBaseDif > windowLen-1 ) {
        avSpacing*= (float) (windowLen-1)/(float) (totBaseDif);
      }
//  compare our local spacing with average spacing to see if we 
//  have a likely hom insertion - use a ratio
//  n.b. unreliable test if we have > 20 bases missing upstream
      float scanRatio=(float) scanDif/avSpacing;
      int len=Math.round(scanRatio)-1;
      int sScan1=sampleScan[i-1];
      int sScan2=sampleScan[i];
      float critScanRatio=1.75f;
      int critHomPC=0;
      if( len == 1 ) {
        int np=sTrace.getNumViablePeaks(sScan1+4,sScan2-4,0);
// if we have 1 peak relax critScanRatio but enforce a strict mut %
// n.b. only do this for examples that would have failed original test
        if( np == 1 && scanRatio < critScanRatio ) {
          critScanRatio=1.6f;
          critHomPC=75;
        }
      }
      if( totBaseDif < 20 && scanRatio > critScanRatio ) {
// check if we have enough unassigned peaks for insertion
// assume peak spacing must be at least 4 scans from neighbours
        int numPeaksInRegion=sTrace.getNumViablePeaks(sScan1+4,sScan2-4,0);
        if( numPeaksInRegion < len ) {
          continue;     // Not enough peaks for insertion
        }
        int insPos=baseNo[i-1];
        String insBases="N";
        if ( len < 11 ) {
          insBases=findInsertBases(sTrace,insPos,len,sScan1,sScan2,critHomPC);
// reject if first base is N, or if len < 6 and some N's returned
          if ( insBases.substring(0,1).equals("N") ) {continue; }
          if ( len < 6 && insBases.indexOf("N") > -1 ) {continue; }
          len=insBases.length();  // reset in case len changed
        } else {
          for( j=1; j<len ; ++j) {
            insBases=insBases.concat("N");
          }
        }
// check if single base ins are really caused by an overloaded peak
// that has been split into 2, each being recognised as amplimer peaks,
// i.e, must have a doublet, one of which must be Overloaded, also the
// ins base must be the same type as the doublet
        if( len == 1 ) {
          if( insBases.equals(refSeq.substring(insPos-1,insPos)) &&
              insBases.equals(refSeq.substring(insPos-2,insPos-1)) &&
              ( sTrace.isBaseOverloaded(insPos,28000) ||
                sTrace.isBaseOverloaded(insPos-1,28000) ) &&
                scanDif <= 20 ) continue;
        }
        int score=50;
        int conf=3;
// check upstream scan diff to check if a small value here caused
// large scan diff at suspected insertion point
        int prevScanDif=sampleScan[i-1]-sampleScan[i-2];
        if( (float) prevScanDif/avSpacing < 0.75f ) {
          score=25;
        } else {
          --conf;
        }
        if( scanRatio > 1.9f ) {
          --conf;
        }
//      if(log.isInfoEnabled()) log.info("HOM INS: Pos: "+insPos+" Length: "+len+" InsSeq: "+insBases+" Ratio: "+scanRatio+" Score: "+score+" Confidence: "+conf);
        // create a Mutation & MutationDetails object
        if( output != null ) {
          int cScan1=controlScan[i-1];
          int cScan2=controlScan[i];
          Mutation mutObj= new Mutation(insPos,insPos+1,Constants.ORA_ID_INSERTION,Constants.ORA_ID_MANUAL_REVIEW,Constants.ORA_ID_HOMOZYGOUS,insBases,"-");
          MutationDetails mutDetailsObj=new MutationDetails(sScan1,sScan2,cScan1,cScan2,score,0.0f,0,false,false);
          output.addMutation(mutObj,mutDetailsObj);
        }
        mutationsPresent=true;
      }
//if(log.isInfoEnabled()) log.info("HOM: "+baseNo[i]+" "+baseDif+" "+scanDif+" "+totBaseDif+" "+avSpacing);
    }
  }
  private float compareMutationSeq(StringBuffer mutationSeq, int refPos, int totalMissing, float critRatio, int maxHetDelLength, int maxHetInsLength) {
// shift parameter is assigned as combo of max possible het ins+del length
  int shift=-maxHetInsLength;
// n.b have to apply toString() method on StringBuffer and
// reinstantiate as memory address's are copied if using straight copy
  StringBuffer origMutationSeq= new StringBuffer(mutationSeq.toString());
  int ndifMin=1000;
  int nLimit=maxHetDelLength+maxHetInsLength+1;
  for( int n=0 ; n < nLimit ; ++n ) {
    int ndif=0;
    int incr=shift+n;
    if( incr == 0 ) continue;  // don't consider zero shift
    if (refPos+incr < 1 ) {
      continue;       // prevents going beyond start of refSeq
    }
    StringBuffer newMutationSeq= new StringBuffer(origMutationSeq.toString());
    for( int i=0; i < mutationSeq.length() ; ++i) {
      String refBase=refSeq.substring(refPos-1+incr,refPos+incr);
      if( ! origMutationSeq.substring(i,i+1).equals(refBase) ) {
        ++ndif;
        newMutationSeq.setCharAt(i,refBase.charAt(0));
      }
      ++incr;
    }
    if( ndif < ndifMin ) {
      ndifMin=ndif;
// n.b need to replace entire contents of str buffer rather than
// re-instantiating as need to ensure same memory address is used
// as we reply on mutationSeq being updated via argument list
      mutationSeq=mutationSeq.replace(0,mutationSeq.length(),newMutationSeq.toString());
    }
int ii=shift+n;
//if(log.isInfoEnabled()) log.info("HET Diffs: shift:"+ii+" "+ndif+" "+ndifMin);
  }
// calculate the ration of corrections we need to make to the refSeq
// but subtract off any trace holes in the mut seq
  int realLen=mutationSeq.length()-totalMissing;
  ndifMin-=totalMissing;
//if(log.isInfoEnabled()) log.info("HET Diffs Min: "+ndifMin+" "+realLen);
  float ratioOfCorrections= (float) ndifMin/ (float) realLen;
  if( ratioOfCorrections > critRatio ) {
    mutationSeq=origMutationSeq;
  }
  return ratioOfCorrections;
  }
  private int getHomOverallConfidence(int pos, int score, int numPossPeaks, SeqTraceAnalysis seqTrace) {
    int confidence=0;      // initial setting to be modified
    double dScore= (double) score;
    int count=0;
    if( seqTrace.getAnalysisIndex(pos-1) == -1 ) { ++count; }
    if( seqTrace.getAnalysisIndex(pos+1) == -1 ) { ++count; }
    confidence= (int) (100.0f*(1.0f-0.1f*(numPossPeaks-1))
                     * (Math.log(2.0*Math.min(5.0,dScore))/Math.log(10.0))
                     * (1.0f-0.25f*count) );
    return confidence;
  }
  private int getHetOverallConfidence(int index, int score, int localMutantRatio, SeqTraceAnalysis cTrace, SeqTraceAnalysis sTrace) {
    int confidence=0;      // initial setting to be incremented
    float peakDrop=scale[index];
    double dScore= (double) score;
//  n.b if localMutantRatio=-1 reset so factor used is 0.9
    if( localMutantRatio == -1 ) { localMutantRatio=45; }
// n.b need to use log10(x) = ln(x) / ln(10) as no log10 in J1.4.2
    confidence= (int) (100.0f*Math.min(1.0f,(11.0f-10.0f*peakDrop)/6.0f)
                    * (Math.log(Math.min(10.0,dScore))/Math.log(10.0))
                     * 0.02f* Math.min(50,localMutantRatio) );
    return confidence;
  }
  private boolean hasBaseVariant(int pos) {
// check if pos is at a position of a HET SNP
    int compIndex=getComparisonIndex(pos);
    String mut=(String) mutations.get(compIndex);
    if ( mut != null ) { return true; }
// check if pos is at a position of a HOM SNP (or hole)
    if ( missingBases != null ) {
      for( int i=0; i < missingBases[0].length; ++i ) {
        if ( missingBases[0][i] == pos ) {
          return true;
        }
      }
    }
    return false;
  }
  private String findInsertBases(SeqTraceAnalysis sTrace, int pos,int len, int sScan1, int sScan2, int stricterMutPC) {
      String baseList[]={"A","C","G","T"};
      int[] sPeakHt = sTrace.getIntensity();
      int[] sScan = sTrace.getScan();
      char[] insBases= new char[len];
      String hom_base;
      float binFactor=1.0f;
      int bestPeak=0,homPeak=0;
      int index=0,bestIndex=0;

// loop over No. of estimated inserted bases
      for( int i=0; i < len; ++i ) {
        insBases[i]= 'N';
// n.b. finding offset below only applicable to mobility corrected data
        int offset=(int) Adjustment.getAdjustmentOffset("N","N",pos-1,pos);
        if( bestIndex > 0 ) {
          sScan1=sScan[bestIndex]+offset;
        } else {
          sScan1+=offset;
        }
        bestPeak =0;  // reset for each trial
        bestIndex=0;
        for( int j = 0 ; j< 4 ; ++j ) {  //  cycle trial bases
          hom_base=baseList[j];

          index=locatePeak(sTrace,sScan1,hom_base,binFactor);
          homPeak=0; // must set
          if( index > 0 ) {
//          if( len == 1 ) {
              if( sScan[index] > sScan2-4 ) { continue; }  // too close
//          }
// check viability of mut peak i.e. is it comparable to locality
// n.b. trial peaks with no like neighbours return -1 in ok[1]
// these are rejected by checkHomMutant
            int[] ok=checkHomMutant(sTrace,pos,index,hom_base);
            if( ok[0] == 1 ) {
              homPeak= sPeakHt[index];
              if( ok[1] > 0 && ok[1] < stricterMutPC ) homPeak=0; // turn off
            } else {
              homPeak=0;      // potential hom peak too small or dye blob
            }
          }
// save largest insert peak for each ins base location
          if( homPeak > bestPeak ) {
            bestPeak = homPeak;
            insBases[i]= hom_base.charAt(0);
            bestIndex=index;
          }
        }
        if( bestIndex > 0 ) {
          sTrace.setRealPeak(bestIndex,AutoCSA.HOM_INS_PEAK);
        }
        ++pos;
      }
      String result=new String(insBases);
// if len > 4 and cannot find last base - assume expected len 1 too big
      if( len > 4 ) {
        if ( result.indexOf("N") > len-2 ) result=result.substring(0,len-1);
      }
      return result;
  }
/**
 * Saves and sets the supplied Het Indel position to a valid value for use in Het Indel detection.
 * The supplied value has 5 subtracted so the real Het position is not passed.
 * Also, the resultant position may be altered so it's a base number that has been analysed.
 * @param ipos the initial position of the suspected Het Indel
 */
  public void setIndelBasePosition(int ipos) {
//  subtract 5 bases off suspected indel pos for safety
//  n.b. ipos-5 always > 5 due to how ipos calculated
    indelBasePosFromQ=(ipos > 0 ? ipos-5 : 0);
    if( indelBasePosFromQ > 0 ) {
// make sure indelBasePosFromQ corresponds to a compared base
      for(int n=indelBasePosFromQ ; n>indelBasePosFromQ-50; --n ) {
        int index=getComparisonIndex(n);
        if( index > -1 ) {
          indelBasePosFromQ=n;
if(log.isInfoEnabled()) log.info("Setting indelBasePosFromQ: "+indelBasePosFromQ);
          break;
        }
      }
    }
  }
  public void finalHetMutationScan(SeqTraceAnalysis cTrace, SeqTraceAnalysis sTrace) {
// firstly test if we have a low crit drop - if so don't use this method
    if( critMutRatio < 0.70f || critMutRatio > 0.81f) return;

    int[] mutation = new int[5];
    String baseList[]={"A","C","G","T"};
    int i,j;
    int resPeak=0,resScan=0,resIntensity;
    String mutationBase=null;
    int mutationScore,curMutationScore;
    boolean mutPeakUnique=false;
    boolean contextChange=false; 
    char uniqueFlag=' ';
    char contextFlag=' ';
    int[] ok = new int[2];

    int[] sPeakHt = sTrace.getIntensity();
    int[] sScan = sTrace.getScan();
//  get array of quality values for tumour (based on peaks only)
    float qualCrit=5.0f;
    float qualCap=2.0f * qualCrit;
    int winForQual=5;
    float[] qual=sTrace.getTraceQualityUsingPeaksOnly(qualCap,winForQual);
    float totQual=0.0f;
//  work over portion of analysed sequence so can use windows

    float qCritForHet;    // depends if we're in an indel region
    int winLen=5;         // window half Length

// n.b. when searching for peaks method locatePeak uses but we must
// use winForQual so reduce using binFactor
    float binFactor=(float) winForQual / peakSearchBin;
    float maxMutRatio;

    for( i = winLen ; i < nPoints-winLen ; ++i ) {

// if we're past a suspected indel pos use larger max drop to try to
// capture as many hets as possible
      if( indelBasePosFromQ == 0 || ( indelBasePosFromQ > 0 && baseNo[i] < indelBasePosFromQ ) ) {
        maxMutRatio=critMutRatio+mutRatioIncrement;
      } else {
        maxMutRatio=critMutRatio+1.5f*mutRatioIncrement;
      }
      if( scale[i] >= critMutRatio && scale[i] < maxMutRatio ) {
// skip base if it or neighbours are Over/Underloaded or dye blobs
        if( sTrace.isBaseOverloaded(baseNo[i-1],28000) ||
            sTrace.isBaseOverloaded(baseNo[i],28000) ||
            sTrace.isBaseOverloaded(baseNo[i+1],28000) ) continue;
        if( sTrace.isBaseUnderloaded(baseNo[i-1]) ||
            sTrace.isBaseUnderloaded(baseNo[i]) ||
            sTrace.isBaseUnderloaded(baseNo[i+1]) ) continue;
        if( sTrace.isBaseDyeBlob(baseNo[i-1]) ||
            sTrace.isBaseDyeBlob(baseNo[i]) ||
            sTrace.isBaseDyeBlob(baseNo[i+1]) ) continue;
// skip base if any trace holes in window of +- winLen bases
        if( baseNo[i+winLen]-baseNo[i-winLen] > 2*winLen ) continue;
// check the actual qual of the current base
// n.b. first define quality cut off for hets
        if( indelBasePosFromQ == 0 || ( indelBasePosFromQ > 0 && baseNo[i] < indelBasePosFromQ ) ) {
          qCritForHet=3.25f;
        } else {
          qCritForHet=3.5f;
        }
        float upperLimForRelaxedQ=critMutRatio+0.5f*mutRatioIncrement;
        float qCritHetChannel;
        if( scale[i] >= critMutRatio && scale[i] < upperLimForRelaxedQ ) {
          qCritHetChannel=3.5f;
        } else {
          qCritHetChannel=4.0f;
        }
        int mutIndex=sTrace.getAnalysisIndex(baseNo[i-winLen]);
        if( qual[mutIndex+winLen] > qCritForHet ) continue;
// skip Quality control if we're past a suspected Indel position
        if( indelBasePosFromQ == 0 || ( indelBasePosFromQ > 0 && baseNo[i] < indelBasePosFromQ ) ) {
          totQual=0.0f;
          for( j = mutIndex ; j <= mutIndex+2*winLen ; ++j ) {
            if( j == mutIndex+winLen ) continue; // skip interrogation base
            totQual+=qual[j];
          }
          totQual/=(float) (2*winLen);
          if( totQual < qualCrit ) continue;  // reject on Av quality
// also insist that at least one neighbour must have q > qualCrit
// n.b. these might be larger than expected as qual only considers peaks
          if( qual[mutIndex+winLen-1] < qualCrit &&
              qual[mutIndex+winLen+1] < qualCrit ) continue;
        }
// now find change, i.e. type, intensity etc.
        //  reset mut data
        for( j = 0 ; j < mutation.length ; ++j ) {
          mutation[j] = 0;
        }
        resPeak=0;
        String curBase =String.valueOf(baseSeq[i]);
        String mut_base=null;
        int index=0;
        for( j = 0 ; j< 4 ; ++j ) {  //  cycle bases
          mut_base=baseList[j];

          if( baseSeq[i] == mut_base.charAt(0) ) continue;
          index=locatePeak(sTrace,sampleScan[i],mut_base,binFactor);
          if( index > 0 ) {
              resPeak= sPeakHt[index];
              resScan= sScan[index];
              ++mutation[4]; // No. of possible mutants
          }
          if( resPeak > mutation[2] ) {
            mutation[0] = index;
            mutation[1] = resScan;
            mutation[2] = resPeak;
            mutationBase= mut_base;
          }
        }
        if ( mutation[0] <= 0 ) continue;   // no mut peak found
// check viability of mut peak i.e. is it comparable to locality
        ok=checkViableMutant(sTrace,baseNo[i],mutation[0],mutationBase,false);
// if mut peak is viable check its validity
        mutationScore=0;
        if( ok[0] == 0 ) continue;   // mut peak not viable
        mutationScore=validateMutantPeak(cTrace,sTrace,i,mutation[0],mutationBase);
        if( mutationScore == 0 ) continue;   // mut peak low score
        curMutationScore=Math.abs(mutationScore);
        if( curMutationScore < 3 ) continue;   // use higher mut score
        mutation[3] = ok[1];
        if( mutationScore < 0 ) {
          mutPeakUnique=true;
          uniqueFlag='U';
        } else {
          mutPeakUnique=false;
          uniqueFlag='-';
        }
//  try to identify a contextual change in neighbourhood of mutation
        float contextInc=(1.0f-scale[i])/3.0f;
        if( i < nPoints-1 && Math.abs(1.0f-scale[i+1]) > contextInc ) {
          contextChange=true;
          contextFlag='C';
        } else {
          contextChange=false;
          contextFlag='-';
        }
        int mutDif=Math.abs(mutation[1]-sampleScan[i]);
// n.b. must set HET_SNP_PEAK before calling getLocalTrueChannelQuality()
        sTrace.setRealPeak(mutation[0],AutoCSA.HET_SNP_PEAK);
// check local channel wise true quality if borderline het
        if ( (indelBasePosFromQ > 0 && baseNo[i] < indelBasePosFromQ)
                  || indelBasePosFromQ == 0 ) {
          float trueQual=sTrace.getLocalTrueChannelQuality(baseNo[i],String.valueOf(mutationBase),3);
          if( trueQual < qCritHetChannel ) {  // cancel mutation
            sTrace.setRealPeak(mutation[0],AutoCSA.NOISE_PEAK);
            continue;
          }
        }
        int conf=75;
        mutations.set(i,mutationBase+"\t"+mutation[1]+"\t"+mutation[2]+"\t"+curMutationScore+"\t"+uniqueFlag+"\t"+contextFlag);
        // create a Mutation & MutationDetails object
        if( output != null ) {
          Mutation mutObj= new Mutation(baseNo[i],baseNo[i],Constants.ORA_ID_SUBSTITUTION,Constants.ORA_ID_MANUAL_REVIEW,getZygosityCode(scale[i]),String.valueOf(mutationBase),curBase);
          MutationDetails mutDetailsObj=new MutationDetails(mutation[1],0,controlScan[i],0,curMutationScore,scale[i],conf,mutPeakUnique,contextChange);
          output.addMutation(mutObj,mutDetailsObj);
        }
      }
    }
  }
  public void controlHetScan(SeqTraceAnalysis cTrace, SeqTraceAnalysis sTrace) {

    int[] mutation = new int[5];
    int[] mutation1 = new int[5];
    String baseList[]={"A","C","G","T"};
    int i,j;
    int resPeak=0,resScan=0,resIntensity;
    String mutationBase=null;
    int mutationScore,curMutationScore;
    boolean mutPeakUnique=false;
    boolean contextChange=false; 
    char uniqueFlag=' ';
    char contextFlag=' ';
    int[] ok = new int[2];

    int[] cPeakHt = cTrace.getIntensity();
    int[] cScan = cTrace.getScan();
    int[] sPeakHt = sTrace.getIntensity();
    int[] sScan = sTrace.getScan();
//  get array of quality values for control (based on peaks only)
    float qualCrit=5.0f;
    float qualCap=2.0f * qualCrit;
    int winForQual=5;
    float[] qual=cTrace.getTraceQualityUsingPeaksOnly(qualCap,winForQual);
    float totQual=0.0f;
//  work over portion of analysed sequence so can use windows

    float qCritForHet=3.0f;
    int winLen=5;         // window half Length
    boolean hetsInControl=true;  // set flag for het peak viability test

// n.b. when searching for peaks method locatePeak uses but we must
// use winForQual so reduce using binFactor
    float binFactor=(float) winForQual / peakSearchBin;

    for( i = winLen ; i < nPoints-winLen ; ++i ) {

// check the actual qual of the current base
        int mutIndex=cTrace.getAnalysisIndex(baseNo[i-winLen]);
        if( qual[mutIndex+winLen] > qCritForHet ) continue;

// skip base if it or neighbours are Over/Underloaded or dye blobs
        if( cTrace.isBaseOverloaded(baseNo[i-1],28000) ||
            cTrace.isBaseOverloaded(baseNo[i],28000) ||
            cTrace.isBaseOverloaded(baseNo[i+1],28000) ) continue;
        if( cTrace.isBaseUnderloaded(baseNo[i-1]) ||
            cTrace.isBaseUnderloaded(baseNo[i]) ||
            cTrace.isBaseUnderloaded(baseNo[i+1]) ) continue;
        if( cTrace.isBaseDyeBlob(baseNo[i-1]) ||
            cTrace.isBaseDyeBlob(baseNo[i]) ||
            cTrace.isBaseDyeBlob(baseNo[i+1]) ) continue;
// skip base if any trace holes in window of +- winLen bases
        if( baseNo[i+winLen]-baseNo[i-winLen] > 2*winLen ) continue;
        totQual=0.0f;
        for( j = mutIndex ; j <= mutIndex+2*winLen ; ++j ) {
          if( j == mutIndex+winLen ) continue; // skip interrogation base
          totQual+=qual[j];
        }
        totQual/=(float) (2*winLen);
        if( totQual < qualCrit ) continue;  // reject on Av quality
// also insist that at least one neighbour must have q > qualCrit
// n.b. these might be larger than expected as qual only considers peaks
//      if( qual[mutIndex+winLen-1] < qualCrit &&
//          qual[mutIndex+winLen+1] < qualCrit ) continue;
        int nbad=0;
        if( qual[mutIndex+winLen-2] < qualCrit ) ++nbad;
        if( qual[mutIndex+winLen-1] < qualCrit ) ++nbad;
        if( qual[mutIndex+winLen+1] < qualCrit ) ++nbad;
        if( qual[mutIndex+winLen+2] < qualCrit ) ++nbad;
        if( nbad > 1 ) continue;
//if(log.isInfoEnabled()) log.info("C-HET: "+baseNo[i]+" "+qual[mutIndex+winLen]+" "+totQual);
//      if(true) continue;
// now find change, i.e. type, intensity etc.
        //  reset mut data
        for( j = 0 ; j < mutation.length ; ++j ) {
          mutation[j] = 0;
          mutation1[j] = 0;
        }
        resPeak=0;
        String curBase =String.valueOf(baseSeq[i]);
        String mut_base=null;
        int index=0;
        for( j = 0 ; j< 4 ; ++j ) {  //  cycle bases
          mut_base=baseList[j];

          if( baseSeq[i] == mut_base.charAt(0) ) continue;
          index=locatePeak(cTrace,controlScan[i],mut_base,binFactor);
          if( index > 0 ) {
              resPeak= cPeakHt[index];
              resScan= cScan[index];
              ++mutation[4]; // No. of possible mutants
          }
          if( resPeak > mutation[2] ) {
            mutation[0] = index;
            mutation[1] = resScan;
            mutation[2] = resPeak;
            mutationBase= mut_base;
          }
        }
        if ( mutation[0] <= 0 ) continue;   // no mut peak found
// check viability of mut peak i.e. is it comparable to locality
        ok=checkViableMutant(cTrace,baseNo[i],mutation[0],mutationBase,hetsInControl);
// if mut peak is viable check its validity
        curMutationScore=0;
        if( ok[0] == 0 ) continue;   // mut peak not viable
        mutation[3] = ok[1];
        mutPeakUnique=false;
        uniqueFlag='-';
        contextChange=false;
        contextFlag='-';
        int mutDif=Math.abs(mutation[1]-controlScan[i]);
// n.b. must set HET_SNP_PEAK before calling getLocalTrueChannelQuality()
        cTrace.setRealPeak(mutation[0],AutoCSA.HET_SNP_PEAK);
// check local channel wise true quality if borderline het
        float trueQual=cTrace.getLocalTrueChannelQuality(baseNo[i],String.valueOf(mutationBase),5);
        if( trueQual < 4.0 ) {  // cancel mutation
          cTrace.setRealPeak(mutation[0],AutoCSA.NOISE_PEAK);
          continue;
        }
        int conf=75;
//if(log.isInfoEnabled()) log.info("HET IN CTRL at Bp: "+baseNo[i]+" ("+curBase+") Change:" +
//               mutationBase+" Scan: "+mutation[1]+" Ht: "+mutation[2]+
//               " Total: "+mutation[4]+" Score: "+curMutationScore+
//               " "+uniqueFlag+" "+contextFlag+" Drop: N/A"+
//               " MutPeak: "+mutation[3]+" MutDiff: "+mutDif+
//               " Confidence: "+conf);
        boolean referenceSub=false;
        boolean germlineSub=false;
        boolean hetSub=false;
        int type=0;
        String hetBase=null;
// n.b. if we're past an indel position in tumour don't look at tumour
        if ( (indelBasePosFromQ > 0 && baseNo[i] < indelBasePosFromQ)
                  || indelBasePosFromQ == 0 ) {
// check if a het at same base has already been found in the tumour
          if( mutations.get(i) == null ) {
// search for het in tumour - look for a 2nd viable het peak
            resPeak=0;
            index=0;
            for( j = 0 ; j< 4 ; ++j ) {  //  cycle bases
              mut_base=baseList[j];

              if( baseSeq[i] == mut_base.charAt(0) ) continue;
              index=locatePeak(sTrace,sampleScan[i],mut_base,binFactor);
              if( index > 0 ) {
                  resPeak= sPeakHt[index];
                  resScan= sScan[index];
                  ++mutation1[4]; // No. of possible mutants
              }
              if( resPeak > mutation1[2] ) {
                mutation1[0] = index;
                mutation1[1] = resScan;
                mutation1[2] = resPeak;
                hetBase= mut_base;
              }
            }
            if ( mutation1[0] > 0 ) {   // mut peak found
// check viability of mut peak i.e. is it comparable to locality
              ok=checkViableMutant(sTrace,baseNo[i],mutation1[0],hetBase,false);
              if( ok[0] == 0 ) hetBase=null; // mut peak not viable
              if( hetBase != null ) {
                mutation1[3] = ok[1];
                int mutDif1=Math.abs(mutation1[1]-sampleScan[i]);
// n.b. must set HET_SNP_PEAK before calling getLocalTrueChannelQuality()
                sTrace.setRealPeak(mutation1[0],AutoCSA.HET_SNP_PEAK);
// check local channel wise true quality if borderline het
                trueQual=sTrace.getLocalTrueChannelQuality(baseNo[i],hetBase,5);
                if( trueQual < 4.0 ) {  // cancel mutation
                  sTrace.setRealPeak(mutation1[0],AutoCSA.NOISE_PEAK);
                  hetBase=null;
                }
              }
            } else {
              hetBase=null;
            }
//if( hetBase != null ) if(log.isInfoEnabled()) log.info("HET IN TUMU at Bp: "+baseNo[i]+" ("+curBase+") Change:" +
//               hetBase+" Scan: "+mutation1[1]+" Ht: "+mutation1[2]+
//               " Total: "+mutation1[4]+" Score: "+curMutationScore+
//               " "+uniqueFlag+" "+contextFlag+" Drop: N/A"+
//               " MutPeak: "+mutation1[3]+" MutDiff: "+mutDif+
//               " Confidence: "+conf);
// test if we have same change (call only ref sub if nothing in tumour)
            if( hetBase == null ) {
              referenceSub=true;
            } else if( hetBase.equals(String.valueOf(mutationBase)) ) {
// call a germline sub but do not call reference sub
              germlineSub=true;
            } else {
// call reference sub and het sub (in tumour)
              referenceSub=true;
              hetSub=true;
              mutations.set(i,hetBase+"\t"+mutation1[1]+"\t"+mutation1[2]+"\t"+curMutationScore+"\t"+uniqueFlag+"\t"+contextFlag);
            }
          } else {
            String[] str= ((String) mutations.get(i)).split("\\t");
//if(log.isInfoEnabled()) log.info("HET already found: "+baseNo[i]+curBase+str[0]);
// test if we have same change
            if( str[0].equals(String.valueOf(mutationBase)) ) {
// alter het sub call to germline sub + do not call reference sub
              if( output != null) output.modifyMutationType(baseNo[i],Constants.ORA_ID_SUBSTITUTION,Constants.ORA_ID_GERMLINE_SUB);
            } else {
// het sub has already been reported so leave, report reference sub
              referenceSub=true;
            }
          }
        }
//if( referenceSub ) if(log.isInfoEnabled()) log.info("CALL: REF_SUB "+baseNo[i]);
//if( hetSub ) if(log.isInfoEnabled()) log.info("CALL: HET_SUB "+baseNo[i]);
//if( germlineSub ) if(log.isInfoEnabled()) log.info("CALL: GERM_SUB "+baseNo[i]);
        // create a Mutation & MutationDetails object
        if( hetSub && output != null ) {
          Mutation mutObj= new Mutation(baseNo[i],baseNo[i],Constants.ORA_ID_SUBSTITUTION,Constants.ORA_ID_MANUAL_REVIEW,Constants.ORA_ID_HETEROZYGOUS,String.valueOf(hetBase),curBase);
          MutationDetails mutDetailsObj=new MutationDetails(mutation1[1],0,controlScan[i],0,curMutationScore,0.0f,conf,mutPeakUnique,contextChange);
          output.addMutation(mutObj,mutDetailsObj);
        }
        if( referenceSub && output != null ) {
          Mutation mutObj= new Mutation(baseNo[i],baseNo[i],Constants.ORA_ID_REFERENCE_SUB,Constants.ORA_ID_MANUAL_REVIEW,Constants.ORA_ID_HETEROZYGOUS,String.valueOf(mutationBase),curBase);
          MutationDetails mutDetailsObj=new MutationDetails(mutation[1],0,controlScan[i],0,curMutationScore,0.0f,conf,mutPeakUnique,contextChange);
          output.addMutation(mutObj,mutDetailsObj);
        }
        if( germlineSub && output != null ) {
          Mutation mutObj= new Mutation(baseNo[i],baseNo[i],Constants.ORA_ID_GERMLINE_SUB,Constants.ORA_ID_MANUAL_REVIEW,Constants.ORA_ID_HETEROZYGOUS,String.valueOf(mutationBase),curBase);
          MutationDetails mutDetailsObj=new MutationDetails(mutation1[1],0,mutation[1],0,curMutationScore,0.0f,conf,mutPeakUnique,contextChange);
          output.addMutation(mutObj,mutDetailsObj);
        }
      }
  }
  private int validateComplexIndelCall(int delPos, String delBases, String insBases, SeqTraceAnalysis sTrace) {
// first check if delBases = insBases n.b. this is possible as peak
// matching algorithm is not the same as the insert finder algorithm
    if( delBases.equals(insBases) ) return -1;
    int delLen=delBases.length();
    int insLen=insBases.length();
// first check if insBases are part of delBases
    if( delLen > 2 && insLen > 2 && (delLen+2) >= insLen ) {
      String newStr=refSeq.substring(delPos-2,delPos-1).concat(delBases).concat(refSeq.substring(delPos-1+delLen,delPos+delLen));
      if( newStr.indexOf(insBases) > -1 ) return 1;
    }
// now check cases of 1 base ins - do we have too many peaks in scan gap
    if( insLen == 1 ) {
      int num=sTrace.getTotalAvailablePeaksInRange(delPos-1,delPos+delLen);
      if ( num > 3 ) return 1;
    }
    return 0;
  }
  private void checkComplexForSimplification(SeqTraceAnalysis cTrace, SeqTraceAnalysis sTrace, int delPos,int len,String insBases,String delBases,int cScan1,int cScan2,int sScan1,int sScan2) {
// this method tests to see if a Complex Indel is really one of;
// sub, ins, del or a del and a sub (only 1 is possible) else
// retains the original compex call
// N.B This can arise as matching algorithm and variant detection
// algorithms use different principles and have different stringencies
// initialise all booleans to false n.b. it is only possible
// for 1 ONLY of the booleans to become true
    boolean delAndSub=false;
    boolean onlySub=false;
    boolean onlyDel=false;
    boolean onlyIns=false;
    int delLen=delBases.length();
    int insLen=insBases.length();
// check if insBases matches at start or end of delBases
    int iMatch=-1;
    if( delLen > insLen ) {
      if( delBases.indexOf(insBases) == 0 ) iMatch=0;
      if( delBases.lastIndexOf(insBases) == delLen-insLen ) iMatch=delLen-insLen;
      if( iMatch > -1 ) onlyDel=true;
    }
// remove final base from both insBases and delBases and test if
// insBases matches on end of delBases, if so we really have a
// deletion of len (len-insLen) followed by insLen-1 wild type bases
// lastly followed by a hom sub
    if( !onlyDel && delLen > insLen && insLen > 1 ) {
      String del=delBases.substring(0,delLen-1);
      String ins=insBases.substring(0,insLen-1);
      int ind=del.lastIndexOf(ins);
      if( ind == del.length()-ins.length() ) delAndSub=true;
    }
// check if delBases matches at start or end of insBases
    if( insLen > delLen ) {
      if( insBases.indexOf(delBases) == 0 ) iMatch=0;
      if( insBases.lastIndexOf(delBases) == insLen-delLen ) iMatch=insLen-delLen;
      if( iMatch > -1 ) onlyIns=true;
    }
//  check if only first or final base differ - really hom sub
    int subPos=0;
    if( delLen == insLen && insLen > 1 ) {
      if( delBases.substring(1,delLen).equals(insBases.substring(1,delLen)) ) subPos=delPos;
      if( delBases.substring(0,delLen-1).equals(insBases.substring(0,delLen-1)) ) subPos=delPos+delLen-1;
      if( subPos > 0 ) onlySub=true;
    }
    if( delAndSub ) {
// save the deletion details
      int realDelLen=delLen-insLen;
      Mutation mutObj= new Mutation(delPos,delPos+realDelLen-1,Constants.ORA_ID_DELETION,Constants.ORA_ID_MANUAL_REVIEW,Constants.ORA_ID_HOMOZYGOUS,"-",delBases.substring(0,realDelLen));
      int newcScan2=cTrace.getScanAtBase(delPos+realDelLen);
      int[] scans=sTrace.getScansOfInsertedBases(insBases,sScan1);
      int newsScan2=scans[0];
      MutationDetails mutDetailsObj=new MutationDetails(sScan1,newsScan2,cScan1,newcScan2,90,0.0f,0,false,false);
      output.addMutation(mutObj,mutDetailsObj);
// n.b. hom subs will NOT have been called yet so check status of 
// corresponding base on normal trace - so can decide on call type
// n.b. if also missing in normal can't change status to missing in both 
// as detection code cannot handle the tumour as bases missed upstream
// also only germline sub if both changes are same - don't know this here
// -> for safety will have to call ref sub (not here) and hom mut (here)
// save the substitution details
      subPos=delPos+delLen-1;
      String wt=delBases.substring(delLen-1,delLen);
      String tt=insBases.substring(insLen-1,insLen);
      mutObj= new Mutation(subPos,subPos,Constants.ORA_ID_SUBSTITUTION,Constants.ORA_ID_MANUAL_REVIEW,Constants.ORA_ID_HOMOZYGOUS,tt,wt);
      int c1=cTrace.getScanAtBase(delPos+delLen-1);
      int s1=scans[insLen-1];
      mutDetailsObj=new MutationDetails(s1,0,c1,0,90,0.0f,0,false,false);
      output.addMutation(mutObj,mutDetailsObj);
    } else if( onlyDel ) {
      int realDelLen=delLen-insLen;
      int realDelPos=0;
      String realDelBases;
      if( iMatch == 0 ) {
        realDelPos=delPos+insLen;
        realDelBases=delBases.substring(insLen-1,insLen-1+realDelLen);
      } else {
        realDelPos=delPos;
        realDelBases=delBases.substring(0,realDelLen);
      }
      Mutation mutObj= new Mutation(realDelPos,realDelPos+realDelLen-1,Constants.ORA_ID_DELETION,Constants.ORA_ID_MANUAL_REVIEW,Constants.ORA_ID_HOMOZYGOUS,"-",realDelBases);
      int newcScan1,newcScan2,newsScan1,newsScan2;
      int[] scans=sTrace.getScansOfInsertedBases(insBases,sScan1);
      if( iMatch == 0 ) {
        newcScan1=cTrace.getScanAtBase(delPos+insLen-1);
        newcScan2=cScan2;
        newsScan1=scans[insLen-1];
        newsScan2=sScan2;
      } else {
        newcScan1=cScan1;
        newcScan2=cTrace.getScanAtBase(delPos+realDelLen);
        newsScan1=sScan1;
        newsScan2=scans[0];
      }
      MutationDetails mutDetailsObj=new MutationDetails(newsScan1,newsScan2,newcScan1,newcScan2,90,0.0f,0,false,false);
      output.addMutation(mutObj,mutDetailsObj);
    } else if( onlyIns ) {
      int realInsLen=insLen-delLen;
      int realInsPos=0;
      String realInsBases;
      if( iMatch == 0 ) {
        realInsPos=delPos+delLen-1;
        realInsBases=insBases.substring(delLen,delLen+realInsLen);
      } else {
        realInsPos=delPos-1;
        realInsBases=insBases.substring(0,realInsLen);
      }
      Mutation mutObj= new Mutation(realInsPos,realInsPos+1,Constants.ORA_ID_INSERTION,Constants.ORA_ID_MANUAL_REVIEW,Constants.ORA_ID_HOMOZYGOUS,realInsBases,"-");
      int newcScan1,newcScan2,newsScan1,newsScan2;
      int[] scans=sTrace.getScansOfInsertedBases(insBases,sScan1);
      if( iMatch == 0 ) {
        newcScan1=cTrace.getScanAtBase(delPos+delLen);
        newcScan2=cScan2;
        newsScan1=scans[delLen-1];
        newsScan2=sScan2;
      } else {
        newcScan1=cScan1;
        newcScan2=cTrace.getScanAtBase(delPos);
        newsScan1=sScan1;
        newsScan2=scans[iMatch];
      }
      MutationDetails mutDetailsObj=new MutationDetails(newsScan1,newsScan2,newcScan1,newcScan2,90,0.0f,0,false,false);
      output.addMutation(mutObj,mutDetailsObj);
    } else if( onlySub ) {
      int ind=0;
      if( subPos == delPos ) { ind=0; } else { ind=delLen-1; }
      String wt=delBases.substring(ind,ind+1);
      String tt=insBases.substring(ind,ind+1);
      Mutation mutObj= new Mutation(subPos,subPos,Constants.ORA_ID_SUBSTITUTION,Constants.ORA_ID_MANUAL_REVIEW,Constants.ORA_ID_HOMOZYGOUS,tt,wt);
      int[] scans=sTrace.getScansOfInsertedBases(insBases,sScan1);
      int c1=cTrace.getScanAtBase(delPos+ind);
      int s1=scans[ind];
      MutationDetails mutDetailsObj=new MutationDetails(s1,0,c1,0,90,0.0f,0,false,false);
      output.addMutation(mutObj,mutDetailsObj);
    } else {
// save the complex details
      Mutation mutObj= new Mutation(delPos,delPos+len-1,Constants.ORA_ID_COMPLEX,Constants.ORA_ID_MANUAL_REVIEW,Constants.ORA_ID_HOMOZYGOUS,insBases,delBases);
      MutationDetails mutDetailsObj=new MutationDetails(sScan1,sScan2,cScan1,cScan2,90,0.0f,0,false,false);
      output.addMutation(mutObj,mutDetailsObj);
    }
  }
}
