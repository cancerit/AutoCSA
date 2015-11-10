package uk.ac.sanger.cgp.autocsa.util;

/**
 *<p> Class provides a set of public static final constants for use throughout autoCSA.</p>
 *
 *Original author:  emd
 *@author $Author$
 *@version $Revision$
 */
public class Constants {

/**
 * null constructor.
 */
  public Constants () {
  }

// trace data parameters
/**
 * STRICT_MIN_PEAK - do not recognise any peaks below this intensity
 */
  public static final int STRICT_MIN_PEAK=250;

/**
 * AVERAGE_BASE_SPACING - nominal value for average base spacing
 */
  public static final int AVERAGE_BASE_SPACING=12;

/**
 * INTENSITY_OVERLOADED - intensity at which peak is deemed Overloaded
 */
  public static final int INTENSITY_OVERLOADED=30000;

/**
 * INTENSITY_UNDERLOADED - intensity at which peak is deemed Underloaded
 */
  public static final int INTENSITY_UNDERLOADED=600;
/**
 * SCAN_SUBTRACTOR_FOR_PEAKS - CSA subtracts this length of trace in
 * scans from both ends to form a range for peak threshold calculations,
 * i.e. this defines a mid-portion of a trace and avoids erroneous
 * scalings at the ends of a trace
 */
  public static final int SCAN_SUBTRACTOR_FOR_PEAKS=1000;
/**
 * CRIT_UNDERLOAD_PEAK_FOR_HET - intensity at which we search below the
 * usual peak intensity threshold for het peaks
 */
  public static final int CRIT_UNDERLOAD_PEAK_FOR_HET=800;
/**
 * CUTOFF_UNDERLOAD_DECREASE - max intensity (over range used to determine
 * peak minima) at which we lower thresholds to recover Underloaded peaks
 */
  public static final int CUTOFF_UNDERLOAD_DECREASE=0;
/**
 * UL_DECREASE_FACTOR - factor by which we lower peak thresholds in order
 * to recover Underloaded peaks (switched on by CUTOFF_UNDERLOAD_DECREASE)
 */
  public static final float UL_DECREASE_FACTOR=2.0f/3.0f;

/**
 * DYE_BLOB_WINDOW_START1 - base No of search region for dye blobs (1)
 */
  public static final int DYE_BLOB_WINDOW_START1=97;

/**
 * DYE_BLOB_WINDOW_START2 - base No of search region for dye blobs (2)
 */
  public static final int DYE_BLOB_WINDOW_START2=137;

/**
 * DYE_BLOB_WINDOW_WIDTH - width in bases of search region for dye blobs
 */
  public static final int DYE_BLOB_WINDOW_WIDTH=16;

/**
 * INPUT_BASE_ORDERING - ordering of data Channels in input file (SCF)
 */
  public static final String INPUT_BASE_ORDERING="GATC";

/**
 * A_CHANNEL_PEAK_RATIO - factor to convert A max intensity to min peak I
 */
  public static final float A_CHANNEL_PEAK_RATIO=0.04f;

/**
 * C_CHANNEL_PEAK_RATIO - factor to convert C max intensity to min peak I
 */
  public static final float C_CHANNEL_PEAK_RATIO=0.04f;

/**
 * G_CHANNEL_PEAK_RATIO - factor to convert G max intensity to min peak I
 */
  public static final float G_CHANNEL_PEAK_RATIO=0.03f;

/**
 * T_CHANNEL_PEAK_RATIO - factor to convert T max intensity to min peak I
 */
  public static final float T_CHANNEL_PEAK_RATIO=0.04f;

/**
 * A_CHANNEL_SHOULDER_RATIO - fraction for min shoulder intensity for A
 */
  public static final float A_CHANNEL_SHOULDER_RATIO=0.20f;
/**
 * A_CHANNEL_SHOULDER_RATIO - fraction for min shoulder intensity for C
 */
  public static final float C_CHANNEL_SHOULDER_RATIO=0.20f;
/**
 * A_CHANNEL_SHOULDER_RATIO - fraction for min shoulder intensity for G
 */
  public static final float G_CHANNEL_SHOULDER_RATIO=0.20f;
/**
 * A_CHANNEL_SHOULDER_RATIO - fraction for min shoulder intensity for T
 */
  public static final float T_CHANNEL_SHOULDER_RATIO=0.20f;
/**
 * SHOULDER_RELAXATION - fraction for reduction in min shoulder intensity
 * n.b. this is applied in the latter part of the trace (after scan=3000)
 */
  public static final float SHOULDER_RELAXATION=0.90f;
/**
 * SCAN_SHOULDER_RELAXATION - scan to start applying SHOULDER_RELAXATION
 */
  public static final int SCAN_SHOULDER_RELAXATION=3000;
// trace analysis/comparison parameters
/**
 * QUALITY_WINDOW_SIZE - width of window (bases) for quality calculation
 */
  public static final int QUALITY_WINDOW_SIZE=5;

/**
 * MIN_ALLOWABLE_BASES - min No of bases required for a comparison
 */
  public static final int MIN_ALLOWABLE_BASES=10;

/**
 * CRIT_OVERLOAD_PC_NORMAL - % limit for NORMAL traces deemed Overloaded
 */
  public static final float CRIT_OVERLOAD_PC_NORMAL=10.0f;

/**
 * CRIT_OVERLOAD_PC_TUMOUR - % limit for TUMOUR traces deemed Overloaded
 */
  public static final float CRIT_OVERLOAD_PC_TUMOUR=20.0f;

/**
 * CRIT_UNDERLOAD_PC_NORMAL - % limit for NORMAL traces deemed Underloaded
 */
  public static final float CRIT_UNDERLOAD_PC_NORMAL=10.0f;

/**
 * CRIT_UNDERLOAD_PC_TUMOUR - % limit for TUMOUR traces deemed Underloaded
 */
  public static final float CRIT_UNDERLOAD_PC_TUMOUR=20.0f;

/**
 * MAX_POSSIBLE_HET_DEL - max possible het DEL characterised by CSA
 */
  public static final int MAX_POSSIBLE_HET_DEL=50;

/**
 * MAX_POSSIBLE_HET_INS - max possible het INS characterised by CSA
 */
  public static final int MAX_POSSIBLE_HET_INS=50;

/**
 * MAX_POSSIBLE_HOM_INS - max possible hom Insertion detected by CSA
 */
  public static final int MAX_POSSIBLE_HOM_INS=25;

/**
 * MAX_HOM_COMPLEX_INS - max possible hom Insertion portion of a
 * Complex detected by CSA
 */
  public static final int MAX_HOM_COMPLEX_INS=10;

// the following are DNA type flags
  public static final int DNA_NORMAL1=0;
  public static final int DNA_NORMAL2=5;
  public static final int DNA_TUMOUR1=1;
  public static final int DNA_TUMOUR2=3;
  public static final int DNA_MRX=4;

// N.B. The remaining constants are all Oracle ID's
// Oracle ID_COMPSTATUS values
  public static final int ORA_ID_DONE_NON_MUTANT=3;
  public static final int ORA_ID_DONE_MUTANT=5;
  public static final int ORA_ID_FAILED_IN_CSA=6;

// Oracle ID_TRACE_STATUS values
  public static final int ORA_ID_RERUN_ANALYSIS=1;
  public static final int ORA_ID_RERUN_UNDERLOADED=6;
  public static final int ORA_ID_RERUN_OVERLOADED=7;
  public static final int ORA_ID_RERUN_MATCH_FAILED=8;
  public static final int ORA_ID_RERUN_NO_ROR_COVERED=10;
  public static final int ORA_ID_RERUN_BAD_TRACE=11;
  public static final int ORA_ID_SUSPECT_INDEL_IN_ROR=12;

// Oracle ID_MUTSTATUS values
  public static final int ORA_ID_MANUAL_REVIEW=3;

// Oracle ID_MUTTYPE values
  public static final int ORA_ID_SUBSTITUTION=1;
  public static final int ORA_ID_INSERTION=2;
  public static final int ORA_ID_DELETION=3;
  public static final int ORA_ID_COMPLEX=4;
  public static final int ORA_ID_REFERENCE_SUB=5;
  public static final int ORA_ID_SPECULATIVE_INDEL=6;
  public static final int ORA_ID_TRACE_HOLE=7;
  public static final int ORA_ID_GERMLINE_SUB=10;

// Oracle ID_ZYGOSITY values
  public static final int ORA_ID_HOMOZYGOUS=1;
  public static final int ORA_ID_HETEROZYGOUS=2;
  public static final int ORA_ID_AMPLIFICATION=3;

}
